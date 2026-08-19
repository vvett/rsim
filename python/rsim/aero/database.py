"""Full-envelope aerodynamic method routing and N-dimensional table output."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from ._mat import write_mat_v5
from .geometry import RocketGeometry
from .grid import AeroGrid
from .methods import (
    BodyGeometryData,
    axisymmetric_body_pressure,
    base_drag_coefficient,
    boattail_pressure_plus_base_drag,
    fin_cp,
    fin_normal_force,
    fin_wave_drag,
    induced_fin_drag,
    jorgensen_body_force,
    preprocess_body,
    roll_damping_derivative,
    skin_friction_coefficient,
    tangent_ogive_wave_drag,
)


FORMAT_VERSION = 2


@dataclass(frozen=True)
class AerodynamicOptions:
    """Configuration choices not represented by database dimensions.

    ``boundary_layer`` is deliberately explicit because the guide forbids a
    silent transition assumption. ``high_mach_switch`` selects Modified
    Newtonian; the guide recommends approximately M>=5 only after validation.
    ``span_efficiency`` is the FDRAG-05 fallback input and should be calibrated.
    """

    boundary_layer: str = "turbulent"
    gamma: float = 1.4
    high_mach_switch: float = 5.0
    span_efficiency: float = 0.9
    body_station_count: int = 121
    body_azimuth_count: int = 48

    def __post_init__(self) -> None:
        if self.boundary_layer not in {"laminar", "turbulent"}:
            raise ValueError("boundary_layer must be 'laminar' or 'turbulent'")
        if self.gamma <= 1.0:
            raise ValueError("gamma must be greater than one")
        if self.high_mach_switch < 5.0:
            raise ValueError("the guide does not recommend Modified Newtonian below Mach 5")
        if not 0.0 < self.span_efficiency <= 1.0:
            raise ValueError("span_efficiency must lie in (0, 1]")
        if self.body_station_count < 21:
            raise ValueError("body_station_count must be at least 21")
        if self.body_azimuth_count < 8:
            raise ValueError("body_azimuth_count must be at least 8")


@dataclass
class AeroCoefficientTable:
    """Regular interpolation table with NumPy ``ij`` mesh ordering.

    Arrays use dimension order ``(mach, alpha_deg, phi_aero_deg, reynolds)``. This is the
    NumPy equivalent of MATLAB ``ndgrid`` and can be passed to regular-grid
    interpolation after selecting the corresponding axes. Components and
    method flags have the same shape as final coefficient arrays.
    """

    axes: dict[str, np.ndarray]
    meshes: dict[str, np.ndarray]
    coefficients: dict[str, np.ndarray]
    components: dict[str, np.ndarray]
    method_flags: dict[str, np.ndarray]
    metadata: dict[str, object]

    def save_npz(self, path: str | Path) -> None:
        """Save plain named arrays and JSON metadata without Python pickles.

        The stable prefixes ``axis__``, ``mesh__``, ``coefficient__``,
        ``component__``, and ``method__`` are intentionally straightforward
        for a future C++ NPZ reader. Existing files are overwritten by NumPy.
        """

        payload: dict[str, np.ndarray] = {
            "format_version": np.array(FORMAT_VERSION, dtype=np.int64),
            "metadata_json": np.array(json.dumps(self.metadata, sort_keys=True)),
        }
        for prefix, values in (
            ("axis", self.axes),
            ("mesh", self.meshes),
            ("coefficient", self.coefficients),
            ("component", self.components),
            ("method", self.method_flags),
        ):
            payload.update({f"{prefix}__{name}": value for name, value in values.items()})
        np.savez_compressed(path, **payload)

    def save_mat(self, path: str | Path) -> None:
        """Write a MATLAB level-5 database with ndgrid-compatible arrays.

        Axes, meshes, coefficients, and numeric components retain the documented
        dimension order. String method flags are stored as uint32 arrays, and
        ``metadata_json_utf8`` contains metadata plus each flag-code dictionary.
        This uses the package's small level-5 writer and does not require SciPy.
        """

        arrays: dict[str, np.ndarray] = {
            "format_version": np.array([[FORMAT_VERSION]], dtype=np.uint32),
        }
        for prefix, values in (
            ("axis", self.axes),
            ("mesh", self.meshes),
            ("coefficient", self.coefficients),
            ("component", self.components),
        ):
            arrays.update({f"{prefix}__{name}": value for name, value in values.items()})

        metadata = dict(self.metadata)
        method_codebooks: dict[str, dict[str, int]] = {}
        for name, flags in self.method_flags.items():
            labels, encoded = np.unique(flags, return_inverse=True)
            # Codes start at one so zero remains available as an unset sentinel.
            codes = (encoded.reshape(flags.shape) + 1).astype(np.uint32)
            arrays[f"method__{name}"] = codes
            method_codebooks[name] = {
                str(label): index + 1 for index, label in enumerate(labels)
            }
        metadata["mat_method_flag_codes"] = method_codebooks
        encoded_metadata = json.dumps(metadata, sort_keys=True).encode("utf-8")
        arrays["metadata_json_utf8"] = np.frombuffer(
            encoded_metadata, dtype=np.uint8
        ).reshape(1, -1)
        write_mat_v5(path, arrays)

    @classmethod
    def load_npz(cls, path: str | Path) -> "AeroCoefficientTable":
        """Load an NPZ produced by :meth:`save_npz` without enabling pickle."""

        groups: dict[str, dict[str, np.ndarray]] = {
            "axis": {},
            "mesh": {},
            "coefficient": {},
            "component": {},
            "method": {},
        }
        with np.load(path, allow_pickle=False) as archive:
            version = int(archive["format_version"])
            if version != FORMAT_VERSION:
                raise ValueError(f"unsupported aerodynamic table format {version}")
            metadata = json.loads(str(archive["metadata_json"]))
            for key in archive.files:
                if "__" not in key:
                    continue
                prefix, name = key.split("__", 1)
                if prefix in groups:
                    groups[prefix][name] = np.array(archive[key], copy=True)
        return cls(
            axes=groups["axis"],
            meshes=groups["mesh"],
            coefficients=groups["coefficient"],
            components=groups["component"],
            method_flags=groups["method"],
            metadata=metadata,
        )


@dataclass(frozen=True)
class _InviscidPoint:
    ca_pressure: float
    ca_base: float
    ca_fin_wave: float
    ca_induced: float
    cn_body: float
    cn_fin: float
    cm_body: float
    cm_fin: float
    fin_cn_alpha: float
    fin_cp_x: float
    cl_p: float
    flags: Mapping[str, str]


@dataclass(frozen=True)
class _RollAwareFinPoint:
    """Individual-fin synthesis resolved into body transverse axes."""

    cy: float
    cz: float
    cmx: float
    cmy: float
    cmz: float
    ca_induced: float
    cn_along_incidence: float
    cn_cross_incidence: float
    cm_along_incidence: float
    cp_x: float
    slope: float
    normal_method: str
    cp_method: str


def _roll_aware_fins(
    mach: float,
    alpha: float,
    phi_aero: float,
    geometry: RocketGeometry,
    options: AerodynamicOptions,
) -> _RollAwareFinPoint:
    """Resolve each fin using the implementation note's low-angle projection.

    Each equivalent incidence is ``alpha*sin(phi_aero-phi_fin)``. The factor
    ``2/N`` makes the small-angle force along the incidence plane match the
    existing aggregated fin-set model for an evenly spaced symmetric set.
    This is a low-angle equivalent-incidence synthesis, not a replacement for
    complete-vehicle CFD or wind-tunnel roll sweeps.
    """

    count = geometry.fins.count
    scale = 2.0 / count
    radial_cp = 0.5 * geometry.reference_diameter + 0.5 * geometry.fins.semi_span

    def resolved(phi: float) -> tuple[np.ndarray, np.ndarray, float, list[float], list[float], str]:
        transverse_force = np.zeros(2)
        transverse_moment = np.zeros(2)
        roll_moment = 0.0
        individual_cn: list[float] = []
        individual_cp: list[float] = []
        normal_method = ""
        for fin_phi in geometry.fins.azimuths_rad:
            # The fin-plane normal is tangent to the body circumference.
            normal = np.array([-np.sin(fin_phi), np.cos(fin_phi)])
            radial = np.array([np.cos(fin_phi), np.sin(fin_phi)])
            alpha_equivalent = alpha * np.sin(phi - fin_phi)
            set_cn, set_slope, normal_method = fin_normal_force(
                mach, float(alpha_equivalent), geometry
            )
            fin_cn = scale * set_cn
            cp_x = fin_cp(mach, float(alpha_equivalent), geometry.fins)
            fin_cm = -fin_cn * (
                cp_x - geometry.moment_reference_x
            ) / geometry.length
            transverse_force += fin_cn * normal
            # Rotate the legacy normal-plane moment convention into body Y/Z.
            transverse_moment += fin_cm * np.array([-normal[1], normal[0]])
            roll_moment += radial_cp / geometry.length * fin_cn * (
                radial[0] * normal[1] - radial[1] * normal[0]
            )
            individual_cn.append(fin_cn)
            individual_cp.append(cp_x)
        return (
            transverse_force,
            transverse_moment,
            roll_moment,
            individual_cn,
            individual_cp,
            normal_method,
        )

    force, moment, cmx, fin_cn_values, fin_cp_values, normal_method = resolved(
        phi_aero
    )
    incidence = np.array([np.cos(phi_aero), np.sin(phi_aero)])
    cross_incidence = np.array([-np.sin(phi_aero), np.cos(phi_aero)])
    cn_along = float(np.dot(force, incidence))
    cn_cross = float(np.dot(force, cross_incidence))
    moment_direction = np.array([-np.sin(phi_aero), np.cos(phi_aero)])
    cm_along = float(np.dot(moment, moment_direction))
    cp_x = (
        geometry.moment_reference_x - geometry.length * cm_along / cn_along
        if abs(cn_along) > 1.0e-12
        else np.nan
    )

    # Scale the existing induced-drag method by individual-fin loading squared.
    aggregate_cn, aggregate_slope, _ = fin_normal_force(mach, alpha, geometry)
    baseline_induced = induced_fin_drag(
        mach, aggregate_cn, geometry.fins,
        geometry.reference_area, options.span_efficiency,
    )
    _, _, _, reference_cn, _, _ = resolved(0.0)
    denominator = sum(value * value for value in reference_cn)
    numerator = sum(value * value for value in fin_cn_values)
    ca_induced = baseline_induced * numerator / denominator if denominator > 1.0e-20 else 0.0
    cp_method = (
        "FALLBACK_BARF_04_NO_ALPHA_MIGRATION"
        if geometry.fins.taper_ratio <= 0.5
        else "CHART_AIAA07_INDIVIDUAL_FIN_CP"
    )
    return _RollAwareFinPoint(
        cy=float(force[0]),
        cz=float(force[1]),
        cmx=float(cmx),
        cmy=float(moment[0]),
        cmz=float(moment[1]),
        ca_induced=float(ca_induced),
        cn_along_incidence=cn_along,
        cn_cross_incidence=cn_cross,
        cm_along_incidence=cm_along,
        cp_x=float(cp_x),
        slope=float(aggregate_slope),
        normal_method=normal_method + "+INDIVIDUAL_FIN_PHI_SYNTHESIS",
        cp_method=cp_method,
    )


def _body_pressure_route(
    mach: float,
    alpha: float,
    geometry: RocketGeometry,
    body: BodyGeometryData,
    options: AerodynamicOptions,
) -> tuple[float, float, float, str]:
    """Route the guide's low/transonic, SOSE, and Modified-Newtonian body methods."""

    if mach < 1.2:
        body_cn, body_cm, _, body_cn_viscous = jorgensen_body_force(
            mach, alpha, geometry, body
        )
        if mach < 0.85:
            ca_pressure = 0.0
            pressure_flag = "ASSEMBLY_ZERO_BELOW_NWL_RANGE"
        elif mach <= 1.15:
            ca_pressure = tangent_ogive_wave_drag(
                mach, geometry.nose.fineness_ratio
            )
            pressure_flag = "CHART_NWL_TR2796_FIG3"
        else:
            chart_endpoint = tangent_ogive_wave_drag(
                1.15, geometry.nose.fineness_ratio
            )
            supersonic = axisymmetric_body_pressure(
                1.2,
                alpha,
                geometry,
                body,
                gamma=options.gamma,
                method="shock_expansion",
                azimuth_count=options.body_azimuth_count,
            )
            fraction = (mach - 1.15) / 0.05
            ca_pressure = (1.0 - fraction) * chart_endpoint + fraction * supersonic.ca
            pressure_flag = "FALLBACK_NWL_TO_SOSE_FAIRING"
        return ca_pressure, body_cn, body_cm, pressure_flag

    _, _, _, body_cn_viscous = jorgensen_body_force(mach, alpha, geometry, body)
    viscous_cm = (
        -body_cn_viscous
        * (body.planform_centroid_x - geometry.moment_reference_x)
        / geometry.length
    )
    if mach >= options.high_mach_switch:
        pressure = axisymmetric_body_pressure(
            mach,
            alpha,
            geometry,
            body,
            gamma=options.gamma,
            method="modified_newtonian",
            azimuth_count=options.body_azimuth_count,
        )
        flag = "FALLBACK_MODIFIED_NEWTONIAN_PANELS"
    else:
        pressure = axisymmetric_body_pressure(
            mach,
            alpha,
            geometry,
            body,
            gamma=options.gamma,
            method="shock_expansion",
            azimuth_count=options.body_azimuth_count,
        )
        flag = (
            "FALLBACK_SOSE_WITH_DETACHED_MN_PANELS"
            if pressure.used_detached_fallback
            else "FALLBACK_SOSE_LOCAL_PANELS"
        )
    return (
        pressure.ca,
        pressure.cn + body_cn_viscous,
        pressure.cm + viscous_cm,
        flag,
    )


def _base_drag_route(mach: float, geometry: RocketGeometry) -> tuple[float, str]:
    """Route power-off cylindrical or conical-afterbody base drag."""

    chart_mach = min(mach, 6.0)
    body_area = geometry.reference_area
    full_base = base_drag_coefficient(chart_mach, 1.0)
    if geometry.boattail and mach >= 1.0:
        body_radius = 0.5 * geometry.reference_diameter
        exit_radius = 0.5 * geometry.boattail.exit_diameter
        slope = np.degrees(
            np.arctan2(body_radius - exit_radius, geometry.boattail.length)
        )
        value = boattail_pressure_plus_base_drag(
            mach,
            float(slope),
            exit_radius / body_radius,
            full_base,
        )
        flag = "CHART_BASE_01_PLUS_EXACT_BASE_02"
    else:
        base_area_ratio = (
            np.pi * geometry.base_diameter**2 / 4.0 / body_area
        )
        value = base_drag_coefficient(chart_mach, float(base_area_ratio))
        flag = "CHART_NSWC_FIG51"
    if mach > 6.0:
        flag += "_HELD_ABOVE_M6"
    return value, flag


def _inviscid_point(
    mach: float,
    alpha: float,
    geometry: RocketGeometry,
    body: BodyGeometryData,
    options: AerodynamicOptions,
) -> _InviscidPoint:
    ca_pressure, cn_body, cm_body, pressure_flag = _body_pressure_route(
        mach, alpha, geometry, body, options
    )
    cn_fin, fin_slope, fin_method = fin_normal_force(mach, alpha, geometry)
    x_fin = fin_cp(mach, alpha, geometry.fins)
    cm_fin = (
        -cn_fin
        * (x_fin - geometry.moment_reference_x)
        / geometry.length
    )
    wave, wave_method = fin_wave_drag(mach, geometry.fins, geometry.reference_area)
    induced = induced_fin_drag(
        mach,
        cn_fin,
        geometry.fins,
        geometry.reference_area,
        options.span_efficiency,
    )
    base, base_method = _base_drag_route(mach, geometry)
    if geometry.fins.taper_ratio <= 0.5:
        fin_cp_method = "FALLBACK_BARF_04_NO_ALPHA_MIGRATION"
    elif mach >= 2.0:
        fin_cp_method = "CHART_AIAA07_FCP_02_03"
    else:
        fin_cp_method = "CHART_AIAA07_FCP_02"
    induced_method = (
        "EXACT_FDRAG_05" if mach <= 1.4 else "ASSEMBLY_ZERO_ABOVE_M1P4"
    )
    return _InviscidPoint(
        ca_pressure=ca_pressure,
        ca_base=base,
        ca_fin_wave=wave,
        ca_induced=induced,
        cn_body=cn_body,
        cn_fin=cn_fin,
        cm_body=cm_body,
        cm_fin=cm_fin,
        fin_cn_alpha=fin_slope,
        fin_cp_x=x_fin,
        cl_p=roll_damping_derivative(fin_slope, geometry),
        flags={
            "body_pressure": pressure_flag,
            "body_normal": (
                "EXACT_JOR_01_08" if mach < 1.2 else pressure_flag + "+EXACT_JOR_VISCOUS"
            ),
            "fin_normal": fin_method,
            "base_drag": base_method,
            "fin_wave_drag": wave_method,
            "fin_induced_drag": induced_method,
            "fin_cp": fin_cp_method,
            "interference": "CHART_AIAA96_3395_KW_ONLY",
        },
    )


def _local_derivative(values: np.ndarray, alpha_rad: np.ndarray) -> np.ndarray:
    if alpha_rad.size == 1:
        return np.full_like(values, np.nan)
    edge_order = 2 if alpha_rad.size >= 3 else 1
    return np.gradient(values, alpha_rad, axis=1, edge_order=edge_order)


def generate_aero_table(
    geometry: RocketGeometry,
    grid: AeroGrid,
    options: AerodynamicOptions | None = None,
) -> AeroCoefficientTable:
    """Generate a full-envelope fixed-geometry aerodynamic coefficient table.

    Method routing follows implementation-guide Secs. 4-10: slender-body plus
    Jorgensen below M=1.2, tangent-ogive NWL transonic wave drag, local-panel
    SOSE from M=1.2, Modified Newtonian from the configured M>=5 switch,
    FDAT/transonic-fairing/Ackeret fin routing, AIAA body-fin upwash and CP
    charts, explicit laminar/turbulent skin friction, power-off base drag, and
    component assembly. The intended range is conventional axisymmetric
    sounding rockets, |alpha|<=10 deg, M>0, and Reynolds numbers compatible
    with the selected attached flat-plate boundary-layer model. Every fallback
    and chart hold is recorded in ``method_flags`` and metadata; validate the
    resulting database before high-consequence use.
    """

    options = options or AerodynamicOptions()
    body = preprocess_body(geometry, options.body_station_count)
    mach_axis = grid.mach.array()
    alpha_axis_deg = grid.alpha_deg.array()
    alpha_axis = np.radians(alpha_axis_deg)
    phi_axis_deg = grid.phi_aero_deg.array()
    phi_axis = np.radians(phi_axis_deg)
    reynolds_axis = grid.reynolds.array()
    shape = (
        mach_axis.size, alpha_axis.size, phi_axis.size, reynolds_axis.size
    )
    meshes_raw = np.meshgrid(
        mach_axis, alpha_axis_deg, phi_axis_deg, reynolds_axis, indexing="ij"
    )
    meshes = {
        "mach": meshes_raw[0],
        "alpha_deg": meshes_raw[1],
        "phi_aero_deg": meshes_raw[2],
        "reynolds": meshes_raw[3],
    }

    component_names = (
        "ca_body_friction",
        "ca_fin_friction",
        "ca_body_pressure",
        "ca_base",
        "ca_fin_wave",
        "ca_fin_induced",
        "cn_body",
        "cn_fin",
        "cm_body",
        "cm_fin",
        "fin_cn_alpha",
        "fin_cp_x",
        "cy_fin",
        "cz_fin",
        "cmx_fin",
        "cmy_fin",
        "cmz_fin",
    )
    components = {name: np.empty(shape) for name in component_names}
    coefficients = {
        name: np.empty(shape)
        for name in (
            "CA", "CN", "CS", "CY", "CZ", "Cm", "Cmx", "Cmy", "Cmz",
            "CD", "CL", "xCP", "Cl_p",
        )
    }
    flag_names = (
        "body_pressure",
        "body_normal",
        "fin_normal",
        "base_drag",
        "fin_wave_drag",
        "fin_induced_drag",
        "fin_cp",
        "interference",
        "skin_friction",
    )
    method_flags = {name: np.empty(shape, dtype="<U64") for name in flag_names}
    fin_wetted_area = 2.0 * geometry.fins.total_area
    fin_reynolds_scale = (
        0.5 * (geometry.fins.root_chord + geometry.fins.tip_chord)
        / geometry.length
    )

    for i, mach in enumerate(mach_axis):
        for j, alpha in enumerate(alpha_axis):
            inviscid = _inviscid_point(
                float(mach), float(alpha), geometry, body, options
            )
            for p, phi_aero in enumerate(phi_axis):
                fins = _roll_aware_fins(
                    float(mach), float(alpha), float(phi_aero), geometry, options
                )
                for key, value in inviscid.flags.items():
                    method_flags[key][i, j, p, :] = value
                method_flags["fin_normal"][i, j, p, :] = fins.normal_method
                method_flags["fin_cp"][i, j, p, :] = fins.cp_method
                incidence_y = np.cos(phi_aero)
                incidence_z = np.sin(phi_aero)
                body_cy = inviscid.cn_body * incidence_y
                body_cz = inviscid.cn_body * incidence_z
                body_cmy = -inviscid.cm_body * incidence_z
                body_cmz = inviscid.cm_body * incidence_y
                cy = body_cy + fins.cy
                cz = body_cz + fins.cz
                cmy = body_cmy + fins.cmy
                cmz = body_cmz + fins.cmz
                cn = cy * incidence_y + cz * incidence_z
                cs = -cy * incidence_z + cz * incidence_y
                cm = -cmy * incidence_z + cmz * incidence_y
                xcp = (
                    geometry.moment_reference_x - geometry.length * cm / cn
                    if abs(cn) > 1.0e-12
                    else np.nan
                )
                for k, reynolds in enumerate(reynolds_axis):
                    cf_body, friction_flag = skin_friction_coefficient(
                        float(reynolds), options.boundary_layer
                    )
                    fin_reynolds = float(reynolds) * fin_reynolds_scale
                    cf_fin, _ = skin_friction_coefficient(
                        fin_reynolds, options.boundary_layer
                    )
                    ca_body_friction = (
                        cf_body * body.wetted_area / geometry.reference_area
                    )
                    ca_fin_friction = (
                        cf_fin * fin_wetted_area / geometry.reference_area
                    )
                    ca = (
                        ca_body_friction + ca_fin_friction
                        + inviscid.ca_pressure + inviscid.ca_base
                        + inviscid.ca_fin_wave + fins.ca_induced
                    )
                    cd = ca * np.cos(alpha) + cn * np.sin(alpha)
                    cl = cn * np.cos(alpha) - ca * np.sin(alpha)
                    index = (i, j, p, k)
                    components["ca_body_friction"][index] = ca_body_friction
                    components["ca_fin_friction"][index] = ca_fin_friction
                    components["ca_body_pressure"][index] = inviscid.ca_pressure
                    components["ca_base"][index] = inviscid.ca_base
                    components["ca_fin_wave"][index] = inviscid.ca_fin_wave
                    components["ca_fin_induced"][index] = fins.ca_induced
                    components["cn_body"][index] = inviscid.cn_body
                    components["cn_fin"][index] = fins.cn_along_incidence
                    components["cm_body"][index] = inviscid.cm_body
                    components["cm_fin"][index] = fins.cm_along_incidence
                    components["fin_cn_alpha"][index] = fins.slope
                    components["fin_cp_x"][index] = fins.cp_x
                    components["cy_fin"][index] = fins.cy
                    components["cz_fin"][index] = fins.cz
                    components["cmx_fin"][index] = fins.cmx
                    components["cmy_fin"][index] = fins.cmy
                    components["cmz_fin"][index] = fins.cmz
                    coefficients["CA"][index] = ca
                    coefficients["CN"][index] = cn
                    coefficients["CS"][index] = cs
                    coefficients["CY"][index] = cy
                    coefficients["CZ"][index] = cz
                    coefficients["Cm"][index] = cm
                    coefficients["Cmx"][index] = fins.cmx
                    coefficients["Cmy"][index] = cmy
                    coefficients["Cmz"][index] = cmz
                    coefficients["CD"][index] = cd
                    coefficients["CL"][index] = cl
                    coefficients["xCP"][index] = xcp
                    coefficients["Cl_p"][index] = inviscid.cl_p
                    method_flags["skin_friction"][index] = friction_flag

    coefficients["CN_alpha"] = _local_derivative(coefficients["CN"], alpha_axis)
    coefficients["Cm_alpha"] = _local_derivative(coefficients["Cm"], alpha_axis)
    body_slope = _local_derivative(components["cn_body"], alpha_axis)
    body_arm = (body.potential_cp_x - geometry.moment_reference_x) / geometry.length
    fin_arm = (
        components["fin_cp_x"] - geometry.moment_reference_x
    ) / geometry.length
    coefficients["Cm_q"] = -2.0 * (
        body_slope * body_arm**2
        + components["fin_cn_alpha"] * fin_arm**2
    )  # DYN-08 quasi-steady fallback

    # CP is indeterminate at exactly zero normal force. Its zero-incidence
    # limit is needed by the runtime stability diagnostic and force lever arm,
    # so interpolate it along the alpha axis from the nearest finite values.
    xcp = coefficients["xCP"]
    for mach_index in range(xcp.shape[0]):
        for phi_index in range(xcp.shape[2]):
            for reynolds_index in range(xcp.shape[3]):
                line = xcp[mach_index, :, phi_index, reynolds_index]
                finite = np.isfinite(line)
                if np.any(finite):
                    line[~finite] = np.interp(
                        alpha_axis[~finite], alpha_axis[finite], line[finite]
                    )
                else:
                    line[:] = geometry.moment_reference_x

    metadata: dict[str, object] = {
        "format": "rsim.aero",
        "format_version": FORMAT_VERSION,
        "dimension_order": ["mach", "alpha_deg", "phi_aero_deg", "reynolds"],
        "angle_units": "degrees on axes; radians in equation evaluation",
        "geometry_units": "SI",
        "reynolds_reference_length": geometry.length,
        "coefficient_reference_area": geometry.reference_area,
        "coefficient_reference_length": geometry.length,
        "stability_reference_length": geometry.reference_diameter,
        "moment_reference_x": geometry.moment_reference_x,
        "phi_aero_period_deg": 360.0 / geometry.fins.count,
        "geometry": asdict(geometry),
        "options": asdict(options),
        "limitations": [
            "Conventional symmetric tangent-ogive/cylinder rocket with one aft fin set.",
            "Low angle of attack only: |alpha| <= 10 deg.",
            "phi_aero dependence uses low-angle individual-fin projected incidence.",
            "Chart and fallback methods are not exact Missile DATCOM reproduction.",
            "Base-pressure chart is held constant above its Mach-6 digitization limit.",
            "Motor-off base drag only; no protuberance, control, plume, or aeroelastic model.",
            "Validate transonic and configuration-specific results against test, CFD, or trusted software.",
        ],
        "sources": [
            "Sounding Rocket Aerodynamic Database Reference implementation guide",
            "Barrowman 1967",
            "Jorgensen NASA TR R-474",
            "USAF Stability and Control DATCOM Section 4",
            "NSWC TR 81-156 digitized charts",
            "NWL TR-2796 Figure 3 digitization",
            "AIAA 96-3395 and AIAA 2007-3937 digitizations",
            "phi_aero_dependence_implementation_note individual-fin synthesis",
            "NASA TR R-100",
        ],
    }
    return AeroCoefficientTable(
        axes={
            "mach": mach_axis,
            "alpha_deg": alpha_axis_deg,
            "phi_aero_deg": phi_axis_deg,
            "reynolds": reynolds_axis,
        },
        meshes=meshes,
        coefficients=coefficients,
        components=components,
        method_flags=method_flags,
        metadata=metadata,
    )
