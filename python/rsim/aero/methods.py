"""Sourced aerodynamic equations used by the baseline method router.

Every public function states its assumptions and validity. Equation tags refer
to ``aero_imp_guide/.../sounding_rocket_aerodynamic_database_reference_updated``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from . import _chart_data as charts
from .geometry import RocketGeometry, TrapezoidalFinSet


def _integrate(values: np.ndarray, coordinates: np.ndarray) -> float:
    """Integrate with the NumPy 2.x name while retaining NumPy 1.x support."""

    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is not None:
        return float(trapezoid(values, coordinates))
    return float(np.trapz(values, coordinates))


@dataclass(frozen=True)
class BodyGeometryData:
    """Numerically preprocessed axisymmetric body geometry."""

    x: np.ndarray
    radius: np.ndarray
    planform_area: float
    planform_centroid_x: float
    wetted_area: float
    potential_cn_alpha: float
    potential_cp_x: float


@dataclass(frozen=True)
class PressureResult:
    """Body pressure coefficients and whether a detached-shock fallback occurred."""

    ca: float
    cn: float
    cm: float
    used_detached_fallback: bool = False


def preprocess_body(
    geometry: RocketGeometry, station_count: int = 121
) -> BodyGeometryData:
    """Precompute body integrals and station-area slender-body quantities.

    Implements GEO-01 through GEO-03 and SB-01 through SB-03 by numerical
    station integration. Assumes a slender axisymmetric tangent-ogive,
    cylinder, and optional conical boattail. ``station_count >= 21`` is
    recommended; convergence should be checked for short/steep transitions.
    Source: Barrowman 1967 and the implementation guide Sections 3 and 4.1.
    """

    if station_count < 5:
        raise ValueError("station_count must be at least 5")
    nose_count = max(3, round(station_count * geometry.nose.length / geometry.length))
    body_count = max(2, round(station_count * geometry.body.length / geometry.length))
    nose_x = np.linspace(0.0, geometry.nose.length, nose_count, endpoint=False)
    nose_r = geometry.nose.radius(nose_x)
    body_start = geometry.nose.length
    body_end = body_start + geometry.body.length
    body_x = np.linspace(body_start, body_end, body_count)
    body_r = np.full_like(body_x, 0.5 * geometry.body.diameter)
    pieces_x = [nose_x, body_x]
    pieces_r = [nose_r, body_r]
    if geometry.boattail:
        tail_count = max(2, station_count - nose_count - body_count + 2)
        tail_x = np.linspace(body_end, geometry.length, tail_count)
        tail_r = np.linspace(
            0.5 * geometry.body.diameter,
            0.5 * geometry.boattail.exit_diameter,
            tail_count,
        )
        pieces_x.append(tail_x[1:])
        pieces_r.append(tail_r[1:])
    x = np.concatenate(pieces_x)
    radius = np.concatenate(pieces_r)

    # GEO-02/GEO-03: projected side area and centroid.
    diameter = 2.0 * radius
    planform_area = _integrate(diameter, x)
    planform_centroid = _integrate(x * diameter, x) / planform_area

    # Surface of revolution frusta used for skin friction and pressure panels.
    dx = np.diff(x)
    dr = np.diff(radius)
    wetted_area = float(np.sum(math.pi * (radius[:-1] + radius[1:]) * np.hypot(dx, dr)))

    # SB-02/SB-03: force increments are proportional to changes in cross-section.
    area = math.pi * radius**2
    delta_area = np.diff(area)
    x_mid = 0.5 * (x[:-1] + x[1:])
    increments = 2.0 * delta_area / geometry.reference_area
    cn_alpha = float(np.sum(increments))
    if abs(cn_alpha) > 1.0e-14:
        cp = float(np.sum(increments * x_mid) / cn_alpha)
    else:
        cp = 0.0
    return BodyGeometryData(
        x=x,
        radius=radius,
        planform_area=planform_area,
        planform_centroid_x=planform_centroid,
        wetted_area=wetted_area,
        potential_cn_alpha=cn_alpha,
        potential_cp_x=cp,
    )


def jorgensen_body_force(
    mach: float,
    alpha_rad: float,
    geometry: RocketGeometry,
    body: BodyGeometryData,
) -> tuple[float, float, float, float]:
    """Return Jorgensen body ``CN, Cm, CN_potential, CN_viscous``.

    Implements JOR-01 through JOR-08 using digitized NSWC TR 81-156 Fig. 70
    ``eta`` and ``C_Dc`` data. Valid for slender circular bodies where the
    crossflow analogy and potential/viscous superposition are reasonable. The
    supplied chart is Mach-dependent but does not resolve Reynolds dependence;
    do not use for asymmetric vortex states, unusual cross-sections, or large
    protuberances. The project database remains limited to |alpha| <= 10 deg.
    Source: Jorgensen NASA TR R-474 and Simon & Blake AIAA 99-4258.
    """

    sin_alpha = math.sin(alpha_rad)
    cn_p = (
        0.5
        * body.potential_cn_alpha
        * math.sin(2.0 * alpha_rad)
        * math.cos(0.5 * alpha_rad)
    )
    if abs(sin_alpha) < 1.0e-15:
        cn_v = 0.0
    else:
        crossflow_mach = mach * abs(sin_alpha)  # JOR-04
        cdc = float(
            np.interp(
                crossflow_mach,
                charts.JORGENSEN_MACH,
                charts.JORGENSEN_CDC,
            )
        )
        fineness = geometry.length / geometry.reference_diameter
        eta = float(
            np.interp(
                fineness,
                charts.JORGENSEN_FINENESS,
                charts.JORGENSEN_ETA,
            )
        )
        # JOR-02; circular cross-section correction is unity.
        cn_v = (
            eta
            * cdc
            * body.planform_area
            / geometry.reference_area
            * sin_alpha
            * abs(sin_alpha)
        )
    # JOR-06 with the guide's x-aft, positive-CN, nose-up-Cm convention.
    cm = (
        -cn_p
        * (body.potential_cp_x - geometry.moment_reference_x)
        / geometry.length
        - cn_v
        * (body.planform_centroid_x - geometry.moment_reference_x)
        / geometry.length
    )
    return cn_p + cn_v, cm, cn_p, cn_v


def datcom_subsonic_fin_slope(
    mach: float, fins: TrapezoidalFinSet, reference_area: float
) -> float:
    """Return subsonic fin-set normal-force slope per radian.

    Implements FDAT-01 through FDAT-03 with ideal thin-section ``kappa=1`` and
    scales the isolated planar slope by total exposed fin area/reference area.
    Valid for 0 <= M < 0.8, thin straight-tapered aft-swept fins, attached
    linear flow, and ordinary aspect ratio. It excludes stall and body
    interference. Source: USAF DATCOM Section 4.1.3.2, Fig. 4.1.3.2-49.
    """

    if not 0.0 <= mach < 0.8:
        raise ValueError("the DATCOM subsonic fin method requires 0 <= Mach < 0.8")
    beta = math.sqrt(1.0 - mach**2)  # FDAT-01
    aspect_ratio = fins.aspect_ratio
    sweep = fins.half_chord_sweep
    denominator = 2.0 + math.sqrt(
        4.0
        + (aspect_ratio * beta) ** 2
        * (1.0 + math.tan(sweep) ** 2 / beta**2)
    )
    planar_slope = 2.0 * math.pi * aspect_ratio / denominator  # FDAT-03
    return planar_slope * fins.total_area / reference_area


def barrowman_body_fin_factor(body_radius: float, semi_span: float) -> float:
    """Return the classical fin-in-presence-of-body factor ``K_FB``.

    Implements BARF-02. It is a low-Mach, low-angle fallback for conventional
    body-mounted fins and must not be combined with the guide's equivalent-angle
    ``K_W/K_B`` synthesis. Source: Barrowman 1967.
    """

    return 1.0 + body_radius / (semi_span + body_radius)


def empirical_upwash_factor(body_radius: float, semi_span: float) -> float:
    """Interpolate the AIAA 96-3395 body-to-fin upwash factor at alpha=0.

    The digitized Fig. 10 data are a chart method for conventional body-mounted
    fins. Input is ``R/s`` and is valid only over the chart range 0 to 1; this
    implementation rejects extrapolation. The factor is used without a
    separate carryover factor, avoiding double counting.
    """

    ratio = body_radius / semi_span
    if ratio < charts.UPWASH_R_OVER_S[0] or ratio > charts.UPWASH_R_OVER_S[-1]:
        raise ValueError("body-radius/fin-span ratio is outside AIAA 96-3395 Fig. 10")
    return float(np.interp(ratio, charts.UPWASH_R_OVER_S, charts.UPWASH_KW))


def ackeret_supersonic_fin_slope(
    mach: float, fins: TrapezoidalFinSet, reference_area: float
) -> float:
    """Return the guide's low-order supersonic fin-set slope check.

    Implements FSUP-02, ``CN_alpha,2D=4/sqrt(M^2-1)``, then applies only
    reference-area bookkeeping for the exposed fins. Valid for M >= 1.4,
    thin fins with attached supersonic leading-edge flow, and small incidence.
    It is explicitly a FALLBACK/check rather than the full finite-fin DATCOM
    chart solution; detached shocks and subsonic leading edges require another
    method. Source: Ackeret linearized thin-airfoil theory and guide Sec. 5.4.
    """

    if mach < 1.4:
        raise ValueError("the baseline supersonic fin route starts at Mach 1.4")
    two_dimensional_slope = 4.0 / math.sqrt(mach**2 - 1.0)
    return two_dimensional_slope * fins.total_area / reference_area


def smoothstep_transonic(
    mach: float, value_at_08: float, value_at_14: float
) -> float:
    """Fair a coefficient between its native Mach 0.8 and 1.4 methods.

    Implements fallback equations TRN-01 through TRN-03. Valid only for
    0.8 <= M <= 1.4 when exact RAS/NSWC chart-native interpolation is not
    available. It is smooth with zero endpoint slopes but is not exact Missile
    DATCOM/RAS reproduction. Source: implementation-guide fallback Sec. 5.3.
    """

    if not 0.8 <= mach <= 1.4:
        raise ValueError("transonic smoothstep requires 0.8 <= Mach <= 1.4")
    tau = (mach - 0.8) / 0.6  # TRN-01
    blend = 3.0 * tau**2 - 2.0 * tau**3  # TRN-02
    return (1.0 - blend) * value_at_08 + blend * value_at_14  # TRN-03


def fin_normal_force(
    mach: float,
    alpha_rad: float,
    geometry: RocketGeometry,
) -> tuple[float, float, str]:
    """Return body-corrected fin ``CN``, slope, and routed method name.

    Routes FDAT-03 below M=0.8, TRN-01..03 from 0.8 through 1.4, and
    the FSUP-02 fallback above 1.4. Uses the digitized AIAA 96-3395 ``K_W``
    factor once. Potential force follows FNL-02's ``sin(2 alpha)/2`` form and
    the digitized AIAA 2007-3937 Fig. 2 nonlinear term is included for grids
    reaching the upper low-angle range. Valid for one conventional fixed fin
    set, |alpha| <= 10 deg, no deflection, and attached flow assumptions of
    the selected native method.
    """

    low = datcom_subsonic_fin_slope(0.8 - 1.0e-9, geometry.fins, geometry.reference_area)
    high = ackeret_supersonic_fin_slope(1.4, geometry.fins, geometry.reference_area)
    if mach < 0.8:
        isolated_slope = datcom_subsonic_fin_slope(mach, geometry.fins, geometry.reference_area)
        method = "EXACT_FDAT_03"
    elif mach < 1.4:
        isolated_slope = smoothstep_transonic(mach, low, high)
        method = "FALLBACK_TRN_03"
    else:
        isolated_slope = ackeret_supersonic_fin_slope(mach, geometry.fins, geometry.reference_area)
        method = "FALLBACK_FSUP_02"
    upwash = empirical_upwash_factor(
        0.5 * geometry.reference_diameter, geometry.fins.semi_span
    )
    slope = isolated_slope * upwash
    potential = 0.5 * slope * math.sin(2.0 * alpha_rad)  # FNL-02

    beta = math.sqrt(abs(1.0 - mach**2))
    chart_coordinate = min(charts.FIN_CNAA_X[-1], beta * abs(math.tan(alpha_rad)))
    nonlinear_ordinate = float(
        np.interp(chart_coordinate, charts.FIN_CNAA_X, charts.FIN_CNAA)
    )
    nonlinear_scale = (
        geometry.fins.total_area / geometry.reference_area * upwash
    )
    nonlinear = (
        nonlinear_ordinate
        * nonlinear_scale
        * math.sin(alpha_rad)
        * abs(math.sin(alpha_rad))
    )
    return potential + nonlinear, slope, method + "+CHART_FNL_01"


def barrowman_fin_cp(fins: TrapezoidalFinSet) -> float:
    """Return trapezoidal-fin CP using corrected guide equation BARF-04.

    Valid as a transparent low-Mach fallback for conventional thin trapezoidal
    fins. It does not reproduce transonic shock movement or unusual planforms.
    Source: Barrowman 1967; the guide corrects the first term to use tip
    leading-edge offset rather than mid-chord-line length.
    """

    cr = fins.root_chord
    ct = fins.tip_chord
    local = (
        fins.tip_leading_edge_offset / 3.0 * (cr + 2.0 * ct) / (cr + ct)
        + (cr + ct - cr * ct / (cr + ct)) / 6.0
    )
    return fins.root_leading_edge_x + local


def fin_cp(mach: float, alpha_rad: float, fins: TrapezoidalFinSet) -> float:
    """Return fin CP with AIAA 2007-3937 migration/forward-shift charts.

    Starts from BARF-04, applies FCP-03's digitized high-Mach forward shift for
    2 <= M <= 4.5 and chart thickness ratios 0.02..0.20, then applies FCP-02
    migration toward the planform centroid for taper ratios greater than 0.5.
    Fig. 12 supports |alpha| to 90 deg, but this package limits the database to
    10 deg. Intermediate 0.5<=AR<=1 is linearly blended between the two
    available conventional families. For unsupported ``taper<=0.5``, migration
    is omitted rather than inventing/extrapolating a chart family. This is a
    CHART route and should not be extrapolated beyond chart geometry/ranges.
    """

    x_ac = barrowman_fin_cp(fins)
    if mach >= charts.XAC_MACH[0]:
        by_thickness = np.array(
            [np.interp(mach, charts.XAC_MACH, row) for row in charts.XAC_SHIFT]
        )
        shift_ratio = float(
            np.interp(
                fins.thickness_ratio,
                charts.XAC_THICKNESS,
                by_thickness,
            )
        )
        x_ac -= shift_ratio * fins.root_chord
    if fins.taper_ratio <= 0.5:
        return x_ac
    absolute_alpha = abs(math.degrees(alpha_rad))
    if mach >= 2.0:
        family = charts.FIN_CP_FACTOR_HIGH_MACH
    elif fins.aspect_ratio <= 0.5:
        family = charts.FIN_CP_FACTOR_AR_LT_05
    elif fins.aspect_ratio >= 1.0:
        family = charts.FIN_CP_FACTOR_AR_GT_1
    else:
        fraction = (fins.aspect_ratio - 0.5) / 0.5
        family = (
            (1.0 - fraction) * charts.FIN_CP_FACTOR_AR_LT_05
            + fraction * charts.FIN_CP_FACTOR_AR_GT_1
        )
    migration = float(np.interp(absolute_alpha, charts.FIN_CP_ALPHA_DEG, family))
    return x_ac + migration * (fins.centroid_x - x_ac)


def skin_friction_coefficient(reynolds: float, state: str) -> tuple[float, str]:
    """Return average smooth flat-plate skin-friction coefficient.

    ``state='laminar'`` implements FRIC-03 (Blasius), valid for an attached
    zero-pressure-gradient laminar boundary layer before transition.
    ``state='turbulent'`` implements fallback FRIC-06, valid approximately for
    smooth attached incompressible turbulent flow at Re > 1e5. Neither route
    resolves transition, roughness, compressible Van Driest corrections,
    curvature, separation, or shock/boundary-layer interaction. Source: USAF
    DATCOM Sec. 4.1.5.1 and the guide Sec. 7.2.
    """

    if reynolds <= 0.0:
        raise ValueError("Reynolds number must be positive")
    if state == "laminar":
        return 1.328 / math.sqrt(reynolds), "EXACT_FRIC_03"
    if state == "turbulent":
        if reynolds <= 1.0:
            raise ValueError("turbulent correlation requires Reynolds number > 1")
        return 0.455 / math.log10(reynolds) ** 2.58, "FALLBACK_FRIC_06"
    raise ValueError("boundary-layer state must be 'laminar' or 'turbulent'")


def tangent_ogive_wave_drag(mach: float, fineness_ratio: float) -> float:
    """Interpolate tangent-ogive transonic ``CA_wave`` from NWL TR-2796 Fig. 3.

    This CHART method is valid only for 0.85 <= M <= 1.15 and digitized nose
    fineness 1.5 <= L/D <= 4.0. It represents compatible tangent ogives and
    should not be extrapolated to blunt or non-ogive noses.
    """

    if not charts.NWL_MACH[0] <= mach <= charts.NWL_MACH[-1]:
        raise ValueError("NWL Fig. 3 is limited to Mach 0.85 through 1.15")
    if not charts.NWL_FINENESS[0] <= fineness_ratio <= charts.NWL_FINENESS[-1]:
        raise ValueError("NWL Fig. 3 is limited to nose fineness 1.5 through 4.0")
    values_at_mach = np.array(
        [np.interp(mach, charts.NWL_MACH, row) for row in charts.NWL_CA_WAVE]
    )
    return float(np.interp(fineness_ratio, charts.NWL_FINENESS, values_at_mach))


def base_drag_coefficient(mach: float, base_area_ratio: float) -> float:
    """Return cylindrical power-off base axial coefficient using BASE-01.

    Interpolates digitized NSWC TR 81-156 Fig. 51 ``-Cp_base`` for 0 <= M <= 6
    and multiplies by ``A_base/Sref``. Assumes a conventional cylindrical
    power-off base at zero/low angle. It is unreliable for motor-on plumes,
    clustered motors, or unusual afterbody/nozzle geometry.
    """

    if not charts.BASE_MACH[0] <= mach <= charts.BASE_MACH[-1]:
        raise ValueError("NSWC Fig. 51 base drag is limited to Mach 0 through 6")
    minus_cp = float(np.interp(mach, charts.BASE_MACH, charts.BASE_MINUS_CP))
    return minus_cp * base_area_ratio


def boattail_pressure_plus_base_drag(
    mach: float,
    slope_deg: float,
    exit_to_body_radius: float,
    cylindrical_base_drag: float,
) -> float:
    """Return conical-afterbody pressure-plus-base drag from BASE-02.

    Implements NASA TR R-100 Eq. (1). Valid for supersonic flow, compatible
    conical converging afterbodies, positive slope magnitude, and preferably
    slope below about 15 deg. ``n=4`` below M=3.5 and ``n=3`` above. The method
    can have large errors for steep slopes and is not a motor-on plume model.
    """

    if mach < 1.0:
        raise ValueError("NASA TR R-100 Eq. (1) is used only for supersonic flow")
    if not 0.0 <= exit_to_body_radius <= 1.0:
        raise ValueError("exit/body radius ratio must lie in [0, 1]")
    exponent = 4 if mach < 3.5 else 3
    pressure = (
        (0.001 * slope_deg + 0.00071 * slope_deg**2)
        / mach
        * (1.0 - exit_to_body_radius**exponent)
    )
    return pressure + cylindrical_base_drag * exit_to_body_radius**3


def prandtl_meyer(mach: float, gamma: float = 1.4) -> float:
    """Return the perfect-gas Prandtl-Meyer angle using SOSE-06.

    Valid only for supersonic Mach and a calorically perfect gas with gamma>1.
    Source: standard expansion-wave relation; guide Sec. 4.6.
    """

    if mach < 1.0 or gamma <= 1.0:
        raise ValueError("Prandtl-Meyer requires Mach >= 1 and gamma > 1")
    root = math.sqrt(mach**2 - 1.0)
    return math.sqrt((gamma + 1.0) / (gamma - 1.0)) * math.atan(
        math.sqrt((gamma - 1.0) / (gamma + 1.0)) * root
    ) - math.atan(root)


def _oblique_shock_angle(mach: float, theta: float, gamma: float) -> float | None:
    mach_angle = math.asin(1.0 / mach)

    def residual(beta: float) -> float:
        numerator = mach**2 * math.sin(beta) ** 2 - 1.0
        denominator = mach**2 * (gamma + math.cos(2.0 * beta)) + 2.0
        rhs = 2.0 / math.tan(beta) * numerator / denominator
        return rhs - math.tan(theta)

    samples = np.linspace(mach_angle + 1.0e-8, math.pi / 2.0 - 1.0e-8, 300)
    values = [residual(float(beta)) for beta in samples]
    for left, right, f_left, f_right in zip(samples[:-1], samples[1:], values[:-1], values[1:]):
        if f_left <= 0.0 <= f_right:
            for _ in range(60):
                middle = 0.5 * (left + right)
                f_middle = residual(float(middle))
                if f_middle >= 0.0:
                    right = middle
                else:
                    left = middle
            return float(0.5 * (left + right))
    return None


def shock_expansion_cp(
    mach: float, turn_angle_rad: float, gamma: float = 1.4
) -> tuple[float, bool]:
    """Return local pressure coefficient for one compression or expansion turn.

    Implements SOSE-01 through SOSE-09 relative to freestream. Positive turns
    select the weak attached oblique-shock root; negative turns use a
    Prandtl-Meyer expansion. Valid for M>1, calorically perfect gas, local 2-D
    turning, and attached compression shocks. Returns ``detached=True`` when no
    weak attached root exists; callers must use another model. It excludes
    viscous and shock-interaction effects. Source: guide Sec. 4.6.
    """

    if mach <= 1.0 or gamma <= 1.0:
        raise ValueError("shock expansion requires Mach > 1 and gamma > 1")
    if abs(turn_angle_rad) < 1.0e-14:
        return 0.0, False
    if turn_angle_rad > 0.0:
        beta = _oblique_shock_angle(mach, turn_angle_rad, gamma)  # SOSE-01
        if beta is None:
            return 0.0, True
        mn1 = mach * math.sin(beta)  # SOSE-02
        pressure_ratio = 1.0 + 2.0 * gamma / (gamma + 1.0) * (mn1**2 - 1.0)  # SOSE-03
    else:
        target = prandtl_meyer(mach, gamma) + abs(turn_angle_rad)  # SOSE-07
        low, high = mach, max(2.0 * mach, mach + 1.0)
        while prandtl_meyer(high, gamma) < target and high < 1.0e4:
            high *= 2.0
        for _ in range(70):
            middle = 0.5 * (low + high)
            if prandtl_meyer(middle, gamma) < target:
                low = middle
            else:
                high = middle
        downstream = 0.5 * (low + high)
        pressure_ratio = (
            (1.0 + 0.5 * (gamma - 1.0) * mach**2)
            / (1.0 + 0.5 * (gamma - 1.0) * downstream**2)
        ) ** (gamma / (gamma - 1.0))  # SOSE-08
    cp = 2.0 / (gamma * mach**2) * (pressure_ratio - 1.0)  # SOSE-09
    return cp, False


def modified_newtonian_cp_max(mach: float, gamma: float = 1.4) -> float:
    """Return stagnation ``Cp,max`` using MN-02 through MN-05.

    Valid for supersonic/high-Mach calorically perfect gas. The surrounding
    modified-Newtonian panel method is recommended around M>=5 only after
    validation and omits real-gas chemistry and viscous coupling.
    Source: NASA TN D-176 and guide Sec. 4.7.
    """

    if mach <= 1.0 or gamma <= 1.0:
        raise ValueError("modified Newtonian Cp,max requires Mach > 1")
    p2_p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (mach**2 - 1.0)  # MN-02
    m2_sq = (1.0 + 0.5 * (gamma - 1.0) * mach**2) / (
        gamma * mach**2 - 0.5 * (gamma - 1.0)
    )  # MN-03
    p02_p2 = (1.0 + 0.5 * (gamma - 1.0) * m2_sq) ** (
        gamma / (gamma - 1.0)
    )  # MN-04
    return (p2_p1 * p02_p2 - 1.0) / (0.5 * gamma * mach**2)  # MN-05


def axisymmetric_body_pressure(
    mach: float,
    alpha_rad: float,
    geometry: RocketGeometry,
    body: BodyGeometryData,
    *,
    gamma: float = 1.4,
    method: str = "shock_expansion",
    azimuth_count: int = 48,
) -> PressureResult:
    """Integrate local pressure over an axisymmetric body surface.

    Implements SOSE-10..12 with independent local 2-D shock/expansion turns,
    or MN-01/MN-06/MN-07 for ``method='modified_newtonian'``. Valid for
    M>=1.2 with attached local shocks in shock-expansion mode; modified
    Newtonian is intended for M>=5 after validation. The panelization is an
    open FALLBACK, not exact Missile DATCOM SOSE bookkeeping. It neglects
    viscous pressure interaction, 3-D shock interaction, and real gas effects.
    """

    if mach < 1.2:
        raise ValueError("body pressure panel integration requires Mach >= 1.2")
    if azimuth_count < 8:
        raise ValueError("azimuth_count must be at least 8")
    flow = np.array([math.cos(alpha_rad), -math.sin(alpha_rad), 0.0])
    dphi = 2.0 * math.pi / azimuth_count
    total_force = np.zeros(3)
    total_pitch_moment = 0.0
    detached = False
    cp_max = modified_newtonian_cp_max(mach, gamma) if method == "modified_newtonian" else None

    for x0, x1, r0, r1 in zip(body.x[:-1], body.x[1:], body.radius[:-1], body.radius[1:]):
        dx = x1 - x0
        dr = r1 - r0
        if dx <= 0.0:
            continue
        slope = dr / dx
        normal_scale = math.sqrt(1.0 + slope**2)
        strip_area = 0.5 * (r0 + r1) * math.hypot(dx, dr) * dphi
        x_mid = 0.5 * (x0 + x1)
        r_mid = 0.5 * (r0 + r1)
        for index in range(azimuth_count):
            phi = (index + 0.5) * dphi
            normal = np.array(
                [-slope, math.cos(phi), math.sin(phi)], dtype=float
            ) / normal_scale
            flow_dot_normal = float(np.dot(flow, normal))
            if method == "modified_newtonian":
                impact = max(0.0, -flow_dot_normal)
                cp = float(cp_max) * impact**2  # MN-01 and MN-06
            elif method == "shock_expansion":
                turn = -math.asin(max(-1.0, min(1.0, flow_dot_normal)))
                cp, panel_detached = shock_expansion_cp(mach, turn, gamma)
                if panel_detached:
                    detached = True
                    impact = max(0.0, -flow_dot_normal)
                    cp = modified_newtonian_cp_max(mach, gamma) * impact**2
            else:
                raise ValueError("pressure method must be shock_expansion or modified_newtonian")
            force = -cp * normal * strip_area  # SOSE-10/MN-07, divided by qbar
            total_force += force
            y_mid = r_mid * math.cos(phi)
            total_pitch_moment += x_mid * force[1] - y_mid * force[0]

    # x points aft; CN is defined opposite +y for the selected positive-alpha flow.
    ca = total_force[0] / geometry.reference_area
    cn = -total_force[1] / geometry.reference_area
    cm = total_pitch_moment / (geometry.reference_area * geometry.length)
    # Shift panel moment from nose datum to the requested reference location.
    cm += cn * geometry.moment_reference_x / geometry.length
    return PressureResult(float(ca), float(cn), float(cm), detached)


def fin_wave_drag(
    mach: float, fins: TrapezoidalFinSet, reference_area: float
) -> tuple[float, str]:
    """Return zero-lift thin double-wedge fin wave drag.

    Implements FDRAG-02 (zero below M=1.05), FDRAG-03 (linear fairing to
    M=1.4), and a transparent Ackeret double-wedge FALLBACK above M=1.4:
    ``Cd_wave=4(t/c)^2/sqrt(M^2-1)`` with exposed-area bookkeeping. Valid for
    thin sharp fins and attached supersonic flow; it is not NWL-TR-3018 or a
    thick/blunt/detached-shock model.
    """

    def supersonic(value_mach: float) -> float:
        local = 4.0 * fins.thickness_ratio**2 / math.sqrt(value_mach**2 - 1.0)
        return local * fins.total_area / reference_area

    if mach < 1.05:
        return 0.0, "EXACT_FDRAG_02"
    if mach < 1.4:
        fraction = (mach - 1.05) / 0.35
        return fraction * supersonic(1.4), "FALLBACK_FDRAG_03"
    return supersonic(mach), "FALLBACK_ACKERET_WAVE"


def induced_fin_drag(
    mach: float,
    fin_cn: float,
    fins: TrapezoidalFinSet,
    reference_area: float,
    efficiency: float,
) -> float:
    """Return sub/transonic fin drag due to normal force using FDRAG-05.

    Valid for M<=1.4, attached lifting flow, untwisted fins, and a supplied
    calibrated span efficiency ``0<e<=1``. Missile DATCOM treats the separate
    induced term as zero above M=1.4 because supersonic pressure/wave loading
    accounts for it elsewhere. Source: USAF DATCOM Eq. 4.1.5.2-h special case.
    """

    if mach > 1.4:
        return 0.0
    if not 0.0 < efficiency <= 1.0:
        raise ValueError("span efficiency must lie in (0, 1]")
    area_ratio = fins.total_area / reference_area
    local_cl = fin_cn / area_ratio if area_ratio > 0.0 else 0.0
    local_cd = local_cl**2 / (math.pi * fins.aspect_ratio * efficiency)
    return local_cd * area_ratio


def roll_damping_derivative(
    fin_slope: float, geometry: RocketGeometry
) -> float:
    """Return Barrowman-type fin roll damping ``Cl_p`` using ROLL-02/03.

    Valid for a rigid rocket, small roll-induced incidence ``|p*r| << V``,
    trapezoidal fins, and the guide normalization ``p_hat=p*Dref/(2V)``.
    It omits unsteady, aeroelastic, and separated-flow effects.
    Source: Barrowman 1967 and guide Sec. 8.3.
    """

    fins = geometry.fins
    r0 = 0.5 * geometry.reference_diameter
    integral = fins.semi_span / 12.0 * (
        (fins.root_chord + 3.0 * fins.tip_chord) * fins.semi_span**2
        + 4.0 * (fins.root_chord + 2.0 * fins.tip_chord) * fins.semi_span * r0
        + 6.0 * (fins.root_chord + fins.tip_chord) * r0**2
    )  # ROLL-03
    panel_slope = fin_slope * geometry.reference_area / fins.total_area
    return (
        -2.0
        * fins.count
        * panel_slope
        / (geometry.reference_area * geometry.reference_diameter**2)
        * integral
    )  # ROLL-02 with Kd=1 for the baseline
