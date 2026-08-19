"""Setup-time helpers for native liquid-propulsion performance tables.

Runtime interpolation and network stepping live in C++.  This module is only
used to generate, validate, and serialize the numeric inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._rsim import CombustionTableData


def _axis(values, name: str) -> np.ndarray:
    result = np.array(values, dtype=np.float64, copy=True).reshape(-1)
    if result.size < 2 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain at least two finite values")
    if np.any(np.diff(result) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _field(values, shape: tuple[int, int], name: str) -> np.ndarray:
    result = np.array(values, dtype=np.float64, order="F", copy=True)
    if result.shape != shape or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return result


@dataclass(frozen=True)
class CombustionTable:
    """Python-side serializable table paired with native interpolation data."""

    chamber_pressure: np.ndarray
    mixture_ratio: np.ndarray
    chamber_temperature: np.ndarray
    cstar: np.ndarray
    cf_vac: np.ndarray
    isp_vac: np.ndarray
    expansion_ratio: float
    oxidizer: str = ""
    fuel: str = ""
    frozen_at_throat: bool = False

    def __post_init__(self) -> None:
        pc = _axis(self.chamber_pressure, "chamber_pressure")
        mr = _axis(self.mixture_ratio, "mixture_ratio")
        shape = (pc.size, mr.size)
        object.__setattr__(self, "chamber_pressure", pc)
        object.__setattr__(self, "mixture_ratio", mr)
        for name in ("chamber_temperature", "cstar", "cf_vac", "isp_vac"):
            object.__setattr__(self, name, _field(getattr(self, name), shape, name))
        if not np.isfinite(self.expansion_ratio) or self.expansion_ratio <= 0.0:
            raise ValueError("expansion_ratio must be positive")

    @property
    def native(self) -> CombustionTableData:
        """Return new native storage; every NumPy value is copied immediately."""

        return CombustionTableData(
            self.chamber_pressure,
            self.mixture_ratio,
            self.chamber_temperature,
            self.cstar,
            self.cf_vac,
            self.isp_vac,
            self.expansion_ratio,
            self.oxidizer,
            self.fuel,
            self.frozen_at_throat,
        )

    @classmethod
    def generate(
        cls,
        chamber_pressure,
        mixture_ratio,
        *,
        expansion_ratio: float,
        oxidizer: str,
        fuel: str,
        frozen_at_throat: bool = False,
    ) -> "CombustionTable":
        """Generate a Pc x MR table with RocketCEA's SI-units wrapper."""

        try:
            from rocketcea.cea_obj_w_units import CEA_Obj
        except ImportError as error:  # pragma: no cover - dependency diagnostic
            raise RuntimeError("RocketCEA 1.2.3 is required to generate tables") from error

        pc = _axis(chamber_pressure, "chamber_pressure")
        mr = _axis(mixture_ratio, "mixture_ratio")
        temperature = np.empty((pc.size, mr.size), order="F")
        cstar = np.empty_like(temperature)
        isp_vac = np.empty_like(temperature)
        cea = CEA_Obj(
            oxName=oxidizer,
            fuelName=fuel,
            pressure_units="Pa",
            cstar_units="m/s",
            temperature_units="K",
            isp_units="sec",
        )
        frozen = int(bool(frozen_at_throat))
        for j, ratio in enumerate(mr):
            for i, pressure in enumerate(pc):
                isp, characteristic_velocity, chamber_t = cea.get_IvacCstrTc(
                    Pc=float(pressure),
                    MR=float(ratio),
                    eps=float(expansion_ratio),
                    frozen=frozen,
                    frozenAtThroat=frozen,
                )
                isp_vac[i, j] = isp
                cstar[i, j] = characteristic_velocity
                temperature[i, j] = chamber_t
        cf_vac = isp_vac * 9.80665 / cstar
        return cls(pc, mr, temperature, cstar, cf_vac, isp_vac,
                   expansion_ratio, oxidizer, fuel, frozen_at_throat)

    def save_mat(self, path: str | Path) -> None:
        """Write an unambiguous MATLAB table using SciPy at setup time."""

        from scipy.io import savemat

        savemat(path, {
            "chamber_pressure": self.chamber_pressure,
            "mixture_ratio": self.mixture_ratio,
            "chamber_temperature": self.chamber_temperature,
            "cstar": self.cstar,
            "cf_vac": self.cf_vac,
            "isp_vac": self.isp_vac,
            "expansion_ratio": float(self.expansion_ratio),
            "oxidizer": self.oxidizer,
            "fuel": self.fuel,
            "frozen_at_throat": int(self.frozen_at_throat),
        })

    @classmethod
    def load_mat(cls, path: str | Path) -> "CombustionTable":
        """Validate a MATLAB file in Python and copy it into normalized arrays."""

        from scipy.io import loadmat

        values = loadmat(path, squeeze_me=True)
        required = {
            "chamber_pressure", "mixture_ratio", "chamber_temperature",
            "cstar", "cf_vac", "isp_vac", "expansion_ratio",
        }
        missing = sorted(required.difference(values))
        if missing:
            raise ValueError("missing combustion-table fields: " + ", ".join(missing))
        return cls(
            values["chamber_pressure"], values["mixture_ratio"],
            values["chamber_temperature"], values["cstar"],
            values["cf_vac"], values["isp_vac"],
            float(values["expansion_ratio"]),
            str(values.get("oxidizer", "")), str(values.get("fuel", "")),
            bool(values.get("frozen_at_throat", False)),
        )


def load_combustion_mat(path: str | Path) -> CombustionTableData:
    """Load, validate, and return a runtime-native combustion table."""

    return CombustionTable.load_mat(path).native
