"""Regenerate an analysis LOX/RP-1 table with optional RocketCEA.

This writes a new file; it does not overwrite the checked-in illustrative
table unless that path is supplied explicitly.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import rsim


HERE = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "data" / "lox_rp1_rocketcea.csv",
    )
    args = parser.parse_args()

    table = rsim.CombustionTable.generate(
        chamber_pressure=np.array([0.15e6, 0.5e6, 1.5e6, 3.0e6, 5.0e6]),
        mixture_ratio=np.array([1.5, 1.8, 2.0, 2.2, 2.5]),
        expansion_ratio=6.0,
        oxidizer="LOX",
        fuel="RP-1",
    )
    rows = []
    for i, pressure in enumerate(table.chamber_pressure):
        for j, ratio in enumerate(table.mixture_ratio):
            rows.append([
                pressure,
                ratio,
                table.chamber_temperature[i, j],
                table.cstar[i, j],
                table.cf_vac[i, j],
                table.isp_vac[i, j],
            ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.output,
        np.asarray(rows),
        delimiter=",",
        header=("chamber_pressure_pa,mixture_ratio,chamber_temperature_k,"
                "cstar_m_per_s,cf_vac,isp_vac_s"),
        comments="",
    )
    print(f"Wrote RocketCEA table to {args.output}")


if __name__ == "__main__":
    main()
