"""Regenerate a combustion CSV with RocketCEA (optional setup utility)."""

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
        default=HERE / "data" / "n2o_ethanol_rocketcea.csv",
    )
    args = parser.parse_args()

    table = rsim.CombustionTable.generate(
        chamber_pressure=np.array([0.4e6, 1.0e6, 2.0e6, 3.0e6, 4.0e6]),
        mixture_ratio=np.array([4.0, 5.0, 6.0, 7.0]),
        expansion_ratio=7.0,
        oxidizer="N2O",
        fuel="Ethanol",
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
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
