"""Small source tables packaged with the aerodynamic implementation.

Values are copied from the verified CSV digitizations under ``aero_imp_guide``
so installed wheels do not depend on the repository-only reference folder.
"""

from __future__ import annotations

import numpy as np


# NSWC TR 81-156, Fig. 51: mean power-off base pressure at alpha=0.
BASE_MACH = np.array([0.0, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0])
BASE_MINUS_CP = np.array([0.13, 0.132, 0.15, 0.2, 0.22, 0.21, 0.185, 0.15, 0.115, 0.09, 0.06, 0.044, 0.033])

# NSWC TR 81-156, Fig. 70: Jorgensen crossflow factors.
JORGENSEN_MACH = np.array([0.0, 0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5])
JORGENSEN_CDC = np.array([1.2, 1.2, 1.25, 1.38, 1.58, 1.78, 1.88, 1.9, 1.87, 1.72, 1.48, 1.28, 1.2])
JORGENSEN_FINENESS = np.array([1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0])
JORGENSEN_ETA = np.array([0.525, 0.56, 0.59, 0.615, 0.645, 0.665, 0.685])

# NWL TR-2796, Fig. 3: tangent-ogive transonic wave axial coefficient.
NWL_MACH = np.array([0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15])
NWL_FINENESS = np.array([1.5, 2.0, 2.5, 3.0, 4.0])
NWL_CA_WAVE = np.array([
    [0.000, 0.006, 0.033, 0.061, 0.092, 0.124, 0.154, 0.181, 0.207, 0.231, 0.252, 0.271, 0.286],
    [0.000, 0.002, 0.018, 0.038, 0.059, 0.080, 0.100, 0.120, 0.139, 0.156, 0.171, 0.184, 0.195],
    [0.000, 0.000, 0.008, 0.020, 0.034, 0.049, 0.064, 0.079, 0.094, 0.107, 0.118, 0.128, 0.136],
    [0.000, 0.000, 0.003, 0.010, 0.020, 0.032, 0.044, 0.057, 0.069, 0.080, 0.089, 0.097, 0.102],
    [0.000, 0.000, 0.000, 0.003, 0.009, 0.016, 0.023, 0.031, 0.038, 0.044, 0.048, 0.050, 0.050],
])

# AIAA 2007-3937, Fig. 2: beta=1 nonlinear fin-force curve.
FIN_CNAA_X = np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16])
FIN_CNAA = np.array([0.0, 0.1, 0.27, 0.5, 0.82, 1.22, 1.72, 2.25, 2.4])

# AIAA 96-3395, Fig. 10: empirical body-to-fin upwash at alpha=0.
UPWASH_R_OVER_S = np.array([0.0, 0.1, 0.2, 0.25, 0.33, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
UPWASH_KW = np.array([1.0, 1.09, 1.21, 1.29, 1.39, 1.46, 1.55, 1.62, 1.70, 1.79, 1.88, 2.0])

# AIAA 2007-3937, Fig. 12: fin CP migration factor for conventional taper>0.5.
FIN_CP_ALPHA_DEG = np.array([0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0])
FIN_CP_FACTOR_AR_LT_05 = np.array([0.0, 0.30, 0.48, 0.65, 0.76, 0.85, 0.91, 0.96, 0.985, 0.997, 1.0])
FIN_CP_FACTOR_AR_GT_1 = np.array([0.0, 0.11, 0.23, 0.44, 0.61, 0.73, 0.84, 0.91, 0.96, 0.99, 1.0])
FIN_CP_FACTOR_HIGH_MACH = np.array([0.0, 0.055, 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 1.0])

# AIAA 2007-3937, Fig. 10: high-Mach aerodynamic-center forward shift.
XAC_MACH = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
XAC_THICKNESS = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.20])
XAC_SHIFT = np.array([
    [0.012, 0.015, 0.018, 0.021, 0.024, 0.027],
    [0.025, 0.030, 0.035, 0.041, 0.046, 0.052],
    [0.037, 0.045, 0.053, 0.061, 0.068, 0.075],
    [0.048, 0.058, 0.068, 0.078, 0.088, 0.098],
    [0.060, 0.071, 0.083, 0.096, 0.107, 0.118],
    [0.080, 0.095, 0.110, 0.125, 0.138, 0.150],
    [0.109, 0.124, 0.141, 0.155, 0.169, 0.180],
])
