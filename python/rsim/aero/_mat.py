"""Minimal MATLAB level-5 writer for numeric aerodynamic table arrays."""

from __future__ import annotations

from pathlib import Path
import struct
from typing import Mapping

import numpy as np


_MI_INT8 = 1
_MI_UINT8 = 2
_MI_INT32 = 5
_MI_UINT32 = 6
_MI_DOUBLE = 9
_MI_MATRIX = 14
_MX_DOUBLE = 6
_MX_UINT8 = 9
_MX_UINT32 = 13


def _element(data_type: int, payload: bytes) -> bytes:
    """Pack one MAT data element and pad its payload to an eight-byte boundary."""
    padding = b"\0" * ((-len(payload)) % 8)
    return struct.pack("<II", data_type, len(payload)) + payload + padding


def _matrix(name: str, values: np.ndarray) -> bytes:
    """Encode one real numeric array using MATLAB's column-major storage order."""
    array = np.asarray(values)
    if array.dtype == np.uint8:
        matrix_class, data_type, stored = _MX_UINT8, _MI_UINT8, array
    elif array.dtype == np.uint32:
        matrix_class, data_type, stored = _MX_UINT32, _MI_UINT32, array.astype("<u4")
    else:
        matrix_class, data_type, stored = _MX_DOUBLE, _MI_DOUBLE, array.astype("<f8")
    dimensions = stored.shape if stored.ndim >= 2 else (stored.size, 1)
    payload = b"".join((
        _element(_MI_UINT32, struct.pack("<II", matrix_class, 0)),
        _element(_MI_INT32, np.asarray(dimensions, dtype="<i4").tobytes()),
        _element(_MI_INT8, name.encode("ascii")),
        _element(data_type, stored.tobytes(order="F")),
    ))
    return _element(_MI_MATRIX, payload)


def write_mat_v5(path: str | Path, arrays: Mapping[str, np.ndarray]) -> None:
    """Write named real arrays in an uncompressed MATLAB level-5 file."""
    description = b"MATLAB 5.0 MAT-file, Platform: rsim, Created by rsim.aero"
    header = description.ljust(116, b" ") + b"\0" * 8 + b"\x00\x01IM"
    with Path(path).open("wb") as output:
        output.write(header)
        for name, values in arrays.items():
            if not name.isascii() or len(name) > 63:
                raise ValueError("MAT variable names must be ASCII and at most 63 characters")
            output.write(_matrix(name, np.asarray(values)))
