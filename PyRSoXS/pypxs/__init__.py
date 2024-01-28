"""
PyRSoXS - Python for Polarized Resonant Soft X-ray Scattering
---

This package includes a set of tools for simulating and analyzing resonant soft x-ray
scattering data. The package allows for modeling with any given sample - detector
geometry, and includes a set of tools for analyzing the resulting data.

Author: Devin Grabner, Harlan Heilman, Brian Collins
"""

from ._rotated_geometry import plot_map, qmap, qmap_histogram, solid_angle_correction
