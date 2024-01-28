Originally written by Devin Grabner, Washington State University, Nov 2023 for use in Wavemetrics Igor Pro

These function where used to replace the NIKA functions for calculating the QMap for the CCD
image and applying an appropriate Solid Angle Correction. This function allows for a rotated
CCD Geometry where the sample doesn't have to be located at the axis of rotation of the CCD.

All of the code was ported to Python by Devin Grabner and Harlan Heilman, Jan 2024

This code is written to model the Q Map, Scattering Angles, Normalized Solid Angle Correction
and Pixel # vs Q histogram for a given CCD geometry and energy. For more information about please see the
the supporting informational PowerPoint in -> /pyrsoxs/docs/CCD Rotation Diagram and Corrections.pptx

If you have any further questions please contact Devin Grabner (devin.grabner@wsu.edu),
Harlan Heilman (harlan.heilman@wsu.edu), or Brian Collins (brian.collins@wsu.edu) for more information.
