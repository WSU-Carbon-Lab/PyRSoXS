"""
Core functions for calculating the qmap for a rotated CCD geometry.
---
This module contains the core functions for calculating the qmap for a rotated CCD
geometry. The qmap is calculated using a complex set of equations that can be described
further in the documentation.

Author: Devin Grabner, Harlan Heilman, Brian Collins
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

from .core.fits import Parameters

# parameters class moved to core/fits/py


def qmap(params: Parameters) -> tuple[np.ndarray, plt.Axes, plt.Axes, plt.Axes]:
    """
    Calculates the qmap for a given CCD geometry and energy.
    ---
    This function calculates the full qmap for a given CCD geometry and energy.
    The qmap is calculated

    Parameters
    ----------
    params : Parameters
        A instanced parameters class containing all the necessary parameters for the
        qmap calculation.

    Returns
    -------
    qmatrix : np.ndarray
        A 2D numpy array containing the calculated qmap.
    ax_rotated_ccd : plt.Axes
        A matplotlib axes object containing the qmap plot.
    ax_correct : plt.Axes
        A matplotlib axes object containing the solid angle correction plot.
    ax_scatter_angle : plt.Axes
        A matplotlib axes object containing the scattering angle plot.
    ax_hist : plt.Axes
        A matplotlib axes object containing the histogram plot.
    """

    # ---------------------------------------------------------------------------
    # Initialization of q-map parameters
    # ---------------------------------------------------------------------------

    # TODO: Auto convert this in the parameters class
    lam = 12398 / params.En  # Wavelength (A)
    phi = params.phi * (np.pi / 180)  # Angle the CCD is rotated to
    c = (
        params.beamCenterY + 0.5
    ) * params.PixelSize  # Distance from the beam center at phi = 0, to bottom of the CCD (mm)

    # Shift phi off poles for tan(phi) calculation
    if phi == np.pi / 2:
        phi += np.finfo(float).eps * 10**7

    # TODO: Add ability to have different Pixel Sizes in Parameters Class
    PixelSizeX = params.PixelSize
    PixelSizeY = params.PixelSize

    # ---------------------------------------------------------------------------
    # Calculation of q-map
    # ---------------------------------------------------------------------------

    qmatrix = np.zeros(
        (params.num_pixels, params.num_pixels)
    )  # Calculated recipicol space q value (A^-1)
    rows, cols = np.meshgrid(np.arange(params.num_pixels), np.arange(params.num_pixels))

    # Use array broadcasting to calculate 'a' and 'b' for all combinations of 'p' and 'q'
    a = (
        (params.d1 + params.d2) * (1 / np.cos(phi))
        - params.d1
        - (((cols + 0.5) * PixelSizeY) + (params.d1 + params.d2) * np.tan(phi) - c)
        * np.sin(phi)
    )
    b = (
        (((cols + 0.5) * PixelSizeY) + (params.d1 + params.d2) * np.tan(phi) - c)
        * np.cos(phi)
    ) ** 2 + ((rows - params.beamCenterX + 0.5) * PixelSizeX) ** 2

    qmatrix = (2 * np.pi * np.sqrt(2) / lam) * np.sqrt(1 - a / np.sqrt(a**2 + b))
    scatter_angle = np.arccos(a / np.sqrt(a**2 + b)) / (np.pi / 180)

    # ---------------------------------------------------------------------------
    # Solid Angle Correction
    # ---------------------------------------------------------------------------
    qmap_correct, ax_correct = solid_angle_correction(np.copy(qmatrix), params)

    # ----------------------------------------------------------------------------
    # Plotting the q-map and related plots
    # ----------------------------------------------------------------------------

    ax_rotated_ccd = plot_map("Rotated_CCD_Contour_Map", qmatrix, "q value (Å$^{-1}$)")
    ax_scatter_angle = plot_map(
        "Scattering_Angle", scatter_angle, "Scattering Angle (degrees)"
    )
    ax_hist = qmap_histogram(np.copy(qmatrix), params.bins)

    print("Maximum Scattering angle:", np.max(scatter_angle), "Degrees")

    return qmap_correct, ax_rotated_ccd, ax_correct, ax_scatter_angle, ax_hist


def solid_angle_correction(
    data_wave: np.ndarray, params: Parameters
) -> tuple[np.ndarray, plt.Axes]:
    """
    Calculates the solid angle correction for a given CCD geometry and energy.
    ---
    This function calculates the solid angle correction for a given CCD geometry and.

    Parameters
    ----------
    data_wave : np.ndarray
        A 2D numpy array containing the calculated qmap.
    params : Parameters
        A instanced parameters class containing all the necessary parameters for the
        qmap calculation.

    Returns
    -------
    solid_angle_correction_map : np.ndarray
        A 2D numpy array containing the calculated solid angle correction.
    ax_solid_angle : plt.Axes
        A matplotlib axes object containing the solid angle correction plot.
    """

    # ---------------------------------------------------------------------------
    # Initialization of q-map parameters
    # ---------------------------------------------------------------------------

    phi = params.phi * (np.pi / 180)  # Angle the CCD is rotated to

    # Shift phi off poles for tan(phi) calculation
    if phi == np.pi / 2:
        phi += np.finfo(float).eps * 10**7

    PixelSizeX = params.PixelSize  # size of a pixel (mm)
    PixelSizeY = params.PixelSize  # size of a pixel (mm)

    c = (
        params.beamCenterY * PixelSizeY
    )  # Distance from the beam center at phi = 0, to bottom of the CCD (mm)

    solid_angle_correction_map = np.zeros_like(data_wave)
    rows, cols = np.meshgrid(np.arange(params.num_pixels), np.arange(params.num_pixels))

    y1 = (
        (cols * PixelSizeY)
        + (params.d1 + params.d2 - (params.d2 + params.d1 * (1 - np.cos(phi))))
        * np.tan(phi)
        - c
    )
    y2 = y1 + PixelSizeY
    x1 = (rows - params.beamCenterX) * PixelSizeX
    x2 = x1 + PixelSizeX
    z1 = np.zeros_like(data_wave) + (params.d2 + params.d1 * (1 - np.cos(phi)))
    # z1 is the distance from the source perpendicular to the Cartesian plane the CCD is in.

    omega = (
        np.arctan(x2 * y2 / (z1 * np.sqrt(x2**2 + y2**2 + z1**2)))
        - np.arctan(x2 * y1 / (z1 * np.sqrt(x2**2 + y1**2 + z1**2)))
        - np.arctan(x1 * y2 / (z1 * np.sqrt(x1**2 + y2**2 + z1**2)))
        + np.arctan(x1 * y1 / (z1 * np.sqrt(x1**2 + y1**2 + z1**2)))
    )  # omega is the solid angle of a pixel

    omega_0 = 4 * np.arctan(
        PixelSizeX
        * PixelSizeY
        / (2 * z1 * np.sqrt(PixelSizeX**2 + PixelSizeY**2 + 4 * z1**2))
    )

    solid_angle_correction_map = omega_0 / omega  # Normalized Solid Angle Correction

    ax_solid_angle = plot_map(
        "Rotated_CCD_SA_correction_Map",
        solid_angle_correction_map,
        "Solid Angle Correction (sr$^2$)",
    )

    return solid_angle_correction_map, ax_solid_angle


def plot_map(window_name: str, matrix: np.ndarray, pixel_value: str) -> plt.Axes:
    """
    Plotting function for the qmap and solid angle correction.

    Parameters
    ----------
    window_name : str
        Name of the figure
    matrix : np.ndarray
        matrix to plot
    pixel_value : str
        colorbar label

    Returns
    -------
    ax : plt.Axes
        matplotlib axes object containing the plot.
    """
    fig, ax = plt.subplots()
    cont = ax.contour(matrix, origin="lower", levels=20, colors="black")
    im = ax.imshow(
        matrix,
        origin="lower",
        cmap="viridis",
        extent=[0, matrix.shape[1], 0, matrix.shape[0]],
    )

    ax.clabel(cont, cont.levels)
    ax.set(
        xlabel="X (pix #)",
        ylabel="Y (pix #)",
        title=f"{window_name} - Contours",
        aspect="auto",
    )
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(pixel_value)
    return ax


def qmap_histogram(input_wave: np.ndarray, num_bins: int) -> plt.Axes:
    """
    Plotting function for the qmap histogram.

    Parameters
    ----------
    input_wave : np.ndarray
        Input qmap
    num_bins : int
        Number of bins for the histogram

    Returns
    -------
    ax: plt.Axes
        matplotlib axes object containing the plot.
    """
    fig, ax = plt.subplots()

    flattened_matrix = input_wave.flatten()
    log_bins = np.logspace(
        np.log10(flattened_matrix.min()), np.log10(flattened_matrix.max()), num_bins
    )

    ax.hist(
        flattened_matrix,
        bins=log_bins,
        log=True,
        color="blue",
        alpha=0.7,
        edgecolor="black",
    )
    ax.set(
        xlabel=r"Q value ($/AA^{-1}$)",
        ylabel="# of pixels",
        title="Qmap Histogram",
        xscale="log",
    )
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_major_formatter().set_useOffset(False)

    return ax
