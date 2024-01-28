from matplotlib.pylab import f


class Constant:
    def __init__(self, value, units=None):
        self.value = value
        if units is None:
            self.units = ""
        else:
            self.units = units

    def __repr__(self):
        return f"{self.value} {self.units}"

    def __str__(self):
        return self.__repr__()


class Parameters:
    """
    Experimental Parameters Class for PyRSoXS

    Parameters
    ----------
    En : float
        Photon Energy in eV
    phi : float
        Angle the CCD is rotated to (Input: degrees)
    d1 : float
        Sample - Detector distance (mm)
    d2 : float
        Distance from axis of rotation to sample (mm)
    beamCenterY : int
        y pixel value of the beam center
    beamCenterX : int
        x pixel value of the beam center
    PixelSize : int
        Size of a pixel (mm)
    num_pixels : int
        Number of pixels along one dimension of the image
    bins : int
        Number of bins used for Histogram

    See Also
    --------
    Constant : Class to hold a constant value and its units

    Notes:
    ------
    This class is used to hold all the experimental parameters needed to calculate
    the qmap and solid angle correction map. The parameters are stored as Constant
    objects which hold the value and units of the parameter. This class also has
    properties for each parameter to make it easier to access the value and units

    Future implementations will include a way to save and load the parameters from
    a experimental data file such a a .fits, .tif, or .h5 file

    Examples:
    ---------
    >>>     Parameters(En = 320, phi = 17.5, d1 = 26, d2 = 23, beamCenterY = 556, beamCenterX = 529, PixelSize = 0.027, num_pixels = 1024, bins = 200)
    >>>     print(Parameters)
    >>>     Parameter        Value   Units
            -----------------------------------
            En               320     eV
            d1               26      mm
            d2               23      mm
            beamCenterY      556     pixels
            beamCenterX      529     pixels
            phi              17.5    degrees
            PixelSize        0.027   mm
            num_pixels       1024    pixels
            bins             200     bins
    """

    # ---------------------------------------------------------------------------
    # Initialize Parameters
    # ---------------------------------------------------------------------------
    def __init__(
        self,
        En: float,
        phi: float,
        d1: float,
        d2: float,
        beamCenterY: int,
        beamCenterX: int,
        PixelSize: float,
        num_pixels: int,
        bins: int,
    ):
        """
        Class to hold all the ccd parameters needed to calculate the qmap and solid angle correction map

        >>>  Parameters(En = 320, phi = 17.5, d1 = 26, d2 = 23, beamCenterY = 556, beamCenterX = 529, PixelSize = 0.027, num_pixels = 1024, bins = 200)
        >>>  print(Parameters)
        >>>



        Args:
            En (float): Photon Energy in eV
            phi (float): Distance from axis of rotation to sample (mm)
            d1 (float): Sample - Detector distance (mm)
            d2 (float): y pixel value of the beam center
            beamCenterY (int): x pixel value of the beam center
            beamCenterX (int): Angle the CCD is rotated to (Input: degrees)
            PixelSize (int): Size of a pixel (mm)
            num_pixels (int): Number of pixels along one dimension of the image
            bins (int): Number of bins used for Histogram
        """
        self._En: Constant = Constant(En, "eV")
        self._d1: Constant = Constant(d1, "mm")
        self._d2: Constant = Constant(d2, "mm")
        self._beamCenterY: Constant = Constant(beamCenterY, "pixels")
        self._beamCenterX: Constant = Constant(beamCenterX, "pixels")
        self._phi: Constant = Constant(phi, "degrees")
        self._PixelSize: Constant = Constant(PixelSize, "mm")
        self._num_pixels: Constant = Constant(num_pixels, "pixels")
        self._bins: Constant = Constant(bins)

    # ---------------------------------------------------------------------------
    # Properties
    # ---------------------------------------------------------------------------
    @property
    def En(self):
        return self._En.value

    @En.setter
    def En(self, value, units="eV"):
        self._En = Constant(value, units)

    @property
    def d1(self):
        return self._d1.value

    @d1.setter
    def d1(self, value, units="mm"):
        self._d1 = Constant(value, units)

    @property
    def d2(self):
        return self._d2.value

    @d2.setter
    def d2(self, value, units="mm"):
        self._d2 = Constant(value, units)

    @property
    def beamCenterY(self):
        return self._beamCenterY.value

    @beamCenterY.setter
    def beamCenterY(self, value, units="pixels"):
        self._beamCenterY = Constant(value, units)

    @property
    def beamCenterX(self):
        return self._beamCenterX.value

    @beamCenterX.setter
    def beamCenterX(self, value, units="pixels"):
        self._beamCenterX = Constant(value, units)

    @property
    def phi(self):
        return self._phi.value

    @phi.setter
    def phi(self, value, units="degrees"):
        self._phi = Constant(value, units)

    @property
    def PixelSize(self):
        return self._PixelSize.value

    @PixelSize.setter
    def PixelSize(self, value, units="mm"):
        self._PixelSize = Constant(value, units)

    @property
    def num_pixels(self):
        return self._num_pixels.value

    @num_pixels.setter
    def num_pixels(self, value, units="pixels"):
        self._num_pixels = Constant(value, units)

    @property
    def bins(self):
        return self._bins.value

    @bins.setter
    def bins(self, value, units="bins"):
        self._bins = Constant(value, units)

    # ---------------------------------------------------------------------------
    # dunders
    # ---------------------------------------------------------------------------

    def __repr__(self):
        s = f"{'Parameter':<15}\t {'Value'}\t {'Units'}\n"
        s += "-" * 35 + "\n"
        for k, v in self.__dict__.items():
            key = k.strip("_")
            s += f"{key:<15}\t {v.value}\t {v.units}\n"
        return s

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if isinstance(other, Parameters):
            return self.__dict__ == other.__dict__
        return False


if __name__ == "__main__":
    test = Parameters(
        En=320,
        phi=17.5,
        d1=26,
        d2=23,
        beamCenterY=556,
        beamCenterX=529,
        PixelSize=0.027,
        num_pixels=1024,
        bins=200,
    )
    print(test)
