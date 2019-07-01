""" Object to process a single raw image"""

import inspect

import numpy as np

from pypeit import msgs
from pypeit.core import procimg
from pypeit.core import flat
from pypeit.images import pypeitimage
from pypeit.par import pypeitpar

from IPython import embed


class ProcessRawImage(pypeitimage.PypeItImage):
    """
    Class to process a raw image

    Args:
        filename (:obj:`str` or None):
            Filename
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.

    Attributes:
        steps (dict):
            Dict describing the steps performed on the image
        _bpm (np.ndarray):
            Holds the bad pixel mask once loaded
        _rawdatasec_img (np.ndarray):
            Holds the rawdatasec_img once loaded
        hdu (fits.HDUList):
            HDUList of the file
    """
    def __init__(self, filename, spectrograph, det, par):

        # Init me
        pypeitimage.PypeItImage.__init__(self, spectrograph, det)
        # Required parameters
        if not isinstance(par, pypeitpar.ProcessImagesPar):
            msgs.error("Bad par input")
        self.par = par  # ProcessImagesPar
        self.filename = filename

        # Attributes
        self._reset_internals()
        self._bpm = None
        self._rawdatasec_img = None
        self.hdu = None

        # Load
        self.load_rawframe()

        # All possible processing steps
        #  Note these have to match the method names below
        self.steps = dict(subtract_bias=False,
                          subtract_overscan=False,
                          subtract_dark=False,
                          trim=False,
                          apply_gain=False,
                          orient=False,
                          flatten=False,
                          )
    @property
    def amps(self):
        """
        Return a list of the amplifier indices, 1-indexed

        Returns:
            list
        """
        return np.unique(self.rawdatasec_img[self.rawdatasec_img > 0]).tolist()

    @property
    def bpm(self):
        """
        Generate and return the bad pixel mask for this image
        Warning:  BPM masks are for processed (e.g. trimmed, rotated) images only!

        Returns:
            np.ndarray:  Bad pixel mask with a bad pixel = 1

        """
        if self._bpm is None:
            self._bpm = self.spectrograph.bpm(shape=self.image.shape,
                                    filename=self.filename,
                                    det=self.det)
        return self._bpm

    # TODO all of these steps below should be consoliated into one method which reads the files. It is silly and extremely
    #  slow to re-read the images every time for each one of these stpes
    @property
    def rawdatasec_img(self):
        """
        Generate and return the datasec image in the Raw reference frame

        Returns:
            np.ndarray

        """
        if self._rawdatasec_img is None:
            self._rawdatasec_img = self.spectrograph.get_rawdatasec_img(self.filename, self.det)
        return self._rawdatasec_img


    @property
    def oscansec_img(self):
        """
        Generate and return the oscansec image

        Returns:
            np.ndarray

        """
        oimg = self.spectrograph.get_oscansec_img(self.filename, self.det)
        return oimg

    def _reset_steps(self):
        """
        Reset all the processing steps to False

        Should consider setting the Image to None too..
        """
        for key in self.steps.keys():
            self.steps[key] = False

    def _reset_internals(self):
        """
        Init or free up memory by resetting the Attributes to None
        """
        self.rawvarframe = None
        self.crmask = None
        self.mask = None
        self.rn2img = None

    def apply_gain(self, force=False):
        """
        Apply the Gain values to self.image

        Args:
            force (bool, optional):

        Returns:
            np.ndarray:  copy of self.image

        """
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Gain was already applied. Returning")
            return self.image.copy()

        gain = np.atleast_1d(self.spectrograph.detector[self.det - 1]['gain']).tolist()
        # Apply
        self.image *= procimg.gain_frame(self.rawdatasec_img, gain, trim=self.steps['trim'])
        self.steps[step] = True
        # Return
        return self.image.copy()

    def process(self, process_steps, pixel_flat=None, illum_flat=None,
                       bias=None, bpm=None, debug=False):
        """
        Process the image

        Note:  The processing steps are currently 'frozen' as is.
          We may choose to allow optional ordering of the steps

        Args:
            process_steps (list):
                List of processing steps
            pixel_flat (np.ndarray, optional):
                Pixel flat image
            illum_flat (np.ndarray, optional):
                Illumination flat
            bias (np.ndarray, optional):
                Bias image
            bpm (np.ndarray, optional):
                Bad pixel mask image

        """
        # Standard order
        #   -- May need to allow for other order some day..
        if 'subtract_overscan' in process_steps:
            self.subtract_overscan()
        if 'trim' in process_steps:
            self.trim()
        if 'subtract_bias' in process_steps: # Bias frame, if it exists, is trimmed
            self.subtract_bias(bias)
        if 'apply_gain' in process_steps:
            self.apply_gain()
        if 'orient' in process_steps:
            self.orient()
        # Flat field
        if 'flatten' in process_steps:
            self.flatten(pixel_flat, illum_flat=illum_flat, bpm=bpm)
        # Return copy of the image
        return self.image.copy()

    def flatten(self, pixel_flat, illum_flat=None, bpm=None, force=False):
        """
        Flat field the image

        Wrapper to flat.flatfield

        Args:
            pixel_flat (np.ndarray):
                Pixel flat image
            illum_flat (np.ndarray, optional):
                Illumination flat image
            bpm (np.ndarray, optional):
                Bad pixel mask image;  if provided, over-rides internal one
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already flat fielded.  Returning the current image")
            return self.image.copy()
        # BPM
        if bpm is None:
            bpm = self.bpm
        # Do it
        self.image = flat.flatfield(self.image, pixel_flat, bpm, illum_flat=illum_flat)
        self.steps[step] = True

    def load_rawframe(self):
        """
        Load a raw image from disk using the Spectrograph method load_raw_frame()

        Also loads up the binning, exposure time, and header of the Primary image
        And the HDUList in self.hdu

        Args:
            filename (str):  Filename

        """
        # Load
        self.image, self.hdu, \
            = self.spectrograph.load_raw_frame(self.filename, det=self.det)
        self.head0 = self.hdu[0].header
        # Shape
        self.orig_shape = self.image.shape
        # Exposure time
        self.exptime = self.spectrograph.get_meta_value(self.filename, 'exptime')
        # Binning
        self.binning = self.spectrograph.get_meta_value(self.filename, 'binning')
        if self.spectrograph.detector[self.det-1]['specaxis'] == 1:
            self.binning_raw = (',').join(self.binning.split(',')[::-1])
        else:
            self.binning_raw = self.binning

    def orient(self, force=False):
        """
        Orient the image in the PypeIt format with spectra running blue (down)
        to red (up).

        Args:
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Orient the image to have blue/red run bottom to top
        # Check if already oriented
        if self.steps[step] and not force:
            msgs.warn("Image was already oriented.  Returning current image")
            return self.image.copy()
        # Orient me
        self.image = self.spectrograph.orient_image(self.image, self.det)
        self.steps[step] = True

    def subtract_bias(self, bias_image, force=False):
        """
        Perform bias subtraction

        Args:
            bias_image (np.ndarray):
                Bias image
            force (bool, optional):
                Force the processing even if the image was already processed
        """
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already bias subtracted.  Returning the current image")
            return self.image.copy()
        # Do it
        self.image -= bias_image
        self.steps[step] = True

    def subtract_overscan(self, force=False):
        """
        Analyze and subtract the overscan from the image

        Args:
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already trimmed")

        temp = procimg.subtract_overscan(self.image, self.rawdatasec_img, self.oscansec_img,
                                         method=self.par['overscan'],
                                         params=self.par['overscan_par'])
        # Fill
        self.steps[step] = True
        self.image = temp

    def trim(self, force=False):
        """
        Trim the image to include only the science data

        Args:
            force (bool, optional):
                Force the processing even if the image was already processed

        """
        step = inspect.stack()[0][3]
        # Check input image matches the original
        if self.orig_shape is not None:
            if self.image.shape != self.orig_shape:
                msgs.warn("Image shape does not match original.  Returning current image")
                return self.image.copy()
        # Check if already trimmed
        if self.steps[step] and (not force):
            msgs.warn("Image was already trimmed.  Returning current image")
            return self.image
        # Do it
        trim_image = procimg.trim_frame(self.image, self.rawdatasec_img < 1)
        # Overwrite
        self.image = trim_image
        self.steps[step] = True

    def __repr__(self):
        return ('<{:s}: file={}, steps={}>'.format(
            self.__class__.__name__, self.filename, self.steps))

