""" Module for Keck/ESI specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit import utils
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import pixels
from pypeit.core import parse


from IPython import embed

class KeckESISpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRES specific code
    """
    def __init__(self):
        # Get it started
        super(KeckESISpectrograph, self).__init__()
        self.spectrograph = 'keck_esi'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'ESI'
        self.numhead = 1
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.1542,
                            darkcurr        = 0.01,
                            saturation      = 64000., # I'm not sure we actually saturate with the DITs???
                            nonlinear       = 1.00,
                            numamplifiers   = 2,
                            gain            = [1.29, 1.29],
                            ronoise         = [2.7, 2.7],
                            datasec         = ['[:,25:1048]', '[:,1049:2072]'],
                            oscansec        = ['[:,2075:2149]', '[:,2155:2229]'],
                            )]
        self.norders = 10
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @property
    def pypeline(self):
        return 'Echelle'

    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_esi'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= [3,4,4,4,4]
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] = 10.0
        #par['calibrations']['tilts']['spat_order'] =  3
        #par['calibrations']['tilts']['spec_order'] =  3

        # Flats
        #par['calibrations']['flatfield']['illumflatten'] = False

        # Extraction
        #par['scienceimage']['bspline_spacing'] = 0.8
        #par['scienceimage']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['method'] = 'skip'

        #par['scienceframe']['process']['sigclip'] = 20.0
        #par['scienceframe']['process']['satpix'] ='nothing'

        # Overscan but not bias
        #  This seems like a kludge of sorts
        #par['calibrations']['biasframe']['useframe'] = 'none'

        # Set the default exposure time ranges for the frame typing
        #par['calibrations']['standardframe']['exprng'] = [None, 20]
        #par['calibrations']['arcframe']['exprng'] = [20, None]
        #par['calibrations']['darkframe']['exprng'] = [20, None]
        #par['scienceframe']['exprng'] = [20, None]
        return par

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA')
        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='TARGNAME')
        meta['decker'] = dict(ext=0, card='SLMSKNAM')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='ELAPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='INSTRUME')
        meta['idname'] = dict(ext=0, card='OBSTYPE')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
        else:
            msgs.error("Not ready for this compound meta")


    def configuration_keys(self):
        return ['decker','dispname']

    def pypeit_file_keys(self):
        pypeit_keys = super(KeckESISpectrograph, self).pypeit_file_keys()
        #pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        pypeit_keys += ['comb_id']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        # TODO: Arcs, tilts, darks?
        if ftype in ['pinhole']:
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['bias']:
            return fitstbl['idname'] == 'Bias'
        if ftype in ['pixelflat', 'trace']:
            return np.any([fitstbl['idname'] == 'DmFlat', fitstbl['idname'] == 'IntFlat'],axis=0)
        if ftype in ['arc', 'tilt']:
            return fitstbl['idname'] == 'Line'
        if ftype in ['science', 'standard']:
            return fitstbl['idname'] == 'Object'
        #
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
        """
        Override parent bpm function with BPM specific to X-Shooter VIS.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        msgs.info("Custom bad pixel mask for NIRES")
        self.empty_bpm(shape=shape, filename=filename, det=det)
        #if det == 1:
        #    self.bpm_img[:, :20] = 1.
        #    self.bpm_img[:, 1000:] = 1.
        return self.bpm_img

    def slit_minmax(self, nfound, binspectral=1):

        # These are the order boundaries determined by eye by JFH. 2025 is used as the maximum as the upper bit is not illuminated
        # Here is the info for all the orders for a good flat
        all_spec_min = np.full(self.norders, -np.inf)
        all_spec_max = np.full(self.norders, np.inf)

        # If the number of slits is less than expected, then take the reddest
        spec_min = all_spec_min[-nfound:]
        spec_max = all_spec_max[-nfound:]

        return spec_min, spec_max

    def slitmask(self, tslits_dict, pad=None):
        """
         Generic routine ton construct a slitmask image from a tslits_dict. Children of this class can
         overload this function to implement instrument specific slitmask behavior, for example setting
         where the orders on an echelle spectrograph end

         Parameters
         -----------
         tslits_dict: dict
            Trace slits dictionary with slit boundary information

         Optional Parameters
         pad: int or float
            Padding of the slit boundaries
         binning: tuple
            Spectrograph binning in spectral and spatial directions

         Returns
         -------
         slitmask: ndarray int
            Image with -1 where there are no slits/orders, and an integer where there are slits/order with the integer
            indicating the slit number going from 0 to nslit-1 from left to right.

         """

        # These lines are always the same
        pad = tslits_dict['pad'] if pad is None else pad
        slitmask = pixels.slit_pixels(tslits_dict['lcen'], tslits_dict['rcen'], tslits_dict['nspat'], pad=pad)

        return slitmask

    def slit2order(self, slit_spat_pos):
        """
        This routine is only for fixed-format echelle spectrographs.
        It returns the order of the input slit based on its slit_pos

        Args:
            slit_spat_pos (float):  Slit position (spatial at 1/2 the way up)

        Returns:
            int: order number

        """
        order_spat_pos = np.array([0.22773035, 0.40613574, 0.56009658,
                                   0.70260714, 0.86335914])
        orders = np.arange(7, 2, -1, dtype=int)
        # Find closest
        iorder = np.argmin(np.abs(slit_spat_pos-order_spat_pos))

        # Check
        if np.abs(order_spat_pos[iorder] - slit_spat_pos) > 0.05:
            msgs.error("Bad echelle input for {:s}".format(self.spectrograph))

        # Return
        return orders[iorder]


    def order_platescale(self, order_vec, binning=None):
        """
        NIRES has no binning

        Args:
            order_vec (np.ndarray):
            binning (optional):

        Returns:
            np.ndarray:

        """
        norders = order_vec.size
        return np.full(norders, 0.1542)

