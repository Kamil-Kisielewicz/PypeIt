""" Module for the ScienceImage class"""
import numpy as np
from pypeit import msgs
from pypeit import utils
from pypeit import ginga
from pypeit.core import coadd2d
from pypeit.core import procimg

from pypeit.images import combinedimage
from pypeit.images import processimage

from IPython import embed

class ScienceImage(object):
    """
    This class will organize and run actions related to
    a Science or Standard star exposure

    ..todo.. Clean this up JFH

    Parameters
    ----------
    file_list : list
      List of raw files to produce the flat field
    spectrograph : pypeit.spectrograph.Spectrograph
    settings : dict-like
    tslits_dict : dict
      dict from TraceSlits class
    tilts : ndarray
      tilts from WaveTilts class
      used for sky subtraction and object finding
    det : int
    sci_bpm : ndarray
      Bad pixel mask for this science image;  can and often does differ
      from the default BPM
    objtype : str
      'science'
      'standard'

    Attributes
    ----------
    frametype : str
      Set to 'science'
    sciframe : ndarray
      Processed 2D frame
    rawvarframe : ndarray
      Variance generated without a sky (or object) model
    modelvarframe : ndarray
      Variance generated with a sky model
    finalvar : ndarray
      Final variance frame
    global_sky : ndarray
      Sky model across the slit/order
    skycorr_box : ndarray
      Local corrections to the sky model
    final_sky : ndarray
      Final sky model; may include 'local' corrections
    obj_model : ndarray
      Model of the object flux
    trcmask : ndarray
      Masks of objects for sky subtraction
    tracelist : list
      List of traces for objects in slits
    inst_name : str
      Short name of the spectrograph, e.g. KASTb
    target_name : str
      Parsed from the Header
    basename : str
      Combination of camera, target, and time
      e.g. J1217p3905_KASTb_2015May20T045733.56
    time : Time
      time object
    specobjs : list
      List of specobjs
    bm: ScienceImageBitMask
      Object used to select bits of a given type
    """

    # Frametype is a class attribute
    frametype = 'science'

    # TODO: Merge into a single parset, one for procing, and one for scienceimage
    def __init__(self, spectrograph, file_list, bg_file_list=[], ir_redux=False, det=1, binning=None, par=None):


        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!
        #self.par = pypeitpar.FrameGroupPar(objtype) if par is None else par
        self.par = spectrograph.default_pypeit_par()['scienceframe'] if par is None else par

        # Start us up
        self.sci_combine = combinedimage.CombinedImage(spectrograph, det, self.par['process'],
                                                       files=file_list, frametype=self.frametype)
        self.bkg_combine = combinedimage.CombinedImage(spectrograph, det, self.par['process'],
                                                       files=bg_file_list, frametype=self.frametype)

        # Instantiation attributes for this object
        self.spectrograph = spectrograph
        # Are we subtracing the sky using background frames? If yes, set ir_redux=True
        self.ir_redux = ir_redux
        if self.ir_redux and (self.bkg_combine.nimages == 0):
            msgs.error('IR reductions require that bg files are specified')
        self.det = det
        self.binning = binning

        # Set some detector parameters that we will need
        self.saturation = self.spectrograph.detector[self.det - 1]['saturation']
        self.mincounts = self.spectrograph.detector[self.det - 1]['mincounts']

        # These attributes will be sert when the image(s) are processed
        self.sci_bpm = None
        self.bias = None
        self.pixel_flat = None
        self.illum_flat = None

        self.steps = []

        # Other bookeeping internals
        self.crmask = None
        self.mask = None


    # JFH TODO This stuff should be eventually moved to processimages?
    def proc(self, bias, pixel_flat, bpm, illum_flat=None, sigrej=None, maxiters=5, show=False):
        """
        Primary wrapper for processing one or more science frames or science frames with bgframes

        Args:
            bias (ndarray, None or str):  Specifies bias subtraction approach and/or provides bias image
            pixel_flat (ndarray):  Pixel flat image
            bpm (ndarray):  Bad pixel mask
            illum_flat (ndarray, optional): Illumination flat
            sigrej (int or float, optional): Rejection threshold for sigma clipping.
                 Code defaults to determining this automatically based on the numberr of images provided.
            maxiters (int, optional):
            show (bool, optional):

        Returns:
            ndarray, ndarray, ndarray, ndarray, ndarray:
              sciimg
              sciivar
              rn2img
              mask
              crmask

        """

        # Process
        self.sci_bpm = bpm
        self.bias = bias
        self.pixel_flat = pixel_flat
        self.illum_flat = illum_flat

        if self.ir_redux:
            self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = self.proc_diff(
                reject_cr=True, sigma_clip=False, sigrej=sigrej, maxiters=maxiters)
        else:
            self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = self.proc_list(
                'sci', reject_cr=True, sigma_clip=False, sigrej=sigrej, maxiters=maxiters)

        # Show the science image if an interactive run, only show the crmask
        if show:
            # Only mask the CRs in this image
            self.show(self.sciimg * (self.crmask == 0), chname='sciimg')

        return self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask

    def proc_list(self, ltype, reject_cr=True, sigma_clip=False, sigrej=None, maxiters=5):
        """
        Process a list of science images

        This includes stacking the images if there is more than 1

        Args:
            ltype (str): Type of images to process ('sci', 'bkg')
            reject_cr (bool, optional):
            sigrej (int or float, optional): Rejection threshold for sigma clipping.
                 Code defaults to determining this automatically based on the numberr of images provided.
            maxiters (int, optional):
            show (bool, optional):

        Returns:
            ndarray, ndarray, ndarray, ndarray, ndarray:
              sciimg
              sciivar
              rn2img
              mask
              crmask

        """
        # Init
        if ltype == 'sci':
            combinedImage = self.sci_combine
        elif ltype == 'bkg':
            combinedImage = self.bkg_combine
        else:
            msgs.error("Bad ltype for proc_list")
        nimg = combinedImage.nimages
        weights = np.ones(nimg)/float(nimg)

        # Load
        #embed(header='198 of sciimg')
        combinedImage.load_images()
        # Process
        process_steps = procimg.init_process_steps(self.bias, self.par['process'])
        process_steps += ['trim', 'apply_gain']
        if (self.pixel_flat is not None) or (self.illum_flat is not None):
            process_steps += ['flatten']
        combinedImage.process_images(process_steps, bias=self.bias, bpm=self.sci_bpm,
                            pixel_flat=self.pixel_flat, illum_flat=self.illum_flat)
        # Load up the stack (and create some internal images0
        img_stack, ivar_stack, rn2img_stack, crmask_stack, mask_stack = combinedImage.build_stack(
            bpm=self.sci_bpm, reject_cr=reject_cr)

        # ToDO The bitmask is not being properly propagated here!
        if nimg > 1:
            img_list = [img_stack]
            var_stack = utils.calc_ivar(ivar_stack)
            var_list = [var_stack, rn2img_stack]
            img_list_out, var_list_out, outmask, nused = coadd2d.weighted_combine(
                weights, img_list, var_list, (mask_stack == 0),
                sigma_clip=sigma_clip, sigma_clip_stack = img_stack, sigrej=sigrej, maxiters=maxiters)
            img = img_list_out[0]
            ivar = utils.calc_ivar(var_list_out[0])
            rn2img = var_list_out[1]
            # assumes everything masked in the outmask is a CR in the individual images
            crmask = np.invert(outmask)
            # Create a mask for this combined image
            processImage = processimage.ProcessImage(None, self.spectrograph, self.det, self.proc_par)
            processImage.image = img
            processImage.rawvarframe = var_list_out[0]
            processImage.crmask = crmask
            mask = processImage.build_mask(bpm=self.sci_bpm, saturation=self.saturation, mincounts=self.mincounts)
        else:
            mask = mask_stack[0, :, :]
            crmask = crmask_stack[0, :, :]
            img = img_stack[0, :, :]
            ivar = ivar_stack[0, :, :]
            rn2img = rn2img_stack[0, :, :]

        return img, ivar, rn2img, mask, crmask


    def proc_diff(self, file_list, bg_file_list, reject_cr=True,
                  sigma_clip=False, sigrej=None, maxiters=5):
        """
        Process a list of science images and their background frames
        Primarily for near-IR reductions

        Wrapper to proc_sci for

        Needed in part to set self.sciframe, although I could kludge it another way..

        Args:
            file_list:
            bg_file_list:
            reject_cr:
            sigma_clip:
            sigrej:
            maxiters:

        Returns:
            tuple: sciimg, sciivar, rn2img, mask, crmask

        """

        sciimg_sci, sciivar_sci, rn2img_sci, mask_sci, crmask_sci = self.proc_list('sci',
            file_list, reject_cr=reject_cr, sigma_clip=sigma_clip, sigrej=sigrej,maxiters=maxiters)
        sciimg_bg, sciivar_bg, rn2img_bg, mask_bg, crmask_bg = self.proc_list('bkg',
            bg_file_list, reject_cr=reject_cr, sigma_clip=sigma_clip, sigrej=sigrej,maxiters=maxiters)

        # Combine the images
        outmask_comb = (mask_sci == 0) & (mask_bg == 0)
        sciimg = sciimg_sci - sciimg_bg
        varcomb = utils.calc_ivar(sciivar_sci) + utils.calc_ivar(sciivar_bg)
        sciivar = utils.calc_ivar(varcomb)*outmask_comb
        rn2img = rn2img_sci + rn2img_bg
        # Let's do some more processing
        processImage = processimage.ProcessImage(None, self.spectrograph, self.det, self.proc_par)
        processImage.image = sciimg
        processImage.rawvarframe = varcomb
        # Now reject CRs again on the differenced image
        crmask_diff = processImage.build_crmask()
        # crmask_eff assumes evertything masked in the outmask_comb is a CR in the individual images
        crmask = crmask_diff | np.invert(outmask_comb)
        # Create a mask for this image now
        mask = processImage.build_mask(bpm=self.sci_bpm, saturation=self.saturation)#, mincounts=self.mincounts)
        #mask = self.build_mask(sciimg, sciivar, crmask, self.sci_bpm, saturation=self.saturation)

        return sciimg, sciivar, rn2img, mask, crmask


    def show(self, image, chname=None):
        """
        Show one of the internal images

        Args:
            image : ndarray, optional
              User supplied image to display

        """

        ch_name = chname if chname is not None else 'image'
        viewer, ch = ginga.show_image(image, chname=ch_name)



    def __repr__(self):
        txt = '<{:s}: nimg={:d}'.format(self.__class__.__name__,
                                        self.nsci)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt



