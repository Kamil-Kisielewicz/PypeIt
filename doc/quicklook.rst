**********
Quick Look
**********

Overview
========

PypeIt provides a set of Quick Look scripts for
quick reductions, presumably at the telescope.
We describe each in turn.

.. _run-calcheck:

pypeit_nires_AB
===============

This script performs a quick A-B reduction of a pair of
Keck/NIRES spectral images.  Currently, the code takes
2min and 20s to process two images with boxcar extractions.
Therefore, there is a first set of nires-output_ in
approximately 1 minute.

Setup
+++++

Before running this script, you will need to

- Download the folder of `NIRES Master calibration frames <https://tinyurl.com/pypeit-nires-masters>`_.
- You may place this folder anywhere.
- Point the Environmental variable *NIRES_MASTERS* at this folder.
   - e.g. export NIRES_MASTERS=/data/Keck_NIRES/Masters_NIRES

Options
+++++++

Here is the usage::

    pypeit_nires_AB /data/Projects/Python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/AB_script/Raw s180604_0089.fits.gz s180604_0090.fits.gz -b 0.5 -h
    usage: pypeit_nires_AB [-h] [-b BOX_RADIUS] full_rawpath fileA fileB

    Script to run PypeIt on a pair of NIRES files (A-B)

    positional arguments:
      full_rawpath          Full path to the raw files
      fileA                 A frame
      fileB                 B frame

    optional arguments:
      -h, --help            show this help message and exit
      -b BOX_RADIUS, --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction


Example
+++++++

Here is an example call::

    pypeit_nires_AB /data/Keck_NIRES/Raw s180604_0089.fits.gz s180604_0090.fits.gz -b 0.5

.. _nires-output:

Output
++++++

If all goes smoothly, the code will generate four spectral
output files, with 2 each with extensions of spec1d and
spec2d.  These can be viewed with :ref:`pypeit-1dspec`
and :ref:`pypeit-2dspec`.  For the 1D spectra, because the
output is only boxcar, you need to use the --extract=BOX option.
