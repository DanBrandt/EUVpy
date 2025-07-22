.. EUVpy documentation master file, created by
   sphinx-quickstart on Tue Jul 22 13:58:54 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EUVpy documentation
===================

.. figure:: images/EUVpy_logo.png
   :align:  center

This Python 3 code contains four major Solar Extreme Ultraviolet (EUV) models for use by atmospheric models primarily
used in the space weather community (such as GITM, TIE-GCM, and WACCM-X).

**Note:** At present, the code has only been tested on Ubuntu. In principle, it should work on Windows and Mac as well.

This repo contains the code to run each model and also generate Solar EUV model outputs in a variety of binning schemes,
most notably those employed by `Solomon and Qian, 2006 <https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2005JA011160>`_,
by `Richards, et al. 1994 <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94ja00518>`_, and those natively used
by the `Aether Model <https://aetherdocumentation.readthedocs.io/en/latest/>`_ developed by Dr. Aaron Ridley at the
University of Michigan.

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   documentation
   methods
   usage
   examples
   EUVpy
