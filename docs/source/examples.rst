Examples
==========

EUVpy can be used to generate EUV irradiance outputs for different models and in different wavelength binning schemes.
A few examples of using its major capabilities are below.

Running Different EUV Models in General
---------------

EUVpy contains 4 different solar EUV models: EUVAC, HEUVAC, HFG, and NEUVAC. Each of these models is parameterized in
terms of the 10.7 cm solar radio flux (F10.7), or a quantity derived from F10.7. As a result, in order to run any of the
aforementioned models, the user must first obtain F10.7 measurements. This can be done as follows:

.. code-block:: python
    # Top level imports
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Qt5Agg')

    # Local imports:
    from EUVpy.tools import processIndices
    from EUVpy.NEUVAC import neuvac
    from EUVpy.empiricalModels.models.EUVAC import euvac
    from EUVpy.empiricalModels.models.HEUVAC import heuvac
    from EUVpy.empiricalModels.models.SOLOMON import solomon

    # Global Plotting Settings:
    import matplotlib.pylab as pylab
    params = {'legend.fontsize': 'large',
              'figure.figsize': (16, 8),
             'axes.labelsize': 'large',
             'axes.titlesize':'x-large',
             'xtick.labelsize':'x-large',
             'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)

    # Get some F10.7 for the entirety of 2018
    f107times, f107, f107a, f107b = processIndices.getCLSF107('2018-01-01', '2018-12-31', truncate=False)

After obtaining F10.7 data, the different models can then easily be executed:

.. code-block:: python
    # Call the models in the EUVAC bins!
    neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')
    euvacFlux, euvacIrr, _, _, _ = euvac.euvac(f107, f107a)
    heuvac_wav, heuvacFlux, heuvacIrr, _, _, _ = heuvac.heuvac(f107, f107a, torr=True)

    # Call the models in the SOLOMON bins!
    neuvacIrrSolomon, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='SOLOMON')
    solomonFluxHFG, solomonIrrHFG = SOLOMON.solomon.solomon(f107, f107a, model='HFG')
    solomonFluxEUVAC, solomonIrrEUVAC = SOLOMON.solomon.solomon(f107, f107a, model='EUVAC')

    # Plotting:
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    # Left plot will be total irradiance in the EUVAC bins
    ax[0].plot(f107times, np.sum(neuvacIrr, axis=-1), label='NEUVAC-37', color='orange')
    ax[0].plot(f107times, np.sum(euvacIrr, axis=-1), label='EUVAC-37', color='green')
    ax[0].plot(f107times, np.sum(heuvacIrr, axis=-1), label='HEUVAC-37', color='red')
    ax[0].set_ylabel('Irradiance (W/m$^2$)')
    ax[0].set_xlabel('Date')
    ax[0].set_title('Total Irradiance (EUVAC Bins)')
    ax[0].legend(loc='best')
    # Right subplot will be total irradiance in the SOLOMON bins
    solomonTable = solomon.solomonTable
    mids = 0.5*(solomonTable[:, 1] + solomonTable[:, 2])*10
    ind2 = 7
    ax[1].plot(f107times, np.sum(neuvacIrrSolomon, axis=-1), label='NEUVAC-22', color='orange')
    ax[1].plot(f107times, np.sum(solomonIrrHFG, axis=-1), label='HFG', color='purple')
    ax[1].plot(f107times, np.sum(solomonIrrEUVAC, axis=-1), label='HEUVAC-22', color='red')
    ax[1].set_ylabel('Irradiance (W/m$^2$)')
    ax[1].set_xlabel('Date')
    ax[1].set_title('Total Irradiance (SOLOMON bins)')
    ax[1].legend(loc='best')
    fig.tight_layout()
    plt.show()

The resulting picture from running the above code should look like the following figure below:

.. figure:: images/Example_1.png
   :align:  center

Time Series of Solar EUV Irradiance
---------------

It can also be useful to look at the evolution of solar EUV in a specific wavelength band over time. As in the preceding
example, we start by obtaining F10.7 data. We will consider a longer stretch of time, and restrict ourselves to
comparing EUVAC and NEUVAC:

.. code-block:: python
    # Get some F10.7 data for the entirety of Solar Cycle 24:
    f107times, f107, f107a, f107b = processIndices.getCLSF107('2008-12-01', '2019-12-31', truncate=False)

    # Call the models:
    neuvacIrradiance, _, _, _ = neuvac.neuvacEUV(f107, f107a)
    euvacFlux, euvacIrr, _, _, _ = euvac.euvac(f107, f107a)

    ind = 11
    mids = 0.5*(euvac.euvacTable[:, 1] + euvac.euvacTable[:, 2])
    plt.figure(figsize=(12, 8))
    plt.plot(f107times, neuvacIrradiance[:, ind], label='NEUVAC', color='tab:orange')
    plt.plot(f107times, euvacIrr[:, ind], label='EUVAC', color='tab:green')
    # plt.plot(f107times, heuvacIrr[:, ind], label='HEUVAC', color='tab:red')
    plt.legend(loc='best')
    plt.xlabel('Date')
    plt.ylabel('Irradiance W/m$^2$')
    plt.title('Solar Irradiance at '+str(mids[ind])+' Angstroms (SC24)')
    plt.tight_layout()
    plt.show()

The result should be the following figure:

.. figure:: images/Example_2.png
   :align:  center

Individual Solar Spectra
---------------

It can also be helpful at times to simply generate the entire spectrum for a particular model, so it may be examined.
To do so, we can simply consider some arbitrary values of F10.7, 81-day averaged F10.7, and 54-day averaged F10.7 in a
backwards-looking window:

.. code-block:: python
    # Sample values for F10.7, F10.7A, and F10.7B
    f107 = 120
    f107a = 85
    f107b = 87

Generally speaking, it's most convenient to view solar spectra in something like a `stair plot'. An example of this can
be found in Figure 8 of `Nishimoto, et al. 2021 <https://link.springer.com/article/10.1186/s40623-021-01402-7>`_. In
order to do that, we need to get the boundaries of the wavelength ranges. We can do that as follows:

.. code-block:: python
    from EUVpy.tools import toolbox
    euvacTable = euvac.euvacTable
    leftsides = euvacTable[:, 1]
    rightsides = euvacTable[:, 2]
    band_indices, band_boundaries = toolbox.band_info(leftsides, rightsides)

Let's compare the NEUVAC, EUVAC, and HEUVAC models:

.. code-block:: python
    neuvacIrr, _, _, _ = neuvac.neuvacEUV(f107, f107b, bands='EUVAC')
    euvacFlux, euvacIrr, _, _, _ = euvac.euvac(f107, f107a)
    heuvac_wav, heuvacFlux, heuvacIrr, _, _, _ = heuvac.heuvac(f107, f107a, torr=True)

    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
    ax.stairs(values=neuvacIrr[0, band_indices], edges=band_boundaries, label='NEUVAC-37', lw=3, color='tab:orange')
    ax.stairs(values=euvacIrr[0, band_indices], edges=band_boundaries, label='EUVAC-37', lw=3, color='tab:green')
    ax.stairs(values=heuvacIrr[0, band_indices], edges=band_boundaries, label='HEUVAC-37', lw=3, color='tab:red')
    ax.set_yscale('log')
    ax.legend(loc='best')
    ax.grid()
    ax.set_xlabel('Wavelength ($\mathrm{\AA}$)')
    ax.set_ylabel('Irradiance (W/m$^2$)')
    ax.set_title('Individual Solar Spectra in EUVAC Bins for (F10.7, F10.7A, F10.7B) = ('+str(f107)+', '+str(f107a)+', '+str(f107b)+')')
    plt.show()

The resulting image should look like the following:

.. figure:: images/Example_3.png
   :align:  center

Irradiance Ensembles
---------------
One of the powerful capabilities EUVpy provides is the ability to generate irradiance ensembles.

Preparing Files for Numerical Models
---------------

=================
GITM
=================

=================
Aether
=================
