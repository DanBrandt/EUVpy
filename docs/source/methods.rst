Methods
==========

The primary purpose of EUVpy is to make the most widely-used Solar EUV models easily available in a single Python package.
These models include `EUVAC <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94ja00518>`_,
`HEUVAC <https://www.sciencedirect.com/science/article/pii/S0273117705008288>`_,
`HFG <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/gl008i011p01147>`_, and
`NEUVAC <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2024SW004043>`_.

The Physics of Solar EUV
----------

The Solar Extreme Ultraviolet Spectrum is receives contributions from a wide variety of sources on the Sun. It strongly
influences space weather on Earth, with 80\% of it being absorbed by the atmosphere above 90 km. This essentially prevents
measurements from being taken by ground-based instruments. Space-based instruments like TIMED/SEE and SORCE/SOLSTICE have
measured the solar spectrum, but the need to process space-based measurements means that solar proxies like F10.7 are
still widely-used for space situational awareness.

The sources of Solar EUV cover a diverse range of altitudes on the Sun, extending from the Photosphere to the Corona.
The emission processes that generate EUV involve multiple charged ions colliding with electrons, which generates
additional ionizations if the energy exceeds a threshold. These ionizations are generally of the form:

.. math::
    \mathrm{X}^{+\mathrm{m}} + \mathrm{e}^{-} \rightarrow \mathrm{X}^{+\mathrm{m}+1} + 2\mathrm{e}^{-}.

The ion left in an excited state then proceeds to a lower energy state via photo emission. Radiative recombination also
plays a role in generating photon emissions through reactions of the following form:

.. math::
    \mathrm{X}^{+\mathrm{m}} + \mathrm{e}^{-} \rightarrow \mathrm{X}^{+\mathrm{m}-1}.

Radiative recombination in particular involves the transfer of a free electron to a captured electron (referred to as a
free-bound process). In the case of free-free emission (i.e. "Bremsstrahlung Radiation"), the electron is accelerated
or decelerated through intraction with the Coloumb potential of an ion. Both free-bound and free-free emissions generate
the Solar EUV continuum, a continuous baseline in the solar spectrum. It is conventional (as in the case of the Flare
Irradiance Spectrum Model) to treat the continuum as being comprised of this background, atop which there are variations
due to variations in the ~11-year Solar Cycle, the rotation of the Sun, as well as due to Solar Flares:

.. math::
    \Phi_{total}(t,\lambda) = \Phi_{continuum}(\lambda) + \Phi_{SC}(t,\lambda) + \Phi_{SR}(t,\lambda) + \Phi_{Flare}(t,\lambda).

Please see `Lilensten, et al. 2008 <https://angeo.copernicus.org/articles/26/269/2008/angeo-26-269-2008.pdf>`_ and
`Chamberlin, et al. 2020 <https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2020sw002588>`_ for more reading on
the Solar EUV spectrum and the FISM model, respectively.

NEUVAC
----------

A challenge in Solar EUV modeling, especially for space weather operations, is the determination of the solar EUV
spectrum quickly and from the ground. Solar proxies, which are correlated strongly with EUV, have historically been used
for this purpose, but their limitations have motivated ways of reconstructing the Solar EUV spectrum instead of relying
on proxies alone.

One such approach is that of the NEUVAC model (`Brandt and Ridley, 2024 <https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2024SW004043>`_),
which reconstructs the solar irradiance of the FISM-2 model from the F10.7 solar proxy alone. The version of NEUVAC
within EUVpy reconstructs FISM-2 specifically because the latter currently represents the community standard of EUV
estimates (since the FISM-2 data is a composite dataset built from measurements from a wide variety of missions), yet
its algorithm employs a 108-day average in its model formulation. It therefore only can provide reliable EUV estimates
54 days into the past, with any estimates more recent than that being provisional. NEUVAC, by comparison, can provide
estimates on the present day.

NEUVAC is able to faithfully reconstruct the FISM-2 irradiance by using a non-linear parametric representation:

.. math::
    \Phi_{total}(t, \lambda) = A_{\lambda}\left[F_{10.7}^{B_{\lambda}}(t)\right] + C_{\lambda}\left[F_{10.7B}^{D_{\lambda}}(t)\right] + E_{\lambda}\left[F_{10.7B}(t) - F_{10.7}(t)\right] + F_{\lambda}

With this representation, NEUVAC is able to capture the variance in EUV in each wavelength band $\lambda$, much better
than linear models which have seen widespread historical use. This is critical especially for wavelength bands near the
soft X-ray, which are incredibly non-stationary and receive contamination from flare emissions. Additionally, NEUVAC
uses Monte Carlo methods to provide quantified uncertainties on irradiance estimates in each wavelength band.



