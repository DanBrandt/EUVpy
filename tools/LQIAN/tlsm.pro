;+
; NAME:
;       tlsm
;
; PURPOSE:
;       To read solar spectra on 1 nm resolution and output to the low-resolution
;       binning scheme for NCAR/TGCM
;
; CATEGORY:
;       Main program.
;
; Notes:
;       The low resolution scheme has logical bins that are grouped
;       based on N2 absorption cross section. Since major species'
;       cross sections have lots of structures in the wave length 
;       range of the logical bins (650 A - 975 A), it is required
;       that input spectra have resolutions that closely matches
;       teh resolution of cross section data for the grouping of logical bins to 
;       make sense. Hinteregger spectrum has resolution that closely
;       matches the resolution of cross section data, therefore, all input spectra
;       will be first scaled to to Hinteregger spectrum resolution for
;       wave length 100 A - 1000 A, using  Hinteregger 
;       reference spectrum. 
;
;       Units of all input specta are converted to the same when 
;       read in: Angstrom for wave length, photon cm^-2 s^-1 for 
;       solar flux.
;
; ROUTINES CALLED:
;       get_spectra, nm_to_hs, rebin, put_spectra
;       
;
; select binning scheme
;
bin_file='./input/low_bins.txt'
;
; Select input solar spectrum
;
solar_file='./input/see__L3_merged_2009316_010.ncdf'
;solar_file='./input/eve2008_min.dat'
;solar_file='./input/FISM_60sec_2005017_v00_00.nc'
;
; read input solar spectra. 
spectra1=get_spectra(solar_file)
;
; scale to 1 nm binning scheme 
spectra2=hi_to_nm(spectra1)
;
; scale to Hinteregger resolution
in_spectra=nm_to_hs(spectra2)
;
; rebin based on selected binning scheme.
rebin, bin_file,in_spectra, out_spectra
;
; output
put_spectra, solar_file, out_spectra

end
