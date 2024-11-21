;+
; NAME:
;	get_spectra
;
; PURPOSE:
;       Given an solar spectra file, read the solar spectra,
;       convert wave length unit to Angstrom and solar flux 
;       unit to photon cm^-2 s^-1 if necessary. 
;
; CATEGORY:
;       utility function called by main program.
;
; CALLING SEQUENCE:
;       out_spectra=get_spectra(solar_file)
;
; INPUTS:
;       solar_file: (='dir_name/file_name')
;
; OUTPUTS:
;       out_spectra: two dimensional array that holds 
;                     output solar spectra. The array 
;                     has at least 3 columns: waves and wavel 
;                     in Angstrom, solar flux in
;                     photon cm^-2 s^-1, any number of
;                     extra flux columns. 
;
; PROCEDURE:
;       some usual/useful conversion of units:
;       1 nm = 10 A
;       1 mW m^-2 = 1 erg cm^-2 s^-1
;       flux(photon cm^-2 s^-1) = (6.242e11/12398)*wave(A) *flux(erg cm^-2 s^-1) 
;
; ROUTINES CALLED:
;       read_gen, read_netcdf
;
function get_spectra, solar_file

  read_netcdf,solar_file,data,attributes
  days=data.date
  n_days=n_elements(days)

  w=data.sp_wave
  f=data.sp_flux
  n_lines=n_elements(w)

  out_spectra=dblarr(n_days+1,n_lines)
  out_spectra[0,*]=double(w)
  for i=1,n_days do begin
     out_spectra[i,*]=double(f[*,i-1])
  endfor
     
  result=size(out_spectra)
  n_columns=result[1]

  ; convert wave unit from nm to Angstrom, solar flux unit
  ; from W M^-2 to photon cm^-2 s^-1 (from w m^-2 to erg cm^-2 s^-1 (1e3),
  ; then from erg cm^-2 s^-1 to photon cm^-2 s^-1)

  out_spectra[0,*]=10*out_spectra[0,*]
  for i=1,n_columns-1 do begin
     out_spectra[i,*]=(6.242e11/12398)*out_spectra[0,*]*out_spectra[i,*] *1e3
  endfor

  ; get bin boundary
  out_spectra=get_boundary(out_spectra)

return, out_spectra
end
