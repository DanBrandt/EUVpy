;+
; NAME:
;       rebin
;
; PURPOSE:
;       To put input solar spectra to a selected binning scheme
;
; CATEGORY:
;       core procedure called by the main program.
;
; CALLING SEQUENCE:
;       rebin,bin_file,in_spectra,out_spectra,sigma_o,sigma_o2,sigma_n,sigma_n2
;
; INPUTS:
;       bin_file:     the selected binning scheme file.
;       in_spectra:   a two dimensional array that holds the input
;                     spectra. It has at least 3 columns:
;                     wave short boundary in A, long boundary in A, 
;                     solar flux in photon cm^-2 s^-1, more fluxes if any.
;                     For example, for Hinteregger solar minimum
;                     spectrum, in_spectra has 5 columns:
;                     wave short, wave long, solar flux, flux from
;                     chromsphere, flux from corona. For Woods 
;                     reference spectrum, it also has 5 columns:
;                     wave short, wave long, solar flux, 
;                     27-day variation, 11-year variation.
;
; OUTPUTS:
;       out_spectra:  a two dimensional array that holds the output
;                     spectra. It has same data structure as in_spectra.
;
; ROUTINES CALLED:
;       getbins, read_gen
;

pro rebin, bin_file,in_spectra,out_spectra

; get binning scheme

getbins, bin_file, wave1,wave2

n_bins=n_elements(wave1)

iflines=(wave1 eq wave2)  ; 1 where bins are lines, 0 for continuum bins

; unpack input spectra

refwvln=0.5*(in_spectra[0,*]+in_spectra[1,*])
n_rows=n_elements(refwvln)
result=size(in_spectra)
n_columns=result[1]
n_cind=n_columns-2
refflux=fltarr(n_cind,n_rows)
refflux[0:n_cind-1,*]=in_spectra[2:n_columns-1,*]

; check lines in the input spectra

lines=wave1*iflines
ifreflines=fltarr(n_rows)
for i=0,n_rows-1 do begin
  for j=0,n_bins-1 do begin
     if ((lines(j) eq refwvln(i)) and (iflines(j) eq 1)) then $
     ifreflines(i)=1
  endfor
endfor

; read cross section data. Logical bins are divided based on N2 absorption 
; cross sections

fenn=read_gen('./input/phfenn.tab' ,1,10,1945)
henken=read_gen('./input/henken.tab',1,3,342)

; Combine Fennely, Henke

ind_h2=where(henken(0,*) le 50)
ind_f2=where(fenn(0,*) gt 50)
ln2=[transpose(henken(0,ind_h2)),transpose(fenn(0,ind_f2))]
an2_ext=[transpose(henken(2,ind_h2)),transpose(fenn(1,ind_f2))]

fennsigan2=interpol( an2_ext,ln2,refwvln)

; Get solar flux into model bins.

modflux=fltarr(n_cind,n_bins)

flag_low=(n_bins lt 45)
count_1=0
count_2=0
count_3=0

for ibin=0,n_bins-1 do begin

  ; the following if statement should handle lines or continua, it does
  ; assume though that the line is one of the Hinteregger lines

  if iflines(ibin) eq 1 then begin

    ind=where((refwvln ge wave1(ibin)) and (refwvln le wave2(ibin)))

  endif else begin

    if flag_low eq 1 then begin

       case wave1[ibin] of

	 650.0: begin
	   case count_1 of
	     0: begin
	       ind=where( (refwvln ge wave1(ibin)) and $
		     (refwvln lt wave2(ibin)) and $
		     (fennsigan2 lt 31) )
	       count_1=count_1+1
             end

             1: ind=where( (refwvln ge wave1(ibin)) and $
                     (refwvln lt wave2(ibin)) and $
                     (fennsigan2 ge 31) )
           endcase
	 end

	 798.0: begin
	   case count_2 of
	     0: begin
	       ind=where( (refwvln ge wave1(ibin)) and $
		     (refwvln lt wave2(ibin)) and $
		     (fennsigan2 lt 4) ) 
               count_2=count_2+1
             end

             1: begin
	       ind=where( (refwvln ge wave1(ibin)) and $
                     (refwvln lt wave2(ibin)) and $
                     (fennsigan2 ge 4)  and (fennsigan2 lt 31) )
               count_2=count_2+1
             end

             2: ind=where( (refwvln ge wave1(ibin)) and $
                     (refwvln lt wave2(ibin)) and $
                     (fennsigan2 ge 31) )
           endcase
	 end

	 913.0: begin
	   case count_3 of
	     0: begin
	       ind=where( (refwvln ge wave1(ibin)) and $
		     (refwvln lt wave2(ibin)) and $
		     (fennsigan2 lt 4) ) 
               count_3=count_3+1
             end

             1: begin
	       ind=where( (refwvln ge wave1(ibin)) and $
                     (refwvln lt wave2(ibin)) and $
                     (fennsigan2 ge 4)  and (fennsigan2 lt 31) )
               count_3=count_3+1
             end
	     2: begin
	       ind=where( (refwvln ge wave1(ibin)) and $
		     (refwvln lt wave2(ibin)) and $
		     (fennsigan2 ge 31) ) 
             end

           endcase
	 end

	 else: begin
           ind=where((refwvln ge wave1(ibin)) and (refwvln lt wave2(ibin)))
	 end

       endcase

    endif else ind=where((refwvln ge wave1(ibin)) and (refwvln lt wave2(ibin)))

    if (ind(0) ne -1) then begin
      indcont=(where(ifreflines(ind) eq 0))
      oldind=ind
      if indcont(0) ne -1 then begin 
	indcont=ind(indcont)
	ind=indcont
      endif
    endif
  endelse
  ;
  if (ind(0) ne -1) then begin
     for cind=0,n_cind-1 do begin
        modflux(cind,ibin)=total(refflux(cind,ind))
     endfor
  endif
endfor

; output 

out_spectra=fltarr(n_columns,n_bins)
out_spectra[0,*]=wave1
out_spectra[1,*]=wave2
out_spectra[2:n_columns-1,*]=modflux[0:n_cind-1,*]

; energy conservation check

ind1=where(in_spectra[0,*] lt wave2[n_bins-1])
ind2=where(refwvln lt wave2[n_bins-1])

if n_columns eq 3 then begin
   print,'in rebin, input:',total(in_spectra[2,ind1])
   print,'in rebin, output:',total(out_spectra[2,*])
endif else if n_columns gt 3 then begin
   print,'in rebin, input:',total(in_spectra[2,ind1]),total(in_spectra[3,ind1])
   print,'in rebin, output:',total(out_spectra[2,*]),total(out_spectra[3,*])
endif

return
end
