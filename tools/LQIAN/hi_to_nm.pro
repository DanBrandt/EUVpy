;+
; NAME:
;       hi_to_nm
;
; PURPOSE:
;       To map 1nm or higher resolution solar spectra 
;       to 1nm binning scheme. It conserves energy.
;
; CATEGORY:
;       utility function called by the main program.
;
; CALLING SEQUENCE:
;       data_out=hi_to_nm(data_in)
;
; INPUTS:
;       data_in:    a two dimension array that holds input 1 nm
;                   or higher resolution spectra. It has 
;                   at least 3 columns (waves and wavel in Angstrom, solar flux, 
;                   additional fluxes if any).
;                   The wave length needs to be in Angstrom.
; end modification
;
; OUTPUTS:
;       data_out:   a two dimensional array that holds output 
;                   spectra on 1 nm binning scheme. It has 
;                   at least three columns:(wave short boundary in A, 
;                   wave long boundary in A, solar flux, 
;                   additional fluxes if any)
;
;
; ROUTINES CALLED:
;       getbins
;

function hi_to_nm,data_in

; get 1 nm bins:

print, 'in hi_to_nm'
print,[data_in[0,*],data_in[1,*],data_in[50,*]]
bin_file='./input/nm_bins.txt'
getbins,bin_file,wave1,wave2
n_bins=n_elements(wave1)

waves=data_in[0,*]
wavel=data_in[1,*]

; construct output spectra data structure

result=size(data_in)
n_columns=result[1]
n_lines=result[2]
data_out=fltarr(n_columns,n_bins)
data_out[0,*]=wave1
data_out[1,*]=wave2
n_cind=n_columns-1

; map input spectra to 1 nm binning scheme
; At this point, it is required that the input spectra
; has at least three columns: wave short, wave long, solar flux,
; additional columns (solar fluxes)

count_ly=0
for i=0,n_bins-1 do begin
  if (wave1[i] eq 1215.67) then begin  ; special treatment for Ly-1
     index=where( (waves le 1215.67) and (wavel ge 1215.67),count_ly)
     if (count_ly ne 0) then begin
	ly_1=index[0]                  ; sometimes the original 
	up_1=index[0]-1                ; data bin that has ly_alpha 
	down_1=index[0]+1              ; solar flux might does not 
	if (up_1 ge 0) then begin      ; include 1215.67A, so search nearby
	; modified to add if statement because TIMED/SEE version 9 include
	; negative data when data is missing
	   if data_in[2,ly_1] gt 0 then begin           
	      temp=[data_in[2,ly_1],data_in[2,up_1]]
           endif else begin
	      temp=[data_in[20,ly_1],data_in[20,up_1]]
           endelse
	   max_val=max(temp,ii)
	   if ii eq 1 then ly_1=up_1
        endif
	if (down_1 le n_lines-1) then begin
	   if data_in[2,ly_1] gt 0 then begin           
	      temp=[data_in[2,ly_1],data_in[2,down_1]]
           endif else begin
	      temp=[data_in[20,ly_1],data_in[20,up_1]]
           endelse
	   max_val=max(temp,ii)
	   if ii eq 1 then ly_1=down_1
        endif
	index[0]=ly_1
	data_out[2:n_cind,i]=data_in[2:n_cind,index[0]]
     endif
  endif else begin
     ind=where( (wavel gt wave1[i]) and (waves lt wave2[i]) )  
     if (ind[0] eq -1) then begin
       data_out[2:n_cind,i]=0
     endif else if (n_elements(ind) eq 1) then begin
       ; I am using linear interpration if data is coarse 
       ; than the binning scheme of nm_bins.txt.ver1 which is very likely at wavelength <3.2nm
       delw_1=wave2[i]-wave1[i]
       delw_2=wavel[ind[0]]-waves[ind[0]]
       if delw_1 le delw_2 then begin
          data_out[2:n_cind,i]=data_in[2:n_cind,ind[0]] * delw_1/delw_2
       endif else begin
          data_out[2:n_cind,i]=data_in[2:n_cind,ind[0]]
       endelse
     endif else begin
       n_ind=n_elements(ind)
       for j=0,n_ind-1 do begin
          ; again, I am using linear interpolation 
          if( (waves[ind[j]] le wave1[i]) and (wavel[ind[j]] gt wave1[i]) ) then begin
	    data_out[2:n_cind,i]=data_out[2:n_cind,i]+(wavel[ind[j]]-wave1[i])/$
	   	       (wavel[ind[j]]-waves[ind[j]])*data_in[2:n_cind,ind[j]]
          endif else if ( (waves[ind[j]] ge wave1[i]) and (wavel[ind[j]] le wave2[i]) ) then begin
	    data_out[2:n_cind,i]=data_out[2:n_cind,i]+data_in[2:n_cind,ind[j]]
          endif else if ( (waves[ind[j]] lt wave2[i]) and (wavel[ind[j]] ge wave2[i]) ) then begin
	    data_out[2:n_cind,i]=data_out[2:n_cind,i]+(wave2[i]-waves[ind[j]])/$
		       (wavel[ind[j]]-waves[ind[j]])*data_in[2:n_cind,ind[j]]
          endif
       endfor
     endelse
  endelse
  if (wave1[i] eq 1200 and count_ly ne 0) then begin
     data_out[2:n_cind,i] = data_out[2:n_cind,i]-data_in[2:n_cind,index[0]]
  endif
endfor

print, 'in hi_to_nm'
print,[data_out[0,*],data_out[1,*],data_out[50,*]]
; energy conservation check

if n_columns eq 3 then begin
  print,'in hi_to_nm, input:',total(data_in[2,where(wavel le 1750)])
  print,'in hi_to_nm, output:',total(data_out[2,where(wave2 le 1750)])
endif else if n_columns gt 3 then begin
  print,'in hi_to_nm, input:',total(data_in[2,where(wavel le 1750)]), $
			      total(data_in[3,where(wavel le 1750)])
  print,'in hi_to_nm, output:',total(data_out[2,*]), $
			       total(data_out[3,*])
endif
return,data_out
end
