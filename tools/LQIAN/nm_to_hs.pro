;+
; NAME:
;       nm_to_hs
;
; PURPOSE:
;       To transfer spectra on 1 nm binning scheme to Hinteregger reference
;       spectrum resolution for wave length 100A-1000 A and
;       add desired lines in FUV.
;
; CATEGORY:
;       utility function called by the main program.
;
; CALLING SEQUENCE:
;       data_out=nm_to_hs(data_in)
;
; INPUTS:
;       data_in: a two dimension array that holds input
;                spectra on 1 nm bins. It has at least 3 
;                columns (short bounday in A, long boundary in A, 
;                solar flux in photon cm^-2 s^-1, additional
;                columns that are flux (photon cm^-2 s^-1) 
;                in nature if any).
;
; OUTPUTS:
;       data_out: a two dimensional array that holds 
;                 high resolution output spectra. It has same
;                 data structure as data_in.
;
; PROCEDURE:
;       The routine uses Hinteregger solar minimum reference 
;       spectrum, calculates flux ratio of each wave length in
;       the 1 nm it belongs to, then apply the ratio to the input 
;       spectra to get spectra in Hinteregger resoltion (only
;       for wave length 100A-1000A and desired lines in FUV). 
;       This process conserves energy. 
;
; ROUTINES CALLED:
;       read_hs  
;

function nm_to_hs,data_in

wave1=data_in[0,*]
wave2=data_in[1,*]
n_bins=n_elements(wave1)
result=size(data_in)
n_columns=result[1]

; read hintergger reference spectrum

hs_file='./input/sc21refw.dat'
solar_activity='low'
a=read_hs(hs_file,solar_activity,1661)
refwvln=a(0,*)
refflux=a(2,*)
n_hs=n_elements(refwvln)

; the following arrays are used to hold the
; mapped data in the Hinteregger wave space
; Since Hinteregger spectrum does not have
; data for the following 4 nm bins (370-380A,
; 380A-390A,420-430A,440-450A), the array are made 
; bigger than Hinteregger array 

n_cind=n_columns-2
ratio=fltarr(n_hs+10)
flux=fltarr(n_cind,n_hs+10)
waves=fltarr(n_hs+10)
wavel=fltarr(n_hs+10)

; Get ratio for each Hinteregger wave length.
; Since Hinteregger spectrum has poor quality 
; in FUV, here HS is only referenced for wave length
; 50A-1050A, and it is enough for our calculation
; purpose for low resolution scheme, which has 
; logical bins in the wave length range 
; 650-975A

i_start=where( ((wave1 le 50) and (wave2 gt 50)), scount)
i_end= where( ((wave1 le 1050) and (wave2 gt 1050)),ecount)
if scount eq 0 then i_start[0]=0
if ecount eq 0 then i_end[0]=n_bins-1

index=0  ; number of Hinteregger spectrum that a ratio is calculated
         ; plus number of 1 nm bins that has no Hinteregger data
for i=i_start[0], i_end[0] do begin
   ind=where( (refwvln ge wave1[i]) and (refwvln lt wave2[i]))
   if ind[0] ne -1 then begin
      tot_flux=total(refflux(ind))
      if tot_flux eq 0 then tot_flux=0.001
      n_ind=n_elements(ind)
      for k=0,n_ind-1 do begin
	ratio[index]=refflux[ind[k]]/tot_flux
	flux[0:n_cind-1,index]=ratio[index]*data_in[2:n_columns-1,i]
	waves[index]=refwvln[ind[k]]
	wavel[index]=waves[index]
	index=index+1
      endfor
   endif else begin
      waves[index]=data_in[0,i]
      wavel[index]=data_in[1,i]
      flux[0:n_cind-1,index]=data_in[2:n_columns-1,i]
      index=index+1
   endelse
endfor

; the data_in is now mapped to Hinteregger wave length 
; space from i_start to i_end
; The output spectra will consists of its original spectra
; for wave length shortward of i_start and longward of
; i_end, and mapped spectra on Hinteregger space from
; i_start to i_end 

if i_end[0] eq n_bins-1 then begin
   n_out=i_start[0]+index 
   data_out=fltarr(n_columns,n_out)
endif else begin
   n_out=i_start[0]+index+(n_bins-1-i_end[0])
   ;print,i_start[0],index,n_out
   data_out=fltarr(n_columns,n_out)
   data_out[*,i_start[0]+index:n_out-1]=data_in[*,i_end[0]+1:n_bins-1]
endelse

if i_start[0] gt 0 then data_out[*,0:i_start[0]-1]=data_in[*,0:i_start[0]-1]

for i=0,index-1 do begin
   data_out[0,i_start[0]+i]=waves[i]
   data_out[1,i_start[0]+i]=wavel[i]
   data_out[2:n_columns-1,i_start[0]+i]=flux[0:n_cind-1,i]
endfor
;print,data_out[0,i_start[0]+index:n_out-1]
;print,data_out[2,i_start[0]+index:n_out-1]
;print,data_in[0,i_end[0]+1:n_bins-1]
;print,data_in[2,i_end[0]+1:n_bins-1]

; GLOW and low resolution scheme have lines. Lines less than 1000A
; have been taken care in the above step. Lines in FUV need to add 
; to the spectra as well.

; energy conservation check

if n_columns eq 3 then begin
  print,'in nm_to_hs, input:',total(data_in[2,*])
  print,'in nm_to_hs, output:',total(data_out[2,*])
endif else if n_columns gt 3 then begin
  print,'in nm_to_hs, input:',total(data_in[2,*]),total(data_in[3,*])
  print,'in nm_to_hs, output:',total(data_out[2,*]),total(data_out[3,*])
endif

;index=where(data_out[0,*] eq 1215.67)
;print,'in nm_to_hs:',data_out[2,index]
print, '      '
return, data_out
end
