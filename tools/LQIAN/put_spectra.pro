;+
; NAME:
;	put_spectra
;
; PURPOSE:
;       To write outputs
;
; CATEGORY:
;       utility procedure called by main program.
;
; CALLING SEQUENCE:
;       put_spectra,solar_file,out_spectra
;
; INPUTS:
;       solar_file:  (='dir_name/file_name'), original input
;                    spectra file name
;       out_spectra: output spectra 
;
;
; PROCEDURE:
;
; ROUTINES CALLED:
;       None.
;
;

pro put_spectra, solar_file,out_spectra

n_bins=n_elements(out_spectra(0,*))

if ((strpos(solar_file,'.ncdf') ne -1) or (strpos(solar_file,'.nc') ne -1)) then begin 
   read_netcdf,solar_file,data_old
   n_bins=n_elements(out_spectra[0,*])
   n_days=n_elements(out_spectra[*,0])-2
   data_new={DATE:lonarr(n_days),UTTIME:dblarr(n_days),YFRAC:dblarr(n_days),WAVE1:dblarr(n_bins), $
             WAVE2:dblarr(n_bins),SP_FLUX:dblarr(n_bins,n_days)}
   if ( (strpos(solar_file,'L3A') ne -1) or (strpos(solar_file,'FISM') ne -1) ) then begin
;      data_new.UTTIME=double(0.5*(data_old.START_TIME+data_old.STOP_TIME)/(3600.))
      data_new.UTTIME=double(data_old.TIME/3600.)
   endif else begin
      data_new.UTTIME=12.
   endelse
   attr_file='./input/see__L3.attr'
   yfrac=dblarr(n_days)
   for i=0,n_days-1 do begin
      iyear=data_old.DATE[i]/1000
      if ( (iyear-iyear/4*4) eq 0) then begin
         yfrac[i]=double(data_old.DATE[i]/1000)+double(((data_old.DATE[i]-(data_old.DATE[i]/1000)*1000)-1+ $
	                   double(data_new.UTTIME[i]/24.))/366.)
      endif else begin
         yfrac[i]=double(data_old.DATE[i]/1000)+double(((data_old.DATE[i]-(data_old.DATE[i]/1000)*1000)-1+ $
	                   double(data_new.UTTIME[i]/24.))/365.)
      endelse
   endfor
   data_new.YFRAC=yfrac
   data_new.DATE=data_old.DATE
   ;
   temp_array=reverse(out_spectra,2)
   data_new.WAVE1=temp_array[0,*]
   data_new.WAVE2=temp_array[1,*]
   for i=0,n_days-1 do data_new.SP_FLUX[*,i]=temp_array[i+2,*]

   result=strsplit(solar_file,'/',/extract)
   basename=result[n_elements(result)-1]
   result=strsplit(basename,'.',/extract)
   fn=result[0]+'_'+strcompress(n_bins,/remove_all)+'.nc'
   file_name=fn
   write_netcdf,data_new,file_name,att_file=attr_file

endif else begin
   file_name='flux_'+strcompress(n_bins,/remove_all)+'.dat'
   openw,1,file_name
   printf,1,'solar spectra, when read into glow, do flux/10^9'
   printf,1,'wave1(A) ', 'wave2',' fluxes(photon cm^-2 s^-1)'
   result=size(out_spectra)
   n_columns=result[1]
;   fmt='$(2f8.2,12(e11.3,:))'
;   fmt='$(2f8.2,4320(e11.3))'    ; fism 3-day flare spectra
;   fmt='$(2f8.2,28800(e11.3))'    ; fism LWS 2010 20day preflare spectra
;   fmt='$(2f8.2,1500(e11.3))'    ; fism LWS 2010 25 hours control spectra
   fmt='$(2f8.2,50(e11.3))'    ; eve2008, 50 day daily data
   for i=0,n_bins-1 do begin
      printf,1,fmt,out_spectra[0,i],out_spectra[1,i],out_spectra[2:n_columns-1,i]
   endfor
   close,1
endelse

return
end
