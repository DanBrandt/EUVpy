;+
; NAME:
;       getbins
;
; PURPOSE:
;       to get the solar flux wavelength bins. 
;
; CATEGORY:
;       utility procedure.
;
; CALLING SEQUENCE:
;       getbins,binfile,wave1,wave2
;
; INPUTS:
;       binfile: binning scheme file name.
;
; OUTPUTS:
;       wave1: short boundary of the bin in Angstrom.
;       wave2: long boundary of the bin in Angstrom.
;
; COMMON BLOCKS:
;       None.
;
; PROCEDURE:
;
; ROUTINES CALLED:
;       None.
;
; MODIFICATION HISTORY:
;       8/94, SMB
;
;+

pro getbins,binfile,wave1,wave2

openr,1,binfile
s=' '
readf,1,s
readf,1,nbins,scale

wave1=fltarr(nbins)
wave2=fltarr(nbins)

for i=0,nbins-1 do begin
	readf,1,b,c
	wave1(i)=b
	wave2(i)=c
endfor
close,1

wave1=wave1*scale
wave2=wave2*scale

return
end
