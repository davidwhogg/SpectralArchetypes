;+
; INPUTS:
;  chisqmax  - chi-squared limit (default 8+2*sqrt(2*8))
; BUGS:
;  - No proper comment header.
;  - Brittle reading of emline file.
;  - Brittle reading of glpsol output.
;  - What does "--mipgap" do?
; LICENSE:
;  Copyright 2009 David W. Hogg (NYU) all rights reserved.
;-
pro lp_emline, chisqmax
if (NOT keyword_set(prefix)) then prefix= 'lp_emline'
if (NOT keyword_set(chisqmax)) then chisqmax= 8.0+2.0*sqrt(2.0*8.0)
if (NOT keyword_set(noisefloor)) then noisefloor= 0.05
if (NOT keyword_set(attenscale)) then attenscale= 0.05
binary= 1
chistr= strtrim(string(chisqmax,format='(F5.2)'),2)
chisqmax= float(chistr)

; read data
indir= '/mount/moon1/ioannis/home/research/projects/sdss/projects/archetype'
filename= indir+'/sdss_emline_matrix.fits.gz'
lindx= [2,3,4,5,6,7,8,9,10]
nline= n_elements(lindx)
emline= (mrdfits(filename))[*,lindx]
emlerr= (mrdfits(filename,1))[*,lindx]

; make invvar and establish noise floor
emivar= 1.0/(emlerr^2+noisefloor^2*emline^2)
emivar[where(emlerr LE 0.0)]= 0.0
emlerr= 0

; trim data to manageable size
halphaindx= nline-1L
halphasn= emline[*,halphaindx]*sqrt(emivar[*,halphaindx])
dindx= (where(halphasn GT 19.5))[0:2999]
emline= emline[dindx,*]
emivar= emivar[dindx,*]
ndata= n_elements(emline[*,halphaindx])
splog, 'got',ndata,' spectra'

; read attenuation and wavelength information
attenu= (mrdfits(filename,2))
wavelength= reform(attenu[*,0],nline)
attenuation= reform(attenu[*,1],nline)

; choose only good candidate archetypes
badindx= where(total((emivar LE 0.0),2) GT 0.5,nbad)
splog, 'threw away',nbad,' candidate archetypes'

; set up matrices for fitting
component0= emline
ones= replicate(1.0,ndata)
component1= emline*(ones#attenuation)

; open file
infilename= prefix+'.lp'
openw, wlun,infilename,/get_lun
printf, wlun,'Minimize'
lp_format_constraint, 0,lindgen(ndata),wlun=wlun,/cost
printf, wlun,''
printf, wlun,'Subject To'

; loop over spectra, finding all candidate archetypes that are okay
;   to represent each spectrum
for dd=0L,ndata-1L do begin

; do all matrix stuff by hand
;   Think x = inv(at w a).(at w).b
   atwb0= (component0*(ones#emivar[dd,*]))#transpose(emline[dd,*])
   atwb1= (component1*(ones#emivar[dd,*]))#transpose(emline[dd,*])
   atwa00= total(component0*(ones#emivar[dd,*])*component0,2)
   atwa01= total(component0*(ones#emivar[dd,*])*component1,2)
   atwa11= total(component1*(ones#emivar[dd,*])*component1,2)
   det= atwa00*atwa11-atwa01*atwa01
   atwainv00= atwa11/det
   atwainv01=-atwa01/det
   atwainv11= atwa00/det
   amp0= atwainv00*atwb0+atwainv01*atwb1
   amp1= atwainv01*atwb0+atwainv11*atwb1
   resid= emline-((amp0#replicate(1.0,nline))*component0+ $
                  (amp1#replicate(1.0,nline))*component1)
   chisq= total(resid*emivar*resid,2)+(amp1/attenscale)^2

; write line if possible
   chisq[badindx]= 2.0*chisqmax
   okay= where(chisq LT chisqmax,nokay)
;   splog, 'spectrum',dd,' nokay',nokay
   if (nokay GT 0) then begin
      lp_format_constraint, dd,okay,wlun=wlun,cost=cost
   endif
endfor

; close file
printf, wlun,''
if keyword_set(binary) then begin
   printf, wlun,'Binary'
   for dd=0L,ndata-1L do printf, wlun,'  a'+lp_format_index(dd)
endif else begin
   printf, wlun,'Bounds'
   lp_format_bounds, lindgen(ndata),wlun=wlun
endelse
printf, wlun,''
printf, wlun,'End'
close, wlun
free_lun, wlun

; run glpsol
outfilename= infilename+'.output'
glpsol= '~/astrometry/projects/archetypes/glpk-4.31/examples/glpsol'
cmd= glpsol+' --cpxlp '+infilename+' -o '+outfilename+' --mipgap 20'
splog, cmd
spawn, cmd

; trim output from glpsol
grepfilename= outfilename+'.grep'
cmd= 'grep "[0-9].a[0-9][0-9][0-9][0-9][0-9][0-9]" '+outfilename $
  +' > '+grepfilename
splog, cmd
spawn, cmd

; read glpsol output
readcol, grepfilename,foo,name,st,activity,format='I,A,A,F'
indx= long(strmid(name,1))
amplitude= fltarr(max(indx)+1)
amplitude[indx]= activity

; write file of indices
use= where(amplitude GT 0.5)
mwrfits, use,prefix+'_archetypes.fits',/create

return
end
