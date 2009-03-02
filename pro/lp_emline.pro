;+
; INPUTS:
;  chisqmax  - chi-squared limit (default 8+2*sqrt(2*8))
; BUGS:
;  - No proper comment header.
;  - Brittle reading of emline file.
;  - Brittle reading of glpsol output.
;  - "--mipgap" may be the fractional excess in the MIP objective
;    over the LP objective.
; LICENSE:
;  Copyright 2009 David W. Hogg (NYU) all rights reserved.
;-
pro lp_emline, chisqmax
if (NOT keyword_set(prefix)) then prefix= 'lp_emline'
if (NOT keyword_set(chisqmax)) then chisqmax= 8.0+2.0*sqrt(2.0*8.0)
if (NOT keyword_set(noisefloor)) then noisefloor= 0.05
if (NOT keyword_set(attenscale)) then attenscale= 0.05
if (NOT keyword_set(num_sub)) then num_sub= 1000
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

; choose only good data
halphaindx= nline-1L
halphasn= emline[*,halphaindx]*sqrt(emivar[*,halphaindx])
seed= -1L
dataindx= (where(halphasn LT 30.0))[0:9999]
dataindx= lindgen(10000)
ndata= n_elements(dataindx)
splog, 'trimmed down to',ndata,' spectra'

; choose only good candidate archetypes
candindx= (where((total((emivar LE 0.0),2) LT 0.5) AND $
                (halphasn GT 19.5)))[0:9999]
ncand= n_elements(candindx)
splog, 'trimmed down to',ncand,' candidate archetypes'

; read attenuation and wavelength information
attenu= (mrdfits(filename,2))
wavelength= reform(attenu[*,0],nline)
attenuation= reform(attenu[*,1],nline)

; set up matrices for fitting
component0= emline[candindx,*]
ones= replicate(1.0,ncand)
component1= emline[candindx,*]*(ones#attenuation)
model= [[[component0]],[[component1]]]

; open file
infilename= prefix+'.lp'
openw, wlun,infilename,/get_lun
printf, wlun,'Minimize'
lp_format_constraint, 0,candindx,wlun=wlun,/cost
printf, wlun,''
printf, wlun,'Subject To'

; loop over spectra, finding all candidate archetypes that are okay
;   to represent each spectrum
for ii=0L,ndata-1L do begin
    dd= dataindx[ii]
    thisflux= reform(emline[dd,*],nline)
    thisivar= reform(emivar[dd,*],nline)
    okay= where(lp_emline_chisq(thisflux,thisivar,model,attenscale) $
                LT chisqmax,nokay)
;    splog, dd,nokay
    if (nokay GT 0) then begin
        lp_format_constraint, dd,candindx[okay],wlun=wlun
    endif
endfor

; close file
printf, wlun,''
if keyword_set(binary) then begin
   printf, wlun,'Binary'
   for jj=0L,ncand-1L do printf, wlun,'  a'+lp_format_index(candindx[jj])
endif else begin
   printf, wlun,'Bounds'
   lp_format_bounds, candindx,wlun=wlun
endelse
printf, wlun,''
printf, wlun,'End'
close, wlun
free_lun, wlun

; run glpsol
splog, 'running glpsol'
outfilename= infilename+'.output'
glpsol= '~/astrometry/projects/archetypes/glpk-4.31/examples/glpsol'
cmd= glpsol+' --cpxlp '+infilename+' -o '+outfilename+' --mipgap 0.05'
splog, cmd
spawn, cmd

; trim output from glpsol
splog, 'trimming glpsol output'
grepfilename= outfilename+'.grep'
cmd= 'grep "[0-9].a[0-9][0-9][0-9][0-9][0-9][0-9]" '+outfilename $
  +' > '+grepfilename
splog, cmd
spawn, cmd

; read glpsol output
splog, 'reading trimmed output'
readcol, grepfilename,foo,name,st,activity,format='I,A,A,F'
indx= long(strmid(name,1))
amplitude= fltarr(max(indx)+1)
amplitude[indx]= activity

; write file of indices
use= candindx[where(amplitude GT 0.5,nuse)]
archfilename= prefix+'_archetypes.fits'
splog, 'writing file of',nuse,' archetype indices to '+archfilename
mwrfits, use,prefix+'_archetypes.fits',/create

return
end
