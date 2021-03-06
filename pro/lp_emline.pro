;+
; INPUTS:
; KEYWORDS:
;  fileonly - if set, make the CPLEX LP file but don't run it.
; BUGS:
;  - No proper comment header.
;  - Brittle reading of emline file.
;  - Brittle reading of glpsol output.
;  - "--mipgap" may be the fractional excess in the MIP objective
;    over the LP objective.
; LICENSE:
;  Copyright 2009 David W. Hogg (NYU) all rights reserved.
;-
pro lp_emline, chisqmax,fileonly=fileonly
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

; put in bad data for unmeasured lines
;   (this permits spectra with unmeasured lines to represent themselves
;   only)
badvalue= -10.0*max(emline)
emline[where(emivar LE 0.0)]= badvalue

; choose subset to use as candidate archetypes
halphaindx= nline-1L
sindx= reverse(sort(emline[*,halphaindx]))
candindx= sindx[0:49999]
ncand= n_elements(candindx)
splog, 'trimmed down to',ncand,' candidate archetypes'

; choose subset to use as data
seed= -1L
rindx= shuffle_indx(n_elements(sindx),seed=seed)
dataindx= rindx[0:49999]
ndata= n_elements(dataindx)
splog, 'trimmed down to',ndata,' spectra'

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
splog, 'writing cplex file '+infilename
spawn, 'date --iso=s'
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
    if ((ii MOD 100) EQ 0) then splog, ii,dd,nokay
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

if keyword_set(fileonly) then begin
    spawn, 'date --iso=s'
    splog, '/FILEONLY : compressing and then stopping'
    cmd= 'gzip -v --best '+infilename
    splog, cmd
    spawn, cmd
    stop
endif

; run glpsol
splog, 'running glpsol'
spawn, 'date --iso=s'
outfilename= infilename+'.output'
glpsol= '~/astrometry/projects/archetypes/glpk-4.31/examples/glpsol'
cmd= glpsol+' --cpxlp '+infilename+' -o '+outfilename+' --mipgap 0.05'
splog, cmd
spawn, cmd
spawn, 'date --iso=s'

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
