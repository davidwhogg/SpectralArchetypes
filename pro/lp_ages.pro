;+
; INPUTS:
;  chisqlim  - delta-chi-squared limit (default 3.0)
; BUGS:
;  - No proper comment header.
;  - Transpose Wong output, but not sure if that is correct.
;  - Brittle reading of glpsol output.
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU) all rights reserved.
;-
pro lp_ages, chisqlim,ioannis=ioannis
if (NOT keyword_set(chisqlim)) then chisqlim= 1.0
chistr= strtrim(string(round(chisqlim)),2)
chisqlim= float(chistr)
if keyword_set(ioannis) then begin
    prefix= 'ages_ioannis.'+chistr
    ages= mrdfits('/mount/moon1/ioannis/home/research/projects/ages/projects/archetype/ages_archetype_chi2grid.fits',1)
; transpose in the next line because John's structure has structure
; elements one for each datum and a row of models in each structure
; element, and when you do john.chi2, you get an array with the models
; in the first index, and we need it the other way.
    yesnogrid= (transpose(ages.chi2) LT chisqlim)
endif else begin
    prefix= 'ages_wong.'+chistr
    ages= mrdfits('~/astrometry/data/archetypes/ages_*.fits.gz',1)
; transpose in the next line because Wong's structure has first index
; looping over models and second over data, and we need the other way.
    yesnogrid= (transpose(ages.chisqgrid) LT chisqlim)
endelse

; make CPLEX LP file
foo= size(yesnogrid,/dimens)
nspectra= foo[0]
nmodel= foo[1]
infilename= prefix+'.lp'
lp_format, yesnogrid,infilename,/binary

; run glpsol
outfilename= infilename+'.output'
cmd= 'glpsol --cpxlp '+infilename+' -o '+outfilename+' --mipgap 0.01'
splog, cmd
spawn, cmd

; trim output from glpsol
grepfilename= outfilename+'.grep'
cmd= 'grep "[0-9].a[0-9][0-9][0-9][0-9][0-9][0-9]" '+outfilename $
  +' > '+grepfilename
splog, cmd
spawn, cmd

; read output
readcol, grepfilename,foo,name,st,activity,format='I,A,A,F'
indx= long(strmid(name,1))
amplitude= fltarr(max(indx)+1)
amplitude[indx]= activity

; sort and pick through templates based on LP amplitude
sindx= reverse(sort(amplitude))
jj= 0L
cover= 0L
use= bytarr(nmodel)
repeat begin
    subgrid= yesnogrid[*,sindx[0:jj]]
    if (jj GT 0) then subgrid= total(subgrid,2)
    oldcover= cover
    cover= round(total(subgrid GT 0))
    dcover= cover-oldcover
    if (dcover GT 0) then use[sindx[jj]]= 1
    jj= jj+1
endrep until ((cover EQ nspectra) or (jj GE nmodel))
useindx= where(use,nuse)

; sort and pick through templates based on responsibility
nonzero= where(amplitude GT 0.0)
resp= round(total(yesnogrid[*,nonzero],1))
sindx= nonzero[reverse(sort(resp))]
jj= 0L
cover= 0L
ruse= bytarr(nmodel)
repeat begin
    subgrid= yesnogrid[*,sindx[0:jj]]
    if (jj GT 0) then subgrid= total(subgrid,2)
    oldcover= cover
    cover= round(total(subgrid GT 0))
    dcover= cover-oldcover
    if (dcover GT 0) then ruse[sindx[jj]]= 1
    jj= jj+1
endrep until ((cover EQ nspectra) or (jj GE nmodel))
ruseindx= where(ruse,nruse)

; choose and save shorter list
splog, 'nuse:',nuse,' nruse:',nruse
if (nruse LT nuse) then begin
    splog, 'responsibility sorting wins'
    use= ruse
endif else begin
    splog, 'amplitude sorting wins'
endelse
useindx= where(use,nuse)
hogg= {ages_id: 0L, responsibility: 0L}
hogg= replicate(hogg,nuse)
hogg.ages_id= (reform((ages.ages_id),nmodel))[useindx]
hogg.responsibility= round(total(yesnogrid[*,useindx],1))

; check
help, yesnogrid[*,useindx]
help, total(yesnogrid[*,useindx],2)
print, minmax(total(yesnogrid[*,useindx],2))
splog, nuse,' templates for delta-chisq < '+chistr+':'
splog, hogg.ages_id

; output
fitsfile= prefix+'.hogg.fits'
mwrfits, hogg,fitsfile,/create
hdr= headfits(fitsfile)
sxaddpar, hdr,'DCHI2LIM',chisqlim
djs_modfits, fitsfile,0,hdr
return
end
