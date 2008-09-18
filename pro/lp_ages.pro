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
pro lp_ages, chisqlim
if (NOT keyword_set(chisqlim)) then chisqlim= 3.0
chistr= string(chisqlim,format='(F3.1)')
chisqlim= float(chistr)
wong= mrdfits('~/astrometry/data/archetypes/ages_*.fits.gz',1)

; make CPLEX LP file
yesnogrid= (transpose(wong.chisqgrid) LT chisqlim)
foo= size(yesnogrid,/dimens)
nspectra= foo[0]
nmodel= foo[1]
infilename= 'ages.lp'
lp_format, yesnogrid,infilename

; run glpsol
outfilename= 'ages.output'
cmd= 'glpsol --cpxlp '+infilename+' -o '+outfilename
splog, cmd
spawn, cmd

; trim output from glpsol
grepfilename= 'ages.output.grep'
cmd= 'grep "[0-9].a[0-9][0-9][0-9][0-9][0-9][0-9]" '+outfilename $
  +' > '+grepfilename
splog, cmd
spawn, cmd

; read output
readcol, grepfilename,foo,name,st,activity,format='I,A,A,F'
indx= long(strmid(name,1))
amplitude= fltarr(max(indx)+1)
amplitude[indx]= activity

; sort and pick through templates
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

; check
help, yesnogrid[*,useindx]
help, total(yesnogrid[*,useindx],2)
print, minmax(total(yesnogrid[*,useindx],2))
splog, nuse,' templates for delta-chisq < '+chistr+':'
splog, wong.ages_id[useindx]

; output
fitsfile= 'agesId.'+chistr+'.fits'
mwrfits, wong.ages_id[useindx],fitsfile,/create
return
end
