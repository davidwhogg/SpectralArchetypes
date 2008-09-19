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
pro lp_ages, chisqlim,binary=binary
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
lp_format, yesnogrid,infilename,/binary

; run glpsol
outfilename= 'ages.output'
cmd= 'glpsol --cpxlp '+infilename+' -o '+outfilename+' --mipgap 0.01'
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
hogg.ages_id= wong.ages_id[useindx]
hogg.responsibility= round(total(yesnogrid[*,useindx],1))

; check
help, yesnogrid[*,useindx]
help, total(yesnogrid[*,useindx],2)
print, minmax(total(yesnogrid[*,useindx],2))
splog, nuse,' templates for delta-chisq < '+chistr+':'
splog, wong.ages_id[useindx]

; output
fitsfile= 'ages_binary_program.'+chistr+'.fits'
mwrfits, hogg,fitsfile,/create
return
end
