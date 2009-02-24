;+
; INPUTS:
;  chisqmax  - chi-squared limit (default 10.0)
; BUGS:
;  - No proper comment header.
;  - Brittle reading of glpsol output.
; LICENSE:
;  Copyright 2009 David W. Hogg (NYU) all rights reserved.
;-
pro lp_emline, chisqmax
if (NOT keyword_set(chisqmax)) then chisqmax= 10.0
if (NOT keyword_set(noisefloor)) then noisefloor= 0.05
chistr= strtrim(string(chisqmax,format='(F5.2)'),2)
chisqmax= float(chistr)

; read data, make invvars, and establish noise floor
indir= '/mount/moon1/ioannis/home/research/projects/sdss/projects/archetype'
filename= indir+'/sdss_emline_matrix.fits.gz'
lindx= [2,3,4,5,6,7,8,9,10]
emline= (mrdfits(filename))[*,lindx]
emlerr= (mrdfits(filename,1))[*,lindx]
attenu= (mrdfits(filename,2))[lindx]
emivar= 1.0/(emlerr^2+noisefloor^2*emline^2)
emivar[where(emlerr LE 0.0)]= 0.0
emlerr= 0

; sort data by h-alpha significance and start with the best
halphaindx= 0
sindx= reverse(sort((emline[*,halphaindx])^2*emivar[*,halphaindx]))
ndata= n_elements(sindx)
splog, 'got',ndata,' spectra'
for dd=0L,ndata-1L do begin
   ii= sindx[dd]
   bad= where(emivar[ii,*] LE 0.0,nbad)
   if (nbad GT 0) then begin
      splog, dd,ii,': at least one bad measurement, skipping as archetype'
   endif else begin
      splog, dd,ii

; set up and solve least-squares problems
   
      stop
   endelse
endfor

; make CPLEX LP file
;; foo= size(yesnogrid,/dimens)
;; nspectra= foo[0]
;; nmodel= foo[1]
;; infilename= prefix+'.lp'
;; lp_format, yesnogrid,infilename,/binary

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
sxaddpar, hdr,'NTEMPLAT',nuse
djs_modfits, fitsfile,0,hdr
return
end
