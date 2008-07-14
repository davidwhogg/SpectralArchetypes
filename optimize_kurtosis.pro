;+
; NAME:
;  optimize_kurtosis
; PURPOSE:
;  search for the direction of maximum kurtosis in a set of points
; COMMENTS:
;  - Doesn't assume zero mean or unit variance; probably should
;    for speed.
; INPUTS:
;  data    - [npoint,ndimen] array of data
; OPTIONAL INPUTS:
;  seed    - for random number generator
;  orth    - make khat orthogonal to this vector or set of vectors
; OUTPUT:
;  khat    - optimal direction
; OPTIONAL OUTPUTS:
;  comp    - data component in khat direction
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU).  All rights reserved.
;-
function optimize_kurtosis, data,comp=comp,seed=seed, $
                            orth=orth
foo= size(data,/dimens)
npoint= foo[0]
ndimen= foo[1]
north= n_elements(orth)/ndimen
if (north GT 0) then orth= reform(orth,ndimen,north)
ntrial= 1000*3^(ndimen-north)
splog, 'starting',ntrial,' trials'
for tt=0L,ntrial-1L do begin
    thiskhat= randomn(seed,3)
    for oo=0L,north-1L do begin
        thisorth= reform(orth[*,oo],ndimen)
        thiskhat= thiskhat-thisorth*((transpose(thisorth)#thiskhat)[0])
    endfor
    norm= sqrt((transpose(thiskhat)#thiskhat)[0])
    thiskhat= thiskhat / norm
    thiscomp= data#thiskhat
    tmp= thiscomp-mean(thiscomp)
    tmp= tmp^2/mean(tmp^2)
    kurtosis= mean(tmp^4)
    if (tt EQ 0) then bestkurtosis= kurtosis+1.0
    if (kurtosis LT bestkurtosis) then begin
        bestkurtosis= kurtosis
        khat= thiskhat
        comp= thiscomp
        splog, tt, khat, bestkurtosis
    endif
endfor
help, comp
return, khat
end
