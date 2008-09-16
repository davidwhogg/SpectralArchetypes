;+
; NAME:
;  optimize_bimodality
; PURPOSE:
;  search for the direction of maximum bimodality in a set of points.
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
function bimodality_skew, data
td= data-mean(data)
var= mean(td^2)
return, mean(td^3)/var^(1.5)
end

function optimize_bimodality, data,comp=comp,seed=seed, $
                              orth=orth
foo= size(data,/dimens)
npoint= foo[0]
ndimen= foo[1]
north= n_elements(orth)/ndimen
if (north GT 0) then orth= reform(orth,ndimen,north)
ntrial= 1L*3L^long(ndimen-north)
khat= dblarr(ndimen)
for ii=0,5 do begin
    amp= 3.0^(1.0-float(ii))
    splog, 'starting',ntrial,' trials with amp=',amp
    for tt=0L,ntrial-1L do begin
        thiskhat= khat+amp*randomn(seed,ndimen)
        for oo=0L,north-1L do begin
            thisorth= reform(orth[*,oo],ndimen)
            thiskhat= thiskhat-thisorth*((transpose(thisorth)#thiskhat)[0])
        endfor
        norm= sqrt((transpose(thiskhat)#thiskhat)[0])
        thiskhat= thiskhat / norm
        thiscomp= data#thiskhat
        scalar= bimodality_scalar(thiscomp)
        if (not keyword_set(bestscalar)) then bestscalar= scalar+1.0
        if (scalar LT bestscalar) then begin
            bestscalar= scalar
            khat= thiskhat
            comp= thiscomp
            splog, tt, khat, bestscalar
        endif
    endfor
endfor
; reverse on skew if necessary
if (bimodality_skew(comp) LT 0.0) then begin
    khat= -1.0*khat
    comp= -1.0*comp
endif
return, khat
end
