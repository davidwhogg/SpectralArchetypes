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

function bimodality_scalar, data,kindx=kindx
td= data-mean(data)
var= mean(td^2)
kmean= dblarr(2)
kmean[0]= mean(td[where(td LE 0.0)])
kmean[1]= mean(td[where(td GE 0.0)])
nzero= 0
repeat begin
    dd2= [[(td-kmean[0])^2],[(td-kmean[1])^2]]
    kindx= fix(dd2[*,1] LT dd2[*,0])
    kvar= mean(dd2[*,0] < dd2[*,1])
    oldkmean= kmean
    oldnzero= nzero
    kmean[0]= mean(td[where(kindx EQ 0,nzero)])
    kmean[1]= mean(td[where(kindx EQ 1,none)])
endrep until (total(abs(oldkmean-kmean)) EQ 0.0)
return, (kvar/var)
end

function optimize_bimodality, data,comp=comp,seed=seed, $
                              orth=orth
foo= size(data,/dimens)
npoint= foo[0]
ndimen= foo[1]
north= n_elements(orth)/ndimen
if (north GT 0) then orth= reform(orth,ndimen,north)
for ii=0,1 do begin
    if (ii EQ 0) then begin
        ntrial= 1000L*3L^long(ndimen-north)
        amp= 1.0
        desc=''
    endif else begin
        ntrial= 10L*3L^long(ndimen-north)
        amp= 0.01
        desc=' refinement'
    endelse
    splog, 'starting',ntrial,desc+' trials'
    for tt=0L,ntrial-1L do begin
        if (ii EQ 0) then begin
            thiskhat= amp*randomn(seed,ndimen)
        endif else begin
            thiskhat= khat+amp*randomn(seed,ndimen)
        endelse
        for oo=0L,north-1L do begin
            thisorth= reform(orth[*,oo],ndimen)
            thiskhat= thiskhat-thisorth*((transpose(thisorth)#thiskhat)[0])
        endfor
        norm= sqrt((transpose(thiskhat)#thiskhat)[0])
        thiskhat= thiskhat / norm
        thiscomp= data#thiskhat
        scalar= bimodality_scalar(thiscomp)
        if (tt EQ 0) then bestscalar= scalar+1.0
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
