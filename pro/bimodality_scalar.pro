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
