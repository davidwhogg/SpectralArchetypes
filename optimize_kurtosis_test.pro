;+
; NAME:
;  optimize_kurtosis_test
; PURPOSE:
;  run optimize_kurtosis on real data
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU).  All rights reserved.
;-
pro optimize_kurtosis_test
seed= -1L
prefix= 'optimize_kurtosis_test'
readcol, '../../data/archetypes/10PCs.dat',x1,x2,x3,x4,x5, $
  format='D,D,D,D'
set_plot, 'ps'
xsize= 7.5 & ysize= 7.5
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
hogg_plot_defaults
psym=6
hogg_plothist, x1, $
  xtitle='PC1'
hogg_plothist, x2, $
  xtitle='PC2'
hogg_plothist, x3, $
  xtitle='PC3'
hogg_plothist, x4, $
  xtitle='PC3'
plot, x1,x2,psym=psym, $
  xtitle='PC1', $
  ytitle='PC2'
plot, x2,x3,psym=psym, $
  xtitle='PC2', $
  ytitle='PC3'
plot, x3,x1,psym=psym, $
  xtitle='PC3', $
  ytitle='PC1'
data= [[x1],[x2],[x3],[x4],[x5]]
khat1= optimize_kurtosis(data,comp=comp1,seed=seed)
orth= khat1
khat2= optimize_kurtosis(data,comp=comp2,seed=seed,orth=orth)
orth= [[orth],[khat2]]
khat3= optimize_kurtosis(data,comp=comp3,seed=seed,orth=orth)
splog, transpose(khat1)#khat1
splog, transpose(khat2)#khat2
splog, transpose(khat3)#khat3
splog, transpose(khat1)#khat2
splog, transpose(khat2)#khat3
splog, transpose(khat3)#khat1
comp3= data # khat3
hogg_plothist, comp1, $
  xtitle='K1'
hogg_plothist, comp2, $
  xtitle='K2'
hogg_plothist, comp3, $
  xtitle='K3'
plot, comp1,comp2,psym=psym, $
  xtitle='K1', $
  ytitle='K2'
plot, comp2,comp3,psym=psym, $
  xtitle='K2', $
  ytitle='K3'
plot, comp3,comp1,psym=psym, $
  xtitle='K3', $
  ytitle='K1'
device, /close
splog, khat1
splog, khat2
splog, khat3
return
end
