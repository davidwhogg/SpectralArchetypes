;+
; NAME:
;  optimize_test
; PURPOSE:
;  run optimize_bimodality on real data
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU).  All rights reserved.
;-
pro optimize_test
seed= -1L
prefix= 'optimize_bimodality_test'
readcol, '../../data/archetypes/10PCs.dat',x1,x2,x3,x4,x5,x6, $
  format='D,D,D,D'
set_plot, 'ps'
xsize= 7.5 & ysize= 7.5
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
hogg_plot_defaults
psym=6
data= [[x1],[x2],[x3],[x4],[x5],[x6]]
khat1= optimize_bimodality(data,comp=comp1,seed=seed)
orth= khat1
khat2= optimize_bimodality(data,comp=comp2,seed=seed,orth=orth)
orth= [[orth],[khat2]]
khat3= optimize_bimodality(data,comp=comp3,seed=seed,orth=orth)
orth= [[orth],[khat2],[khat3]]
khat4= optimize_bimodality(data,comp=comp4,seed=seed,orth=orth)
splog, transpose(khat1)#khat1
splog, transpose(khat2)#khat2
splog, transpose(khat3)#khat3
splog, transpose(khat4)#khat4
splog, transpose(khat1)#khat2
splog, transpose(khat1)#khat3
splog, transpose(khat1)#khat4
splog, transpose(khat2)#khat3
splog, transpose(khat2)#khat4
splog, transpose(khat3)#khat4
comp3= data # khat3
hogg_plothist, comp1,xtitle='K1'
hogg_plothist, comp2,xtitle='K2'
hogg_plothist, comp3,xtitle='K3'
hogg_plothist, comp4,xtitle='K4'
plot, comp1,comp2,psym=psym, $
  xtitle='K1', $
  ytitle='K2'
plot, comp1,comp3,psym=psym, $
  xtitle='K1', $
  ytitle='K3'
plot, comp1,comp4,psym=psym, $
  xtitle='K1', $
  ytitle='K4'
plot, comp2,comp3,psym=psym, $
  xtitle='K2', $
  ytitle='K3'
plot, comp2,comp4,psym=psym, $
  xtitle='K2', $
  ytitle='K4'
plot, comp3,comp4,psym=psym, $
  xtitle='K3', $
  ytitle='K4'
device, /close
splog, khat1
splog, khat2
splog, khat3
splog, khat4
return
end
