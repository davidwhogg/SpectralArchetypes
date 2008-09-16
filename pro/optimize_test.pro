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
readcol, '../../data/archetypes/10PCs.dat', $
  x1,x2,x3,x4,x5,x6,x7,x8,x9, $
  format='D,D,D,D,D,D,D,D,D'
set_plot, 'ps'
xsize= 7.5 & ysize= 7.5
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
hogg_plot_defaults
psym=6
data= [[x1],[x2],[x3],[x4],[x5],[x6],[x7],[x8],[x9]]
khat1= optimize_bimodality(data,comp=comp1,seed=seed)
foo= bimodality_scalar(comp1,kindx=kindx)
data0= data[where(kindx EQ 0),*]
data1= data[where(kindx EQ 1),*]
orth= khat1
khat2= optimize_bimodality(data0,comp=comp2,seed=seed,orth=orth)
khat3= optimize_bimodality(data1,comp=comp3,seed=seed,orth=orth)
splog, transpose(khat1)#khat1
splog, transpose(khat2)#khat2
splog, transpose(khat3)#khat3
; splog, transpose(khat4)#khat4
splog, transpose(khat1)#khat2
splog, transpose(khat1)#khat3
; splog, transpose(khat1)#khat4
splog, transpose(khat2)#khat3
; splog, transpose(khat2)#khat4
; splog, transpose(khat3)#khat4
hogg_plothist, data#khat1,xtitle='K1',title='all'
hogg_plothist, data0#khat2,xtitle='K2',title='red'
hogg_plothist, data#khat2,xtitle='K2',title='all'
hogg_plothist, data1#khat3,xtitle='K3',title='blue'
hogg_plothist, data#khat3,xtitle='K3',title='all'
; hogg_plothist, data#khat4,xtitle='K4'
plot, data0#khat1,data0#khat2,psym=psym, $
  xtitle='K1',ytitle='K2',title='red'
plot, data#khat1,data#khat2,psym=psym, $
  xtitle='K1',ytitle='K2',title='all'
plot, data1#khat1,data1#khat3,psym=psym, $
  xtitle='K1',ytitle='K3',title='blue'
plot, data#khat1,data#khat3,psym=psym, $
  xtitle='K1',ytitle='K3',title='all'
device, /close
splog, khat1
splog, khat2
splog, khat3
return
end
