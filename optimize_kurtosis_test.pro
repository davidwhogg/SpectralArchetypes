;+
; NAME:
;  optimize_kurtosis_test
; PURPOSE:
;  run optimize_kurtosis on real data
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU).  All rights reserved.
;-
pro optimize_kurtosis_test
prefix= 'optimize_kurtosis_test'
read_ascii, '../../data/archetypes/3PC.dat',x1,x2,x3, $
            format='D,D,D'
set_plot, 'ps'
xsize= 7.5 & ysize= 7.5
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
hogg_plot_defaults
psym=6
hogg_plothist, x1, $
               xlabel='x_1'
hogg_plothist, x2, $
               xlabel='x_2'
hogg_plothist, x3, $
               xlabel='x_3'
plot, x1,x2,psym=psym, $
      xlabel='x_1', $
      ylabel='x_2'
plot, x2,x3,psym=psym, $
      xlabel='x_2', $
      ylabel='x_3'
plot, x3,x1,psym=psym, $
      xlabel='x_3', $
      ylabel='x_1'
data= [[x1],[x2],[x3]]
device, /close
return
end
