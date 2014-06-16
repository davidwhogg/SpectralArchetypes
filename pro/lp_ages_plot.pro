;+
; BUGS:
;  - No proper comment header.
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU) all rights reserved.
;-
pro lp_ages_plot
prefix= 'ages_ntemplate_dchi2lim'
file= findfile('./ages_ioannis.*.hogg.fits',count=nfile)
dchi2lim= fltarr(nfile)
ntemplate= lonarr(nfile)
for ff=0,nfile-1 do begin
    hdr= headfits(file[ff])
    dchi2lim[ff]= sxpar(hdr,'DCHI2LIM')
    ntemplate[ff]= sxpar(hdr,'NTEMPLAT')
endfor
set_plot, 'ps'
xsize= 7.5 & ysize= 7.5
device, file=prefix+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
hogg_plot_defaults
hogg_usersym, 16,/fill
plot, alog10(dchi2lim),ntemplate,psym=8, $
  xrange=[-0.2,5.0],xtitle='log!d10!n( delta-chi-squared )', $
  yrange=[0,1461],ytitle='number of archetypes required', $
  title= 'Moustakas AGES sample'
xyouts, alog10(dchi2lim),ntemplate,' '+strtrim(string(ntemplate),2), $
  orientation=45.0
device,/close
return
end
