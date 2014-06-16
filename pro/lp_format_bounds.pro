;+
; BUGS:
;  - No proper comment header.
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU) all rights reserved.
;-
pro lp_format_bounds, aindex,wlun=wlun
if NOT keyword_set(wlun) then wlun= -1
term= ' 0 <= a'+lp_format_index(aindex)+' <= 1'
nterm= n_elements(term)
for jj=0L,nterm-1L,1L do begin
    printf, wlun,term[jj]
endfor
return
end
