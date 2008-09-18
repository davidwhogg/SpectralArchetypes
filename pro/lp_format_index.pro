;+
; BUGS:
;  - No proper comment header.
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU) all rights reserved.
;-
function lp_format_index, index
return, string(index,format='(I07)')
end
