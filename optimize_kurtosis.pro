;+
; NAME:
;  optimize_kurtosis
; PURPOSE:
;  search for the direction of maximum kurtosis in a set of points
; COMMENTS:
;  - Doesn't assume zero mean or unit variance; probably should
;    for speed.
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU).  All rights reserved.
;-
pro optimize_kurtosis, data
return khat
end
