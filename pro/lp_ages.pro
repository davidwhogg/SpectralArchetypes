;+
; BUGS:
;  - No proper comment header.
;  - Transpose Wong output, but not sure if that is correct.
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU) all rights reserved.
;-
pro lp_ages
foo= mrdfits('~/astrometry/data/archetypes/ages_*.fits.gz',1)
help, foo,/struct
yesnogrid= (transpose(foo.chisqgrid) LT 3.0)
foo= 0
help, yesnogrid
lp_format, yesnogrid,'ages.lp'
return
end
