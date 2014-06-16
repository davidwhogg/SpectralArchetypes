;+
; NAME:
;  lp_format
; PURPOSE:
;  write CPLEX LP format file for archetypes problem
; INPUTS:
;  matrix    - integer or binary matrix, with first index going over
;              spectra and second index going over models, and
;              non-zero entries where the model represents the
;              spectrum
;  filename  - for output, default to ./idl.lp
; KEYWORDS:
;  binary    - make a binary programming .lp file instead of the
;              linear programming default
; OUTPUTS:
;  [cplex lp file]
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU) all rights reserved.
;-
pro lp_format, matrix,filename,binary=binary
if (NOT keyword_set(filename)) then filename= './idl.lp'
if (NOT keyword_set(matrix)) then begin
    seed= -1L
    matrix= (randomu(seed,10000,10000) LT 0.001)
    help, matrix
endif
foo= size(matrix,/dimens)
nspectra= foo[0]
nmodel= foo[1]
splog, 'found',nspectra,' spectra and',nmodel,' models'
openw, wlun,filename,/get_lun
printf, wlun,'Minimize'
lp_format_constraint, 0,lindgen(nmodel),wlun=wlun,/cost
printf, wlun,''
printf, wlun,'Subject To'
for cc=0L,nspectra-1L do begin
    row= reform(matrix[cc,*],nmodel)
    if total(row GT 0) then begin
        lp_format_constraint, cc,where(row GT 0),wlun=wlun
    endif
endfor
printf, wlun,''
if keyword_set(binary) then begin
    printf, wlun,'Binary'
    for ii=0L,nmodel-1L do printf, wlun,'  a'+lp_format_index(ii)
endif else begin
    printf, wlun,'Bounds'
    lp_format_bounds, lindgen(nmodel),wlun=wlun
endelse
printf, wlun,''
printf, wlun,'End'
close, wlun
free_lun, wlun
;;cmd= 'gzip -fv --best '+filename
;;splog, cmd
;;spawn, cmd
return
end

