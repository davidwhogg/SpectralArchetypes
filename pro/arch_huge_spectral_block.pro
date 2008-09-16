pro arch_huge_spectral_block

; read in some vetted SDSS catalog
sample= 'dr7'
letter= 'bsafe'
post= '0'
lssfile= '/global/data/sdss/lss/'+sample+'/'+letter+'/'+post $
  +'/post_catalog.'+sample+letter+post+'.fits*'
lsshdr= headfits(lssfile)
lssindx= (mrdfits(lssfile,1)).object_position
specfile= sxpar(lsshdr,'VAGC_RED')+'/object_sdss_spectro.fits*'
spec= mrdfits(specfile,1,rows=temporary(lssindx))
sindx= sort(spec.plate)
spec= spec[sindx]
okay= where(spec.plate GT 0)
spec= spec[okay]
plates= spec[uniq(spec.plate)].plate

; make lambda 
loglam= alog10(3900.0)+1e-4*findgen(3600)
flux= fltarr(n_elements(loglam),n_elements(spec))
invvar= flux

; read in each plate's spectra and resample
for pp=0L,n_elements(plates)-1L do begin
    plate= plates[pp]
    splog, 'working on plate',plate
    pindx= where(spec.plate EQ plate)
    thisspec= spec[pindx]
    readspec, thisspec.plate,thisspec.fiberid,mjd=thisspec.mjd, $
      flux=thisflux,invvar=thisinvvar,loglam=thisloglam
    for jj=0L,n_elements(thisspec)-1L do begin
        jjloglam= thisloglam[*,jj]-alog10(1.0+thisspec[jj].z)
        combine1fiber, jjloglam,thisflux[*,jj],thisinvvar[*,jj], $
          newloglam=loglam,newflux=newflux,newivar=newinvvar
        flux[*,pindx[jj]]= newflux
        invvar[*,pindx[jj]]= newinvvar
    endfor
endfor

; write out huge spectral block
fluxname= 'arch_flux.fits'
mwrfits, flux,fluxname,/create
invvarname= 'arch_invvar.fits'
mwrfits, temporary(invvar),invvarname,/create

; make color block
nlam= n_elements(loglam)
colorvec= fltarr(3,nlam)
colorvec[0,0:nlam/3-1]= 1.0
colorvec[1,nlam/3:2*nlam/3-1]= 1.0
colorvec[2,2*nlam/3:nlam-1]= 1.0
colorname= 'arch_color.fits'
mwrfits, flux##colorvec,colorname,/create
return
end
