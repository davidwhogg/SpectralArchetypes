pro arch_huge_spectral_block.pro

; read in some vetted SDSS catalog
sample= 'dr7'
letter= 'bsafe'
post= '0'
lssfile= '/global/data/sdss/lss/'+sample+'/'+letter+'/'+post $
  +'/post_catalog.'+sample+letter+post_'.fits*'
lsshdr= headfits(lssfile)
lssindx= (mrdfits(lssfile,1)).object_position
specfile= sxpar(lsshdr,'VAGC_RED')+'object_sdss_spectro.fits'
spec= mrdfits(specfile,1,rows=lssindx)
help
stop

; set up standard wavelength grid

; read in each spectrum and resample

; write out huge spectral block

return
end
