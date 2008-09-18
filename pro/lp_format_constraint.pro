pro lp_format_constraint, cnumber,aindex,wlun=wlun,cost=cost
if NOT keyword_set(wlun) then wlun= -1
term= ' + 1 a'+lp_format_index(aindex)
nterm= n_elements(term)
if keyword_set(cost) then begin
    term[0]= ' cost:    '+term[0]
endif else begin
  term[0]= ' c'+lp_format_index(cnumber)+':'+term[0]
  term[nterm-1L]= term[nterm-1L]+' >= 1'
endelse
linestart= ''
for jj=0L,nterm-1L,5L do begin
    printf, wlun,linestart+strjoin(term[jj:((jj+4L)<(nterm-1L))],'')
    linestart= '          '
endfor
return
end
