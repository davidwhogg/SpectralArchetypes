function lp_emline_chisq, emline,emivar,model,attenscale
component0= model[*,*,0]
component1= model[*,*,1]
longones= replicate(1.0,n_elements(component0[*,0]))
shortones= replicate(1.0,n_elements(component0[0,*]))
; Think x = inv(at w a).(at w).b
atwb0= (component0*(longones#emivar))#emline
atwb1= (component1*(longones#emivar))#emline
atwa00= total(component0*(longones#emivar)*component0,2)
atwa01= total(component0*(longones#emivar)*component1,2)
atwa11= total(component1*(longones#emivar)*component1,2)
det= atwa00*atwa11-atwa01*atwa01
atwainv00= atwa11/det
atwainv01=-atwa01/det
atwainv11= atwa00/det
amp0= atwainv00*atwb0+atwainv01*atwb1
amp1= atwainv01*atwb0+atwainv11*atwb1
resid= longones#emline $
  -((amp0#shortones)*component0+ $
    (amp1#shortones)*component1)
return, total(resid*(longones#emivar)*resid,2)+(amp1/attenscale)^2
end
