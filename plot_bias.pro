function load_colors

  restore,filename="./data/color_tables/color_table_fig_bias.dat"
  return,[[r],[g],[b]]
end


function calculate_bias, d0d1Ls_obs=d0d1Ls_obs, d0d1Ls_SR=d0d1Ls_SR, Nr=Nr

   ;region/years/realization/shipclass/variable(d0,d1,Ls)/observational ice thickness
  bias =fltarr(4,40,Nr,4,3,3)
  for ireg=0,3 do begin
     for sk=0,3 do begin
        for stat=0,2 do begin         
           ;denormalize
           d0o=rebin(reform(( d0d1Ls_obs[ireg,11:50,sk,stat,0]-(findgen(40)+1979))*365.-91.),40,Nr)
           d1o=rebin(reform(( d0d1Ls_obs[ireg,11:50,sk,stat,1]-(findgen(40)+1979))*365.-91.),40,Nr)
           Lso=rebin(reform(( d0d1Ls_obs[ireg,11:50,sk,stat,2])*365.),40,40)

           d0m=reform(d0d1Ls_SR[ireg,19:58,*,sk,0,0])
           d1m=reform(d0d1Ls_SR[ireg,19:58,*,sk,1,0])
           Lsm=reform(d0d1Ls_SR[ireg,19:58,*,sk,2,0])

           bias[ireg,*,*,sk,0,stat]= d0m-d0o
           bias[ireg,*,*,sk,1,stat]= d1m-d1o
           bias[ireg,*,*,sk,2,stat]= Lsm-Lso
        endfor 
     endfor 
  endfor 

  return,bias
end


pi=4.*atan(1.)
mydevice = !d.name
!P.MULTI = [0, 3,2,2]

SET_PLOT, 'ps'
DEVICE , filename='./figs/fig_bias.eps',  /encapsulated,/helvetica, decomposed=0,BITS_PER_PIXEL=8, COLOR=1,xsize=6.0*2.54,ysize=7.0*2.54, scale=1
device, SET_FONT='Helvetica',/tt_font


!p.thick=2
!x.thick=2
!y.thick=2


rgb=load_colors()
TVLCT,rgb[*,0],rgb[*,1],rgb[*,2]

pos3=fltarr(4,3)
pos3[*,2]=[0.08,0.72,0.95,0.99]
pos3[*,1]=[0.08,0.38,0.95,0.65]
pos3[*,0]=[0.08,0.04,0.95,0.31]

Nr=2
if(N_elements(d0d1Ls_SR) EQ 0) then begin
   restore,filename=save_dir+"d0d1Ls_SR.dat",/verbose
   d0d1Ls_SR[*,140,*,*,*,*]=d0d1Ls_SR[*,139,*,*,*,*]
   d0d1Ls_SR(where(d0d1Ls_SR EQ -1))=!Values.F_nan
endif 
if(N_elements(d0d1Ls_obs) EQ 0) then begin
   restore,filename=save_dir+"d0d1Ls_obs.dat",/verbose
   d0d1Ls_obs(where(d0d1Ls_obs EQ -1))=!Values.F_nan
endif

bias=calculate_bias(d0d1Ls_obs=d0d1Ls_obs,d0d1Ls_SR=d0d1Ls_SR,Nr=Nr)


strsk=['PC3','PC7','1B','NIS']
strreg=['Arctic Bridge','NWP south','NWP north','Beaufort']
strvar=['First Day Bias (days)','Last Day Bias (days)','Season Length Bias (days)']
stryref1=['Model melts late' ,'Model freezes late' ,'Model long']
stryref2=['Model melts early','Model freezes early','Model short']

;y range
yran=fltarr(2,3)
yran[*,0]=[-120,100]
yran[*,1]=[-120,100]
yran[*,2]=[-200,100]
yoff1=[-133,-133,-220]
yoff2=[-149,-149,-243]

for var=0,2 do begin
   CASE var of
      0: ytickvX=[-100,-50,0,50,100]
      1: ytickvX=[-100,-50,0,50,100]
      2: ytickvX=[-200,-150,-100,-50,0,50,100]
   ENDCASE 
   plot,fltarr(5),xrange=[0,4],position=pos3[*,var],charsize=1.5,font=1, linestyle=1, $
        xstyle=1,ystyle=1,xtickformat='(A1)',xminor=5, $
        yrange=yran[*,var],yticks=5,ytickv=ytickvX,yminor=5

   dx=0.025
   for ireg=0,3 do begin
      x0=ireg
      x=(findgen(4)-1.5)*0.20+0.5
      x2=(findgen(4)-1.5)*0.20+0.5+0.1
      for sk=0,3 do begin

         polyfill,[-1,-1,1,1,-1]*dx+x0+x[sk],[-1,1,1,-1,-1]*2*stddev(mean(bias[ireg,*,*,sk,var,0],dim=2,/nan),/nan),color=4
         polyfill,[-1,-1,1,1,-1]*dx+x0+x[sk],[-1,1,1,-1,-1]*1*stddev(mean(bias[ireg,*,*,sk,var,0],dim=2,/nan),/nan),color=3

         y0=mean(mean(bias[ireg,*,*,sk,var,0],dim=3,/nan),/nan) ;thick ice RIO (most restrictive to shipping)
         y1=mean(mean(bias[ireg,*,*,sk,var,1],dim=3,/nan),/nan) ;mean ice RIO
         y2=mean(mean(bias[ireg,*,*,sk,var,2],dim=3,/nan),/nan) ;thin ice RIO (least restrictive to shipping)
         plots,[-2,2]*dx+x0+x[sk],[y0,y0],color=0,thick=6
         plots,[-1,1]*dx+x0+x[sk],[y1,y1],color=0,thick=3
         plots,[ 0,0]*dx   +x0+x[sk],[y2,y2],color=0,thick=0
         plots,[0,0]+x0+x[sk],[y0,y2],color=0
         xyouts,x0+x[sk],yoff1[var],strsk[sk],font=1,alignment=0.5,charsize=0.9
         xyouts,4.15, 50,stryref1[var],font=1,alignment=0.5,charsize=0.9,orientation=90
         xyouts,4.15,-65,stryref2[var],font=1,alignment=0.5,charsize=0.9,orientation=90
      endfor
      xyouts,x0+0.5,yoff2[var],strreg[ireg],font=1,alignment=0.5,charsize=1.2
   endfor 
   xyouts,0.035,0.5*(pos3[1,var]+pos3[3,var]),strvar[var],font=1,alignment=0.5,charsize=1.1,orientation=90,/normal
   
endfor 


DEVICE, /close
!P.MULTI = [0, 1, 1]
SET_PLOT, mydevice



end

