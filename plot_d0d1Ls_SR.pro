function load_colors

  restore,filename="./data/color_tables/color_table_fig2.dat"
  return,[[r],[g],[b]]
end 

function load_plot_positions, reg=reg

  ppos=fltarr(4,3,4)

  ;top row, columns=r=0,1,2,3
  ppos[*,0,0]=[0.05,0.56,0.24,0.87]
  ppos[*,0,1]=[0.30,0.56,0.49,0.87]
  ppos[*,0,2]=[0.55,0.56,0.74,0.87]
  ppos[*,0,3]=[0.80,0.56,0.99,0.87]

  ;middle row, columns=r=0,1,2,3
  ppos[*,1,0]=[0.05,0.23,0.24,0.54]
  ppos[*,1,1]=[0.30,0.23,0.49,0.54]
  ppos[*,1,2]=[0.55,0.23,0.74,0.54]
  ppos[*,1,3]=[0.80,0.23,0.99,0.54]

  ;bottom rpw, columns=r=0,1,2,3  
  ppos[*,2,0]=[0.05,0.06,0.24,0.21]
  ppos[*,2,1]=[0.30,0.06,0.49,0.21]
  ppos[*,2,2]=[0.55,0.06,0.74,0.21]
  ppos[*,2,3]=[0.80,0.06,0.99,0.21]
      
  return,ppos[*,*,reg]
end 


function load_TSyear

  ;annual mean GMST calculated from each realization of LENS
  ;PI based on GMST 200year average of LENS control run 1500-1699, GMSTaa=287.122

  restore,filename="./data/auxilary/TSgmaa.dat"
  dt10_lens=fltarr(40,251)
  ;1850-2100
  dT10_lens[0,*]=smooth(TSgmaa[0,*],10,/edge_mirror)-287.122
  for r=1,39 do begin
     ;1920-2100
     dT10_lens[r,70:250]=smooth(TSgmaa[r,70:250],10,/edge_mirror)-287.122
  endfor
  
  TdateNN_lens=fltarr(40,40)
  for i=0,39 do begin
     level=i*0.1+0.1
     for r=0,39 do begin
        TdateNN_lens[i,r]= min(where(dT10_lens[r,*] gt level))
     endfor
  endfor
  
  return,TdateNN_lens
end

pro plot_row1,ppos=ppos, dx=dx, Lsm_risk=Lsm_risk,reg=reg, $
              t1=t1,t2=t2,t4=t4,tdateNN_lens=tdateNN_lens,offset=offset,th=th,lc=lc

  regstr=['Arctic Bridge','NWP Southern Route','NWP Northern Route','Beaufort Region']
  plot,dx,fltarr(141),yrange=[0,365],charsize=1.4,xrange=[1960,2100], $
       xstyle=9,ystyle=9,font=1,position=ppos[*,0],xtickformat='(A1)',xminor=2, yminor=4, xticklen=0.04
        
  plots,[0,0]+t1,[0,365],color=lc,linestyle=ls,thick=th
  plots,[0,0]+t2,[0,365],color=lc,linestyle=ls,thick=th
  plots,[0,0]+t4,[0,365],color=lc,linestyle=ls,thick=th        
  
  skc=[4,6,8,10]
  for y=0,139 do begin
     for sk=0,3 do begin 
       ;Ls1 regional mean time series
        ts=reform(Lsm_risk[reg,y,sk,*])
        yx=ts[1]                ;RIO=-10
        yn=ts[0]                ;RIO=0
        polyfill,[0,0,1,1,0]+y+1960,[yn,yx,yx,yn,yn],color=skc[sk],noclip=1
     endfor 
  endfor
  skc=[3,5,7,9]
  ;add heavy lines for sk=0 when at 0/364
  x=where(Lsm_risk[reg,0:139,0,0] EQ 364,nX) 
  for i=0,nX-1 do plots,dx(x(i)),Lsm_risk[reg,x(i),0,0]-2,color=skc[0],psym=6,symsize=0.1,thick=3

  ss=[1,1]
  cobs=[40,60,150,190]*0.
  
        ;dark line for RIO=0 (added risk to RIO=-10 shaded above)
  for sk=3,0,-1 do oplot,dx,Lsm_risk[reg,*,sk,0],color=skc[sk]        
        
  xyouts,ppos[0,0]-0.030,0.5*(ppos[1,0]+ppos[3,0]),'Season Length',font=1,alignment=0.5,orientation=90,charsize=0.90,/normal
        
  axis,yaxis=0,xrange=[0,365],color=0,ytickformat='(A1)',yminor=4,ystyle=1
  axis,yaxis=1,xrange=[0,365],color=0,ytickformat='(A1)',yminor=4,ystyle=1
  axis,xaxis=0,xrange=[1960,2100],color=0,xtickformat='(A1)',xminor=2,xticks=7,xstyle=1, xticklen=0.04
  axis,xaxis=1,xrange=[1960,2100],color=0,xtickformat='(A1)',xminor=1,xticks=1,xstyle=1, xticklen=0.04
  for i=0,7 do plots,rebin(mean(tdateNN_lens[i*5+4,*],dim=2)+offset,2),[355-(i mod 2)*10,365],color=0,linestyle=0
  tkstr=strtrim(string(findgen(4)*1.0+1.0,format='(F3.1)'),2)
  for i=0,3 do xyouts,rebin(mean(tdateNN_lens[i*10+9,*],dim=2)+offset,2),375,tkstr[i],font=1,charsize=0.75,alignment=0.5
  
  strlabel='Warming Level (!Z(00B0)C)'
  xyouts,0.5*(ppos[0,0]+ppos[2,0]),0.91,strlabel,font=1,alignment=0.5,charsize=0.85,/normal
  xyouts,0.5*(ppos[0,0]+ppos[2,0]),0.97,regstr[reg],font=1,alignment=0.5,charsize=1.1,/normal
end 


pro plot_row2, ppos=ppos, dx=dx, d0m_risk=d0m_risk,d1m_risk=d1m_risk, reg=reg, $
               t1=t1,t2=t2,t4=t4,tdateNN_lens=tdateNN_lens,offset=offset,th=th,lc=lc

  plot,dx,fltarr(141),charsize=1.2, font=1,position=ppos[*,1], $
       xrange=[1960,2100],xstyle=9,xtickformat='(A1)',xminor=2, xticklen=0.04, $
       yrange=[0,365]    ,ystyle=9,ytickformat='(A1)',yminor=1,yticks=12

  plots,[0,0]+t1,[0,365],color=lc,linestyle=ls,thick=th
  plots,[0,0]+t2,[0,365],color=lc,linestyle=ls,thick=th
  plots,[0,0]+t4,[0,365],color=lc,linestyle=ls,thick=th
  
  skc=[4,6,8,10]
  for y=0,139 do begin
     for sk=0,3 do begin
        col=skc[sk]
        ;d0 regional mean time series
        ts=reform(d0m_risk[reg,y,sk,*])
        yn=ts[1]                ;RIO=-10
        yx=ts[0]                ;RIO=0
        polyfill,[0,0,1,1,0]+y+1960,[yn,yx,yx,yn,yn],color=col,noclip=1
        ;d1 regional mean time series
        ts=reform(d1m_risk[reg,y,sk,*])
        yx=ts[1]                ;RIO=-10
        yn=ts[0]                ;RIO=0
        polyfill,[0,0,1,1,0]+y+1960,[yn,yx,yx,yn,yn],color=col,noclip=1
     endfor
     ;sk 2/1 overlap
     ts2=reform(d0m_risk[reg,y,2,*])
     ts1=reform(d0m_risk[reg,y,1,*])
     yn=ts2[1]
     yx=ts1[0]
     if (yx gt yn) then polyfill,[0,0,1,1,0]+y+1960,[yn,yx,yx,yn,yn],color=12,noclip=0
     ts2=reform(d1m_risk[reg,y,2,*])
     ts1=reform(d1m_risk[reg,y,1,*])
     yn=ts1[0]
     yx=ts2[1]
     if (yx gt yn) then polyfill,[0,0,1,1,0]+y+1960,[yn,yx,yx,yn,yn],color=12,noclip=0           
     ;sk 3/2 overlap
     ts3=reform(d0m_risk[reg,y,3,*])
     ts2=reform(d0m_risk[reg,y,2,*])
     yn=ts3[1]
     yx=ts2[0]
     if (yx gt yn) then polyfill,[0,0,1,1,0]+y+1960,[yn,yx,yx,yn,yn],color=11,noclip=0
     ts3=reform(d1m_risk[reg,y,3,*])
     ts2=reform(d1m_risk[reg,y,2,*])
     yn=ts2[0]
     yx=ts3[1]
     if (yx gt yn) then polyfill,[0,0,1,1,0]+y+1960,[yn,yx,yx,yn,yn],color=11,noclip=0           
  endfor
  skc=[3,5,7,9]
  ;add heavy lines for sk=0 when at 0/364
  x=where(d1m_risk[reg,0:139,0,0] EQ 364,nX) & for i=0,nX-1 do plots,dx(x(i)),d1m_risk[reg,x(i),0,0]-1,color=skc[0],psym=6,symsize=0.1
  x=where(d0m_risk[reg,0:139,0,0] EQ   0,nX) & for i=0,nX-1 do plots,dx(x(i)),d0m_risk[reg,x(i),0,0]+1,color=skc[0],psym=6,symsize=0.1


  ;dark line for RIO=0 (added risk to RIO=-10 shaded above)
  for sk=3,0,-1 do oplot,dx,d0m_risk[reg,*,sk,0],color=skc[sk]
  for sk=3,0,-1 do oplot,dx,d1m_risk[reg,*,sk,0],color=skc[sk]

  xyouts,ppos[0,1]-0.020,0.5*(ppos[1,1]+ppos[3,1]),'Time of Year',font=1,alignment=0.5,orientation=90,charsize=0.90,/normal
  mstr=['A','M','J','J','A','S','O','N','D','J','F','M']
  for m=0,11 do xyouts,1954,(m+0.25)/12.*365.,mstr[m],font=1,alignment=0.5,charsize=0.75

  axis,yaxis=0,xrange=[0,365],color=0,ytickformat='(A1)',yminor=1,ystyle=1,yticks=12
  axis,yaxis=1,xrange=[0,365],color=0,ytickformat='(A1)',yminor=1,ystyle=1,yticks=12
  axis,xaxis=0,xrange=[1960,2100],color=0,xtickformat='(A1)',xminor=2,xticks=7,xstyle=1, xticklen=0.04
  axis,xaxis=1,xrange=[1960,2100],color=0,xtickformat='(A1)',xminor=1,xticks=1,xstyle=1, xticklen=0.04
  for i=0,7 do plots,rebin(mean(tdateNN_lens[i*5+4,*],dim=2)+offset,2),[355-(i mod 2)*10,365],color=0,linestyle=0
end 


pro plot_row3, ppos=ppos, dx=dx, d1_risk=d1_risk, reg=reg, $
               t1=t1,t2=t2,t4=t4,tdateNN_lens=tdateNN_lens,offset=offset,th=th,lc=lc,Nr=Nr

  plot,dx,fltarr(141)+!Values.F_nan,yrange=[-0.2,1.2],charsize=1.4,xrange=[1960,2100], $
       xstyle=9,font=1,position=ppos[*,2],xtickformat='(A1)',linestyle=1, $
       ystyle=1,yticks=2,ytickv=[0,0.5,1.],yminor=5,xticks=7,xminor=2,xticklen=0.08

  plots,[0,0]+t1,[-0.2,1.2],color=lc,linestyle=ls,thick=th
  plots,[0,0]+t2,[-0.2,1.2],color=lc,linestyle=ls,thick=th
  plots,[0,0]+t4,[-0.2,1.2],color=lc,linestyle=ls,thick=th

  skc=[3,5,7,9]
  for sk=3,0,-1 do begin
     if(Nr GT 1) then Fsampled=total(Finite(d1_risk[reg,*,0:Nr-1,sk,0]),3)/float(Nr) else Fsampled=Finite(d1_risk[reg,*,0:Nr-1,sk,0])
     oplot,dx,Fsampled,color=skc[sk]
  endfor 
  xyouts,ppos[0,2]-0.0390,0.5*(ppos[1,2]+ppos[3,2]),'Navigation',font=1,alignment=0.5,orientation=90,charsize=0.90,/normal
  xyouts,ppos[0,2]-0.0240,0.5*(ppos[1,2]+ppos[3,2]),'Probability',font=1,alignment=0.5,orientation=90,charsize=0.90,/normal

  axis,xaxis=1,xrange=[1960,2100],color=0,xtickformat='(A1)',xminor=1,xticks=1,xstyle=1
  for i=0,7 do plots,rebin(mean(tdateNN_lens[i*5+4,*],dim=2)+offset,2),[1.15-(i mod 2)*0.05,1.2],color=0,linestyle=0

  for y=1960,2100,20 do xyouts,y-1,-0.43,strtrim(string(y,format='(i4)'),2),font=1,alignment=0.5,charsize=0.70
  xyouts,0.5*(ppos[0,0]+ppos[2,0]),0.005,'Year',font=1,alignment=0.5,charsize=0.85,/normal

end 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAIN


mydevice = !d.name
!P.MULTI = [0, 3,2,2]
SET_PLOT, 'ps' 
DEVICE , filename='./figs/fig_SR.eps',  /encapsulated,/helvetica, decomposed=0,BITS_PER_PIXEL=8, COLOR=1,xsize=7.5*2.54,ysize=3.75*2.54, scale=1
device, SET_FONT='Helvetica',/tt_font


!p.thick=2
!x.thick=2
!y.thick=2
th=5 ;linethickness

rgb=load_colors()
TVLCT,rgb[*,0],rgb[*,1],rgb[*,2]
lc=2

Nr=2 ;number of realizations processed
if(N_elements(d0d1Ls_SR) EQ 0) then begin
   restore,filename=save_dir+"d0d1Ls_SR.dat",/verbose
endif 

;missing final year
d0d1Ls_SR[*,140,*,*,*,*]=d0d1Ls_SR[*,139,*,*,*,*]
d0d1Ls_SR(where(d0d1Ls_SR EQ -1))=!Values.F_nan        
;combine RIO cutoffs        
d0_risk=fltarr(4,141,Nr,4,2)
d1_risk=fltarr(4,141,Nr,4,2)
Ls_risk=fltarr(4,141,Nr,4,2)
d0_risk[*,*,0:Nr-1,*,*]=d0d1Ls_SR[*,*,*,*,0,*]
d1_risk[*,*,0:Nr-1,*,*]=d0d1Ls_SR[*,*,*,*,1,*]
Ls_risk[*,*,0:Nr-1,*,*]=d0d1Ls_SR[*,*,*,*,2,*]
d0_risk(where(d0_risk EQ -1))=!Values.F_nan
d1_risk(where(d1_risk EQ -1))=!Values.F_nan
Ls_risk(where(Ls_risk EQ -1))=!Values.F_nan
;ensemble mean
d0m_risk=mean(d0_risk,dim=3,/nan)
d1m_risk=mean(d1_risk,dim=3,/nan)
Lsm_risk=mean(Ls_risk,dim=3,/nan)
        

dx=findgen(140)+1960
for reg=0,3 do begin
       
   ppos=load_plot_positions(reg=reg)
   
   ;find year each realization's 10yr running average > given temp
   TdateNN_lens=load_TSyear()   ;[wl,realization], wl=[0.1,0.2,0.3,0.4....4.0]
   
   ;+1,+2,+4K ranges and means
   offset=1850
   t1=mean(TdateNN_lens[ 9,*]+offset)+1850
   t2=mean(TdateNN_lens[19,*]+offset)+1850
   t4=mean(TdateNN_lens[39,*]+offset)+1850
        
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Top Row
   plot_row1, ppos=ppos, dx=dx, Lsm_risk=Lsm_risk,reg=reg, $
              t1=t1,t2=t2,t4=t4,tdateNN_lens=tdateNN_lens,offset=offset,th=th,lc=lc
        
   ;;;;;;;;;;;;;;;;;;;
   ;Second row (d0/d1)
   plot_row2, ppos=ppos, dx=dx, d0m_risk=d0m_risk,d1m_risk=d1m_risk, reg=reg, $
              t1=t1,t2=t2,t4=t4,tdateNN_lens=tdateNN_lens,offset=offset,th=th,lc=lc

   ;;;;;;;;;;;;;;
   ;Third row (fraction of realizations with accessible route)
   plot_row3, ppos=ppos, dx=dx,d1_risk=d1_risk,reg=reg, $
              t1=t1,t2=t2,t4=t4,tdateNN_lens=tdateNN_lens,offset=offset,th=th,lc=lc, Nr=Nr
endfor 
     
DEVICE, /close
!P.MULTI = [0, 1, 1]
SET_PLOT, mydevice 

end


