function is_inside, x, y, px, py
   ;  x, y  - x,y coordinate of the point.
   ; px,py - The x,y coordinates of the polygon vertices
   ; return value = 1 if the point is inside the polygon, otherwise 0

  sx = Size(px)
  sy = Size(py)
  IF (sx[0] EQ 1) THEN NX=sx[1] ELSE RETURN, -1 ; Error if px not a vector
  IF (sy[0] EQ 1) THEN NY=sy[1] ELSE RETURN, -1 ; Error if py not a vector
  IF (NX EQ NY) THEN N = NX ELSE RETURN, -1     ; Incompatible dimensions

  polysides = nx

  oddNodes=0
  j=polysides-1

  for i=0L,polysides-1 do begin
     if (((py[i] LT y && py[j]  GE y)||(py[j] LT y && py[i] GE y))) then begin
        if (px[i]+(y-py[i])/(py[j]-py[i])*(px[j]-px[i]) LT x) then oddNodes=~oddNodes
     endif
     j=i
  endfor

  IF (oddNodes) THEN RETURN, 1 ELSE RETURN, 0
end    

pro plot_stipple,di=di,dj=dj,path=path,xy=xy,posp=posp,ss=ss,color=col

  xmax=21590.
  ymax=27940.
  
  xc=cos(findgen(21)/20.*2*!pi)
  xs=sin(findgen(21)/20.*2*!pi)
  usersym,xc,xs,/fill
  
  for i=posp[0],posp[2],di/xmax do begin 
     for j=posp[1],posp[3],dj/ymax do begin 
        polyin_count=0 
        for p=0L, (N_ELEMENTS(path) - 1 ) do begin 
           S=[indgen(path[p].N), 0] 
           in=is_inside(i,j,reform(xy[0,PATH[p].OFFSET + S]),reform(xy[1,PATH[p].OFFSET + S])) 
           if (in) then polyin_count=polyin_count+1 
        endfor 
        if(polyin_count EQ 1) then plots,i,j,psym=8,symsize=ss,color=col,/normal
     endfor
  endfor   
end

function load_model_aux, name=name

  filename="./data/auxilary/model_aux.nc"
  nCDFObject = Obj_New('NCDF_DATA', filename)
  var=ncdfObject -> ReadVariable(name)

  return,var
end

function load_colors

  restore,filename="./data/color_tables/color_table_fig1.dat"
  return,[[r],[g],[b]]  
end 

function load_plot_positions

  ppos=fltarr(4,12)
  ppos[*,0] =[0.0400,0.645,0.2900,0.895]
  ppos[*,1] =[0.2766,0.645,0.5266,0.895]
  ppos[*,2] =[0.5133,0.645,0.7633,0.895]
  ppos[*,3] =[0.7500,0.645,1.0000,0.895]
  ppos[*,4] =[0.0400,0.37,0.2900,0.62]
  ppos[*,5] =[0.2766,0.37,0.5266,0.62]
  ppos[*,6] =[0.5133,0.37,0.7633,0.62]
  ppos[*,7] =[0.7500,0.37,1.0000,0.62]        
  ppos[*,8] =[0.0400,0.095,0.2900,0.345]        
  ppos[*,9] =[0.2766,0.095,0.5266,0.345]
  ppos[*,10]=[0.5133,0.095,0.7633,0.345]
  ppos[*,11]=[0.7500,0.095,1.0000,0.345]

  return,ppos
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


function regrid_d0d1Ls_to_lonlat, var_in=var_in,lon_out=lon_out,lat_out=lat_out

  tlon = load_model_aux(name='TLON')
  tlat = load_model_aux(name='TLAT')
  tlon(where(tlon GT 500.))=0.
  tlat(where(tlat GT 500.))=0.
   
  donm=[0,31,59,90,120,151,181,212,243,273,304,334,365]

  Nx=N_elements(lon_out)
  Ny=N_elements(lat_out)

  var_out=fltarr(Nx,Ny,4,3) ;[x,y,sk,WL]

  triangulate, tlon,tlat,tri
  for sk=0,3 do begin
     print,sk
     for i=0,2 do begin 
        var_out[*,*,sk,i]=griddata(tlon,tlat,var_in[*,*,sk,i],/nearest_neighbor,triangles=tri,xout=lon_out,yout=lat_out,/grid)
     endfor 
  endfor 

  return,var_out
end 


pro plot_colorbar1

  ;colorbar1
  dx=0.03 & x0=0.095 & dy=0.007 & y0=0.925 ;above plot1
  color=[27,28,29,30]
  label=['364','273','182','91','1']
  for i=0,3 do begin
     polyfill,[0,dx,dx,0,0]+x0+i*dx,[-1,-1,1,1,-1]*dy+y0,color=color[i],/normal
  endfor

  dx1=[-0.005,0.005,0.005,-0.005,-0.005]
  plots,[0.,dx,dx,0.,0.]*4.+x0+dx1,[-1,-1,1,1,-1]*dy+y0,color=0,/normal
  polyfill,[0.,1,1,0.,0.]*0.005+x0+4*dx,[-1,-1,1,1,-1]*dy+y0,color=0,/normal
  for i=0,4 do begin
     xyouts,x0+0.00+dx*i,y0-0.025,label[i],/normal,charsize=1.2,font=1,alignment=0.5
  endfor
  xyouts,0.155,0.94,'Season Length [days]',font=1,charsize=1.2,alignment=0.5,/normal

end 


pro plot_colorbar2
  
   ;colorbar2
  dx=0.028571 & x0=0.225 & dy=0.01 & y0=0.04
  color=[15,16,17,18,19,20,21,22,23,24,25,26]
  label=['0','15','30','45','60','90','120','150','180','240','300','365']
  for i=0,3 do begin
     polyfill,[0,dx,dx,0,0]+x0+i*dx,[-1,-1,1,1,-1]*dy+y0,color=color[i],/normal
  endfor
  for i=4,7 do begin
     polyfill,[0,2*dx,2*dx,0,0]+x0+4*dx+(i-4)*dx*2.,[-1,-1,1,1,-1]*dy+y0,color=color[i],/normal
  endfor
  for i=8,10 do begin
     polyfill,[0,3*dx,3*dx,0,0]+x0+12*dx+(i-8)*dx*3.,[-1,-1,1,1,-1]*dy+y0,color=color[i],/normal
  endfor
  ;2mo,6mo highlighting contours        
  plots,[ 4, 4]*dx+x0,[-1,1]*dy+y0,color=3,linestyle=0,thick=10,/normal
  plots,[12,12]*dx+x0,[-1,1]*dy+y0,color=4,linestyle=0,thick=10,/normal
  for i=0,4 do begin
     xyouts,x0-0.01+dx*i,y0-0.035,label[i],/normal,font=1,charsize=1.5
  endfor
  for i=5,8 do begin
     xyouts,x0-0.01+4*dx+(i-4)*dx*2.,y0-0.035,label[i],/normal,font=1,charsize=1.5
  endfor
  for i=9,11 do begin
     xyouts,x0-0.01+12*dx+(i-8)*dx*3.,y0-0.035,label[i],charsize=1.5,/normal,font=1
  endfor        
  xyouts,0.525,0.06,'Additional Season Length [days]',font=1,charsize=1.5,alignment=0.5,/normal
end 

;;;;;;;;;;;;;;;;;;;;
;main

mydevice = !d.name
!P.MULTI = [0, 3,2,2]
SET_PLOT, 'ps' 
DEVICE , filename='./figs/fig_map.eps',  /encapsulated,/helvetica, decomposed=0,BITS_PER_PIXEL=8, COLOR=1,xsize=11.*2.54,ysize=7.*2.54, scale=1
device, SET_FONT='Helvetica',/tt_font


!p.thick=2
!x.thick=2
!y.thick=2

rgb=load_colors()
TVLCT,rgb[*,0],rgb[*,1],rgb[*,2]
land_color=1
coast_outline=2
color_2mon=3
color_6mon=4

        
ppos=load_plot_positions()

clevs4=[-5,0.5,91.,183.,274.,363.5,366]
carray_hcl4=[0,30,29,28,27,255]        
dclevs=[-10,-1,1,15,30,45,60,90,120,150,180,240,300,365]
dcarray_hcl=[255,255,15,16,17,18,19,20,21,22,23,24,25,26]

lon_out=findgen(621)*0.25-190.+360.
lat_out=findgen(221)*0.25+35.   

;regrid model data
if(N_elements(d0d1Ls_WL) EQ 0) then begin 
   restore,filename=save_dir+"d0d1Ls_WL.dat",/verbose
endif
if(N_elements(Lsm_out) EQ 0) then begin
   ;only regrid season length
   Lsm_out=regrid_d0d1Ls_to_lonlat(var_in=d0d1Ls_WL[*,*,*,*,2],lon_out=lon_out,lat_out=lat_out)
endif 
        

;;;;;;;;
;plotting        

;Ref ship class: GTC at +1K
aref=Lsm_out[*,*,3,0]
sk=[3,2,1,0]  ;GTC,1B,PC7,PC3
for wl=0,2 do begin 
   for k=0,3 do begin
      skk=sk[k]
      p=k+4*wl
      MAP_SET, /CONIC, 55, -75, -25,STANDARD_PARALLELS=[50,80], /ISOTROPIC, LIMIT=[50, 210, 90, 320], pos=ppos[*,p],/noerase,/noborder,clip=1
      daref=Lsm_out[*,*,skk,wl]-aref
      ;upper left plot only
      if((k EQ 0)&&(wl EQ 0)) then contour, aref,lon_out,lat_out,cell_fill=1, /overplot, levels=clevs4,c_colors=carray_hcl4
      ;all other positions
      if((k GT 0)||(wl GT 0)) then begin
         ;difference
         contour,daref,lon_out,lat_out,cell_fill=1, /overplot, levels=dclevs,c_colors=dcarray_hcl
         contour,daref,lon_out,lat_out,cell_fill=0, /overplot, levels=[60,180],c_colors=[color_2mon,color_6mon],c_thick=3
              
         ;stipple where raw array > 0.
         posp=ppos[*,p]
         posp[2]=posp[2]-0.03
         posp[3]=posp[3]-0.005
         a=(aref LT 364.)*(Lsm_out[*,*,skk,wl] GT 362.);*(1.-mask)
         contour,a,lon_out,lat_out,PATH_INFO=path, PATH_XY=xy, xstyle=1,ystyle=1,levels=0.2,/overplot  
         plot_stipple, di=100,dj=200,path=path,xy=xy,posp=posp,ss=0.25,color=0
      endif
      map_continents,/fill_continents,color=land_color
      map_continents,/coasts,color=coast_outline,limit=[40,180,88,330]

   endfor 
endfor 

SKstr=['Not Ice Strengthened','Category C: 1B','Category A: PC7','Category A: PC3']
for i=0,3 do xyouts,0.5*(ppos[0,i]+ppos[2,i])-0.02,0.975,SKstr[i],alignment=0.5,font=1,charsize=1.5,/normal
WLstr=['+1!Z(00B0)C','+2!Z(00B0)C','+4!Z(00B0)C']
for i=0,2 do xyouts,0.015,0.5*(ppos[1,i*4]+ppos[3,i*4]),WLstr[i],font=1,alignment=0.5,orientation=90,charsize=1.5,/normal


plot_colorbar1

plot_colorbar2
       
        
        DEVICE, /close
        !P.MULTI = [0, 1, 1]
        SET_PLOT, mydevice 

end


