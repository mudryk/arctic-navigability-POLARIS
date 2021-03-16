function load_colors

  restore,filename="./data/color_tables/color_table_fig3.dat"
  return,[[r],[g],[b]]
end

function load_TSyear

  ;annual mean GMST calculated from each realization of LENS
  ;PI based on GMST 200year average of
  ;LENS control run 1500-1699, GMSTaa=287.122

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

function load_community_info, lonlat=lonlat,names=names
                                ;real world locations and names

  missing_lonlat=(N_elements(lonlat) EQ 0)
  missing_names =(N_elements(names) EQ 0)
  if (missing_names)  then names=0
  if (missing_lonlat) then lonlat=0
  if (missing_lonlat && missing_names) then begin
     lonlat=1
     names=0
  endif

  Clonlat=fltarr(50,2)
  Cnames=strarr(50)

  line=''
  file='./data/auxilary/northern_communities.txt'
  OPENR, lun,file,/get_lun
  READF, lun, line
  for i=0,49 do begin
     READF, lun, line
     linex=strsplit(line,/extract)
     N=n_elements(linex)
     Clonlat[i,0]=-float(linex[N-1])
     Clonlat[i,1]= float(linex[N-2])
     Cnames[i]=strjoin(linex[2:N-3],' ',/single)
  endfor
  free_lun,lun
  close,lun

  if(lonlat) then return, Clonlat
  if(names)  then return, Cnames

end

pro plot_table,Nadx=Nadx,Nd1OW=Nd1OW,NdOW=NdOW,plotOW=plotOW,carray=carray,sk=sk,nCR=nCR,pp=pp,dC=dC

  Ndmon=[31,28,31,30,31,30,31,31,30,31,30,31,31] ;Jan-Dec
  ms=findgen(12)+0.5
  
  plot,ms,fltarr(6),yrange=[0,50.],xticks=12,charsize=1.7,xrange=[0,12], xminor=1, ystyle=1,$ 
       Xtickformat='(A1)',font=1,position=pp,ytickformat='(A1)',yticks=6
  for c=0,49 do begin
     for m=0,11 do begin
        polyfill,[-0.5,-0.5,0.5,0.5,-0.5]+ms[m],[0,1,1,0,0]+49-c,color=carray[floor(Nadx[c,m,sk]*(Nadx[c,m,sk] GT 0))]
     endfor
  endfor
  if(plotOW) then begin 
     for c=0,49 do begin
        for m=0,11 do begin
           ;if (NdOW[c,m,sk] EQ Ndmon[m]) then polyfill,[-0.5,-0.5,0.5,0.5,-0.5]+ms[m],[0,1,1,0,0]+49-c,color=1,/line_fill,orientation=30
           if (NdOW[c,m,sk] GE 27.5) then polyfill,[-0.5,-0.5,0.5,0.5,-0.5]+ms[m],[0,1,1,0,0]+49-c,color=1,/line_fill,orientation=30
           if (Nd1OW[c,m,sk] GE 27.5) then polyfill,[-0.5,-0.5,0.5,0.5,-0.5]+ms[m],[0,1,1,0,0]+49-c,color=1,/line_fill,orientation=150
        endfor
     endfor
  endif 
  for x=0,12 do plots,[0,0]+x,[0,50.],thick=3,color=0
  for i=0,4 do plots,[0,12],[0,0]+nCR[i] ,thick=3,color=0 ;region1
  nCRm=0.5*(nCR+shift(nCR,1))
  monstr=['J','F','M','A','M','J','J','A','S','O','N','D']
  for m=0,11 do xyouts,ms[m],-1.5,monstr[m],font=0,alignment=0.5,charsize=1.0

  retSymbol = '!9' + String("260B) + '!X' ;" ;this is just to fake-close the double quote
  xyouts,0.5*(pp[0]+pp[2]),0.91,dC+retSymbol+'C Warming Relative to Present Day',font=0,color=0,/normal,charsize=0.8,alignment=0.5
end 


pro plot_legend

  ;regions 1-4 along right
  x0=0.95 & y0=0.88
  xyouts,x0,y0,'1',font=0,color=0,charsize=0.9,alignment=0.5,/normal
  x0=0.95 & y0=0.79
  xyouts,x0,y0,'2',font=0,color=0,charsize=0.9,alignment=0.5,/normal
  x0=0.95 & y0=0.55
  xyouts,x0,y0,'3',font=0,color=0,charsize=0.9,alignment=0.5,/normal
  x0=0.95 & y0=0.13
  xyouts,x0,y0,'4',font=0,color=0,charsize=0.9,alignment=0.5,/normal


  dx=0.03 & dy=0.02
  ;bottom  
  ;<3 days
  x0=0.05 & y0=0.02
  plots,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,linestyle=0,color=0,/normal,thick=3
  xyouts,x0+dx+0.01,y0+0.2*dy,' < 3 days',font=0,charsize=0.8,/normal,alignment=0.,color=0
  ;3-7 days
  x0=0.21 & y0=0.02 
  polyfill,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,color=6,/normal
  plots,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,linestyle=0,color=0,/normal,thick=3
  xyouts,x0+dx+0.01,y0+0.2*dy,' 3- 7 days',font=0,charsize=0.8,/normal,alignment=0.,color=0        
  ;7-14 days
  x0=0.37 & y0=0.02
  polyfill,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,color=7,/normal
  plots,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,linestyle=0,color=0,/normal,thick=3
  xyouts,x0+dx+0.01,y0+0.2*dy,' 7-14 days',font=0,charsize=0.8,/normal,alignment=0.,color=0        
  ;14-21 days
  x0=0.53 & y0=0.02
  polyfill,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,color=8,/normal
  plots,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,linestyle=0,color=0,/normal,thick=3
  xyouts,x0+dx+0.01,y0+0.2*dy,'14-21 days',font=0,charsize=0.8,/normal,alignment=0.,color=0
  ;21-28 days
  x0=0.69 & y0=0.02
  polyfill,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,color=9,/normal
  plots,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,linestyle=0,color=0,/normal,thick=3
  xyouts,x0+dx+0.01,y0+0.2*dy,'21-28 days',font=0,charsize=0.8,/normal,alignment=0.,color=0
  ;>28days
  x0=0.85 & y0=0.02
  polyfill,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,color=10,/normal
  plots,x0+[0,1,1,0,0]*dx,y0+[1,1,0,0,1]*dy,linestyle=0,color=0,/normal,thick=3
  xyouts,x0+dx+0.01,y0+0.2*dy,'>28 days',font=0,charsize=0.8,/normal,alignment=0.,color=0        

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;MAIN



Clonlat=load_community_info(/lonlat)
Cnames=load_community_info(/names)

;number of realizations processed
Nr=40                          
if(N_elements(Nad) EQ 0) then begin
   restore,filename=save_dir+"Nad.dat",/verbose
endif
if(N_elements(Nad2v1) EQ 0) then begin
   Nad2v1=fltarr(50,12,4)       ;2K vs 1K
   Nad4v1=fltarr(50,12,4)       ;4K vs 1K
   Nd1OW =fltarr(50,12,4)
   Nd2OW =fltarr(50,12,4)
   Nd4OW =fltarr(50,12,4)
   
   y1=reform(TdateNN_lens[ 9,*]-(1960-1850))
   y2=reform(TdateNN_lens[19,*]-(1960-1850))
   y4=reform(TdateNN_lens[39,*]-(1960-1850))
   d21=fltarr(Nr)
   d41=fltarr(Nr)
   d1OW=fltarr(Nr)
   d2OW=fltarr(Nr)
   d4OW=fltarr(Nr)
   for k=0,49 do begin
      for m=0,11 do begin
         mi=(m+9)mod 12
         for sk=0,3 do begin                 
            for r=0,Nr-1 do begin
               d21[r]=(Nad[k,mi,y2[r],r,sk]-Nad[k,mi,y1[r],r,sk])
               d41[r]=(Nad[k,mi,y4[r],r,sk]-Nad[k,mi,y1[r],r,sk])
               d1OW[r]=Nad[k,mi,y1[r],r,sk]
               d2OW[r]=Nad[k,mi,y2[r],r,sk]
               d4OW[r]=Nad[k,mi,y4[r],r,sk]
            endfor 
            Nad2v1[k,m,sk]=mean(d21)
            Nad4v1[k,m,sk]=mean(d41)
            Nd1OW[k,m,sk]=mean(d1OW)
            Nd2OW[k,m,sk]=mean(d2OW)
            Nd4OW[k,m,sk]=mean(d4OW)
         endfor 
      endfor 
   endfor 
endif 

   
mydevice = !d.name
!P.MULTI = [0, 1,2,2]
SET_PLOT, 'ps' 
DEVICE , filename='./figs/fig_Nad.eps',  /encapsulated,/helvetica, decomposed=0,BITS_PER_PIXEL=8, COLOR=1,xsize=7.5*2.54,ysize=7.*2.54, scale=1
device, SET_FONT='Helvetica',/tt_font


!p.thick=2
!x.thick=2
!y.thick=2


rgb=load_colors()
TVLCT,rgb[*,0],rgb[*,1],rgb[*,2]


p1=[0.23,0.10,0.57,0.90]
p2=[0.59,0.10,0.93,0.90]


carray=fltarr(32)
carray[ 0: 2]  =255             ;<3days
carray[ 3: 6]=  6               ;<3-7days
carray[ 7:13]=  7               ;1-2weeks
carray[14:20]=  8               ;2-3weeks
carray[21:27]=  9               ;3+weeks
carray[28:31]= 10               ;4+weeks


nCR=[50,45,30,4,0]

;ship class
sk=3

;2K - 1K
plot_table,Nadx=Nad2v1,Nd1OW=Nd1OW,NdOW=Nd2OW,plotOW=1,carray=carray,sk=sk,nCR=nCR,pp=p1,dC='+1' 
;4K - 1K
plot_table,Nadx=Nad4v1,Nd1OW=Nd1OW,NdOW=Nd4OW,plotOW=1,carray=carray,sk=sk,nCR=nCR,pp=p2,dC='+3' 

;row labels
device, SET_FONT='Times',/tt_font
for k=0,49 do xyouts,-13.,49.2-k,Cnames[k],font=0,color=0,alignment=1,charsize=0.65
device, SET_FONT='Helvetica',/tt_font

xyouts,0.94,0.91,'Region'  ,font=0,color=0,/normal,charsize=0.8,alignment=0.,orientation=0
skstr=[' PC3',' Supply Ship ',' Supply Ship ',' Pleasure Craft ']
xyouts,0.5*(p1[2]+p2[0]),0.965,'Northern Community Increases in'+skstr[sk]+ 'Navigability',font=0,color=0,/normal,charsize=1.2,alignment=0.5,orientation=0        

plot_legend


DEVICE, /close
!P.MULTI = [0, 1, 1]
SET_PLOT, mydevice 

end


