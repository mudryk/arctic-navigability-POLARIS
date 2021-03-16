function load_year_full, year=year, exp=exp, var=var

  dir="./data/model/"
  case exp of
     '20TR': base='b.e11.B20TRC5CNBDRD.f09_g16'
     'RCP' : base='b.e11.BRCP85C5CNBDRD.f09_g16'
  endcase 

  yearstr=strtrim(string(year),2)

  v =fltarr(320,104,365,40)
  rl=['r001-010','r011-020','r021-030','r031-040']
  for set=0,3 do begin 
     rs=indgen(10)+set*10
     rli=rl[set]

     file=dir+rli+'/'+base+'.'+var+'.'+rli+'.NNN.'+yearstr+'.nc'
     nCDFObject = Obj_New('NCDF_DATA', file)
     v10=ncdfObject -> ReadVariable(var)
     v10=v10*(v10 LT 500)
     v[*,*,*,rs] = v10
  endfor
  
  return,v
end 

function load_year_r2, year=year, exp=exp, var=var

  dir="./data/model/"
  case exp of
     '20TR': base='b.e11.B20TRC5CNBDRD.f09_g16'
     'RCP' : base='b.e11.BRCP85C5CNBDRD.f09_g16'
  endcase 

  yearstr=strtrim(string(year),2)

  rli='r001-002'

  file=dir+rli+'/'+base+'.'+var+'.'+rli+'.NNN.'+yearstr+'.nc'
  nCDFObject = Obj_New('NCDF_DATA', file)
  v=ncdfObject -> ReadVariable(var)
  v=v*(v LT 500)
  
  return,v
end 

function load_model_aux, name=name

  filename="./data/aux/model_aux.nc"
  nCDFObject = Obj_New('NCDF_DATA', filename)
  var=ncdfObject -> ReadVariable(name)
  
  return,var
end

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

function convert_hice_rio, hice=hice, sk=sk

  ;RV table: 12 ice types/thickness ranges, 4 ship classes
  ; 0,+0-10cm, 10-15,15-30, 30-50,50-70,70-95,95-120, 120-200, 200-250,250-300,>300
  RV=fltarr(12,4)              
  ;RV[*,0]= [+3.0,+3.0,+3.0,+3.0,+2.0,+2.0,+2.0,+2.0,+2.0,+1.0, 0.0,-1.0] ;PC3
  ;RV[*,1]= [+3.0,+2.0,+2.0,+2.0,+1.0,+1.0, 0.0,-1.0,-2.0,-3.0,-3.0,-3.0] ;PC7
  ;RV[*,2]= [+3.0,+2.0,+2.0,+1.0, 0.0,-1.0,-2.0,-3.0,-4.0,-5.0,-6.0,-6.0] ;1b
  ;RV[*,3]= [+3.0,+1.0, 0.0,-1.0,-2.0,-3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-8.0] ;GTC
  ;originally used values
  RV[*,0]= [+3.0,+3.0,+3.0,+3.0,+2.0,+2.0,+2.0,+2.0,+2.0,+1.0, 0.0,-1.0] ;PC3
  RV[*,1]= [+3.0,+2.0,+2.0,+2.0,+1.0,+1.0, 0.0,-1.0,-2.0,-3.0,-3.0,-3.0] ;PC7
  RV[*,2]= [+3.0,+2.0,+2.0,+1.0, 0.0,-1.0,-2.0,-3.0,-3.0,-4.0,-5.0,-5.0] ;1b
  RV[*,3]= [+3.0,+1.0, 0.0,-1.0,-2.0,-3.0,-3.0,-3.0,-4.0,-5.0,-6.0,-6.0] ;GTC
  
  
  ; convert thickness to RIO
  ; RIO=C_i*RV_i, but only 1 ice type per grid cell (average thickness)
  ; multiplied *10 to put units into "tenths"
  RIO =(RV[ 0,sk]*(hice eq 0.00)             + $
        RV[ 1,sk]*(hice gt 0.00)*(hice lt 0.10) + $
        RV[ 2,sk]*(hice ge 0.10)*(hice lt 0.15) + $
        RV[ 3,sk]*(hice ge 0.15)*(hice lt 0.30) + $
        RV[ 4,sk]*(hice ge 0.30)*(hice lt 0.50) + $
        RV[ 5,sk]*(hice ge 0.50)*(hice lt 0.70) + $              
        RV[ 6,sk]*(hice ge 0.70)*(hice lt 0.95) + $
        RV[ 7,sk]*(hice ge 0.95)*(hice lt 1.20) + $               
        RV[ 8,sk]*(hice ge 1.20)*(hice lt 2.00) + $
        RV[ 9,sk]*(hice ge 2.00)*(hice lt 2.50) + $
        RV[10,sk]*(hice ge 2.50)*(hice lt 3.00) + $               
        RV[11,sk]*(hice ge 3.00))*10.

  return,RIO
end 

function calculate_d0d1Ls,ts=ts,threshold=threshold
  ;for 1 time series (can be for composite region RIO or grid cell RIO)

  dm= mean(where(-ts GE 0.))    ;negative of ts was sent for rest of routine
  x1=where(ts LE 0.)
  nX1=N_elements(x1)
  if (nX1 EQ 1) then begin
     d0=dm
     d1=dm
     Ls=0
  endif else begin
     N=N_elements(ts)-1
     d0=max(where(ts[ 0: dm] GT threshold))
     if (d0 EQ -1) then d0=0
     d1=min(where(ts[dm:N] GT threshold))
     if (d1 EQ -1) then d1=N else d1=d1-1+floor(dm)
     Ls=d1-d0
  endelse

  return, [d0,d1,Ls]
end


function calculate_d0d1Ls_map, RIO=RIO, Nr=Nr
  ;based on RIO of each grid cell

  N=320L*104
  d0d1Ls_map=fltarr(320L*104,Nr,3)
  
  for r=0,Nr-1 do begin
     for i=0,N-1 do begin
        tsn=-RIO[i,*,r]           
        d0d1Ls=calculate_d0d1Ls(ts=tsn,th=0.)
        d0d1Ls_map[i,r,*]=d0d1Ls
     endfor 
  endfor 
    
  return, reform(d0d1Ls_map,320,104,Nr,3)
end   

function load_shipping_routes
   ;binary mask for all mode grid cells within all shipping route
  
   tlon = load_model_aux(name='TLON')
   tlat = load_model_aux(name='TLAT')
   tlon(where(tlon GT 500.))=0.
   tlat(where(tlat GT 500.))=0.

   shipping_routes=fltarr(320,104,4)
   files=strarr(4)
   ;list of vertices forming polygons for region
   files[0]="./data/aux/arcticbridge_region.csv"
   files[1]="./data/aux/nwp_s_region.csv"
   files[2]="./data/aux/nwp_n_region.csv"
   files[3]="./data/aux/beaufort_region.csv"

   for route=0,3 do begin       
      data=read_csv(files[route])
      xp=data.field1+360.
      yp=data.field2
      for i=0,319 do begin
         for j=0,103 do begin
            ;set to 1 if within polygon
            shipping_routes[i,j,route]= is_inside(tlon[i,j],tlat[i,j],xp,yp)
         endfor
      endfor
   endfor 

   return,shipping_routes
end 

function load_SRxw, SR=SR,route=route
  ;model row/cols and grid cell areas with a specific route
  
  tmask = load_model_aux(name='tmask')
  tmask=tmask[*,*,0]
  tmask(where(tmask GT 2.,nX))=0.

  tarea = load_model_aux(name='tarea')
  tarea=tarea[*,*,0]

  x=where (SR[*,*,route]*tmask GT 0.)
  w=tarea(x)

  ;locations (row/column) and grid cell areas for all
  ;for all model grid cells within shipping route polygon
  return,[[x],[w]]

end

function calculate_d0d1Ls_SR, RIO=RIO, Nr=Nr
  ;based on area-weighted RIO of composite shipping routes (four considered)
  
  d0d1Ls_SR=fltarr(4,Nr,3,2)
  SR=load_shipping_routes()
  for route=0,3 do begin
     xw=load_SRxw(SR=SR,route=route)
     ;grid cell locations for route
     x=xw[*,0] 
     N=N_elements(x)     
     ;grid cell areas
     w=rebin(xw[*,1]/1.e9,N,365)
     atot=total(w[*,0])
     ci=w/atot

     for r=0,Nr-1 do begin
        RIOx=     RIO[x,*,r] 
        tsn=-total(RIOx*ci,1)      
        ;RIO=0 threshold
        d0d1Ls=calculate_d0d1Ls(ts=tsn,threshold=0)
        d0d1Ls_SR[route,r,*,0]=d0d1Ls
        ;RIO=-10 threshold (put in positive bc ts=tsn)
        d0d1Ls=calculate_d0d1Ls(ts=tsn,threshold=10)
        d0d1Ls_SR[route,r,*,1]=d0d1Ls        
     endfor
  endfor
  
  return,d0d1Ls_SR
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
  file='./data/aux/northern_communities.txt'
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

function get_community_locations
  ;find nearest ocean grid cell in model
  
  ;auxilary model info
  tmask = load_model_aux(name='tmask')
  tmask=tmask[*,*,0]
  tmask(where(tmask GT 2.))=0.
  tlon = load_model_aux(name='TLON')
  tlat = load_model_aux(name='TLAT')

  ;lon/lat of real-world locations
  Clonlat=load_community_info(/lonlat)

  ;nearest model gridcell over water
  Mlonlat=fltarr(50,2)
  for k=0,49 do begin
     lon0=Clonlat[k,0]+360
     lat0=Clonlat[k,1]
     ;nearest model gridcell (probably over land)
     dx=tlon-lon0
     dy=tlat-lat0
     dr=sqrt(dx*dx+dy*dy)
     ij=where(dr EQ min(dr))
     jj=ij/320
     ii=ij- ij/320 * 320
     ;fan out around
     dr=fltarr(7,7)+!Values.F_nan
     for j=-3,3 do begin
        for i=-3,3 do begin
           ;only look over water
           if(tmask[ii+i,jj+j] EQ 1) then begin
              dx=tlon[ii+i,jj+j]-tlon[ii,jj]
              dy=tlat[ii+i,jj+j]-tlat[ii,jj]
              dr[i+3,j+3]=sqrt(dx*dx+dy*dy)
           endif
        endfor
     endfor
     ij=where(dr EQ min(dr))
     ;ocean cell = original cell + offsets
     Mlonlat[k,1]=jj + (ij/7 -3)
     Mlonlat[k,0]=ii + (ij - ij/7*7 -3)
  endfor 
  return,Mlonlat  
end 

function calculate_N_access_days,RIO=RIO,Nr=Nr

  ;starting day of each month for ice year (AMJJASONDJFM)
  donm=[0,30,61,91,122,153,183,214,244,275,306,334,365]
  
  Nad_year=fltarr(50,12,Nr)
  locs=get_community_locations()  
  for k=0,49 do begin
     i=locs[k,1]*320+locs[k,0]
     for r=0,Nr-1 do begin
        ;daily accessibility for community location
        ts=RIO[i,*,r]
        for m=0,11 do begin
           Nad_year[k,m,r]=total(ts[donm[m]:donm[m+1]-1] GT 0.)
        endfor 
     endfor 
  endfor 
  
  return,Nad_year
  
end 

function load_TSyear  
  ;annual mean GMST calculated from each realization of LENS
  ;PI based on GMST 200year average of LENS control run 1500-1699,
  ;GMSTaa=287.122

  restore,filename="./data/aux/TSgmaa.dat"
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


function calculate_d0d1Ls_WL, d0d1Ls_in=d0d1Ls_in,year0=year0,Nr=Nr

  tlon = load_model_aux(name='TLON')
  tlat = load_model_aux(name='TLAT')
  tlon(where(tlon GT 500.))=0.
  tlat(where(tlat GT 500.))=0.
   
  d0d1Ls_WL=fltarr(320,104,Nr,4,3,3) ;[x,y,r,sk,WL,var]

  ;find year each realization's 10yr running average > given temp
  TdateNN_lens=load_TSyear() ;[wl,realization], wl=[0.1,0.2,0.3,0.4....4.0]

  ;gather navigability stats at appropriate year for each realization
  wl=[9,19,39] ;+1.0,+2.0,+4.0
  for sk=0,3 do begin
     print,sk
     for r=0,Nr-1 do begin
        for i=0,2 do begin 
           wli=wl[i]
           y1=TdateNN_lens[wli,r]-(year0-1850)
           d0d1Ls_WL[*,*,r,sk,i,*]=mean(d0d1Ls_in[*,*,y1-1:y1+1,r,sk,*],dim=3)
        endfor 
     endfor 
  endfor 

  ;average
  d0d1Ls_WL(where(d0d1Ls_WL EQ -1))=!Values.F_nan
  d0d1Ls_WL=mean(d0d1Ls_WL,dim=3,/nan)
  
  return, d0d1Ls_WL
end 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAIN

   full=0 ;process all 40 CESM-LE realizations (~100GB - data not provided)
   r2=1   ;process 2 CESM-LE realization (r=0,1) data linked(?)
   save_data=0
   save_dir="./data/out/"  

   Nr=40
   if(r2) then Nr=2
   Nt=141
   
   ;output arrays
   ;First day of transit/last day of transit, season length
   ;for all model grid cells
   d0d1Ls_map =fltarr(320,104,Nt,Nr,4,3)
   ;and for 4 Canadian shipping routes (composite regions)
   d0d1Ls_SR =fltarr(4,Nt,Nr,4,3,2)
   ;Number of monthly ship-access days
   Nad=fltarr(50,12,Nt,Nr,4)

   year0=1960
   hice =fltarr(320,104,365,Nr)
   ;apr-dec of previous calendar year           
   if(full) then hi=load_year_full(year=year0,exp='20TR',var='hi_d')
   if(r2  ) then hi=load_year_r2(year=year0,exp='20TR',var='hi_d')
   hice[*,*,0:274,*] = hi[*,*,90:364,*]

   for y=0,139 do begin   
   ;for y=0,1 do begin
      print,y+year0+1
      exp='20TR'
      if (y+year0+1 GE 2006) then exp='RCP'

      ;jan-mar of current calendar year      
      if(full) then hi=load_year_full(year=y+year0+1,exp=exp,var='hi_d')
      if(r2  ) then hi=load_year_r2(year=y+year0+1,exp=exp,var='hi_d')
      hice[*,*,275:364,*]=hi[*,*,  0:89,*]

      
      if(((y+year0+1) EQ 1969)&& full) then begin
         ;missing output from 1969, realization #35 (r=34 in code - Cindexing)
         ;see https://www.cesm.ucar.edu/projects/community-projects/LENS/known-issues.html         
         ;only days 333,363,364 = 0 for hi_d
         hice[*,*,333,34]=0.5*(hice[*,*,332,34]+hice[*,*,334,34])
         hice[*,*,363:364,34]=hice[*,*,[362,362],34]
      endif 
      
      ;perform a range of calculations for each ship class
      for sk=0,3 do begin 
         ;convert hice to RIO for ship class sk
         RIO=convert_hice_rio(hice=hice,sk=sk)
         RIO=reform(RIO,320L*104,365,Nr)

         ;now do calculations
         ;calc 1: (full d0,d1,Ls maps) for 1981-2100
         d0d1Ls_map[*,*,y,*,sk,*]=calculate_d0d1Ls_map(RIO=RIO,Nr=Nr)
         
         ;calc 2: shipping routes
         d0d1Ls_SR[*,y,0:Nr-1,sk,*,*]=calculate_d0d1Ls_SR(RIO=RIO,Nr=Nr)

         ;calc 3: monthly ship-access days for 50 communities
         Nad[*,*,y,*,sk]=calculate_N_access_days(RIO=RIO,Nr=Nr)
         
         ;calc 3 community accessible days
      endfor 
      ;apr-dec added to following ice year
      hice[*,*,0:274,*]=hi[*,*,90:364,*]
   endfor 

   ;convert years to +1,+2,+4K above PI
   d0d1Ls_WL = calculate_d0d1Ls_WL(d0d1Ls_in=d0d1Ls_map,year0=year0,Nr=Nr)
   
   ;option to save and reload various parts of the output   
   if(save_data) then begin
      if(file_test(save_dir)) then begin
         if(full) then begin 
            ;save,filename=save_dir+"d0d1Ls_map.dat",d0d1Ls_map
            save,filename=save_dir+"d0d1Ls_WL.dat",d0d1Ls_WL
            save,filename=save_dir+"d0d1Ls_SR.dat",d0d1Ls_SR
            save,filename=save_dir+"Nad.dat",Nad
         endif
         if(r2) then begin 
            save,filename=save_dir+"d0d1Lsr2_map.dat",d0d1Ls_map
            save,filename=save_dir+"d0d1Lsr2_SR.dat",d0d1Ls_SR
            save,filename=save_dir+"Nadr2.dat",Nad
         endif         
      endif 
   endif 


end 



