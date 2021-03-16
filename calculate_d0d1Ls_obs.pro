function find_d0_d1,ts=ts,th=th

  dm=mean(where(-ts ge 0.))
  x1=where(ts LE 0.)
  nX1=N_elements(x1)
  if (nX1 EQ 1) then begin
     d0=dm
     d1=dm
  endif else begin
     ;threshold has reverse logic because using -ts
     N=N_elements(ts)-1
     d0=max(where(ts[ 0: dm] GT th))
     if (d0 EQ -1) then d0=0
     d1=min(where(ts[dm:N] GT th))
     if (d1 EQ -1) then d1=N else d1=d1-1+floor(dm)
  endelse                     

  return,[d0,d1]
end 

function calculate_d0d1LsO,ts=ts,time=time,y=y
  ;routine for observations is different than for model

  i0=min(where(floor(time) EQ y))
  i1=min(where(time GE (y+1.25)) )

  tsn=-ts[i0:i1]
  ds=find_d0_d1(ts=tsn,th=0)
  if(ds[1] LT ds[0]) then begin ;multiple open periods, use last
     i0=i0+ds[0]
     tsn=-ts[i0:i1]
     ds=find_d0_d1(ts=tsn,th=0)
  endif 
  if((ds[1]-ds[0]) eq (i1-i0)) then begin ;open all year
     d0=(y+  91/365.)
     d1=(y+1+91/365.)
     Ls=1.
  endif else if ((ds[1]-ds[0]) eq 0.) then begin ;never open
     d0=-1.
     d1=-1.
     Ls= 0.
  endif else begin 
                                ;melt
     if((ds[0]+i0) gt i0) then begin 
        x0=time(ds[0]+i0)
        y0=ts[ds[0]+i0]
        m=(y0-ts[ds[0]+1+i0])/(x0-time[ds[0]+1+i0])           
        day0=x0-y0/m
     endif else day0=float(y+91./365.)
                                ;freeze
     if((ds[1]+i0) lt i1) then begin
        x0=time(ds[1]+i0)
        y0=ts[ds[1]+i0]
        m=(y0-ts[ds[1]+1+i0])/(x0-time[ds[1]+1+i0])           
        day1=x0-y0/m
     endif else day1=float(y+1.+91./365.)
     d0=day0
     d1=day1
     Ls=day1-day0
  endelse
  return,[d0,d1,Ls]
end


function calculate_d0d1LS_obsR,SIA=SIA,time=time,reg=reg
  

  RV=fltarr(4,5,3) ;|sk,ice type,stat|, stat=min,avg,max RV choice
  RV[0,*,0]=[-1,     2, 3,  3,3] ;Minimum
  RV[1,*,0]=[-3,    -2, 2,  2,3]
  RV[2,*,0]=[-5,    -3, 1,  2,3]
  RV[3,*,0]=[-6,    -4,-1,  1,3]
  RV[0,*,1]=[ 0  ,   2, 3.0,3,3] ;Average
  RV[1,*,1]=[-3  ,-0.5, 2.0,2,3]
  RV[2,*,1]=[-4.5,-1.5, 1.5,2,3]
  RV[3,*,1]=[-5.5,-3.0,-0.5,1,3]
  RV[0,*,2]=[ 1,   2,   3,  3,3] ;Max
  RV[1,*,2]=[-3,   1,   2,  2,3]
  RV[2,*,2]=[-4,   0,   1,  2,3]
  RV[3,*,2]=[-5,  -2,   0,  1,3]

  
  d0d1Ls_obsR=fltarr(51,4,3,3)

  N=N_elements(time)
  for sk=0,3 do begin 
     for stat=0,2 do begin
        ts=total(SIA*rebin(RV[sk,*,stat],1,5,N),1)*10.

        y0=0
        if(reg EQ 0) then y0 = 3
        for yy=y0,50 do begin
           y=yy+1968
           
           d0d1Ls_obsR[yy,sk,stat,*]=calculate_d0d1LsO(ts=ts,time=time,y=y)
        endfor
     endfor
  endfor
  return,d0d1Ls_obsR
end 

function read_time,file=file

  donm=[0,31,59,90,120,151,181,212,243,273,304,334,365]
  
  data = READ_CSV(file, N_TABLE_HEADER=1)
  N=N_elements(data.field01)
  
  date=intarr(N,3)
  for i=0,N-1 do date[i,*]=strsplit(data.field01[i],/extract,'/')
  
  yyyy=date[*,2]
  mm=date[*,1]
  dd=date[*,0]
  time=yyyy+(dd+donm[mm-1])/365.

  return,time
end 


function read_SIA,file=file

  data = READ_CSV(file, N_TABLE_HEADER=1)
  N=N_elements(data.field01)
  
  ;MYI/FYI/YNGI/NEWI/NOI
  SIA=fltarr(5,N)
  SIA[0,*]=data.field08/data.field06
  SIA[1,*]=data.field09/data.field06
  SIA[2,*]=data.field10/data.field06
  SIA[3,*]=data.field11/data.field06
  SIA[4,*]=(data.field06-data.field07)/data.field06

  return,SIA
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAIN

files=strarr(4)
files[0]="./data/observations/Observed_Sea_Ice_Data_AB.csv"
files[1]="./data/observations/Observed_Sea_Ice_Data_NWPs.csv"
files[2]="./data/observations/Observed_Sea_Ice_Data_NWPn.csv"
files[3]="./data/observations/Observed_Sea_Ice_Data_BF.csv"

save_data=0
save_dir="./data/out/"

;regions,years (1968-2018),sk,stat
d0d1Ls_obs =fltarr(4,51,4,3,3)-1

for reg=0,3 do begin 

   time=read_time(file=files[reg])
   SIA=read_SIA(file=files[reg])
   
   d0d1Ls_obs[reg,*,*,*,*]=calculate_d0d1Ls_obsR(SIA=SIA,time=time,reg=reg)

endfor

if(save_data) then begin
   if(file_test(save_dir)) then save,filename=save_dir+"d0d1Ls_obs.dat",d0d1Ls_obs
endif 

end
