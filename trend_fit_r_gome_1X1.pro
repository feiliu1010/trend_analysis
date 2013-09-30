;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;define the function of the form F(x)=a+b*x+c*sin(d*x+e)
;a and d is known
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro gfunct, X, A, F, pder

ax = sin( 12*!DPI/6 * X + A[2] )
bx = cos( 12*!DPI/6 * X + A[2] )
F = A[0] * X + A[1] * ax

;If the procedure is called with 4 parameters, calculate the partial derivatives
if n_params() ge 4 then $
  pder= [[X], [ax], [A[1] * bx]]

end

pro trend_fit_r_gome_1x1
;input data
nlon_g=360
nlat_g=180
;global density
global = fltarr(nlon_g,nlat_g)

nlon = 64
nlat = 40
;china density,limit=[15,72,55,136]
no2 = fltarr(nlon,nlat)
;star point of grid for China in global map
slon = (72+180)/1
slat = (90-55)/1

flag = 0U
flag2 = 0U
For year = 1996, 2005 do begin
For month = 1,12 do begin

Yr4  = string(year,format='(i4.4)')
Mon2 = string(month,format='(i2.2)')
nymd = year * 10000L + month * 100L + 1 * 1L

if nymd eq 19960101 then continue
if nymd eq 19960201 then continue
if nymd eq 19960301 then continue
if nymd eq 19980101 then continue
if nymd eq 20050301 then continue
if nymd eq 20050401 then continue
if nymd eq 20050501 then continue
if nymd eq 20050601 then continue
if nymd eq 20050701 then continue
if nymd eq 20050801 then continue
if nymd eq 20050901 then continue
if nymd eq 20051001 then continue
if nymd eq 20051101 then continue
if nymd eq 20051201 then continue
header = strarr(7,1)
if (year lt 2003) or ((year eq 2003) and (month lt 4)) then begin
filename ='/home/liufei/Data/Satellite/NO2/GOME/no2_'+Yr4+Mon2+'_1x1.asc'
endif else begin
filename ='/home/liufei/Data/Satellite/NO2/SCIAMACHY/no2_'+Yr4+Mon2+'_1x1.asc'
endelse 
openr,lun,filename,/get_lun
readf,lun,header,global
no2 = global[slon:slon+nlon-1,slat:slat+nlat-1]/100

if flag  then begin
;m is used for counting the total number of months
m = m+1
;define no2_month to save no2 for all months
no2_month =[ [no2_month],[no2] ]
endif else begin
m = 1U
no2_month = no2
flag = 1U
endelse
;print,Yr4,Mon2
close,/all

endfor
endfor
print,'m',m
print,'size of no2_month',size(no2_month)



;convert no2_month to 3-D array no2_data
no2_data =fltarr (nlon,nlat,m)
For num = 0,m-1 do begin
no2_data[*,*,num] = no2_month[0:nlon-1,nlat*num:nlat*(num+1)-1]
endfor
;find no2_data with reasonable value( filter value < 0)
loc = where(no2_data lt 0)
;print,'loc',loc
no2_data[loc] = -999.0
undefine,loc
print,'size of no2_data',size(no2_data)

;********************
;;maximum=fltarr(nlon,nlat)
;deduce the fisrt reasonal data of January( constant in fit function)
;;For J = 0,nlat-1 do begin
;;FOR I = 0,nlon-1 do begin
;For J = 10,10 do begin
;FOR I = 18,18 do begin
;;	loc = where(no2_data[I,J,*] gt 0)
	;if (array_equal(loc,[-1])) then begin
	;no2_data[I,J,*] = -999.0
	;endif else begin
;;	first = no2_data[loc[0]]
	;print,'first',first
	;print,'no2_data[I,J,2]_before',no2_data[I,J,2]
;;	no2_data[I,J,*]=no2_data[I,J,*]-first
;;		if maximum[I,J] then begin
;;			maximum[I,J]=[maximum[I,J],first]
;;		endif else begin
;;		maximum[I,J]=[first]
;;		endelse
	;print,'no2_data[I,J,2]_after',no2_data[I,J,2]
	;endelse
;;	undefine,loc
;;	undefine,first
;;endfor
;;endfor
;*************************

;regard the data of 199604
nodata = 0
For J = 0,nlat-1 do begin
FOR I = 0,nlon-1 do begin
	if no2_data[I,J,0] ge 0 then begin
		first = no2_data[I,J,0] 
		no2_data[I,J,*]=no2_data[I,J,*]-first
	endif else begin
	no2_data[I,J,*]=-999.0
	nodata=nodata+1
	endelse
	if ( I eq 44) and (J eq 27) then begin
		print,'first,A[44,27]',first
	endif
endfor
endfor
;*********************************

;regard the data of 1996 mean
;;nodata = 0
;;For J = 0,nlat-1 do begin
;;FOR I = 0,nlon-1 do begin
;;	if not array_equal(where(no2_data[I,J,0:8] ge 0.0),[-1]) then begin
;;		first=mean(no2_data[where(no2_data[I,J,0:8] ge 0.0)]) 
;;       	no2_data[I,J,*]=no2_data[I,J,*]-first
;;	endif else begin
;;		 no2_data[I,J,*]=-999.0
;;		nodata=nodata+1
;;	endelse 
;;endfor
;;endfor
;*********************************
print,'nodata in 199604 ',nodata


print,'size of no2_data',size(no2_data)
print,'MAX OF no2_data',MAX(no2_data[*,*,*]),'MIN OF no2_data',MIN(no2_data[*,*,*])
;print,'MAX OF FIRST',MAX(maximum)
;print,'loc of -999.0',where(no2_data eq -999.0)
;Compute the fit to the function we have just defined. First, define the independent and dependent variables
;define para to save fit parameters
para = fltarr(nlon,nlat)

;deltab=deltan/n^(3/2)*((1+cor)/(1-cor))^1/2
;deltan is the variance in the remainder
;n is the number of months with availble data
cor=fltarr(nlon,nlat)
deltan=fltarr(nlon,nlat)
deltab=fltarr(nlon,nlat)


nodata_para=0
For J = 0,nlat-1 do begin
FOR I = 0,nlon-1 do begin

;For J = 0,20 do begin
;FOR I = 0,24 do begin
    flag2 = 0U
    nodata_no2 =0
    For num = 0,m-1 do begin
	;filter the value < 0
;*******************
;;	  if no2_data[I,J,num] ge -MAX(maximum) then begin
;******************
	 if no2_data[I,J,num] gt -999.0 then begin
		if flag2 then begin
			Y = [Y,no2_data[I,J,num] ]
			x = [x,num+1]
		endif else begin
			Y = no2_data[I,J,num] 
			x = [num+1]
			flag2 = 1U
		endelse
          endif else begin
		nodata_no2 = nodata_no2+1
 	  endelse
     endfor
	
  if nodata_no2 eq m then begin
		para[I,J]=-999.0
		nodata_para=nodata_para+1
  endif else begin	

     
   ;print, I,J
   ;print,Y,x

    if N_elements(x) le 3 then begin
		para[I,J]=-999.0
		deltab[I,J]=1.0
		print,'Elements <3 -999.0',I,J
    endif else begin
	weights = make_array(N_elements(x),value=1.0) 
	x =float( x)/12
	;Provide an initial guess of the function's parameters.
	A = [0.025,3,1.5708]
;;**************filter the data without complement time series
		yfit = CURVEFIT(x,y,weights,A,SIGMA,FUNCTION_NAME='gfunct',ITMAX=50)
		para[I,J]= A[0]
		;if para[I,J] le 0.0 then begin
			;print,'I',I,'J',J,'PARA <0',para[I,J]
		;endif

		ax = sin( !DPI/6 * X + A[2] )
		F = A[0] * X + A[1] * ax
		remain = F-Y
;;	if N_elements(x) ge 12*5 then begin
		;print,'remain',remain
		;print,'x',x
		lag_cor = [1]
		cor = A_CORRELATE(remain,lag_cor,/double)
;		cor = FFT(remain,-1)
;		print,'cor',cor
;		deltan[I,J]= VARIANCE(remain)
		deltan[I,J]=STDDEV(remain)
;		n = N_elements(x)
		n=9
;		print,deltan[I,J]
		deltab[I,J] = deltan[I,J]/(n^(1.5))*(((1+cor[0])/(1-cor[0]))^(0.5))
;		print,deltab[I,J]
;;	endif else begin
;;		para[I,J]=-999.0	
;;		deltab[I,J]=1.0
		;print,'999 value','I',I,'J',J,para[I,J]
;;        endelse
    endelse

  endelse

  if ( I eq 44) and (J eq 27) then begin
	print,'x',x
	print,'y',y
	print,'B',A[0]
	print,'C', A[1]
	print,'E', A[2]
  endif
endfor
endfor

print,'nodata_para',nodata_para


;select abs(B/deltab)>2
para2=para
loc = where(abs(para/deltab) le 2)
para2[loc]= -999.0
undefine,loc
print,'para2 <0.2',N_elements(para2[where(para2 lt -0.2)])

para[where (para ne -999.0)]=para[where (para ne -999.0)] 
para2[where (para2 ne -999.0)]=para2[where (para2 ne -999.0)]

header_output = [['ncols 64'],['nrows 40'],['xllcorner 72'],['yllcorner 15'],['cellsize 1.0'],['nodata_value -999.0'],['version 2.0']]
outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit_gome_1X1.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para
close,lun

outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit2_gome_1X1.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para2
close,lun
end
