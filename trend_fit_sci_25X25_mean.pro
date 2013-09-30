;define x = t/12, F(x)-C = B*x+k11*sin(2Pi*x)+k12*sin(4Pi*x)+k13*sin(6Pi*x)+k14*sin(8Pi*x)+k21*cos(2Pi*x)+k22*cos(4Pi*x)+k23*cos(6Pi*x)+k24*cos(8Pi*x)
;define a procedure to return F(x) and the partial derivatives,given x. Note that A is an array containing the values B,k11,k12,k13,k14,k21,k22,k23,k24,C
;the partial derivatives are computed:
;dF/dB = x
;dF/dk11 = sin(!DPI*2*x)
;dF/dk12 = sin(!DPI*4*x)
;dF/dk13 = sin(!DPI*6*x)
;dF/dk14 = sin(!DPI*8*x)
;dF/dk21 = cos(!DPI*2*x)
;dF/dk22 = cos(!DPI*4*x)
;dF/dk23 = cos(!DPI*6*x)
;dF/dk24 = cos(!DPI*8*x)
pro gfunct,x,A,F,pder
        k11x = A[1]*sin(!DPI/6*x)
        k12x = A[2]*sin(!DPI/3*x)
        k13x = A[3]*sin(!DPI/2*x)
        k14x = A[4]*sin(!DPI*(2.0/3.0)*x)
        k21x = A[5]*cos(!DPI/6*x)
        k22x = A[6]*cos(!DPI/3*x)
        k23x = A[7]*cos(!DPI/2*x)
        k24x = A[8]*cos(!DPI*(2.0/3.0)*x)
        F =(1.0/12.0)* A[0]*x+k11x+k12x+k13x+k14x+k21x+k22x+k23x+k24x
; If the procedure is called with four parameters, calculate the partial derivatives.
        if N_PARAMS() ge 4 then $
                pder = [[(1.0/12.0)*x],[sin(!DPI/6*x)],[sin(!DPI/3*x)],[ sin(!DPI/2*x)],[sin(!DPI*(2.0/3.0)*x)],[cos(!DPI/6*x)],[cos(!DPI/3*x)],[cos(!DPI/2*x)],[cos(!DPI*(2.0/3.0)*x)]]
end



pro trend_fit_sci_25x25_mean
;input data
nlon_g=1440
nlat_g=720
;global density
global = fltarr(nlon_g,nlat_g)

nlon = 256
nlat = 160
;china density,limit=[15,72,55,136]
no2 = fltarr(nlon,nlat)
first = fltarr(nlon,nlat)
average = fltarr(nlon,nlat)
number = fltarr(nlon,nlat)
;star point of grid for China in global map
slon = (72+180)/0.25
slat = (90-55)/0.25

flag = 0U
flag2 = 0U
For year = 2002, 2011 do begin
For month = 1,12 do begin

	Yr4  = string(year,format='(i4.4)')
	Mon2 = string(month,format='(i2.2)')
	nymd = year * 10000L + month * 100L + 1 * 1L

	if nymd eq 20020101 then continue
	if nymd eq 20020201 then continue
	if nymd eq 20020301 then continue
	if nymd eq 20020401 then continue
	if nymd eq 20020501 then continue
	if nymd eq 20020601 then continue
	if nymd eq 20020701 then continue

	if nymd eq 20110901 then continue
	if nymd eq 20111001 then continue
	if nymd eq 20111101 then continue
	if nymd eq 20111201 then continue


	header = strarr(7,1)
	filename ='/z6/satellite/SCIAMACHY/no2/KNMI_v2.0/'+Yr4+'/no2_'+Yr4+Mon2+'.grd'
	openr,lun,filename,/get_lun
	readf,lun,header,global
	no2 = global[slon:slon+nlon-1,slat:slat+nlat-1]/100

	;calculate long-term avearage
	For J=0, nlat-1 do begin
	For I = 0,nlon-1 do begin
		if no2[I,J] ge 0.0 then begin
			average[I,J]+= no2[I,J]
			number[I,J] +=1
		endif
	endfor
	endfor


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

For J=0, nlat-1 do begin
        For I = 0,nlon-1 do begin
                if number[I,J] gt 0.0 then begin
                        average[I,J] = average[I,J]/number[I,J]   
;filter the data 9-year average lt 1*10^15 molecules/cm2
			if average[I,J] lt 1.0 then begin
				average[I,J] = -999.0
			endif
                endif else begin
			average[I,J] = -999.0
		endelse
        endfor
endfor


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
;;first = fltarr(nlon,nlat)
;deduce the fisrt reasonal data ( constant in fit function)
;;For J = 0,nlat-1 do begin
;;FOR I = 0,nlon-1 do begin
;;	loc = where(no2_data[I,J,*] ge 0)
;;	if (array_equal(loc,[-1])) then begin
;;		no2_data[I,J,*] = -999.0
;;	endif else begin
;;		first[I,J] = no2_data[I,J,loc[0]]
;;		no2_data[I,J,*]=no2_data[I,J,*]-first[I,J]
;;			if maximum[I,J] then begin
;;				maximum[I,J]=[maximum[I,J],first[I,J]]
;;			endif else begin
;;				maximum[I,J]=[first[I,J]]
;;			endelse
;;	endelse
;;	undefine,loc
;;endfor
;;endfor

;*************************

;regard the data of 200208
;nodata = 0
;For J = 0,nlat-1 do begin
;FOR I = 0,nlon-1 do begin
;	if no2_data[I,J,13] ge 0 then begin
;		first[I,J] = no2_data[I,J,13] 
;		no2_data[I,J,*]=no2_data[I,J,*]-first[I,J] 
;	endif else begin
;	no2_data[I,J,*]=-999.0
;	nodata=nodata+1
;	endelse
;	if ( I eq 174) and (J eq 110) then begin
;		print,'first,A[174,110]',first
;	endif
;endfor
;endfor
;*********************************

;regard the data of 2002 mean
;;nodata = 0
;;first = fltarr(nlon,nlat)
;;For J = 0,nlat-1 do begin
;;FOR I = 0,nlon-1 do begin
;;	loc = where(no2_data[I,J,0:4] ge 0.0)
;;	number_loc = n_elements(loc)
;;	if not array_equal(loc,[-1]) then begin
;;	   For n = 0,number_loc-1 do begin
;;		first[I,J]+=no2_data[I,J,loc[n]]
;;	   endfor
;;         first[I,J] = first[I,J]/(number_loc-1) 
;;	   no2_data[I,J,*]=no2_data[I,J,*]-first[I,J]
;;	endif else begin
;;           no2_data[I,J,*]=-999.0
;;           nodata=nodata+1
;;	endelse 
;;	undefine,loc
;;	undefine,number_loc 
;;endfor
;;endfor
;*********************************

;regard the data of 2002-2012 mean
nodata = 0
first = fltarr(nlon,nlat)
For J = 0,nlat-1 do begin
FOR I = 0,nlon-1 do begin
      if average[I,J] ne -999.0 then begin
         first[I,J] = average[I,J]
         no2_data[I,J,*]=no2_data[I,J,*]-first[I,J]
      endif else begin
           no2_data[I,J,*]=-999.0
           nodata=nodata+1
      endelse
endfor
endfor
;*********************************



print,'first[163,79]',first[163,79]
print,'nodata in 2002 ',nodata

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

nodata_delta_min = 0
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
;	value = -max(maximum)  
;	if no2_data[I,J,num] ge value then begin
;******************
	 if no2_data[I,J,num] gt -999.0 then begin
		if flag2 then begin
			Y = [Y,no2_data[I,J,num] ]
			x = [x,num+7]
		endif else begin
			Y = no2_data[I,J,num] 
			x = [num+7]
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

	if N_elements(x) le 9 then begin
		para[I,J]=-999.0
		deltab[I,J]=1.0
		print,'Elements <9 -999.0',I,J
	endif else begin
		weights = make_array(N_elements(x),value=1.0) 
		;Provide an initial guess of the function's parameters.
		A = [0.02,2.666,4.475,2.725,3.785,2.003,4.498,6.22,5.0345]
		yfit = CURVEFIT(x,y,weights,A,SIGMA,FUNCTION_NAME='gfunct',ITMAX=100)
		para[I,J]= A[0]
		;if para[I,J] le 0.0 then begin
			;print,'I',I,'J',J,'PARA <0',para[I,J]
		;endif

		k11x = A[1]*sin(!DPI/6*x)
	        k12x = A[2]*sin(!DPI/3*x)
	       	k13x = A[3]*sin(!DPI/2*x)
	       	k14x = A[4]*sin(!DPI*(2.0/3.0)*x)
	        k21x = A[5]*cos(!DPI/6*x)
	        k22x = A[6]*cos(!DPI/3*x)
	        k23x = A[7]*cos(!DPI/2*x)
	        k24x = A[8]*cos(!DPI*(2.0/3.0)*x)
	        F = A[0]*(1.0/12.0)*x+k11x+k12x+k13x+k14x+k21x+k22x+k23x+k24x
		remain = F-Y
		;print,'remain',remain
		lag_cor = [1]
		cor = A_CORRELATE(remain,lag_cor,/double)
		;print,'cor',cor
		deltan[I,J]=STDDEV(remain)
		;n is total year
		;n=9
		n=N_elements(x)/12
		deltab[I,J] = deltan[I,J]/(n^(1.5))*(((1+cor[0])/(1-cor[0]))^(0.5))
		;filter deltab lt delta_min
		delta_min = 0.65 - first[I,J] + 0.3*mean(y)
		if  deltab[I,J] lt delta_min then begin
			deltab[I,J] = delta_min
			nodata_delta_min +=1
		endif
	endelse

    endelse

    if ( I eq 163) and (J eq 79) then begin
;	print,'x',x
;	print,'y',y
	print,'B',A[0]
;	print,'C', A[1]
;	print,'E', A[2]
	plot,x,y+first[I,J],psym=2,$
        xrange=[0,108],yrange=[min(y+first[I,J]),max(y+first[I,J])],$
	title='fit_0.25X0.25',xtitle='number of months',ytitle='no2/(10^15moles/cm2)'
 	oplot,x,F+first[I,J],psym=0
        image = tvrd(true =1)
        write_jpeg,'no2_fit_0.25X0.25.jpg',image,true=1
   endif
endfor
endfor

print,'nodata_para',nodata_para
print,'nodata_delta_min',nodata_delta_min
;filter the data 9-year average lt 1*10^15 molecules/cm2
loc = where(average  eq -999.0)
print,'number of background',n_elements(loc)
para[loc]=-999.0
undefine,loc
;tw is the value of  the student's t-distribution for a significance level of 0.05 and the degrees of free given for the time series
;select abs(B/deltab)>2
tw = T_CVF(0.05,N_elements(x)-1)
para2=para
loc = where(abs(para/deltab) le tw)
para2[loc]= -999.0
print,'nodata_para2',N_elements(loc)
undefine,loc

loc = where( (para2 gt -999.0) and (average gt 0.0) )
para3 = fltarr(nlon,nlat)
For J=0, nlat-1 do begin
para3[*,J] = replicate(-999.0,nlon)
endfor
para3[loc]= para2[loc]/average[loc]
undefine,loc

loc = where( (para gt -999.0) and (average gt 0.0) )
para4 = fltarr(nlon,nlat)
For J=0, nlat-1 do begin
para4[*,J] = replicate(-999.0,nlon)
endfor
para4[loc]= para[loc]/average[loc]


header_output = [['ncols 256'],['nrows 160'],['xllcorner 72'],['yllcorner 15'],['cellsize 0.25'],['nodata_value -999.0']]
outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit_sci_0.25X0.25.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para
close,lun

outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit2_sci_0.25X0.25.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para2
close,lun

outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit3_sci_0.25X0.25.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para3
close,lun

outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit4_sci_0.25X0.25.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para4
close,lun
end
