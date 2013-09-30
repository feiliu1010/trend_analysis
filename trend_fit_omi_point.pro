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
;dF/dC = 1
pro gfunct,x,A,F,pder
        k11x = A[1]*sin(!DPI/6*x)
        k12x = A[2]*sin(!DPI/3*x)
        k13x = A[3]*sin(!DPI/2*x)
        k14x = A[4]*sin(!DPI*(2.0/3.0)*x)
        k21x = A[5]*cos(!DPI/6*x)
        k22x = A[6]*cos(!DPI/3*x)
        k23x = A[7]*cos(!DPI/2*x)
        k24x = A[8]*cos(!DPI*(2.0/3.0)*x)
        F =(1.0/12.0)* A[0]*x+k11x+k12x+k13x+k14x+k21x+k22x+k23x+k24x+A[9]
; If the procedure is called with four parameters, calculate the partial derivatives.
        if N_PARAMS() ge 4 then $
                pder = [[(1.0/12.0)*x],[sin(!DPI/6*x)],[sin(!DPI/3*x)],[ sin(!DPI/2*x)],[sin(!DPI*(2.0/3.0)*x)],[cos(!DPI/6*x)],[cos(!DPI/3*x)],[cos(!DPI/2*x)],[cos(!DPI*(2.0/3.0)*x)],[replicate(1.0,N_ELEMENTS(x))]]
end



pro trend_fit_omi_point
;input data
nlon_g=2880
nlat_g=1440
;global density
global = fltarr(nlon_g,nlat_g)

nlon = 256*2
nlat = 160*2
;china density,limit=[15,72,55,136]
no2 = fltarr(nlon,nlat)
average = fltarr(nlon,nlat)
number = fltarr(nlon,nlat)
;star point of grid for China in global map
slon = (72+180)/0.125
slat = (90-55)/0.125

flag = 0U
flag2 = 0U
For year = 2004, 2012 do begin
For month = 1,12 do begin

	Yr4  = string(year,format='(i4.4)')
	Mon2 = string(month,format='(i2.2)')
	nymd = year * 10000L + month * 100L + 1 * 1L
        
	if nymd eq 20040101 then continue
        if nymd eq 20040201 then continue
        if nymd eq 20040301 then continue
        if nymd eq 20040401 then continue
        if nymd eq 20040501 then continue
        if nymd eq 20040601 then continue
        if nymd eq 20040701 then continue
        if nymd eq 20040801 then continue
        if nymd eq 20040901 then continue

	if nymd eq 20121101 then continue
	if nymd eq 20121201 then continue


	header = strarr(7,1)
	filename='/z6/satellite/OMI/no2/KNMI_L3/v2.0/no2_'+Yr4+Mon2+'.grd'
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
print,'MAX OF no2_data',MAX(no2_data[*,*,*]),'MIN OF no2_data',MIN(no2_data[*,*,*])


;scenario_info=[scenario no.,absolute trend, background value, relative trend]
scenario_info=fltarr(13,4)
For scenario= 1,13 do begin

;Compute the fit to the function we have just defined. First, define the independent and dependent variables
;define para to save fit parameters
para = fltarr(nlon,nlat)
;deltab=deltan/n^(3/2)*((1+cor)/(1-cor))^1/2
;deltan is the variance in the remainder
cor=fltarr(nlon,nlat)
deltan=fltarr(nlon,nlat)
deltab=fltarr(nlon,nlat)

no2_point=fltarr(nlon,nlat,m)
no2_point_num=fltarr(nlon,nlat,m)
no2_point_average=fltarr(nlon,nlat)
;power plant:[112.5738889,35.46555556]
I = 325
J = 156
point_I=I
point_J=J
case scenario of 
;0.125x0.125
1: I=I
;0.25x0.25_left_down
2: I=I-1
;0.25x0.25_left_up
3: I=I-1
;0.25x0.25_right_down
4: I=I
;0.25x0.25_right_up
5: I=I
;0.5x0.5_left_down
6: I=I-2
;0.5x0.5_left_up
7: I=I-2
;0.5x0.5_righdown
8: I=I-1
;0.5x0.5_right_up
9: I=I-1
;1x1_left_down
10:I=I-5
;1x1_left_up
11:I=I-5
;1x1_right_down
12:I=I-2
;1x1_right_up
13:I=I-2
else:print,'wrong scenario'
endcase

case scenario of
1: J=J
2: J=J
3: J=J-1
4: J=J
5: J=J-1
6: J=J-2
7: J=J-3
8: J=J-2
9: J=J-3
10:J=J-6
11:J=J-7
12:J=J-6
13:J=J-7
else:print,'wrong scenario'
endcase

case scenario of
1: inter =0
2: inter =1
3: inter =1
4: inter =1
5: inter =1
6: inter =3
7: inter =3
8: inter =3
9: inter =3
10:inter =7
11:inter =7
12:inter =7
13:inter =7
else:print,'wrong scenario'
endcase


For num = 0,m-1 do begin
	For xdim = I, I+inter do begin
		For ydim = J, J+inter do begin
			if no2_data[xdim,ydim,num] gt -999.0 then begin
				no2_point[point_I,point_J,num]+=no2_data[xdim,ydim,num]
				no2_point_num[point_I,point_J,num]+=1
			endif
		endfor
	endfor
endfor
print,'no2_point[point_I,point_J,num]',no2_point[point_I,point_J,num-1]
print,'no2_point_num[point_I,point_J,num]',no2_point_num[point_I,point_J,num-1]

For  num = 0,m-1 do begin
	 no2_point[point_I,point_J,num]= no2_point[point_I,point_J,num]/ no2_point_num[point_I,point_J,num]
endfor

loc= where( no2_point[point_I,point_J,*] gt -999.0)
no2_point_average[point_I,point_J] = mean( no2_point[point_I,point_J,loc])
undefine,loc

flag2 = 0U
For num = 0,m-1 do begin
	 if no2_point[point_I,point_J,num] gt -999.0 then begin
		if flag2 then begin
			Y = [Y,no2_point[point_I,point_J,num] ]
			x = [x,num+8]
		endif else begin
			Y = no2_point[point_I,point_J,num] 
			x = [num+8]
			flag2 = 1U
		endelse
          endif 
endfor
	

if N_elements(x) le 10 then begin
		para[point_I,point_J]=-999.0
		deltab[point_I,point_J]=1.0
		print,'Elements <10 -999.0',point_I,point_J
endif else begin
		weights = make_array(N_elements(x),value=1.0) 
		;Provide an initial guess of the function's parameters.
		A = [0.02,2.666,4.475,2.725,3.785,2.003,4.498,6.22,5.0345,5]
		yfit = CURVEFIT(x,y,weights,A,SIGMA,FUNCTION_NAME='gfunct',ITMAX=100)
		para[point_I,point_J]= A[0]

		k11x = A[1]*sin(!DPI/6*x)
	        k12x = A[2]*sin(!DPI/3*x)
	       	k13x = A[3]*sin(!DPI/2*x)
	       	k14x = A[4]*sin(!DPI*(2.0/3.0)*x)
	        k21x = A[5]*cos(!DPI/6*x)
	        k22x = A[6]*cos(!DPI/3*x)
	        k23x = A[7]*cos(!DPI/2*x)
	        k24x = A[8]*cos(!DPI*(2.0/3.0)*x)
	        F = A[0]*(1.0/12.0)*x+k11x+k12x+k13x+k14x+k21x+k22x+k23x+k24x+A[9]
		remain = F-Y
		;print,'remain',remain
		lag_cor = [1]
		cor = A_CORRELATE(remain,lag_cor,/double)
		;print,'cor',cor
		deltan[point_I,point_J]=STDDEV(remain)
		;n is total year
		;n=9
		n=N_elements(x)/12
		deltab[point_I,point_J] = deltan[point_I,point_J]/(n^(1.5))*(((1+cor[0])/(1-cor[0]))^(0.5))
		;filter deltab lt delta_min
;		delta_min = 0.65 - A[9] + 0.3*mean(y)
;		if  deltab[I,J] lt delta_min then begin
;			deltab[I,J] = delta_min
;			nodata_delta_min +=1
;		endif
	
endelse


	print,'scenario',scenario
	print,'B',A[0]
	print,'background'+string(point_I)+string(point_J),no2_point_average[point_I,point_J]
	print,'trend',A[0]/no2_point_average[point_I,point_J]
	plot,x,y,psym=2,$
        yrange=[0,max(y)],$
	title='fit_scenario_'+string(scenario),xtitle='number of months',ytitle='no2/(10^15moles/cm2)'
 	oplot,x,F,psym=0
        image = tvrd(true =1)
        write_jpeg,'no2_fit_scenario'+string(scenario)+'.jpg',image,true=1
	scenario_info[scenario-1,*]=[scenario,A[0],no2_point_average[point_I,point_J],A[0]/no2_point_average[point_I,point_J]]

endfor

header_output=['scenario_NO.','absolute trend', 'background value', 'relative trend']
outfile = '/home/liufei/Data/Satellite/NO2/trend/scenario_info.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output, scenario_info
close,lun


end
