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



pro region_trend_fit_25x25
filter=1
sensor='SCIAMACHY'
;sensor='GOME'
;sensor='OMI'
case sensor of 
'OMI': extra=10
'SCIAMACHY': extra=0
'GOME': extra=0
endcase
ncol=12

;OMI:2004/10~2011/10
;nrow=85
;SCIAMACHY:2003/05~2011/12
;nrow=104
;GOME:1996/04~2003/04(exclude 199801)
;nrow=84
case sensor of
'OMI': nrow=85
'SCIAMACHY': nrow=104
'GOME': nrow=84
endcase


header=strarr(2,1)
month=fltarr(ncol,nrow)
if filter then begin
	inputfile ='/home/liufei/Data/Trend/filter_region_monthly_mean_'+sensor+'_0.25X0.25.asc'
endif else begin
	inputfile ='/home/liufei/Data/Trend/region_monthly_mean_'+sensor+'_0.25X0.25.asc' 
endelse
openr,lun,inputfile,/get_lun
readf,lun,header,month
close,/all
print,'no data in month', n_elements(where( month le 0))-1

pp=fltarr(ncol,nrow)
if filter then begin
	inputfile ='/home/liufei/Data/Trend/filter_region_monthly_mean_PP_'+sensor+'_0.25X0.25.asc'
endif else begin
	inputfile ='/home/liufei/Data/Trend/region_monthly_mean_PP_'+sensor+'_0.25X0.25.asc'
endelse
openr,lun,inputfile,/get_lun
readf,lun,header,pp
close,/all
free_lun,lun
print,'no data in pp',n_elements( where( pp le 0))-1

no_pp=fltarr(ncol,nrow)
if filter then begin
        inputfile ='/home/liufei/Data/Trend/filter_region_monthly_mean_no_PP_'+sensor+'_0.25X0.25.asc'
endif else begin
        inputfile ='/home/liufei/Data/Trend/region_monthly_mean_no_PP_'+sensor+'_0.25X0.25.asc'
endelse
openr,lun,inputfile,/get_lun
readf,lun,header,no_pp
close,/all
free_lun,lun
print,'no data in no_pp',n_elements( where( no_pp le 0))-1


;Compute the fit to the function we have just defined. First, define the independent and dependent variables
;define para to save fit parameters
para = fltarr(ncol,3)
;repara is the relative trend
repara= fltarr(ncol,3)
;deltab=deltan/n^(3/2)*((1+cor)/(1-cor))^1/2
;deltan is the variance in the remainder
;cor=fltarr(ncol,2)
deltan=fltarr(ncol,3)
deltab=fltarr(ncol,3)

For control=0,2 do begin
temp=fltarr(ncol,nrow)
;control=0: month, control=1: pp
    case control of 
	0: temp=month
	1: temp=pp
	2: temp=no_pp
    endcase
    
    For i=0,ncol-1 do begin
	flag=0
	nodata_num=0
	For j=0,nrow-1 do begin
	    if temp[i,j] gt 0 then begin
		if flag then begin
			x=[x,j+extra]
			y=[y,temp[i,j]]
		endif else begin
			x=[j+extra]
			y=[temp[i,j]]
			flag=1
		endelse
	    endif else begin
		nodata_num+=1
	    endelse
	endfor
	if nodata_num gt (nrow-24) then begin
		para[i,control]=-999
	endif else begin
	weights = make_array(N_elements(x),value=1.0) 
	;Provide an initial guess of the function's parameters.
	A = [0.02,2.666,4.475,2.725,3.785,2.003,4.498,6.22,5.0345,5]
	yfit = CURVEFIT(x,y,weights,A,SIGMA,FUNCTION_NAME='gfunct',ITMAX=100)
	para[i,control]= A[0]
	repara[i,control]=para[i,control]/(total(y)/n_elements(y))

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
	deltan[i,control]=STDDEV(remain)
	;n is total year
	n=N_elements(x)/12
	deltab[i,control] = deltan[i,control]/(n^(1.5))*(((1+cor[0])/(1-cor[0]))^(0.5))
	;filter deltab lt delta_min
;		delta_min = 0.65 - A[9] + 0.3*mean(y)
;		if  deltab[I,J] lt delta_min then begin
;			deltab[I,J] = delta_min
;			nodata_delta_min +=1
;		endif
	;tw is the value of  the student's t-distribution for a significance level of 0.05 and the degrees of free given for the time series
	;select abs(B/deltab)>tw
	tw = T_CVF(0.05,N_elements(x)-1)
	if ( abs(para[i,control]/deltab[i,control]) le tw) then begin
        	print,'not through tw',i,control
	endif

	plot,x,y,psym=2
	oplot,x,f,psym=0
	
	endelse
    endfor
endfor

print,'para_month',para[*,0]
print,'para_pp',para[*,1]
print,'para_no_pp',para[*,2]
print,'repara_month',repara[*,0]
print,'repara_pp',repara[*,1]
print,'repara_no_pp',repara[*,2]

end
