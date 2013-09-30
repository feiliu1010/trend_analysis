pro region_monthly_data
;global density
nlon_g=1440
nlat_g=720
global = fltarr(nlon_g,nlat_g)
;china density,limit=[10,70,60,150]
nlon = 320
nlat = 200
no2 = fltarr(nlon,nlat)
average = fltarr(nlon,nlat)
number = fltarr(nlon,nlat)
;star point of grid for China in global map
cellsize=0.25
slon = (70+180)/cellsize
slat = (90-60)/cellsize
;region stands for mask of 3 qu 9 qun
header=strarr(6,1)
region=fltarr(nlon,nlat)
filename_region = '/home/liufei/Data/Trend/Input/regions_25.asc'
openr,lun,filename_region,/get_lun
readf,lun,header,region
close,/all
;no is NO. of regions
;region_num is the total number of target regions
no = region(UNIQ(region, SORT(region)))
region_num=N_elements(no)
print,'NO.',no
print,'region_num',region_num
;pp stands for mask of pp control emissions grid
pp=fltarr(nlon,nlat)
filename_pp ='/home/liufei/Data/Trend/PP_Control_2010_annual_0.25.asc'
openr,lun,filename_pp,/get_lun
readf,lun,header,pp
close,/all
free_lun,lun

flag = 0U
;filter: 1=use filter data;0 = no use filter data
filter =1
;don't forget to modify time period as well
sensor='SCIAMACHY'
;sensor='GOME'
;sensor='OMI'
;sig stands for grid with significant change
sig=fltarr(nlon,nlat)
if filter then begin
	filename_sig='/home/liufei/Data/Trend/filter_relative_trend_'+sensor+'_0.25X0.25.asc'
endif else begin
	filename_sig='/home/liufei/Data/Trend/relative_trend_'+sensor+'_0.25X0.25.asc'
endelse
openr,lun,filename_sig,/get_lun
readf,lun,header,sig
close,/all
free_lun,lun
;GOME(1996/04-2003/06):1996/04-2003/04
;For year = 1996, 2003 do begin
;SCIAMACHY(2002/08-2012/03):2003/05-2012/03
For year = 2003, 2011 do begin
;OMI(2004/10~)
;For year = 2004, 2011 do begin
    For month = 1,12 do begin

	Yr4  = string(year,format='(i4.4)')
	Mon2 = string(month,format='(i2.2)')
	nymd = year * 10000L + month * 100L + 1 * 1L

	if STRCMP(sensor,'SCIAMACHY') then begin 
		 if total(nymd eq [20030101,20030201,20030301,20030401]) eq 1.0 then continue
	endif else if STRCMP(sensor,'GOME') then  begin
		if total(nymd eq [19960101,19960201,19960301,19980101,$
			20030501,20030601,20030701,20030801,20030901,20031001,20031101,20031201]) eq 1.0 then continue
	endif else begin
		if total(nymd eq [20040101,20040201,20040301,20040401,20040501,20040601,20040701,20040801,20040901,$
			20111101,20111201]) eq 1.0 then continue
	endelse


	header = strarr(7,1)
	if STRCMP(sensor,'SCIAMACHY') then begin
		filename ='/z6/satellite/SCIAMACHY/no2/KNMI_v2.0/'+Yr4+'/no2_'+Yr4+Mon2+'.grd'
	endif else if  STRCMP(sensor,'GOME') then begin
		filename='/z6/satellite/GOME/no2/KNMI_v2.0/'+Yr4+'/no2_'+Yr4+Mon2+'.grd'
	endif else begin
		filename='/home/liufei/Data/Satellite/NO2/OMI/no2_'+Yr4+Mon2+'_0.25x0.25.asc'
		header = strarr(6,1)
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
no2_data[loc] = -999.0
undefine,loc
print,'size of no2_data',size(no2_data)
print,'MAX OF no2_data',MAX(no2_data[*,*,*]),'MIN OF no2_data',MIN(no2_data[*,*,*])

;calculate average of region with significant change
region_month=fltarr(region_num-1,m)
region_pp=fltarr(region_num-1,m)
region_no_pp=fltarr(region_num-1,m)
For num = 0,m-1 do begin
	For i= 1,region_num-1 do begin
		temp_data=fltarr (nlon,nlat)
		temp_data=no2_data[*,*,num]
		temp_grid=where( (temp_data gt 0) and (region eq no[i]) and (sig ne -999.0) )
		temp_num= N_elements(temp_grid)
		if array_equal(temp_grid,[-1]) then begin
			print,'no valid data in region',i,'day number',num
			region_month[i-1,num]= -999
		endif else begin
			region_month[i-1,num]= total(temp_data[temp_grid])/temp_num
		endelse
	
		temp_pp= where( (temp_data gt 0) and (region eq no[i]) and (sig gt -999.0) and (pp gt 0) )
		temp_num_pp=N_elements(temp_pp)
		if array_equal(temp_pp,[-1]) then begin
			region_pp[i-1,num]=-999
			print,'no_pp',i,num
		endif else begin
		region_pp[i-1,num]= total(temp_data[temp_pp])/temp_num_pp
		endelse
		
		temp_no_pp= where( (temp_data gt 0) and (region eq no[i]) and (sig gt -999.0) and (pp le 0) )
                temp_num_no_pp=N_elements(temp_no_pp)
                if array_equal(temp_no_pp,[-1]) then begin
                        region_no_pp[i-1,num]=-999
                        print,'all_pp',i,num
                endif else begin
                region_no_pp[i-1,num]= total(temp_data[temp_no_pp])/temp_num_no_pp
                endelse
		
		undefine,temp_data
		undefine,temp_grid
		undefine,temp_num
		undefine,temp_pp
		undefine,temp_num_pp
		undefine,temp_no_pp
                undefine,temp_num_no_pp
	endfor
endfor

header_output = [['ncols'+string(region_num-1)],['nrows'+string(m)]]
if filter then begin
	outfile = '/home/liufei/Data/Trend/filter_region_monthly_mean_'+sensor+'_0.25X0.25.asc'
endif else begin
	outfile = '/home/liufei/Data/Trend/region_monthly_mean_'+sensor+'_0.25X0.25.asc'
endelse
openw,lun,outfile,/get_lun
printf,lun,header_output,region_month
close,lun

if filter then begin
	outfile = '/home/liufei/Data/Trend/filter_region_monthly_mean_PP_'+sensor+'_0.25X0.25.asc'
endif else begin
	outfile = '/home/liufei/Data/Trend/region_monthly_mean_PP_'+sensor+'_0.25X0.25.asc'
endelse
openw,lun,outfile,/get_lun
printf,lun,header_output,region_pp
close,lun

if filter then begin
        outfile = '/home/liufei/Data/Trend/filter_region_monthly_mean_no_PP_'+sensor+'_0.25X0.25.asc'
endif else begin
        outfile = '/home/liufei/Data/Trend/region_monthly_mean_no_PP_'+sensor+'_0.25X0.25.asc'
endelse
openw,lun,outfile,/get_lun
printf,lun,header_output,region_no_pp
close,lun


end
