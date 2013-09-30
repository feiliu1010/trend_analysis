pro regrid_KNMI_L3
;input data
;***********sci
;nlon_g=1440
;nlat_g=720
;***********omi
nlon_g=1440*2
nlat_g=720*2
;global density
global = fltarr(nlon_g,nlat_g)

;output data,1X1
nlon = 360
nlat = 180
global_2 = fltarr(nlon,nlat)

flag = 0U
flag2 = 0U
;For year = 1996, 2012 do begin
For year = 2004, 2012 do begin
For month = 1,12 do begin
Yr4  = string(year,format='(i4.4)')
Mon2 = string(month,format='(i2.2)')
nymd = year * 10000L + month * 100L + 1 * 1L

if nymd eq 19960101 then continue
if nymd eq 19960201 then continue
if nymd eq 19960301 then continue
if nymd eq 19980101 then continue

if nymd eq 20121101 then continue
if nymd eq 20121201 then continue

;***********sci no data
if nymd eq 20020101 then continue
if nymd eq 20020201 then continue
if nymd eq 20020301 then continue
if nymd eq 20020401 then continue
if nymd eq 20020501 then continue
if nymd eq 20020601 then continue
if nymd eq 20020701 then continue
;***********omi no data
if nymd eq 20040101 then continue
if nymd eq 20040201 then continue
if nymd eq 20040301 then continue
if nymd eq 20040401 then continue
if nymd eq 20040501 then continue
if nymd eq 20040601 then continue
if nymd eq 20040701 then continue
if nymd eq 20040801 then continue
if nymd eq 20040901 then continue
header = strarr(7,1)
;if (year lt 2003) or ((year eq 2003) and (month lt 4)) then begin
;filename = '/z6/satellite/GOME/no2/KNMI_v2.0/'+Yr4+'/no2_'+Yr4+Mon2+'.grd'
;endif else begin
;filename ='/z6/satellite/SCIAMACHY/no2/KNMI_v2.0/'+Yr4+'/no2_'+Yr4+Mon2+'.grd'
;endelse
filename ='/z6/satellite/OMI/no2/KNMI_L3/v2.0/no2_'+Yr4+Mon2+'.grd'
openr,lun,filename,/get_lun
readf,lun,header,global
close,/all
free_lun,lun


xj=0
;**********sci
;For J=0,nlat_g-4 do begin
;xi=0
;For I=0,nlon_g-4 do begin
;sub = global[I:I+3,J:J+3]
;loc = where(sub ge 0.0)
;if Array_equal(loc,[-1]) then begin
;global_2[xi,xj]= -999.0
;endif else begin
;global_2[xi,xj]=mean( sub[loc])
;endelse
;undefine,sub
;undefine,loc
;xi = xi+1
;I = I+3
;endfor
;J = J+3
;xj = xj +1
;endfor

;**********omi
For J=0,nlat_g-8 do begin
xi=0
For I=0,nlon_g-8 do begin
sub = global[I:I+7,J:J+7]
loc = where(sub ge 0.0)
if Array_equal(loc,[-1]) then begin
global_2[xi,xj]= -999.0
endif else begin
global_2[xi,xj]=mean( sub[loc])
endelse
undefine,sub
undefine,loc
xi = xi+1
I = I+7
endfor
J = J+7
xj = xj +1
endfor

;print,'xi',xi,'xj',xj
print,nymd,'mean of global',mean(global[where (global ge 0.0)])
print,nymd,'number of global',N_elements(global[where (global ge 0.0)])
print,nymd,'mean of global2',mean(global_2[where (global_2 ge 0.0)] )
print,nymd,'number of global2',N_elements(global_2[where (global_2 ge 0.0)])

;if (year lt 2003) or ((year eq 2003) and (month lt 4)) then begin
;outfile = '/home/liufei/Data/Satellite/NO2/GOME/no2_'+Yr4+Mon2+'_1x1.asc'
;endif else begin
;outfile = '/home/liufei/Data/Satellite/NO2/SCIAMACHY/no2_'+Yr4+Mon2+'_1x1.asc'
;endelse
outfile = '/home/liufei/Data/Satellite/NO2/OMI/no2_'+Yr4+Mon2+'_1x1.asc'
header_output =[['ncols 360'],['nrows  180'],['xllcorner -180'],['yllcorner -90'],['cellsize 1.0'],['nodata_value -999.0']]
openw,lun,outfile,/get_lun
printf,lun,header_output
printf,lun,global_2
close,/all
free_lun,lun
endfor
endfor

end
