pro power_control_grid
	
	;limit=[10,70,60,150]
	nlon = 320
	nlat = 200
	;control_pp is grid with pp emission accounts for >0.6
	control_pp=fltarr(nlon,nlat)
	;nox is 2010 total NO2 emssion from MEIC
	nox=fltarr(nlon,nlat)
	nox_power= fltarr(nlon,nlat)

	season='annual'
	year = 2010
	header = strarr(6,1)
	Yr4  = string(year,format='(i4.4)')

	;Emissions from power plants account for >60% emssions in 2010
        filename_all = '/home/liufei/Data/Decline/Input/2010__all__NOx_25.asc'
        header2 = strarr(6,1)
        openr,lun,filename_all,/get_lun
        readf,lun,header2,nox
        close,/all
        free_lun,lun

	filename_power = '/home/liufei/Data/Decline/Input/2010__power__NOx_25.asc'
	header2 = strarr(6,1)
	openr,lun,filename_power,/get_lun
        readf,lun,header2,nox_power
	close,/all
	free_lun,lun

	For I=0,nlon-1 do begin
		For J=0,nlat-1 do begin
			if (nox_power[I,J] gt 0) and (nox_power[I,J]/nox[I,J] gt 0.6) then begin
				control_pp[I,J]=nox_power[I,J]
			endif else begin
				control_pp[I,J]=-999.0
			endelse
		endfor
	endfor			
	
	header_output = [['ncols '+string(nlon)],['nrows '+string(nlat)],['xllcorner 70'],['yllcorner 10'],['cellsize 0.25'],['nodata_value -999.0']]
	outfile='/home/liufei/Data/Trend/PP_Control_'+Yr4+'_'+Season+'_0.25.asc'
	openw,lun,outfile,/get_lun
	printf,lun,header_output,control_pp
	close,/all
	free_lun,lun

	
	;Grids with power plants >=600 MW in 2011
	;filename_LOC = '/home/liufei/Data/Decline/Input/pp_600_loc.asc'
        ;header2 = strarr(6,1)
        ;openr,lun,filename_LOC,/get_lun
        ;readf,lun,header2,PP_Loc
        ;close,/all
        ;free_lun,lun
	
	end	
