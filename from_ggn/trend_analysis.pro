;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;define the function of the form F(x)=a+b*x+c*sin(d*x+e)
;a and d is known                                       
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro gfunct, X, A, F, pder

ax = sin( 0.5236 * X + A[2] )
bx = cos( 0.5236 * X + A[2] )
F = A[0] * X + A[1] * ax

;If the procedure is called with 4 parameters, calculate the partial derivatives
if n_params() ge 4 then $
  pder= [[X], [ax], [A[1] * bx]]

end



;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;Compute the fit to the function we have just defined
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro trend_analysis

FORWARD_FUNCTION CTM_Grid, CTM_Type, CTM_Get_Data, CTM_Make_DataInfo, Nymd2Tau

InType = CTM_Type( 'GENERIC', Resolution=[1d0,1d0])
InGrid = CTM_Grid( InType )

xmid = InGrid.xmid
ymid = InGrid.ymid

limit=[15,70,55,136]

i1_index = where(xmid ge limit[1] and xmid le limit[3])
j1_index = where(ymid ge limit[0] and ymid le limit[2])
I1 = min(i1_index, max = I2)
J1 = min(j1_index, max = J2)
I0 = ( 70 + 180 ) / 1
J0 = ( 15 + 90 ) / 1
print,I1,I2,J1,J2


;prepare the data needed
;filename1 = '/z3/wangsiwen/Satellite/no2/SCIAMACHY_Bremen_v1.5/scia_no2_v1-5_200301.bpch'
filename1 = '/z3/gengguannan/satellite/no2/GOME_Bremen/gome_no2_v3-0_199604_1x1.bpch'

ctm_get_data,datainfo1,filename = filename1,tracer=1
a=*(datainfo1[0].data)


xx = I2 - I1 + 1
yy = J2 - J1 + 1

no2 = fltarr(xx,yy,177)
data = fltarr(xx,yy,177)
B = fltarr(xx,yy)
C = fltarr(xx,yy)
E = fltarr(xx,yy)
count = fltarr(xx,yy)
status = fltarr(xx,yy)


for y = 1996,2010 do begin
for m = 1,12 do begin

Yr4  = String( y, format = '(i4.4)' )
Mon2 = String( m, format = '(i2.2)' )

nymd = y * 10000L + m * 100L + 1 * 1L
Tau0 = nymd2tau(NYMD)
print,nymd


if nymd eq 19960101 then continue
if nymd eq 19960201 then continue
if nymd eq 19960301 then continue
if nymd eq 19960401 then continue


if y lt 2003 $
  then filename2 = '/z3/gengguannan/satellite/no2/GOME_Bremen/gome_no2_v3-0_'+ Yr4 + Mon2 +'_1x1.bpch' $
  else filename2 = '/z3/gengguannan/satellite/no2/SCIAMACHY_Bremen_v1.5/scia_no2_v1-5_'+ Yr4 + Mon2 +'_1x1.bpch'


ctm_get_data,datainfo2,filename = filename2,tau0=nymd2tau(NYMD),tracer=1
data18=*(datainfo2[0].data)

k = m+(y-1996)*12-5

for I = I1,I2 do begin
  for J = J1,J2 do begin
    no2[I-I0,J-J0,k] = data18[I,J]
    data[I-I0,J-J0,k] = data18[I,J]-a[I,J]
  endfor
endfor

CTM_CLEANUP

endfor
endfor


;fit
for I = 0,xx-1 do begin
  for J = 0,yy-1 do begin
    X0 = indgen(177)+1
    X = make_array(1)
    Y = make_array(1)
    Y0 = make_array(1)
    flag = 1
    for p = 0,n_elements(X0)-1 do begin
      if no2[I,J,p] gt 0 then begin
        if flag eq 1 then begin
          Y = [data[I,J,p]]
          Y0 = [no2[I,J,p]]
          X = [X0[p]]
        endif else begin
          Y = [Y,data[I,J,p]]
          Y0 = [Y0,no2[I,J,p]]
          X = [X,X0[p]]
        endelse
      endif
      flag = 0
    endfor

    print,n_elements(X)

    weights = make_array(177,value=1.0)

    ; Provide an initial guess of the function's parameters
    A = [0.03,3,1.5708]

    ; Compute the parameters.
    yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct', ITMAX=50)

    ; Print the parameters returned in A.
    PRINT, 'Function parameters: ', A
;    PRINT, STATUS

    B[I,J]=A[0]
    C[I,J]=A[1]
    E[I,J]=A[2]
    count[I,J]=n_elements(Y)
;    status[I,J]=STATUS

    CTM_CLEANUP

  endfor
endfor

outfile1 = '/home/gengguannan/result/trend_analysis_result_1x1.hdf'
outfile2 = '/home/gengguannan/result/trend_analysis_result_1x1.bpch'

IF (HDF_EXISTS() eq 0) then message, 'HDF not supported'

;open the HDF file
FID = HDF_SD_Start(Outfile1,/RDWR,/Create)

HDF_SETSD, FID, B, 'B',             $
           Longname='monthly trend',$
           Unit='unitless',         $
           FILL=-999.0
HDF_SETSD, FID, C, 'C',             $
           Longname='amplitude',    $
           Unit='unitless',         $
           FILL=-999.0
HDF_SETSD, FID, E, 'E',             $
           Longname='phase shift',  $
           Unit='unitless',         $
           FILL=-999.0
HDF_SETSD, FID, count, 'count',     $
           Longname='count',        $
           Unit='unitless',         $
           FILL=-999.0
HDF_SETSD, FID, status, 'status',   $
           Longname='status',       $
           Unit='unitless',         $
           FILL=-999.0

HDF_SD_End, FID


  ; Make a DATAINFO structure
   success = CTM_Make_DataInfo( B,                       $
                                ThisDataInfo,            $
                                ThisFileInfo,            $
                                ModelInfo=InType,        $
                                GridInfo=InGrid,         $
                                DiagN='IJ-AVG-$',        $
                                Tracer=1,                $
                                Tau0= nymd2tau(19960101),$
                                Unit='unitless',         $
                                Dim=[xx, yy, 0, 0],      $
                                First=[1L, 1L, 1L],      $
                                /No_vertical )

   CTM_WriteBpch, ThisDataInfo, FileName = OutFile2

end
