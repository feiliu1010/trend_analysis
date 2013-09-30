pro trend_analysis_filter1

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
 

max_season = fltarr(InGrid.IMX,InGrid.JMX)
max_mask = fltarr(InGrid.IMX,InGrid.JMX)


for y = 1997,2010 do begin

no2 = fltarr(InGrid.IMX,InGrid.JMX,4)
max = fltarr(InGrid.IMX,InGrid.JMX)

for k = 0,3 do begin

Yr4  = String( y, format = '(i4.4)' )

nymd = y * 10000L + 1 * 100L + 1 * 1L
Tau0 = nymd2tau(NYMD)
print,nymd

season = ['JJA','SON','DJF','MAM']

if y lt 2003 $
  then filename = '/z3/gengguannan/satellite/no2/GOME_Bremen/gome_no2_v3-0_'+ Yr4 +'_'+ season[k] +'.1x1.bpch' $
  else filename = '/z3/gengguannan/satellite/no2/SCIAMACHY_Bremen_v1.5/scia_no2_v1-5_'+ Yr4 +'_'+ season[k] +'.1x1.bpch'

ctm_get_data,datainfo,filename = filename,tau0=nymd2tau(NYMD),tracer=1
data18=*(datainfo[0].data)


for I = I1,I2 do begin
  for J = J1,J2 do begin
    no2[I,J,k] = data18[I,J]
  endfor
endfor

endfor

CTM_Cleanup

;max season
for I = I1,I2 do begin
  for J = J1,J2 do begin
    max[I,J] = where(no2[I,J,*] eq max(no2[I,J,*]))
  endfor
endfor

for I = I1,I2 do begin
  for J = J1,J2 do begin
    if (max[I,J] gt 0) then max_season[I,J] = max_season[I,J] + 1
  endfor
endfor

endfor

for I = I1,I2 do begin
  for J = J1,J2 do begin
    if max_season[I,J] eq max(max_season) then max_mask[I,J] = 1
  endfor
endfor

print,total(max_mask)

OutFile = '/home/gengguannan/result/trend_analysis_mask1_1x1.bpch'

  ; Make a DATAINFO structure
   success = CTM_Make_DataInfo( max_mask,                $
                                ThisDataInfo,            $
                                ThisFileInfo,            $
                                ModelInfo=InType,        $
                                GridInfo=InGrid,         $
                                DiagN='LANDMAP',         $
                                Tracer=802,              $
                                Tau0= nymd2tau(19850101),$
                                Unit='unitless',         $
                                Dim=[InGrid.IMX,InGrid.JMX, 0, 0],      $
                                First=[1L, 1L, 1L],      $
                                /No_vertical )

   CTM_WriteBpch, ThisDataInfo, FileName = OutFile



end
