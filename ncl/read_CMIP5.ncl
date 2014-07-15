load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

;---Directory
indir = asciiread("indir_regrid&read.dat",-1,"string")
outdir = "maps/"
outdir2 = asciiread("indir_reanalysis.dat",-1,"string")
system("if ! test -d " + indir+outdir + "; then mkdir " + indir+outdir +"; fi")

;---Models
;---CMIP5:
   model = "BCC"
;   model = "CCSM4"
;   model = "GISS_121"
;   model = "GISS_122"  
;   model = "GISS_123"
;   model = "GISS_124"
;   model = "GISS_125"
;   model = "GISS_126"
;   model = "GISS_127"
;   model = "GISS_128"
;   model = "GISS_1221"
;   model = "IPSL" ; don't use, 'cause file for 105001-124912 time period is missing
;   model = "MIROC"
;   model = "MPI"
;   model = "MRI" 
;---PMIP3
;   model = "CSIRO"
;   model = "HadCM3"

;---Input files: NetCDF surface runoff
var = "evspsbl"
if (var .eq. "mrros")
 varlong = "surface runoff"
end if
if (var .eq. "evspsbl")
 varlong = "evaporation"
end if
if (var .eq. "pr")
 varlong = "precipitation"
end if
cmd = "find "+indir+"regrid/"+model+"_reg_"+var+"* -type f -print"
f = systemfunc(cmd)
;print(f)
files = addfiles(f, "r")
ListSetType(files,"cat")
fnc_name = model+"_reg_"+var+"_Amon_bcc-csm1-1_past1000_r1i1p1_085001-185012.nc"
fnc = addfile(indir+"regrid/"+fnc_name,"r")

;---Input files: watershed coordinates
fwat1_name = "basin_515.dat" ; robust
fwat2_name = "basin_899.dat" ; fine
nps = 899 ; number of points
if (nps .eq. 515) then
 wat = "Volgabasin_robust"
else
 wat = "Volgabasin_fine"
end if
watershed = asciiread(fwat2_name,(/nps,3/),"float")
nlat = watershed(:,2)
nlon = watershed(:,1)

;---Input files: weights
fwgt1_name = "weights_Volga_polygon.dat"
fwgt2_name = "weights_Volga_rect.dat"
fwgt3_name = "weights_CaspianSea_rect.dat"
nmps1 = 50 ; number of mesh points
nmps2 = 135
nmps3 = 35
weights1 = asciiread(fwgt1_name,(/nmps1,3/),"float")
weights2 = asciiread(fwgt2_name,(/nmps2/),"float")
weights3 = asciiread(fwgt3_name,(/nmps3/),"float")

;---Variables
time = files[:]->time
lat = fnc->lat
lon = fnc->lon
evspsbl = files[:]->evspsbl*2592000 ; kg/m2*c(in month) -> [kg/m2]/month = mm/month

evspsbl!0 = "time"
evspsbl!1 = "lat"
evspsbl!2 = "lon"
evspsbl&lat = lat
evspsbl&lon = lon
evspsbl&lat@units = "degrees_north"
evspsbl&lon@units = "degrees_east"

evspsbl_yr = month_to_annual(evspsbl,0)
;printVarSummary(evspsbl_yr)

;---Map
maparea = 2 ; 1 - Volgabasin, 2 - CaspianSea
title = model+"_"+var+"_"+maparea+"_m_"+wat
wkstype = "png"
wkstype@wkWidth = 2000
wkstype@wkHeight = 1800
wks1 = gsn_open_wks(wkstype,indir+outdir+title)
gsn_define_colormap(wks1,"WhiteBlue")

evspsblres = True
evspsblres@gsnDraw = False
evspsblres@gsnFrame = False
mpres = evspsblres

evspsblres@gsnAddCyclic = False
evspsblres@cnFillOn = True
evspsblres@cnLinesOn = False
evspsblres@cnLineColor = "Black"
evspsblres@lbOrientation = "Vertical"

evspsblres@mpLimitMode = "LatLon"
if (maparea .eq. 1) then
 evspsblres@mpMinLatF = 40.
 evspsblres@mpMaxLatF = 65.
 evspsblres@mpMinLonF = 30.
 evspsblres@mpMaxLonF = 65.
 evspsblres@tiMainFontHeightF = 0.02
else
 evspsblres@mpMinLatF = 35.
 evspsblres@mpMaxLatF = 50.
 evspsblres@mpMinLonF = 45.
 evspsblres@mpMaxLonF = 57.
 evspsblres@tiMainFontHeightF = 0.01
end if
evspsblres@mpGridAndLimbDrawOrder = "Predraw"
evspsblres@tiMainString = model+" "+varlong
evspsblres@gsnLeftString = ""
evspsblres@gsnRightString = "mm/month"
evspsblres@mpOutlineBoundarySets = "Geophysical"
;mpres@mpDataBaseVersion = "HighRes"

evspsblid = gsn_csm_contour_map(wks1,evspsbl(7803,:,:),evspsblres) ; March 1500 year
;mpid = gsn_csm_map(wks1,mpres)

;---Polygon: draw a watershed on the map
dum = new(nps,graphic)
do i = 0,nps-2
 watres = True
 watres@gsLineColor = "black"
 watres@gsLineThicknessF = 3
 dum(i) = gsn_add_polyline(wks1,evspsblid,nlon(i:i+1),nlat(i:i+1),watres)
end do

;---Polygon: rewrite lat, lon into 1D array
plat = new(dimsizes(lon)*dimsizes(lat),double)
plon = new(dimsizes(lon)*dimsizes(lat),double)
LON = new(dimsizes(lon)*dimsizes(lat),integer)
LAT = new(dimsizes(lon)*dimsizes(lat),integer)
k = 0
do i = 0, dimsizes(lon)-1
 do j = 0, dimsizes(lat)-1
  LON(k) = i
  LAT(k) = j
  plat(k) = lat(j) ; 1D array
  plon(k) = lon(i) ; 1D array
  k = k+1
 end do
end do

;---Polygon: find indices of mesh points (j_platin,i_plonin) inside polygon
inout = gc_inout(plat,plon,nlat,nlon)
inoutres = True
inoutres@gsMarkerIndex = 2
inoutres@gsMarkerColor = "black"
inoutres@gsMarkerSizeF = 0.007
inoutres@gsMarkerThicknessF = 2
platin = where(inout,plat,plat@_FillValue)
plonin = where(inout,plon,plon@_FillValue)
platout = where(inout,plat@_FillValue,plat)
plonout = where(inout,plon@_FillValue,plon)

j_platin = ind(.not.ismissing(platin))
i_plonin = ind(.not.ismissing(plonin))
j_platout = ind(.not.ismissing(platout))
j_plonout = ind(.not.ismissing(plonout))
if (maparea .eq. 1) then
 inoutiin = gsn_add_polymarker(wks1,evspsblid,plonin,platin,inoutres)
 inoutidout1 = gsn_add_polymarker(wks1,evspsblid,plonout(1194:1202),platout(1194:1202),inoutres)
 inoutidout2 = gsn_add_polymarker(wks1,evspsblid,plonout(1265:1269),platout(1265:1269),inoutres)
 inoutidout3 = gsn_add_polymarker(wks1,evspsblid,plonout(1271:1273),platout(1271:1273),inoutres)
 inoutidout4 = gsn_add_polymarker(wks1,evspsblid,plonout(1336:1339),platout(1336:1339),inoutres)
 inoutidout5 = gsn_add_polymarker(wks1,evspsblid,plonout(1344),platout(1344),inoutres)
 inoutidout6 = gsn_add_polymarker(wks1,evspsblid,plonout(1407:1410),platout(1407:1410),inoutres)
 inoutidout7 = gsn_add_polymarker(wks1,evspsblid,plonout(1415),platout(1415),inoutres)
 inoutidout8 = gsn_add_polymarker(wks1,evspsblid,plonout(1478:1481),platout(1478:1481),inoutres)
 inoutidout9 = gsn_add_polymarker(wks1,evspsblid,plonout(1485:1486),platout(1485:1486),inoutres)
 inoutidout10 = gsn_add_polymarker(wks1,evspsblid,plonout(1549:1552),platout(1549:1552),inoutres)
 inoutidout11 = gsn_add_polymarker(wks1,evspsblid,plonout(1556:1557),platout(1556:1557),inoutres)
 inoutidout12 = gsn_add_polymarker(wks1,evspsblid,plonout(1620:1623),platout(1620:1623),inoutres)
 inoutidout13 = gsn_add_polymarker(wks1,evspsblid,plonout(1627:1628),platout(1627:1628),inoutres)
 inoutidout14 = gsn_add_polymarker(wks1,evspsblid,plonout(1691),platout(1691),inoutres)
 inoutidout15 = gsn_add_polymarker(wks1,evspsblid,plonout(1693),platout(1693),inoutres)
 inoutidout16 = gsn_add_polymarker(wks1,evspsblid,plonout(1698:1699),platout(1698:1699),inoutres)
 inoutidout17 = gsn_add_polymarker(wks1,evspsblid,plonout(1763:1764),platout(1763:1764),inoutres)
 inoutidout18 = gsn_add_polymarker(wks1,evspsblid,plonout(1769:1770),platout(1769:1770),inoutres)
 inoutidout19 = gsn_add_polymarker(wks1,evspsblid,plonout(1833:1835),platout(1833:1835),inoutres)
 inoutidout20 = gsn_add_polymarker(wks1,evspsblid,plonout(1840:1841),platout(1840:1841),inoutres)
 inoutidout21 = gsn_add_polymarker(wks1,evspsblid,plonout(1904:1907),platout(1904:1907),inoutres)
 inoutidout22 = gsn_add_polymarker(wks1,evspsblid,plonout(1912),platout(1912),inoutres)
 inoutidout23 = gsn_add_polymarker(wks1,evspsblid,plonout(1975:1977),platout(1975:1977),inoutres)
 inoutidout24 = gsn_add_polymarker(wks1,evspsblid,plonout(1983),platout(1983),inoutres)
 inoutidout25 = gsn_add_polymarker(wks1,evspsblid,plonout(2046:2049),platout(2046:2049),inoutres)
 inoutidout26 = gsn_add_polymarker(wks1,evspsblid,plonout(2054),platout(2054),inoutres)
 inoutidout27 = gsn_add_polymarker(wks1,evspsblid,plonout(2117:2120),platout(2117:2120),inoutres)
 inoutidout28 = gsn_add_polymarker(wks1,evspsblid,plonout(2125),platout(2125),inoutres)
 inoutidout29 = gsn_add_polymarker(wks1,evspsblid,plonout(2188:2192),platout(2188:2192),inoutres)
 inoutidout30 = gsn_add_polymarker(wks1,evspsblid,plonout(2194:2196),platout(2194:2196),inoutres)
else
 inoutidin = gsn_add_polymarker(wks1,evspsblid,plonin,platin,inoutres)
 inoutidin1 = gsn_add_polymarker(wks1,evspsblid,plonin(1692),platin(1692),inoutres)
 inoutidout2 = gsn_add_polymarker(wks1,evspsblid,plonout(1686:1691),platout(1686:1691),inoutres)
 inoutidout3 = gsn_add_polymarker(wks1,evspsblid,plonout(1763),platout(1763),inoutres)
 inoutidin4 = gsn_add_polymarker(wks1,evspsblid,plonin(1762),platin(1762),inoutres)
 inoutidout5 = gsn_add_polymarker(wks1,evspsblid,plonout(1757:1761),platout(1757:1761),inoutres)
 inoutidout6 = gsn_add_polymarker(wks1,evspsblid,plonout(1828:1834),platout(1828:1834),inoutres)
 inoutidout7 = gsn_add_polymarker(wks1,evspsblid,plonout(1899:1905),platout(1899:1905),inoutres)
 inoutidout8 = gsn_add_polymarker(wks1,evspsblid,plonout(1970:1976),platout(1970:1976),inoutres)
end if

;---Map: overlay data and polygon on the map
;overlay(mpid,evspsblid)
draw(evspsblid)
;draw(mpid)
frame(wks1)

;---------------Annual sums of surface runoff from the catchment area------------------
;---Intro
ndim = dimsizes(evspsbl_yr)
nyr = ndim(0)
summ = new(nyr,float)
summ3 = new(nyr,float)
summ = 0.
summ3 = 0.
accano = new(nyr,float)
accano = 0.
wgt1 = new((/dimsizes(lat),dimsizes(lon)/),float) ; 50 values
wgt2 = new((/nmps2/),float) ; 135 values
wgt3 = new((/nmps3/),float); 35 values
wgt1 = 0.
wgt2 = 0.
wgt3 = 0.

;---Constants and limits
pi = 4.*atan(1.)
dx = 2.*pi/180. ; 2 degrees grid step
dy = 2.*pi/180. ; 2 degrees grid step

if (maparea .eq. 1)
 y1 = LAT(1194) ; left bottom corner
 y2 = LAT(1202) ; left upper corner
 y3 = LAT(2196) ; right upper corner
 ;print("y1(1194) "+y1+" y2(1202) "+y2)
 x1 = LON(1194)
 x2 = LON(2187) ; right bottom corner
 x3 = LON(2196)
 ;print("x1(1194) "+x1+" x2(2187) "+x2)
else
 y1 = LAT(1686) ; left bottom corner
 y2 = LAT(1692) ; left upper corner
 y3 = LAT(1976) ; right upper corner
; print("y1(1686) "+y1+" y2(1692) "+y2)
 x1 = LON(1686)
 x2 = LON(1970) ; right bottom corner
 x3 = LON(1976)
; print("x1(1686) "+x1+" x2(1970) "+x2)
end if
;---Types of summation
SUMpol = False ; sum1 from area inside polygon (50 mesh points)
SUMrec = True ; sum2 from rectangular area (135 mesh points)
SUMcel = False ; sum3 from one grid cell

;---Summation
if (SUMpol) then
do n = 0, nmps-1
 wgt1(LAT(j_platin(n)),LON(i_plonin(n))) = weights1(n,2)
end do

do yr = 0, nyr-1
summ(yr) = 0
summ3(yr) = 0
 do n = 0, nmps-1
  area1 = doubletofloat(4.*12387.69*10.^6*cos(platin(j_platin(n))*4.*atan(1.)/180.)) ; 2*2degrees*111.3^2 10^6 m*cos(latitude)
  summ(yr) = summ(yr) + evspsbl_yr(yr,LAT(j_platin(n)),LON(i_plonin(n)))*area1*wgt1(LAT(j_platin(n)),LON(i_plonin(n))) ; [kg]
 end do
end do
end if

if (SUMrec .and. maparea .eq. 1) then
do yr = 0, nyr-1
summ(yr) = 0
bas = 0.
m = 0
 do lo = x1,x2 ; 46:62 LAT
  do la = y1,y2 ; 32:60 LON
   if ismissing(evspsbl_yr(yr,la,lo)) then
    evspsbl_yr(yr,la,lo) = 0.
   end if
   area2 = doubletofloat(6371000.^2.*cos(la*pi/180.)*dx*2*sin(dy/2)) ; grid cell area
   summ(yr) = summ(yr) + evspsbl_yr(yr,la,lo)*area2*weights2(m) ; [kg]
   summ3(yr) = summ3(yr) + evspsbl_yr(yr,63,23)*area2*weights2(69) ; [kg]
;   bas = bas + area2*weights2(m) ; [m2]
   m = m+1
  end do
 end do
end do
end if

if (SUMrec .and. maparea .eq. 2) then
 do yr = 0, nyr-1
 summ(yr) = 0
 bas = 0.
 m = 0
  do lo = x1,x2 ; 36:48 LAT
   do la = y1,y2 ; 46:54 LON
    if ismissing(evspsbl_yr(yr,la,lo)) then
     evspsbl_yr(yr,la,lo) = 0.
    end if
    area2 = doubletofloat(6371000.^2.*cos(la*pi/180.)*dx*2*sin(dy/2)) ; grid cell area
    summ(yr) = summ(yr) + evspsbl_yr(yr,la,lo)*area2*weights3(m) ; [kg]
;    summ3(yr) = summ3(yr) + evspsbl_yr(yr,63,23)*area2*weights3(69) ; [kg]
    bas = bas + area2*weights2(m) ; [m2]
    m = m+1
   end do
  end do
 end do
end if

print("basin [km2]: "+bas/10.^6) ; [km2]

;---Volume of runoff
vol = summ/(10.^12) ; /10.^3 kg/m3 * 10.^9 m3 -> [km3]

;---Accumulated anomalies
avevol = avg(vol) ; average volume
;print("average volume: "+avevol)
N = num(.not.ismissing(vol))
;print("number of missed values: "+N)

if (var .eq. "evspsbl") then
 accano = cumsum((vol-avevol)/avevol,2) ; accumulated amomalies
 fout1 = outdir2+"f_"+var+"_accano_"+maparea+"_"+model+".dat"
 fout2 = outdir2+"f_"+var+"_volume_"+maparea+"_"+model+".dat"
else
 accano = cumsum(vol,2)
 fout2 = outdir2+"f_"+var+"_volume_"+maparea+"_"+model+".dat"
end if
asciiwrite(fout1,sprintf("%f",accano))
asciiwrite(fout2,sprintf("%f",vol))

runavevol = runave(vol,10,0)

;---Graphs
g2_name = model+"_"+var+"_"+maparea+"_g_vol"
g3_name = model+"_"+var+"_"+maparea+"_g_runavevol"
g4_name = model+"_"+var+"_"+maparea+"_g_accano"
wks2 = gsn_open_wks("png",indir+outdir+g2_name)
wks3 = gsn_open_wks("png",indir+outdir+g3_name)
wks4 = gsn_open_wks("png",indir+outdir+g4_name)

xyres2 = True
xyres3 = xyres2
xyres4 = xyres2
xyres2@tiMainString = model+" Volume of "+varlong
xyres3@tiMainString = model+" 10 years running average"
xyres4@tiMainString = model+" Accumulated anomalies"

plot2 = gsn_csm_xy(wks2,ispan(0,dimsizes(vol)-1,1),vol,xyres2)
plot3 = gsn_csm_xy(wks3,ispan(0,dimsizes(runavevol)-1,1),runavevol,xyres3)
plot4 = gsn_csm_xy(wks4,ispan(0,dimsizes(accano)-1,1),accano,xyres4)

end
