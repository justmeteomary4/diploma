load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

;---Directory
indir = asciiread("indir.dat",-1,"string")
outdir = "maps/"
system("if ! test -d " + indir+outdir + "; then mkdir " + indir+outdir +"; fi")

;---Models
;model = "BCC"
;model = "CCSM4"
;model = "CSIRO"
model = "GISS_121"
;model = "HadCM3"
;model = "IPSL"
;model = "MIROC"
;model = "MPI"

;---Input files: NetCDF surface runoff
cmd = "find "+indir+"regrid/"+model+"_reg* -type f -print"
f = systemfunc(cmd)
;print(f)
files = addfiles(f, "r")
ListSetType(files,"cat")
fnc_name = model+"_reg_mrros_Lmon_GISS-E2-R_past1000_r1i1p121_085001-089912.nc"
fnc = addfile(indir+"regrid/"+fnc_name,"r")

;---Input files: watershed coordinates
fwat1_name = "basin_515.dat" ; robust
fwat2_name = "basin_899.dat" ; fine
nps = 899 ; number of points
if (nps .eq. 515) then
 wat = "robust"
else
 wat = "fine"
end if
watershed = asciiread(fwat2_name,(/nps,3/),"float")
nlat = watershed(:,2)
nlon = watershed(:,1)

;---Input files: weights
fwgt1_name = "weights.dat"
fwgt2_name = "weights_rec.dat"
nmps1 = 50 ; number of mesh points
nmps2 = 135
weights1 = asciiread(fwgt1_name,(/nmps1,3/),"float")
weights2 = asciiread(fwgt2_name,(/nmps2/),"float")

;---Variables
time = files[:]->time
lat = fnc->lat
lon = fnc->lon
mrros = files[:]->mrros*2592000 ; kg/m2*c -> [kg/m2] = mm
;print(time(0:20))

mrros!0 = "time"
mrros!1 = "lat"
mrros!2 = "lon"
mrros&lat = lat
mrros&lon = lon
mrros&lat@units = "degrees_north"
mrros&lon@units = "degrees_east"

mrros_yr = month_to_annual(mrros,0)
;printVarSummary(mrros_yr)

;---Map
title = model+"_m_"+wat+"_test"
wkstype = "png"
wkstype@wkWidth = 2000
wkstype@wkHeight = 1800
wks1 = gsn_open_wks(wkstype,indir+outdir+title)
gsn_define_colormap(wks1,"WhiteBlue")

mrrosres = True
mrrosres@gsnDraw = False
mrrosres@gsnFrame = False
txres = mrrosres
mpres = mrrosres

mrrosres@gsnAddCyclic = False
mrrosres@cnFillOn = True
mrrosres@cnLinesOn = False
mrrosres@cnLineColor = "Black"
;mrrosres@cnLevelSpacingF = 5.0
;mrrosres@cnLineThicknessF = 2.0
;mrrosres@lbLabelAutoStride = True
mrrosres@lbLabelPosition = "Right"

mpres@mpLimitMode = "LatLon"
mpres@mpMinLatF = 40.
mpres@mpMaxLatF = 65.
mpres@mpMinLonF = 30.
mpres@mpMaxLonF = 65.
mpres@mpGridAndLimbDrawOrder = "Predraw"
mpres@tiMainString = model+" Surface runoff"
mpres@gsnLeftString = ""
mpres@gsnRightString = "mm/year"

;txres@txFontHeightF = 0.01
;txres@txJust = "BottomLeft"

mrrosid = gsn_csm_contour(wks1,mrros(252,:,:),mrrosres)
mpid = gsn_csm_map(wks1,mpres)

;---Polygon: draw a watershed on the map
dum = new(nps,graphic)
do i = 0,nps-2
 watres = True
 watres@gsMarkerIndex = 16
 watres@gsMarkerSizeF = 200
 dum(i) = gsn_add_polyline(wks1,mrrosid,nlon(i:i+1),nlat(i:i+1),watres)
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
inoutres2 = True
inoutres2@gsMarkerIndex = 20
inoutres2@gsMarkerColor = "red"
inoutres2@gsMarkerSizeF = 0.007
platin = where(inout,plat,plat@_FillValue)
plonin = where(inout,plon,plon@_FillValue)
j_platin = ind(.not.ismissing(platin))
i_plonin = ind(.not.ismissing(plonin))
inoutid = gsn_add_polymarker(wks1,mrrosid,plonin,platin,inoutres) 
inoutid2 = gsn_add_polymarker(wks1,mrrosid,plonin(1341),platin(1341),inoutres2)

;---Map: overlay data and polygon on the map
overlay(mpid,mrrosid)
draw(mpid)
frame(wks1)

;---------------Annual sums of surface runoff from the catchment area------------------
;---Intro
ndim = dimsizes(mrros_yr)
nyr = ndim(0)
summ = new(nyr,float)
summ3 = new(nyr,float)
summ = 0.
summ3 = 0.
accano = new(nyr,float)
accano = 0.
wgt1 = new((/dimsizes(lat),dimsizes(lon)/),float) ; 50 values
wgt2 = new((/nmps2/),float) ; 135 values
wgt1 = 0.
wgt2 = 0.

;---Constants and limits
pi = 4.*atan(1.)
dx = 2.*pi/180. ; 2 degrees grid step
dy = 2.*pi/180. ; 2 degrees grid step

y1 = LAT(1194) ; left bottom corner
y2 = LAT(1202) ; left upper corner
y3 = LAT(2196) ; right upper corner
print("y1(1194) "+y1+" y2(1202) "+y2)
x1 = LON(1194)
x2 = LON(2187) ; right bottom corner
x3 = LON(2196)
print("x1(1194) "+x1+" x2(2187) "+x2)

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
  summ(yr) = summ(yr) + mrros_yr(yr,LAT(j_platin(n)),LON(i_plonin(n)))*area1*wgt1(LAT(j_platin(n)),LON(i_plonin(n))) ; [kg]
 end do
end do
end if

if (SUMrec) then
do yr = 0, nyr-1
summ(yr) = 0
bas = 0
m = 0
 do lo = x1,x2 ; 46:62 LAT
  do la = y1,y2 ; 32:60 LON
   if ismissing(mrros_yr(yr,la,lo)) then
    mrros_yr(yr,la,lo) = 0.
   end if
   area2 = doubletofloat(6371000.^2.*cos(la*pi/180.)*dx*2*sin(dy/2)) ; grid cell area
   summ(yr) = summ(yr) + mrros_yr(yr,la,lo)*area2*weights2(m) ; [kg]
   summ3(yr) = summ3(yr) + mrros_yr(yr,63,23)*area2*weights2(69) ; [kg]
;   bas = bas + area2*weights2(m) ; [m2]
   m = m+1
  end do
 end do
end do
end if

;print("basin [km2]: "+bas/10.^6) ; [km2]

;---Volume of runoff
vol = summ/(10.^12) ; /10.^3 kg/m3 * 10.^9 m3 -> [km3]

fout = indir+outdir+model+"_f_volume.dat"
asciiwrite(fout,sprintf("%f",vol))

;---Accumulated anomalies
avevol = avg(vol)
print("average volume: "+avevol)
N = num(.not.ismissing(vol))
print(""+N)

;accano = cumsum(summ,0)./avesum.

runavevol = runave(vol,10,0)

;---Graphs
g2_name = model+"_g_vol"
g3_name = model+"_g_runavevol"
g4_name = model+"_g_accano"
wks2 = gsn_open_wks("png",indir+outdir+g2_name)
wks3 = gsn_open_wks("png",indir+outdir+g3_name)
wks4 = gsn_open_wks("png",indir+outdir+g4_name)

xyres2 = True
xyres3 = xyres2
xyres4 = xyres2
xyres2@tiMainString = model+" Volume of runoff"
xyres3@tiMainString = model+" 10 years running average"
xyres4@tiMainString = model+" Accumulated anomalies"

plot2 = gsn_csm_xy(wks2,ispan(0,dimsizes(vol)-1,1),vol,xyres2)
plot3 = gsn_csm_xy(wks3,ispan(0,dimsizes(runavevol)-1,1),runavevol,xyres3)
;plot4 = gsn_csm_xy(wks4,ispan(0,dimsizes(accano)-1,1),accano,xyres4)

end
