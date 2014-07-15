load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

;---Directory
indir = asciiread("indir_reanalysis.dat",-1,"string")
outdir1 = "/graphs_ncl/"
system("if ! test -d " + indir+outdir1 + "; then mkdir " + indir+outdir1 +"; fi")

;---List of models
model = (/"BCC","CCSM4","CSIRO","GISS_121","GISS_122","GISS_123","GISS_124",\
          "GISS_125","GISS_126","GISS_127","GISS_128","GISS_1221","HadCM3",\
          "MIROC","MPI","MRI"/)

;model = (/"BCC","CCSM4","GISS_121","GISS_122","GISS_123","GISS_124","GISS_125",\
;          "GISS_126","GISS_127","GISS_128","GISS_1221","MIROC","MPI","MRI"/)
nm = dimsizes(model)

;---Input files:
nyr = 1000 ; number of years
vol = new((/nm,nyr/),float)
fmodel_name = new(nm,string)
var = "mrros_volume"
maparea = "1"

do m = 0, nm-1
 fmodel_name(m) = "f_"+var+"_"+maparea+"_"+model(m)+".dat"
 vol(m,:) = asciiread(indir+fmodel_name(m),(/nyr/),"float")
end do

;---Input files: paleo reconstractions
rec_name = "all"
frec1_name = "paleo_Combined_withFillValue.dat"
frec2_name = "paleo_Klige_withFillValue.dat"
frec3_name = "paleo_Kroonenberg_withFillValue.dat"
frec4_name = "paleo_Rychagov_withFillValue.dat"
reco1_raw = asciiread(indir+frec1_name,(/nyr,2/),"float")
reco2_raw = asciiread(indir+frec2_name,(/nyr,2/),"float")
reco3_raw = asciiread(indir+frec3_name,(/nyr,2/),"float")
reco4_raw = asciiread(indir+frec4_name,(/nyr,2/),"float")
reco1_raw@_FillValue = -9999.
reco2_raw@_FillValue = -9999.
reco3_raw@_FillValue = -9999.
reco4_raw@_FillValue = -9999.
reco1 = linmsg(reco1_raw(:,1),0)
reco2 = linmsg(reco2_raw(:,1),0)
reco3 = linmsg(reco3_raw(:,1),0)
reco4 = linmsg(reco4_raw(:,1),0)

;---Accumulated anomalies
avevol = new(nm,float)
accano = new((/nm,nyr/),float)
fmodel_ave = new(nm,string)
fmodel_accano = new(nm,string)
do m = 0, nm-1
 avevol(m) = avg(vol(m,:)) ; average volume
 accano(m,:) = cumsum((vol(m,:)-avevol(m))/avevol(m),2) ; accumulated amomalies
end do

;---Graphs: specify X axis
tyr = new((/1,nyr/),"integer")
t = ispan(0,nyr-1,1)
t_offset = 850
tyr(0,:) = t + t_offset

;---Graphs: specify Y axis
ymin1 = min(accano)
ymin2 = min(reco3)
ymax1 = max(accano)
ymax2 = max(reco3)
print("accanomin: "+ymin1+" accanomax: "+ymax1+" recmin: "+ymin2+" recmax: "+ymax2)

;---Graph 1: each model vs all paleo reconstractions
mod1 = True
mod1@gsnFrame = False
mod1@gsnMaximize = True

mod1@xyLineColor = "blue"
mod1@xyDashPattern = 0 ; solid
mod1@xyLineThicknessF = 2.5

mod1@trXMinF = t_offset
mod1@trXMaxF = t_offset + nyr
mod1@trYMinF = ymin1
mod1@trYMaxF = ymax1

mod1@tiMainFontHeightF = 0.05
mod1@tiXAxisString = "Time, years"
mod1@tiYAxisString = "m/year"
mod1@tiXAxisFontHeightF = 0.02
mod1@tiYAxisFontHeightF = 0.02

rec0 = True
rec1 = True
rec2 = True
rec3 = True
rec4 = True
rec1@gsnFrame = False

rec1@xyLineColor = "firebrick2"
rec1@xyDashPattern = 9 ; dashed
rec1@xyLineThicknessF = 4.0

rec1@trXMinF = t_offset
rec1@trXMaxF = t_offset + nyr
rec1@trYMinF = ymin1
rec1@trYMaxF = ymax1

rec2 = rec1
rec2@xyLineColor = "darkviolet"
rec2@xyDashPattern = 3
rec3 = rec1
rec3@xyLineColor = "darkslategrey"
rec3@xyDashPattern = 12
rec4 = rec1
rec4@xyLineColor = "forestgreen"
rec4@xyDashPattern = 4

;rec4@gsnCenterString = "Caspian Sea level paleoreconstractions (850-1850)"
;rec4@tiMainFontHeightF = 0.05
;rec4@tiXAxisString = "Time, years"
;rec4@tiYAxisString = "m/year"
;rec4@tiXAxisFontHeightF = 0.02
;rec4@tiYAxisFontHeightF = 0.02

do m = 0, nm-1
; g1_name = rec_name
 g1_name = "g_accano_"+model(m)+"_vs_"+rec_name+"_withCSIROHadCM3"
 wks1 = gsn_open_wks("png",indir+outdir1+g1_name)
 mod1@gsnCenterString=model(m)+" accumulated anomalies vs "+rec_name+" reconstractions"
 plot1 = gsn_csm_xy2(wks1,tyr,accano(m,:),reco1,mod1,rec1)
 plot2 = gsn_csm_xy2(wks1,tyr,accano(m,:),reco2,mod1,rec2)
 plot3 = gsn_csm_xy2(wks1,tyr,accano(m,:),reco3,mod1,rec3)
 plot4 = gsn_csm_xy2(wks1,tyr,accano(m,:),reco4,mod1,rec4)
; plot1 = gsn_csm_xy2(wks1,tyr,reco1,reco2,rec1,rec2)
; plot2 = gsn_csm_xy2(wks1,tyr,reco3,reco4,rec3,rec4)
 overlay(plot1,plot2)
 overlay(plot1,plot3)
 overlay(plot1,plot4)
 frame(wks1)
end do

;---Graph 2: emsemble
ensem = new(nyr,float)
do yr = 0, nyr-1
ensem(yr) = avg(accano(:,yr))
end do
ensmin = min(ensem)
ensmax = max(ensem)
print("ensmin: "+ensmin+" ensmax: "+ensmax)

ens = True
ens@gsnFrame = False
ens@gsnMaximize = True

ens@xyLineColor = "black"
ens@xyDashPattern = 0 ; solid
ens@xyLineThicknessF = 3.0

ens@trXMinF = t_offset
ens@trXMaxF = t_offset + nyr
ens@trYMinF = ymin1;ymin2
ens@trYMaxF = ymax1;ymax2

ens@tiXAxisString = "Time, years"
ens@tiYAxisString = "m/year"
ens@tiXAxisFontHeightF = 0.02
ens@tiYAxisFontHeightF = 0.02

g2_name = "g_enseble_vs_"+rec_name+"_withCSIROHadCM3"
wks2 = gsn_open_wks("png",indir+outdir1+g2_name)
ens@gsnCenterString = "Ensemble accumulated anomalies vs "+rec_name+" reconstractions"
plot5 = gsn_csm_xy2(wks2,tyr,ensem,reco1,ens,rec1)
plot6 = gsn_csm_xy2(wks2,tyr,ensem,reco2,ens,rec2)
plot7 = gsn_csm_xy2(wks2,tyr,ensem,reco3,ens,rec3)
plot8 = gsn_csm_xy2(wks2,tyr,ensem,reco4,ens,rec4)
overlay(plot5,plot6)
overlay(plot5,plot7)
overlay(plot5,plot8)
frame(wks2)

end
