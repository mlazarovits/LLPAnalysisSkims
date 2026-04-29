from FileProcessor import FileProcessor
from ROOT import RDataFrame, TFile, TGraph
import os
from array import array
import cmsstyle as CMS

def GetTotalMassSplit(filename):
    mgl = filename[filename.find("mGl-")+4:]
    mgl = float(mgl[:mgl.find("_")])
    mn1 = filename[filename.find("mN1-")+4:]
    mn1 = float(mn1[:mn1.find("_")])
    return mgl - mn1

#import kerberos credentials to conda env if not already there
kerb = os.getenv("KRB5CCNAME")
if(kerb is None):
    print("Setting kerebos credentials")
    os.environ["KRB5CCNAME"] = "API:"


fileprocessor = FileProcessor()
files = fileprocessor.GetFiles("SMS_gogoGZ")
files_10cm = [file for file in files if "ct0p1" in file]
files_50cm = [file for file in files if "ct0p5" in file]

#(mass_split, ms_max)
dms, mspeaks = array('d'), array('d')
for file in files_10cm:
    dm = GetTotalMassSplit(file)
    if(dm < 200):
        continue
    rdf = RDataFrame("kuSkimTree",file)
    rdf = rdf.Define("Ms","rjr_Ms[0]")
    #do met cut + pts cut
    rdf_presel = rdf.Filter("rjrPTS[0] < 150 && selCMet > 150") 
    #do N-1 Rs cut
    rdf0 = rdf_presel.Filter("rjr_Rs[0] > 0.15")
    msmean = rdf0.Mean("Ms").GetValue()
    dms.append(dm)
    mspeaks.append(msmean)
gr_10cm = TGraph(len(dms),dms,mspeaks)
gr_10cm.SetTitle("10cm")
gr_10cm.SetName("gr_10cm")
gr_10cm.SetMarkerColor(CMS.p10.kBlue)
gr_10cm.SetMarkerStyle(20)

dms, mspeaks = array('d'), array('d')
for file in files_50cm:
    dm = GetTotalMassSplit(file)
    if(dm < 200):
        continue
    rdf = RDataFrame("kuSkimTree",file)
    rdf = rdf.Define("Ms","rjr_Ms[0]")
    #do met cut + pts cut
    rdf_presel = rdf.Filter("rjrPTS[0] < 150 && selCMet > 150") 
    #do N-1 Rs cut
    rdf0 = rdf_presel.Filter("rjr_Rs[0] > 0.15")
    msmean = rdf0.Mean("Ms").GetValue()
    dms.append(dm)
    mspeaks.append(msmean)
gr_50cm = TGraph(len(dms),dms,mspeaks)
gr_50cm.SetTitle("50cm")
gr_50cm.SetName("gr_50cm")
gr_50cm.SetMarkerColor(CMS.p10.kOrange)
gr_50cm.SetMarkerStyle(20)


leg = CMS.cmsLeg(0.6,0.6,0.9,0.7)
leg.AddEntry(gr_10cm)
leg.AddEntry(gr_50cm)
leg.SetTextSize(0.04)

name = "ms_peak"
x_min = min(dms)-50
x_max = max(dms)+50
y_min = min(mspeaks)-100
y_max = max(mspeaks)+100
x_label = "total mass splitting [GeV]"
y_label = "mean Ms [GeV]"
CMS.SetLumi(None,run = "Run 3")
CMS.SetEnergy(13.6)
cv = CMS.cmsCanvas(name, x_min, x_max, y_min, y_max, x_label, y_label,square=False, extraSpace=0.05, iPos=0, with_z_axis=False, yTitOffset=1.43)
CMS.GetCmsCanvasHist(cv).GetYaxis().SetTitleSize(0.04)
CMS.GetCmsCanvasHist(cv).GetYaxis().SetLabelSize(0.04)
CMS.GetCmsCanvasHist(cv).GetXaxis().SetTitleSize(0.04)
CMS.GetCmsCanvasHist(cv).GetXaxis().SetLabelSize(0.04)
cv.cd()
gr_10cm.Draw("P")
gr_50cm.Draw("sameP")
leg.Draw("same")

ofilename = "msPeaks.root"
ofile = TFile.Open(ofilename,"RECREATE")
ofile.cd()
cv.Write()
ofile.Close()
print("wrote to",ofilename)

 
	
	 
