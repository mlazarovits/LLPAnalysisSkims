import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input","-i",help="file with canvases to draw",required=True)
parser.add_argument("--odir","-o",help="output directory for plots",required=True)
parser.add_argument("--format","-f",help="format to save plots",default="png")
args = parser.parse_args()

file = ROOT.TFile.Open(args.input,"READ")
for key in file.GetListOfKeys():
	can = file.Get(key.GetName())
	can.Draw()
	can.SaveAs(args.odir+"/"+can.GetName()+"."+args.format)
file.Close()
