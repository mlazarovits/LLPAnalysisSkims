import ROOT
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input","-i",help="file with canvases to draw",required=True)
parser.add_argument("--odir","-o",help="output directory for plots",required=True)
parser.add_argument("--format","-f",help="format to save plots",default="png")
parser.add_argument("--remake",help="remake existing files",default=False,action='store_true')
args = parser.parse_args()

if not os.path.isdir(args.odir):
	print("Creating directory",args.odir)
	os.mkdir(args.odir)

file = ROOT.TFile.Open(args.input,"READ")
for key in file.GetListOfKeys():
	can = file.Get(key.GetName())
	can.Draw()
	if os.path.exists(args.odir+"/"+can.GetName()+"."+args.format) and not args.remake:
		continue
	can.SaveAs(args.odir+"/"+can.GetName()+"."+args.format)
file.Close()
