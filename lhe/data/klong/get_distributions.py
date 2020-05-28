#
# Run `python get_distributions.py` to convert klong_pythia events into kinematics.root
#
from klong_pythia import data
from math import sqrt
import ROOT

# histo helpers
def add_bin(h, x1, x2):
    h.SetBinContent(x1, h.GetBinContent(x1) + h.GetBinContent(x2))
    h.SetBinError(x1, sqrt(pow(h.GetBinError(x1),2) + pow(h.GetBinError(x2),2)))
    h.SetBinContent(x2, 0)
    h.SetBinError(x2, 0)
def add_bin_2d(h, x1, y1, x2, y2):
    h.SetBinContent(x1,y1, h.GetBinContent(x1,y1) + h.GetBinContent(x2,y2))
    h.SetBinError(x1,y1, sqrt(pow(h.GetBinError(x1,y1),2) + pow(h.GetBinError(x2,y2),2)))
    h.SetBinContent(x2,y2, 0)
    h.SetBinError(x2,y2, 0)
def remove_overflow(h):
    add_bin(h, 1, 0)
    add_bin(h, h.GetNbinsX(), h.GetNbinsX()+1)
def remove_overflow_2d(h):
    for binx in range(h.GetNbinsX()+2):
        add_bin_2d(h, binx, h.GetNbinsY(), binx, h.GetNbinsY()+1)
	add_bin_2d(h, binx, 1, binx, 0)
    for biny in range(h.GetNbinsY()+2):
        add_bin_2d(h, h.GetNbinsX(), biny, h.GetNbinsX()+1, biny)
        add_bin_2d(h, 1, biny, 0, biny)

#
f = ROOT.TFile("kinematics.root","recreate")
h = ROOT.TH2F("klong","; K_{L} p_{T} [GeV]; K_{L} p_{z} [GeV]",5,0,1,15,0,30)

n=0
for ks in data:
    for k in ks:
        if len(k)!=5:
            print 'bad val', k
            continue
        x,y,z,m,e = k
        h.Fill(sqrt(x*x+y*y),z)
        n+=1
        
remove_overflow_2d(h)
h.Scale(1./h.Integral())
h.Write()
f.Close()

