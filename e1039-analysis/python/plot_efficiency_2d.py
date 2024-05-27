from ROOT import *
import glob
#from DQacceptance import *
#print(Nevts.keys())

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(55)

#infiles = glob.glob("/cms/data/seaquest/users/dsperka/DarkQuest/15May2024/Eta*.root")
#infiles = glob.glob("/cms/data/seaquest/users/dsperka/DarkQuest/16May2024/Brem*.root")
infiles = glob.glob("/cms/data/seaquest/users/dsperka/DarkQuest/16May2024/Eta*.root")

obs = "truthdimuon_mass"
h_total      = {}
h_trig       = {}
h_reco       = {}
h_recotrig   = {}
eff_trig     = {}
eff_reco     = {}
eff_recotrig = {}
c = {}

#eta
xbins = 35; xlow = 0.2; xhigh = 0.55;
#brem
#xbins = 50; xlow = 0.2; xhigh = 1.2;

h_total["total"] = TH2F("total_"+"total","total_"+"total",xbins,xlow,xhigh,38,-7.8,-4.0)
h_trig["total"] = TH2F("trig_"+"total","trig_"+"total",xbins,xlow,xhigh,38,-7.8,-4.0)
h_reco["total"] = TH2F("reco_"+"total","reco_"+"total",xbins,xlow,xhigh,38,-7.8,-4.0)
h_recotrig["total"] = TH2F("recotrig_"+"total","recotrig_"+"total",xbins,xlow,xhigh,38,-7.8,-4.0)

h_total["total"].Sumw2()
h_trig["total"].Sumw2()
h_reco["total"].Sumw2()
h_recotrig["total"].Sumw2()

for f in infiles:
    ch = TChain("Events")
    ch.Add(f)

    mass = f.split("/")[-1].split("_")[1]
    logeps = f.split("/")[-1].split("_")[5].replace(".root","")
    eps = '{:0.2e}'.format(pow(10.,float(logeps)))
    massf = '{:0.3f}'.format(float(mass))
    #print(eps,mass)
    #print(Nevts[eps][massf])


    h_total[mass+"_"+logeps] = TH2F("total_"+mass+"_"+logeps,"total_"+mass+"_"+logeps,xbins,xlow,xhigh,38,-7.8,-4.0)
    h_trig[mass+"_"+logeps] = TH2F("trig_"+mass+"_"+logeps,"trig_"+mass+"_"+logeps,xbins,xlow,xhigh,38,-7.8,-4.0)
    h_reco[mass+"_"+logeps] = TH2F("reco_"+mass+"_"+logeps,"reco_"+mass+"_"+logeps,xbins,xlow,xhigh,38,-7.8,-4.0)
    h_recotrig[mass+"_"+logeps] = TH2F("recotrig_"+mass+"_"+logeps,"recotrig_"+mass+"_"+logeps,xbins,xlow,xhigh,38,-7.8,-4.0)
    
    h_total[mass+"_"+logeps].Sumw2()
    h_trig[mass+"_"+logeps].Sumw2()
    h_reco[mass+"_"+logeps].Sumw2()
    h_recotrig[mass+"_"+logeps].Sumw2()

    #ch.Draw(logeps+":"+obs+">>total_"+mass+"_"+logeps,"weight*(n_truthdimuons>0)","goff")
    #ch.Draw(logeps+":"+obs+">>trig_"+mass+"_"+logeps,"weight*(n_truthdimuons>0 && nim_trigger[3]>0)","goff")
    #ch.Draw(logeps+":"+obs+">>reco_"+mass+"_"+logeps,"weight*(n_truthdimuons>0 && n_dimuons>0)","goff")
    #ch.Draw(logeps+":"+obs+">>recotrig_"+mass+"_"+logeps,"weight*(n_truthdimuons>0 && n_dimuons>0 && nim_trigger[3]>0)","goff")

    ch.Draw(logeps+":"+obs+">>total_"+mass+"_"+logeps,"(n_truthdimuons>0)","goff")
    ch.Draw(logeps+":"+obs+">>trig_"+mass+"_"+logeps,"(n_truthdimuons>0 && nim_trigger[3]>0)","goff")
    ch.Draw(logeps+":"+obs+">>reco_"+mass+"_"+logeps,"(n_truthdimuons>0 && n_dimuons>0)","goff")
    ch.Draw(logeps+":"+obs+">>recotrig_"+mass+"_"+logeps,"(n_truthdimuons>0 && n_dimuons>0 && nim_trigger[3]>0)","goff")

    efftrig = 0.0
    effreco = 0.0
    effrecotrig = 0.0
    if (h_total[mass+"_"+logeps].Integral()>0):
        efftrig = h_trig[mass+"_"+logeps].Integral()/h_total[mass+"_"+logeps].Integral()
        effreco = h_reco[mass+"_"+logeps].Integral()/h_total[mass+"_"+logeps].Integral()
        effrecotrig = h_recotrig[mass+"_"+logeps].Integral()/h_total[mass+"_"+logeps].Integral()
    
        print(f.split("/")[-1],efftrig,effreco,effrecotrig)
    
    h_total["total"].Add(h_total[mass+"_"+logeps])
    h_trig["total"].Add(h_trig[mass+"_"+logeps])
    h_reco["total"].Add(h_reco[mass+"_"+logeps])
    h_recotrig["total"].Add(h_recotrig[mass+"_"+logeps])

obs="total"

#eff_trig[obs] = TEfficiency(h_trig[obs],h_total[obs])
#eff_trig[obs].SetStatisticOption(TEfficiency.kBJeffrey);
#eff_reco[obs] = TEfficiency(h_reco[obs],h_total[obs]);
#eff_reco[obs].SetStatisticOption(TEfficiency.kBJeffrey);    
#eff_recotrig[obs] = TEfficiency(h_recotrig[obs],h_total[obs]);    
#eff_recotrig[obs].SetStatisticOption(TEfficiency.kBJeffrey);

eff_trig[obs] = h_trig[obs].Clone("eff_trig")
eff_trig[obs].Divide(h_total[obs])

eff_reco[obs] = h_reco[obs].Clone("eff_reco")
eff_reco[obs].Divide(h_total[obs])

eff_recotrig[obs] = h_recotrig[obs].Clone("eff_recotrig")
eff_recotrig[obs].Divide(h_total[obs])

c["reco"] = TCanvas("c_"+obs,"c_"+obs,800,600)
c["reco"].SetLeftMargin(0.18)
c["reco"].SetRightMargin(0.18)
c["reco"].SetBottomMargin(0.12)
c["reco"].cd()
c["reco"].SetLogz()
eff_reco[obs].GetXaxis().SetTitle("mass (GeV)")
eff_reco[obs].GetYaxis().SetTitle("log(#varepsilon)")
eff_reco[obs].GetZaxis().SetTitle("eff(reco.)")
eff_reco[obs].GetZaxis().SetRangeUser(0.01,1.0)
eff_reco[obs].Draw("colz")
#c["reco"].SaveAs("plots/eff_reco_mass_vs_logeps.pdf")
    

c["trig"] = TCanvas("c_"+obs,"c_"+obs,800,600)
c["trig"].SetLeftMargin(0.18)
c["trig"].SetRightMargin(0.18)
c["trig"].SetBottomMargin(0.12)
c["trig"].cd()
c["trig"].SetLogz()
eff_trig[obs].GetXaxis().SetTitle("mass (GeV)")
eff_trig[obs].GetYaxis().SetTitle("log(#varepsilon)")
eff_trig[obs].GetZaxis().SetTitle("eff(trig.)")
eff_trig[obs].GetZaxis().SetRangeUser(0.01,1.0)
eff_trig[obs].Draw("colz")
#c["trig"].SaveAs("plots/eff_trig_mass_vs_logeps.pdf")
    
c["recotrig"] = TCanvas("c_"+obs,"c_"+obs,800,600)
c["recotrig"].SetLeftMargin(0.18)
c["recotrig"].SetRightMargin(0.18)
c["recotrig"].SetBottomMargin(0.12)
c["recotrig"].cd()
c["recotrig"].SetLogz()
eff_recotrig[obs].GetXaxis().SetTitle("mass (GeV)")
eff_recotrig[obs].GetYaxis().SetTitle("log(#varepsilon)")
eff_recotrig[obs].GetZaxis().SetTitle("eff(reco)#times eff(trig.)")
eff_recotrig[obs].GetZaxis().SetTitleOffset(1.2)
eff_recotrig[obs].GetZaxis().SetRangeUser(0.001,0.1)
eff_recotrig[obs].Draw("colz")
#c["recotrig"].SaveAs("plots/eff_recotrig_mass_vs_logeps.pdf")
    
