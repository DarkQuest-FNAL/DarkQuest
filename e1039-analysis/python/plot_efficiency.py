from ROOT import *
gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

ch = TChain("Events")
#ch.Add("/cms/data/seaquest/users/dsperka/DarkQuest/16May2024/Brem_0.85*.root")
#ch.Add("/cms/data/seaquest/users/dsperka/DarkQuest/16May2024/Brem_*.root")
ch.Add("/home/dsperka/DarkQuest-FNAL/DarkQuest/e1039-analysis/SimHits/macro/output.root")

observables  = ["truthdimuon_pz[0]"]
h_total      = {}
h_trig       = {}
h_reco       = {}
h_recotrig   = {}
eff_trig     = {}
eff_reco     = {}
eff_recoandtrig = {}
eff_recotrig = {}
eff_trigreco = {}
c = {}

for obs in observables:

    if obs=="truthdimuon_pz[0]":
        h_total[obs] = TH1F("total_"+obs,"total_"+obs,20,0,80)
        h_trig[obs] = TH1F("trig_"+obs,"trig_"+obs,20,0,80)
        h_reco[obs] = TH1F("reco_"+obs,"reco_"+obs,20,0,80)
        h_recotrig[obs] = TH1F("recotrig_"+obs,"recotrig_"+obs,20,0,80)
    if obs=="truthdimuon_mass":
        h_total[obs] = TH1F("total_"+obs,"total_"+obs,55,0,0.55)
        h_trig[obs] = TH1F("trig_"+obs,"trig_"+obs,55,0,0.55)
        h_reco[obs] = TH1F("reco_"+obs,"reco_"+obs,55,0,0.55)
        h_recotrig[obs] = TH1F("recotrig_"+obs,"recotrig_"+obs,55,0,0.55)


    h_total[obs].Sumw2()
    h_trig[obs].Sumw2()
    h_reco[obs].Sumw2()
    h_recotrig[obs].Sumw2()
    
    #ch.Draw(obs+">>total_"+obs,"weight*(n_truthdimuons>0)","goff")
    #ch.Draw(obs+">>trig_"+obs,"weight*(n_truthdimuons>0 && nim_trigger[3]>0)","goff")
    #ch.Draw(obs+">>reco_"+obs,"weight*(n_truthdimuons>0 && n_dimuons>0)","goff")
    #ch.Draw(obs+">>recotrig_"+obs,"weight*(n_truthdimuons>0 && n_dimuons>0 && nim_trigger[3]>0)","goff")

    #den = "(truthtrack_rectrack_id[0]>-1) && (truthtrack_rectrack_id[1]>-1)"
    den = "1==1"
    ch.Draw(obs+">>total_"+obs,"weight*("+den+")","goff")
    ch.Draw(obs+">>trig_"+obs,"weight*("+den+" && nim_trigger[3]>0)","goff")
    ch.Draw(obs+">>reco_"+obs,"weight*("+den+" && truthdimuon_recoed[0]>0)","goff")
    ch.Draw(obs+">>recotrig_"+obs,"weight*("+den+" && truthdimuon_recoed[0]>0 && nim_trigger[3]>0)","goff")

    #ch.Draw(obs+">>total_"+obs,"(n_truthdimuons>0)","goff")
    #ch.Draw(obs+">>trig_"+obs,"(n_truthdimuons>0 && nim_trigger[3]>0)","goff")
    #ch.Draw(obs+">>reco_"+obs,"(n_truthdimuons>0 && n_dimuons>0)","goff")
    #ch.Draw(obs+">>recotrig_"+obs,"(n_truthdimuons>0 && n_dimuons>0 && nim_trigger[3]>0)","goff")        

    eff_trig[obs] = TEfficiency(h_trig[obs],h_total[obs])
    eff_trig[obs].SetStatisticOption(TEfficiency.kFCP);
    eff_reco[obs] = TEfficiency(h_reco[obs],h_total[obs]);
    eff_reco[obs].SetStatisticOption(TEfficiency.kFCP);    
    eff_recoandtrig[obs] = TEfficiency(h_recotrig[obs],h_total[obs]);    
    eff_recoandtrig[obs].SetStatisticOption(TEfficiency.kFCP);

    eff_recotrig[obs] = TEfficiency(h_recotrig[obs],h_trig[obs]);    
    eff_recotrig[obs].SetStatisticOption(TEfficiency.kFCP);

    eff_trigreco[obs] = TEfficiency(h_recotrig[obs],h_reco[obs]);    
    eff_trigreco[obs].SetStatisticOption(TEfficiency.kFCP);

    c[obs] = TCanvas("c_"+obs,"c_"+obs,800,600)
    c[obs].cd()
    c[obs].SetLogy()
    h_total[obs].GetYaxis().SetRangeUser(0.001,9.0)
    h_total[obs].GetXaxis().SetTitle(obs)
    h_total[obs].SetLineWidth(0)
    h_total[obs].SetMarkerSize(0)
    h_total[obs].Draw("h")
    eff_reco[obs].SetLineColor(2)
    eff_reco[obs].SetMarkerColor(2)
    eff_reco[obs].Draw("Psame")
    eff_trig[obs].SetLineColor(4)
    eff_trig[obs].SetMarkerColor(4)
    eff_trig[obs].Draw("PSame")
    eff_recoandtrig[obs].SetLineColor(1)
    eff_recoandtrig[obs].SetMarkerColor(1)
    eff_recoandtrig[obs].Draw("Psame")
    eff_recotrig[obs].SetLineColor(6)
    eff_recotrig[obs].SetMarkerColor(6)
    eff_recotrig[obs].Draw("Psame")
    eff_trigreco[obs].SetLineColor(8)
    eff_trigreco[obs].SetMarkerColor(8)
    eff_trigreco[obs].Draw("Psame")
    
    leg = TLegend(0.18,0.74,0.48,0.89)
    leg.AddEntry(eff_reco[obs],"eff(reco.)","lp")
    leg.AddEntry(eff_trig[obs],"eff(trig.)","lp")
    leg.AddEntry(eff_recoandtrig[obs],"eff(reco. and trig.)","lp")
    leg.AddEntry(eff_recotrig[obs],"eff(reco.|trig.)","lp")
    leg.AddEntry(eff_trigreco[obs],"eff(trig.|reco.)","lp")
    leg.Draw("same")
    
    c[obs].SaveAs("plots/efficiency_"+obs+".pdf")
