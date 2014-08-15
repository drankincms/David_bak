from ROOT import *
gStyle.SetOptStat(0)

files={}
files["elec"]="RootFiles_For2DLimits_23Sep/electron_BestJetJet2W_M_WprimeRight_Histos-twobtags_85_140_dr03_lep50.root"
files["mu"]="RootFiles_For2DLimits_23Sep/muon_BestJetJet2W_M_WprimeRight_Histos-twobtags_85_140_dr03_lep50.root"

c=TCanvas()
#c.SetLogy()
#c.SetBottomMargin(0.2)
c.Draw()

yDiv=0.35
uPad=TPad("uPad","",0,yDiv,1,1) #for actual plots
uPad.SetLogy()
uPad.SetBottomMargin(0.0)
uPad.SetRightMargin(.05)
uPad.SetLeftMargin(.18)
uPad.Draw()

lPad=TPad("lPad","",0,0,1,yDiv) #for runner
lPad.SetTopMargin(0.0)
lPad.SetBottomMargin(0.4)
lPad.SetRightMargin(.05)
lPad.SetLeftMargin(.18)
lPad.SetGridy()
lPad.Draw()

for channel in ["elec","mu"]:
    file=TFile(files[channel])
    w=file.Get(channel+"_invmass_twobtags__wjets")
    tt=file.Get(channel+"_invmass_twobtags__ttbar")

    b=w.Clone("background")
    b.SetTitle("")
    b.SetLineColor(kBlack)
    b.GetYaxis().SetLabelSize(0.08)
    b.Add(tt)
    
    bp=tt.Clone("background"); bm=tt.Clone("background");
    bp.SetLineColor(kRed); bm.SetLineColor(kBlue)
    for sys in ["jes","jer","btag","PU","matching","q2scale","zerotagshape","hf",'shape']:
        bp.Reset(); bm.Reset();
        
        wp=file.Get(channel+"_invmass_twobtags__wjets__"+sys+"__plus")
        wm=file.Get(channel+"_invmass_twobtags__wjets__"+sys+"__minus")
        if not wp:
            wp=file.Get(channel+"_invmass_twobtags__wjets")
            wm=file.Get(channel+"_invmass_twobtags__wjets")
        bp.Add(wp)
        bm.Add(wm)
        
        ttp=file.Get(channel+"_invmass_twobtags__ttbar__"+sys+"__plus")
        ttm=file.Get(channel+"_invmass_twobtags__ttbar__"+sys+"__minus")
        if not ttp:
            ttp=file.Get(channel+"_invmass_twobtags__ttbar")
            ttm=file.Get(channel+"_invmass_twobtags__ttbar")
        bp.Add(ttp)
        bm.Add(ttm)

        l=TLegend(0.35, 0.65, 0.89, 0.89)
        l.SetFillColor(kWhite)
        l.SetBorderSize(0)
        
        l.AddEntry(b,"nominal")
        l.AddEntry(bp,sys+" plus")
        if sys!='zerotagshape' and sys!='PU': l.AddEntry(bm,sys+" minus")

        uPad.cd()
        b.SetMaximum(10**4)
        b.SetMinimum(1.5)
        b.Draw()
        b.GetYaxis().SetTitleSize(0.12)
        b.GetYaxis().SetTitleOffset(0.6)        
        b.GetYaxis().SetTitle('Events/bin')
        bp.Draw("SAME")
        if sys!='zerotagshape' and sys!='PU': bm.Draw("SAME")
        l.Draw("SAME")

        latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.06)
        latex2.SetTextAlign(31) # align right
        if channel == 'elec':
            latex2.DrawLatex(0.89, 0.93, "CMS Preliminary, 19.5 fb^{-1} at #sqrt{s} = 8 TeV");
        if channel == 'mu':
            latex2.DrawLatex(0.89, 0.93, "CMS Preliminary, 19.5 fb^{-1} at #sqrt{s} = 8 TeV");
            
        latex3 = TLatex()
        latex3.SetNDC()
        latex3.SetTextSize(0.06)
        latex3.SetTextAlign(31) # align right
        if (channel == 'elec'):
            latex3.DrawLatex(0.9, 0.80, "e+jets N_{b-tags} = 2");
        if (channel == 'mu'):
            latex3.DrawLatex(0.9, 0.80, "#mu+jets N_{b-tags} = 2");
                                                                                                                    

        lPad.cd(0)
        rp=bp.Clone("rp"); rm=bm.Clone("rm")
        rp.GetXaxis().SetTitleSize(0.16)
        rp.GetXaxis().SetLabelSize(0.12)
        rp.GetYaxis().SetTitleSize(0.10)
        rp.GetYaxis().SetTitleOffset(0.45)
        rp.GetYaxis().SetLabelSize(0.08)
        rp.GetYaxis().CenterTitle()
        rp.SetTitle(";M(tb) [GeV];#frac{sys-nom}{nom}")
        rp.GetYaxis().SetNdivisions(5,2,0)
        rp.Add(b,-1); rm.Add(b,-1)
        rp.Divide(b); rm.Divide(b)
        rp.SetMaximum(.22); rp.SetMinimum(-.22)
        rp.Draw()
        if sys!='zerotagshape' and sys!='PU': rm.Draw("SAME")

        c.SaveAs(channel+"_"+sys+"_2btags.pdf")
