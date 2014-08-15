from ROOT import *
gStyle.SetOptStat(0)

def setTDRStyle():

    # For the canvas:
    ROOT . gStyle . SetCanvasBorderMode(0)
    ROOT . gStyle . SetCanvasColor(0) # must be kWhite but I dunno how to do that in PyROOT
    ROOT . gStyle . SetCanvasDefH(600) #Height of canvas
    ROOT . gStyle . SetCanvasDefW(600) #Width of canvas
    ROOT . gStyle . SetCanvasDefX(0)   #POsition on screen
    ROOT . gStyle . SetCanvasDefY(0)


# For the Pad:
    ROOT . gStyle . SetPadBorderMode(0);
    # ROOT . gStyle . SetPadBorderSize(Width_t size = 1);
    ROOT . gStyle . SetPadColor(0); # kWhite
    ROOT . gStyle . SetPadGridX(0); #false
    ROOT . gStyle . SetPadGridY(0); #false
    ROOT . gStyle . SetGridColor(0);
    ROOT . gStyle . SetGridStyle(3);
    ROOT . gStyle . SetGridWidth(1);

# For the frame:
    ROOT . gStyle . SetFrameBorderMode(0);
    ROOT . gStyle . SetFrameBorderSize(1);
    ROOT . gStyle . SetFrameFillColor(0);
    ROOT . gStyle . SetFrameFillStyle(0);
    ROOT . gStyle . SetFrameLineColor(1);
    ROOT . gStyle . SetFrameLineStyle(1);
    ROOT . gStyle . SetFrameLineWidth(1);

# For the histo:
    # ROOT . gStyle . SetHistFillColor(1);
    # ROOT . gStyle . SetHistFillStyle(0);
    ROOT . gStyle . SetHistLineColor(1);
    ROOT . gStyle . SetHistLineStyle(0);
    ROOT . gStyle . SetHistLineWidth(1);
    # ROOT . gStyle . SetLegoInnerR(Float_t rad = 0.5);
    # ROOT . gStyle . SetNumberContours(Int_t number = 20);

    ROOT . gStyle . SetEndErrorSize(2);
    #ROOT . gStyle . SetErrorMarker(20);   #/ I COMMENTED THIS OUT
    #ROOT . gStyle . SetErrorX(0.);

    #ROOT . gStyle . SetMarkerStyle(20);


#For the fit/function:
    ROOT . gStyle . SetOptFit(1011);
    ROOT . gStyle . SetFitFormat("5.4g");
    ROOT . gStyle . SetFuncColor(2);
    ROOT . gStyle . SetFuncStyle(1);
    ROOT . gStyle . SetFuncWidth(1);

#For the date:
    ROOT . gStyle . SetOptDate(0);
    # ROOT . gStyle . SetDateX(Float_t x = 0.01);
    # ROOT . gStyle . SetDateY(Float_t y = 0.01);

# For the statistics box:
    ROOT . gStyle . SetOptFile(0);
    ROOT . gStyle . SetOptStat(0); # To display the mean and RMS:   SetOptStat("mr");
    ROOT . gStyle . SetStatColor(0); # kWhite
    ROOT . gStyle . SetStatFont(42);
    #ROOT . gStyle . SetStatFontSize(0.025);
    ROOT . gStyle . SetStatFontSize(0.04);
    ROOT . gStyle . SetStatTextColor(1);
    ROOT . gStyle . SetStatFormat("6.4g");
    ROOT . gStyle . SetStatBorderSize(1);
    ROOT . gStyle . SetStatH(0.1);
    ROOT . gStyle . SetStatW(0.15);
    # ROOT . gStyle . SetStatStyle(Style_t style = 1001);
    # ROOT . gStyle . SetStatX(Float_t x = 0);
    # ROOT . gStyle . SetStatY(Float_t y = 0);

# Margins:
    ROOT . gStyle . SetPadTopMargin(0.07);
    ROOT . gStyle . SetPadBottomMargin(0.13);
    ROOT . gStyle . SetPadLeftMargin(0.16);
    #ROOT . gStyle . SetPadRightMargin(0.12);
    ROOT . gStyle . SetPadRightMargin(0.03);

# For the Global title:

    ROOT . gStyle . SetOptTitle(0);
    ROOT . gStyle . SetTitleFont(42);
    ROOT . gStyle . SetTitleColor(1);
    ROOT . gStyle . SetTitleTextColor(1);
    ROOT . gStyle . SetTitleFillColor(10);
    ROOT . gStyle . SetTitleFontSize(0.05);
    # ROOT . gStyle . SetTitleH(0); # Set the height of the title box
    # ROOT . gStyle . SetTitleW(0); # Set the width of the title box
    # ROOT . gStyle . SetTitleX(0); # Set the position of the title box
    # ROOT . gStyle . SetTitleY(0.985); # Set the position of the title box
    # ROOT . gStyle . SetTitleStyle(Style_t style = 1001);
    # ROOT . gStyle . SetTitleBorderSize(2);

# For the axis titles:

    ROOT . gStyle . SetTitleColor(1, "XYZ");
    ROOT . gStyle . SetTitleFont(42, "XYZ");
    ROOT . gStyle . SetTitleSize(0.06, "XYZ");
    # ROOT . gStyle . SetTitleXSize(Float_t size = 0.02); # Another way to set the size?
    # ROOT . gStyle . SetTitleYSize(Float_t size = 0.02);
    ROOT . gStyle . SetTitleXOffset(0.9);
    ROOT . gStyle . SetTitleYOffset(1.25);
    # ROOT . gStyle . SetTitleOffset(1.1, "Y"); # Another way to set the Offset

# For the axis labels:

    ROOT . gStyle . SetLabelColor(1, "XYZ");
    ROOT . gStyle . SetLabelFont(42, "XYZ");
    ROOT . gStyle . SetLabelOffset(0.007, "XYZ");
    ROOT . gStyle . SetLabelSize(0.05, "XYZ");

# For the axis:

    ROOT . gStyle . SetAxisColor(1, "XYZ");
    ROOT . gStyle . SetStripDecimals(1); # kTRUE
    ROOT . gStyle . SetTickLength(0.03, "XYZ");
    ROOT . gStyle . SetNdivisions(510, "XYZ");
    ROOT . gStyle . SetPadTickX(1);  # To get tick marks on the opposite side of the frame
    ROOT . gStyle . SetPadTickY(1);

# Change for log plots:
    ROOT . gStyle . SetOptLogx(0);
    ROOT . gStyle . SetOptLogy(0);
    ROOT . gStyle . SetOptLogz(0);

# Postscript options:
    ROOT . gStyle . SetPaperSize(20.,20.);
    # ROOT . gStyle . SetLineScalePS(Float_t scale = 3);
    # ROOT . gStyle . SetLineStyleString(Int_t i, const char* text);
    # ROOT . gStyle . SetHeaderPS(const char* header);
    # ROOT . gStyle . SetTitlePS(const char* pstitle);

    # ROOT . gStyle . SetBarOffset(Float_t baroff = 0.5);
    # ROOT . gStyle . SetBarWidth(Float_t barwidth = 0.5);
    # ROOT . gStyle . SetPaintTextFormat(const char* format = "g");
    # ROOT . gStyle . SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    # ROOT . gStyle . SetTimeOffset(Double_t toffset);
    # ROOT . gStyle . SetHistMinimumZero(kTRUE);

setTDRStyle()

def plotIt(channel,mass,btags,inputDir):

    
    #############################################################

    if channel=='mu': 
        inputFile=inputDir+'/muon_BDT_wp'+mass+'R_Histos-'+btags+'.root'
        #inputFile=inputDir+'/muon_BDT_wp'+mass+'R_Histos-'+btags+'_BDT0.0.root'
    elif channel=='elec': 
        inputFile=inputDir+'/electron_BDT_wp'+mass+'R_Histos-'+btags+'.root'
        #inputFile=inputDir+'/electron_BDT_wp'+mass+'R_Histos-'+btags+'_BDT0.0.root'

    print inputFile

    systematics=['jes__plus',
                 'jes__minus',
                 'jer__plus',
                 'jer__minus',
                 'btag__plus',
                 'btag__minus',
                 'hf__plus',
                 'hf__minus',
                 'shape__plus',
                 'shape__minus']

    #if (channel=='elec'): samples=['wjets','ttbar','qcd']
    if (channel=='elec'): samples=['wjets','ttbar']
    if (channel=='mu'): samples=['wjets','ttbar']
    #############################################################

    input=TFile(inputFile)

    data=input.Get(channel+'_BDT_'+btags+'__DATA').Clone('data')
    
    
    wjets=input.Get(channel+'_BDT_'+btags+'__wjets').Clone("wjets")
    wjets.SetLineColor(1)
    wjets.SetLineWidth(2)
    wjets.SetFillColor(kGreen-3)
    ttbar=input.Get(channel+'_BDT_'+btags+'__ttbar').Clone("ttbar")
    ttbar.SetLineColor(1)
    ttbar.SetLineWidth(2)
    ttbar.SetFillColor(kRed-7)

    #if (channel=='elec'):
    #    qcd=input.Get(channel+'_BDT_'+btags+'__qcd').Clone("qcd")
    #    qcd.SetLineColor(1)
    #    qcd.SetLineWidth(2)
    #    qcd.SetFillColor(kYellow)
                
    background=input.Get(channel+'_BDT_'+btags+'__ttbar').Clone("background")
    background.Add(input.Get(channel+'_BDT_'+btags+'__wjets'))
    #if (channel == 'elec'): background.Add(qcd)
    uncertainty=background.Clone("uncertainty")
    pull=background.Clone('pull')

    stack = THStack('stack', 'stack')
    #if (channel == 'elec'): stack.Add(qcd)
    stack.Add(wjets)
    stack.Add(ttbar)

    signal=input.Get(channel+'_BDT_'+btags+'__wp'+mass+'R').Clone("signal")
    signal.SetLineColor(1)
    signal.SetLineWidth(2)
    signal.SetFillColor(kAzure-2)
    stack.Add(signal)

    for sample in samples:
        nominal=input.Get(channel+'_BDT_'+btags+'__'+sample)
        for systematic in systematics:
            diff=nominal.Clone("diff")

            name=channel+'_BDT_'+btags+'__'+sample+'__'+systematic
            hist=input.Get(name)
            if hist:
                diff.Add(input.Get(name),-1)
                print systematic
                for bin in range(1,uncertainty.GetNbinsX()+1):
                    if (sample != 'qcd'): uncertainty.SetBinError(bin,sqrt(uncertainty.GetBinError(bin)**2+diff.GetBinContent(bin)**2))
                    if (channel=='elec' and sample=='qcd' and systematic=='btag__minus'): uncertainty.SetBinError(bin,sqrt(uncertainty.GetBinError(bin)**2+1.0*qcd.GetBinContent(bin)**2))
                    if (uncertainty.GetBinError(bin)>0):
                        pull.SetBinContent(bin,(data.GetBinContent(bin)-background.GetBinContent(bin))/uncertainty.GetBinError(bin))
                    else:
                        pull.SetBinContent(bin,0)
                    #pull.SetBinContent(bin,0) ### don't show pulls
                        
    data.SetLineColor(kBlack)
    data.SetMarkerStyle(20)
    background.SetFillColor(kRed-7)
    background.SetLineColor(kRed-7)
    background.SetTitle('')
    uncertainty.SetFillStyle(3344)
    uncertainty.SetFillColor(kBlack)
    uncertainty.SetLineColor(kBlack)
    uncertainty.SetMarkerSize(0)
    gStyle.SetHatchesLineWidth(1)
    pull.SetTitle("")
    pull.GetYaxis().SetTitle('#sigma(data-bkgd)')

    background.GetXaxis().SetLabelSize(0)
    background.GetYaxis().CenterTitle()
    background.GetYaxis().SetLabelSize(0.08)
    background.GetYaxis().SetTitleSize(0.12)
    background.GetYaxis().SetTitleOffset(.75)

    pull.SetLineColor(kBlack)
    pull.SetFillColor(kRed)
    pull.SetLineWidth(2)
    pull.GetXaxis().SetLabelSize(.15)
    pull.GetXaxis().SetTitleSize(0.18)
    pull.GetXaxis().SetTitleOffset(0.95)
    pull.GetXaxis().SetTitle('BDT Disciminant')

    pull.GetYaxis().SetLabelSize(0.125)
    pull.GetYaxis().SetTitleSize(0.1)
    pull.GetYaxis().SetTitleOffset(.55)
    pull.GetYaxis().SetNdivisions(5);
    pull.SetMaximum(3)
    pull.SetMinimum(-3)

    prelimTex=TLatex()
    prelimTex.SetNDC()
    prelimTex.SetTextSize(0.04)
    prelimTex.SetTextAlign(31) # align right
    lumi=19.6
    prelimTex.DrawLatex(0.9, 0.95, "CMS Preliminary, "+str(lumi)+" fb^{-1} at #sqrt{s} = 8 TeV");

    channelTex = TLatex()
    channelTex.SetNDC()
    channelTex.SetTextSize(0.06)
    channelTex.SetTextAlign(31)
    if channel=='elec': text='e'
    elif channel=='mu': text='#mu'
    if (btags == 'ge1btags'): text+="+jets N_{b-tags} #geq 1"
    if (btags == 'onebtags'): text+="+jets N_{b-tags} = 1"
    if (btags == 'twobtags'): text+="+jets N_{b-tags} = 2"
    channelTex.DrawLatex(0.5, 0.83, text);
                                         
    legend=TLegend(0.55,0.70,0.90,0.90)
    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.AddEntry(data,"Data")
    legend.AddEntry(signal,"W'_{R} m = "+mass+" GeV, #sigma = 1.0 pb", "F")
    legend.AddEntry(ttbar,"t#bar{t} + Single-Top","F")
    legend.AddEntry(wjets,"W#rightarrowl#nu + Z/#gamma*#rightarrowl^{+}l^{-} + VV","F")
    #if (channel == 'elec'): legend.AddEntry(qcd,"QCD","F")
    legend.AddEntry(uncertainty,"Stat #oplus Sys Uncertainty","F")
        
    c=TCanvas("","",1000,800)

    yDiv=0.35
    uPad=TPad("uPad","",0,yDiv,1,1) #for actual plots
    uPad.SetTopMargin(0.07)
    uPad.SetBottomMargin(0)
    uPad.SetRightMargin(.05)
    uPad.SetLeftMargin(.18)
    uPad.SetLogy()
    uPad.Draw()

    lPad=TPad("lPad","",0,0,1,yDiv) #for sigma runner
    lPad.SetTopMargin(0)
    lPad.SetBottomMargin(.4)
    lPad.SetRightMargin(.05)
    lPad.SetLeftMargin(.18)
    lPad.SetGridy()
    lPad.Draw()

    uPad.cd()
    stack.SetMaximum( 20*stack.GetMaximum() )
    stack.SetMinimum( 0.5 )    
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("Events")
    data.Draw('SAME')
    uncertainty.Draw('SAME E2')
    legend.Draw("SAME")
    prelimTex.DrawLatex(0.9, 0.95, "CMS Preliminary, "+str(lumi)+" fb^{-1} at #sqrt{s} = 8 TeV");
    channelTex.DrawLatex(0.5, 0.83, text);
    #c.SaveAs(channel+'_BDT_wp'+mass+'R_'+btags+'_shapeSys.pdf')

    lPad.cd()
    pull.Draw("HIST")

    c.SaveAs(inputDir+'/'+channel+'_BDT_wp'+mass+'R_'+btags+'_shapeSys.pdf')
    #c.SaveAs(inputDir+'/'+channel+'_BDT_wp'+mass+'R_'+btags+'_shapeSys_bdtCut0pt0.pdf')
    



inputDir='/uscms_data/d2/dsperka/8TeV/MakeTBntuples/29Jan/CMSSW_5_3_6/src/BDT/wprimetb_BDT/RootFiles_ForBDTLimits_1and2_ScaleWjetsForQCD_AllVars/'    

if __name__=='__main__':
    #plotIt('elec','800','ge1btags',inputDir)
    #plotIt('elec','1000','ge1btags',inputDir)
    #plotIt('elec','1200','ge1btags',inputDir)
    #plotIt('elec','1400','ge1btags',inputDir)
    #plotIt('elec','1600','ge1btags',inputDir)
    #plotIt('elec','1800','ge1btags',inputDir)
    #plotIt('elec','2000','ge1btags',inputDir)
    #plotIt('elec','2200','ge1btags',inputDir)
    #plotIt('elec','2400','ge1btags',inputDir)
    #plotIt('elec','2600','ge1btags',inputDir)
    #plotIt('elec','2800','ge1btags',inputDir)
    #plotIt('elec','3000','ge1btags',inputDir)
             
    #plotIt('mu','800','ge1btags',inputDir)
    #plotIt('mu','1000','ge1btags',inputDir)
    #plotIt('mu','1200','ge1btags',inputDir)
    #plotIt('mu','1400','ge1btags',inputDir)
    #plotIt('mu','1600','ge1btags',inputDir)
    #plotIt('mu','1800','ge1btags',inputDir)
    #plotIt('mu','2000','ge1btags',inputDir)
    #plotIt('mu','2200','ge1btags',inputDir)
    #plotIt('mu','2400','ge1btags',inputDir)
    #plotIt('mu','2600','ge1btags',inputDir)
    #plotIt('mu','2800','ge1btags',inputDir)
    #plotIt('mu','3000','ge1btags',inputDir)            

    plotIt('elec','800','onebtags',inputDir)
    plotIt('elec','1000','onebtags',inputDir)
    plotIt('elec','1200','onebtags',inputDir)
    plotIt('elec','1400','onebtags',inputDir)
    plotIt('elec','1600','onebtags',inputDir)
    plotIt('elec','1800','onebtags',inputDir)
    plotIt('elec','2000','onebtags',inputDir)
    plotIt('elec','2200','onebtags',inputDir)
    plotIt('elec','2400','onebtags',inputDir)
    plotIt('elec','2600','onebtags',inputDir)
    plotIt('elec','2800','onebtags',inputDir)
    plotIt('elec','3000','onebtags',inputDir)
            
    plotIt('mu','800','onebtags',inputDir)
    plotIt('mu','1000','onebtags',inputDir)
    plotIt('mu','1200','onebtags',inputDir)
    plotIt('mu','1400','onebtags',inputDir)
    plotIt('mu','1600','onebtags',inputDir)
    plotIt('mu','1800','onebtags',inputDir)
    plotIt('mu','2000','onebtags',inputDir)
    plotIt('mu','2200','onebtags',inputDir)
    plotIt('mu','2400','onebtags',inputDir)    
    plotIt('mu','2600','onebtags',inputDir)
    plotIt('mu','2800','onebtags',inputDir)
    plotIt('mu','3000','onebtags',inputDir)
           
    plotIt('elec','800','twobtags',inputDir)
    plotIt('elec','1000','twobtags',inputDir)
    plotIt('elec','1200','twobtags',inputDir)
    plotIt('elec','1400','twobtags',inputDir)
    plotIt('elec','1600','twobtags',inputDir)
    plotIt('elec','1800','twobtags',inputDir)
    plotIt('elec','2000','twobtags',inputDir)
    plotIt('elec','2200','twobtags',inputDir)
    plotIt('elec','2400','twobtags',inputDir)
    plotIt('elec','2600','twobtags',inputDir)
    plotIt('elec','2800','twobtags',inputDir)
    plotIt('elec','3000','twobtags',inputDir)            

    plotIt('mu','800','twobtags',inputDir)
    plotIt('mu','1000','twobtags',inputDir)
    plotIt('mu','1200','twobtags',inputDir)
    plotIt('mu','1400','twobtags',inputDir)
    plotIt('mu','1600','twobtags',inputDir)
    plotIt('mu','1800','twobtags',inputDir)
    plotIt('mu','2000','twobtags',inputDir)
    plotIt('mu','2200','twobtags',inputDir)
    plotIt('mu','2400','twobtags',inputDir)
    plotIt('mu','2600','twobtags',inputDir)
    plotIt('mu','2800','twobtags',inputDir)
    plotIt('mu','3000','twobtags',inputDir)
           
           
