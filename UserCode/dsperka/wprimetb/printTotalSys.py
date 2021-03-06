from ROOT import *
gStyle.SetOptStat(0)

def plotIt(channel,mass,btags,inputDir):

    uncertainty = 0

    #############################################################
    if channel == 'elec':
        inputFile=inputDir+'electron_BestJetJet2W_M_WprimeRight_Histos-'+btags+'_85_140_dr03_lep50.root'
    if channel == 'mu':
        inputFile=inputDir+'muon_BestJetJet2W_M_WprimeRight_Histos-'+btags+'_85_140_dr03_lep50.root'

    systematics=['jes',
                 'jer',
                 'btag',
                 'hf',
                 'shape',
                 'PU',
                 'matching,'
                 'q2scale,'
                 'zerotagshape,'
                 ]

    #############################################################

    input=TFile(inputFile)

    data=input.Get(channel+'_invmass_'+btags+'__DATA').Clone('data')
    #wjets=input.Get(channel+'_invmass_'+btags+'__wjets').Clone("wjets")
    wjets=input.Get(channel+'_invmass_'+btags+'__wjets__btag__plus').Clone("wjets")
    #ttbar=input.Get(channel+'_invmass_'+btags+'__ttbar').Clone("ttbar")
    ttbar=input.Get(channel+'_invmass_'+btags+'__ttbar__btag__plus').Clone("ttbar")
    
    background=input.Get(channel+'_invmass_'+btags+'__ttbar').Clone("background")
    background.Add(input.Get(channel+'_invmass_'+btags+'__wjets'))
    
    signal=input.Get(channel+'_invmass_'+btags+'__wp'+mass).Clone("signal")

    uncertainty = 0

    uncertainty += background.GetEntries() # mc statistical uncertainty
    uncertainty += (0.026*background.Integral())**2 # luminosity
    uncertainty += (0.08*ttbar.Integral())**2 # ttbar rate
    if channel == 'elec':
        uncertainty += (0.03*background.Integral())**2 # id. + trig
    if channel == 'mu':
        uncertainty += (0.02*background.Integral())**2 # id. + trig


    for systematic in systematics:
                
        namewjetsplus=channel+'_invmass_'+btags+'__wjets__'+systematic+'__plus'
        namewjetsminus=channel+'_invmass_'+btags+'__wjets__'+systematic+'__minus'
 
        namettbarplus=channel+'_invmass_'+btags+'__ttbar__'+systematic+'__plus'
        namettbarminus=channel+'_invmass_'+btags+'__ttbar__'+systematic+'__minus'
      
        histwjetsplus=input.Get(namewjetsplus)
        histwjetsminus=input.Get(namewjetsminus)
        if not (histwjetsplus and histwjetsminus):
            histwjetsplus = input.Get(channel+'_invmass_'+btags+'__wjets')
            histwjetsminus = input.Get(channel+'_invmass_'+btags+'__wjets')
 
        histttbarplus=input.Get(namettbarplus)
        histttbarminus=input.Get(namettbarminus)
        if not (histttbarplus and histttbarminus):
            histttbarplus = input.Get(channel+'_invmass_'+btags+'__ttbar')
            histttbarminus = input.Get(channel+'_invmass_'+btags+'__ttbar')
            
        totalplus=histwjetsplus.Clone("totalplus")
        totalplus.Add(histttbarplus)

        totalminus=histwjetsplus.Clone("totalminus")
        totalminus.Add(histttbarminus)

        diff = totalplus.Integral() - totalminus.Integral()
        if systematic == 'PU' or systematic == 'zerotagshape' or systematic == 'shape':
            uncertainty += (diff)**2
        else: 
            uncertainty += (0.5*diff)**2
             
            
    uncertainty = sqrt(uncertainty)

    print channel+' '+btags
    print 'data:      '+str(data.Integral())
    print 'wjets:     '+str(wjets.Integral())
    print 'ttbar:     '+str(ttbar.Integral())
    print 'tot. bkg.: '+str(background.Integral())+' +/- '+str(uncertainty)
    print 'bkg./data: '+str((wjets.Integral()+ttbar.Integral())/data.Integral())

inputDir='/uscms_data/d2/dsperka/8TeV/MakeTBntuples/29Jan/CMSSW_5_3_6/src/UserCode/dsperka/wprimetb/RootFiles_For2DLimits_23Sep_TEST/'

if __name__=='__main__':
    plotIt('elec','800','onebtagspre',inputDir)
    plotIt('elec','800','twobtagspre',inputDir)

    plotIt('elec','800','onebtags',inputDir)
    plotIt('elec','800','twobtags',inputDir)
    
    plotIt('mu','800','onebtagspre',inputDir)
    plotIt('mu','800','twobtagspre',inputDir)
    
    plotIt('mu','800','onebtags',inputDir)
    plotIt('mu','800','twobtags',inputDir)
   
               
