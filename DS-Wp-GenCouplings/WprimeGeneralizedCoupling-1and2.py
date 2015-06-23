import ROOT, sys, os, re, string
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2F,TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad



#xsec_RT = {'800':2.3022,'900':1.3818,'1000':0.85538,'1100':0.54325,'1200':0.35203,'1300':0.23219,'1400':0.15547,'1500':0.10518,'1600':0.072012,'1700':0.049683,'1800':0.034576,'1900':0.024249,'2000':0.017124,'2100':0.012176,'2200':0.0087191,'2300':0.0062918,'2400':0.0045757,'2500':0.0033568,'2600':0.0024870,'2700':0.0018624,'2800':0.0014102,'2900':0.0010818,'3000':0.00084115}


#xsec_LT = {'800':3.1080,'900':2.2731,'1000':1.8087,'1100':1.547,'1200':1.3870,'1300':1.2945,'1400':1.2390,'1500':1.2061,'1600':1.1869,'1700':1.1761,'1800':1.1705,'1900':1.1678,'2000':1.1673,'2100':1.1680,'2200':1.1692,'2300':1.1711,'2400':1.1727,'2500':1.1746,'2600':1.1763,'2700':1.1780,'2800':1.1797,'2900':1.1810,'3000':1.1825}

#xsec_MT = {'800':5.4166,'900':3.6684,'1000':2.6815,'1100':2.1031,'1200':1.7539,'1300':1.5389,'1400':1.4043,'1500':1.3194,'1600':1.2650,'1700':1.2305,'1800':1.2090,'1900':1.1954,'2000':1.1872,'2100':1.1824,'2200':1.1798,'2300':1.1787,'2400':1.1784,'2500':1.1791,'2600':1.1792,'2700':1.1803,'2800':1.1813,'2900':1.1825,'3000':1.1835}

xsec_RT = {'800':2.3022,'900':1.3818,'1000':0.85538,'1100':0.54325,'1200':0.35203,'1300':0.23219,'1400':0.15547,'1500':0.10518,'1600':0.072012,'1700':0.049683,'1800':0.034576,'1900':0.024249,'2000':0.017124,'2100':0.012176,'2200':0.0087191,'2300':0.0062918,'2400':0.0045757,'2500':0.0033568,'2700':0.0018624,'2800':0.0014102,'2900':0.0010818}

xsec_LT = {'800':3.1080,'900':2.2731,'1000':1.8087,'1100':1.547,'1200':1.3870,'1300':1.2945,'1400':1.2390,'1500':1.2061,'1600':1.1869,'1700':1.1761,'1800':1.1705,'1900':1.1678,'2000':1.1673,'2100':1.1680,'2200':1.1692,'2300':1.1711,'2400':1.1727,'2500':1.1746,'2700':1.1780,'2800':1.1797,'2900':1.1810}


xsec_MT = {'800':5.4166,'900':3.6684,'1000':2.6815,'1100':2.1031,'1200':1.7539,'1300':1.5389,'1400':1.4043,'1500':1.3194,'1600':1.2650,'1700':1.2305,'1800':1.2090,'1900':1.1954,'2000':1.1872,'2100':1.1824,'2200':1.1798,'2300':1.1787,'2400':1.1784,'2500':1.1791,'2700':1.1803,'2800':1.1813,'2900':1.1825}


xsec_SM = 5.55
#Lumi=3556
Lumi=1

ListWprimeMass = ['wp800','wp900','wp1000','wp1100','wp1200','wp1300','wp1400','wp1500','wp1600','wp1700','wp1800','wp1900','wp2000','wp2100','wp2200','wp2300','wp2400','wp2500','wp2700','wp2800','wp2900']

RootFiles = {}

inDir = "/uscms_data/d2/dsperka/8TeV/MakeTBntuples/29Jan/CMSSW_5_3_6/src/UserCode/dsperka/wprimetb/RootFiles_For2DLimits_23Sep/"

RootFiles['Right_elec'] = TFile(inDir+"electron_BestJetJet2W_M_WprimeModRight_Histos-onetwobtags_85_140_dr03_lep50.root")
RootFiles['Left_elec'] =  TFile(inDir+"electron_BestJetJet2W_M_WprimeLeft_Histos-onetwobtags_85_140_dr03_lep50.root")
RootFiles['Mix_elec'] =   TFile(inDir+"electron_BestJetJet2W_M_WprimeMix_Histos-onetwobtags_85_140_dr03_lep50.root")

RootFiles['Right_mu'] = TFile(inDir+"muon_BestJetJet2W_M_WprimeModRight_Histos-onetwobtags_85_140_dr03_lep50.root")
RootFiles['Left_mu'] =  TFile(inDir+"muon_BestJetJet2W_M_WprimeLeft_Histos-onetwobtags_85_140_dr03_lep50.root")
RootFiles['Mix_mu'] =   TFile(inDir+"muon_BestJetJet2W_M_WprimeMix_Histos-onetwobtags_85_140_dr03_lep50.root")



def makeHistos_al_ar(als, ars, lepton):
       
    CombinedHist = {}
    keys = []

    al_step=als
    ar_step=ars
    al= 0.10000*al_step
    ar= 0.10000*ar_step
    print  ' al =', al, ' ar =',ar
    print '----'

    prefix = 'Wprime_Histos_'
    suffix = '.root'
    mid = '-'
    fileName= prefix + lepton + mid + `al_step` + mid  +`ar_step` + suffix
    print ' output to FILE', fileName
    f = TFile(fileName, "RECREATE")

    rt_wgt=(-al*al + al*al*al*al + ar*ar*ar*ar- al*al*ar*ar)
    lt_wgt=(al*al - al*al*ar*ar)
    st_wgt=(1-al*al) 
    mt_wgt=(al*al*ar*ar)
    print ' ar,al,mix,sm weights ', rt_wgt, ' ', lt_wgt, ' ', mt_wgt, ' ', st_wgt
    #if (rt_wgt <0) : print ' NEGATIVE WEIGHT R', al,' ', ar
    #if (lt_wgt <0) : print ' NEGATIVE WEIGHT L', al,' ', ar
    #if (mt_wgt <0) : print ' NEGATIVE WEIGHT M', al,' ', ar
    #if (st_wgt <0) : print ' NEGATIVE WEIGHT S', al,' ', ar
 


    for key in RootFiles['Mix_'+lepton].GetListOfKeys():

        if keys.count( str(key.GetName()) ) == 1:
            print ' already used ', str(key.GetName()) 
            continue
        keys . append( str(key.GetName()) )
        print 'Name / Title = ',key.GetName(),' / ',key.GetTitle()

        hist = {}
        name = str(key.GetName())
        print  ' Name ', name
        
        hist['Right'] = RootFiles['Right_'+lepton].Get( str(key.GetName()) )
        hist['Left'] =  RootFiles['Left_'+lepton].Get(  str(key.GetName()) )
        hist['Mix'] =   RootFiles['Mix_'+lepton].Get(   str(key.GetName()) )

        isOneBtag = name.find("onebtags")

        if (name.find("jes__plus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__jes__plus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__jes__plus" )
        elif  (name.find("jes__minus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__jes__minus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__jes__minus" )
        elif  (name.find("btag__plus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__btag__plus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__btag__plus" )
        elif  (name.find("btag__minus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__btag__minus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__btag__minus" )
        elif  (name.find("jer__plus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__jer__plus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__jer__plus" )
        elif  (name.find("jer__minus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__jer__minus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__jer__minus" )
        elif  (name.find("PU__plus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__PU__plus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__PU__plus" )
        elif  (name.find("PU__minus") > 0) :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb__PU__minus" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb__PU__minus" )
        else :
            if (isOneBtag>0): hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_onebtags__topstb" )
            else: hist['SM']=     RootFiles['Right_'+lepton].Get(lepton+"_invmass_twobtags__topstb" )
            
        isThereWprime = name.find("wp")
        #print ' wp is there? ',isThere


        for Mass in ListWprimeMass:
            wp = name.find(Mass)
            #print ' wp mass ', wp
            if (wp>0) :
                rt_xsec = xsec_RT[Mass.lstrip('wp')]
                lt_xsec = xsec_LT[Mass.lstrip('wp')]
                mt_xsec = xsec_MT[Mass.lstrip('wp')]
                st_xsec = xsec_SM
                print ' xsec ', Mass ,' ', rt_xsec ,' ', lt_xsec ,' ', mt_xsec ,' ', st_xsec
                
        if ( isThereWprime >0 ) :
            #print 'Name / Title = ',key.GetName(),' / ',key.GetTitle()
            #print 'Hist = ', hist['SM1b'].Integral(), ' ',   hist['Right'].Integral(), ' ',  hist['Left'].Integral(), ' ', hist['Mix'].Integral()
            NewHist = hist['SM']
            print ' NewHist = ', NewHist.Integral()
            #ast_wgt=(st_wgt*st_xsec*Lumi)-1
            ast_wgt=st_wgt-1
            #ast_wgt=st_wgt
            NewHist.Add ( hist['SM'], ast_wgt )
            #print ' Hist = ', NewHist.Integral()
            NewHist.Add ( hist['Right'], rt_wgt*Lumi*rt_xsec)
            #print ' Hist = ', NewHist.Integral()
            NewHist.Add ( hist['Mix'], mt_wgt*Lumi*mt_xsec )
            #print ' Hist = ', NewHist.Integral()
            NewHist.Add ( hist['Left'], lt_wgt*Lumi*lt_xsec  )
            #print ' Hist = ', NewHist.Integral()
            #NewHist.SetName( str(key.GetTitle()) )
            NewHist.SetName( str(hist['Mix'].GetTitle()) )
            NewHist.SetTitle( str(hist['Mix'].GetTitle()) )

            print ' ar,al,mix,sm weights ', rt_wgt, ' ', lt_wgt, ' ', mt_wgt, ' ', st_wgt
            print 'Combined Hist = ',NewHist.Integral(), ' ',hist['SM'].Integral(), ' ',hist['Right'].Integral(),' ',hist['Left'].Integral(),' ', hist['Mix'].Integral()

            nbins = NewHist.GetXaxis().GetNbins()
            #print "Number of bins    ", str(nbins)
            ntotzero=0
            anegContent=0.
            allContent=0.
#            print ' inits check ',ntotzero, ' ', anegContent, ' ', allContent
            allContent=(NewHist.Integral())
            for j in range(1,nbins+1):
                binc = NewHist.GetBinContent(j)
                if (binc < 0.): 
                    #print "Negative bin content  ",str(j)," ",binc, ' ', hist['Right'].GetBinContent(j), ' ', hist['Left'].GetBinContent(j), ' ', hist['Mix'].GetBinContent(j), ' ',  rt_wgt*hist['Right'].GetBinContent(j), ' ',lt_wgt* hist['Left'].GetBinContent(j), ' ', mt_wgt*hist['Mix'].GetBinContent(j)
                    ntotzero+=1
                    anegContent+=binc
                    xa=0.0000001
                    NewHist.SetBinContent(j,xa)
                    binc = NewHist.GetBinContent(j)
                    #print "set bin content to zero ",str(j)," ",binc
                    #if (binc !=0.0000001) : print 'something went wrong in setting the bin content to zero'
            if (ntotzero>0):
                print  ' Name ', name
                print 'numbins negative ',ntotzero, '  ', anegContent, '  ', allContent
                zerofrac=anegContent/(allContent)
                if (zerofrac < -0.02 ) :
                    print 'zerofrac ', zerofrac

            CombinedHist[ str(key.GetTitle()) ] = NewHist
            
            del hist
            del NewHist
            
        else :
            CombinedHist[ str(key.GetTitle()) ] =  hist['Mix']
            del hist
            #print ''

            
        CombinedHist[ str(key.GetTitle()) ].Write()
        #print 'Wrote ',str(key.GetTitle())
al=0.
ar=1.
al_step=0
ar_step=0
lepton="elec"
for ar_step in range(0,11):
    for al_step in range(0 , 11):
        print 'al_step = ',al_step, ' ar_step =',ar_step
        makeHistos_al_ar(al_step, ar_step, lepton)
lepton="mu"
for ar_step in range(0,11):
    for al_step in range(0 , 11):
        print 'al_step = ',al_step, ' ar_step =',ar_step
        makeHistos_al_ar(al_step, ar_step, lepton)
                        
