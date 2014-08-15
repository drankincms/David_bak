//
// PileUp Utility
//
// Zaixing Mao  
// 

//#include <cmath>
#include "TVector3.h"
#include "LumiReweightingStandAlone.h"

namespace PU{

//**************************************************
//*******            USER OPTIONS          *********
//**************************************************
    
// Distribution used for Summer2012 MC.
// https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
float MCDist_Summer2012_S10[60] = {
    2.560E-06,
    5.239E-06,
    1.420E-05,
    5.005E-05,
    1.001E-04,
    2.705E-04,
    1.999E-03,
    6.097E-03,
    1.046E-02,
    1.383E-02,
    1.685E-02,
    2.055E-02,
    2.572E-02,
    3.262E-02,
    4.121E-02,
    4.977E-02,
    5.539E-02,
    5.725E-02,
    5.607E-02,
    5.312E-02,
    5.008E-02,
    4.763E-02,
    4.558E-02,
    4.363E-02,
    4.159E-02,
    3.933E-02,
    3.681E-02,
    3.406E-02,
    3.116E-02,
    2.818E-02,
    2.519E-02,
    2.226E-02,
    1.946E-02,
    1.682E-02,
    1.437E-02,
    1.215E-02,
    1.016E-02,
    8.400E-03,
    6.873E-03,
    5.564E-03,
    4.457E-03,
    3.533E-03,
    2.772E-03,
    2.154E-03,
    1.656E-03,
    1.261E-03,
    9.513E-04,
    7.107E-04,
    5.259E-04,
    3.856E-04,
    2.801E-04,
    2.017E-04,
    1.439E-04,
    1.017E-04,
    7.126E-05,
    4.948E-05,
    3.405E-05,
    2.322E-05,
    1.570E-05,
    5.005E-06
};

// User generated Data distribution
float DataDist_Oct2012[60] = {
    12238.2,
    32262.2,
    88488.8,
    225526,
    487946,
    2.47713e+06,
    1.4766e+07,
    4.44375e+07,
    1.0279e+08,
    1.95543e+08,
    3.33172e+08,
    5.09762e+08,
    6.45795e+08,
    7.16998e+08,
    7.62148e+08,
    8.00239e+08,
    8.22388e+08,
    8.21958e+08,
    8.0435e+08,
    7.75997e+08,
    7.41057e+08,
    7.02468e+08,
    6.61859e+08,
    6.17413e+08,
    5.64036e+08,
    4.97929e+08,
    4.20604e+08,
    3.37939e+08,
    2.56828e+08,
    1.83743e+08,
    1.23728e+08,
    7.88409e+07,
    4.78733e+07,
    2.77934e+07,
    1.53873e+07,
    8.07166e+06,
    3.98747e+06,
    1.85096e+06,
    810013,
    337106,
    135102,
    52853.3,
    20418,
    7846.44,
    3007.4,
    1148.47,
    435.702,
    163.72,
    60.8208,
    22.3317,
    8.11182,
    2.91896,
    1.0414,
    0.368285,
    0.128918,
    0.0445702,
    0.0151809,
    0.00508254,
    0.00166966,
    0.00053755
};

//**************************************************

//___Define Variables
std::vector< float > DataDist;
std::vector< float > MCDist;

reweight::LumiReWeighting lw;

    
//___Set the LumiWeights___
void SetDist(){
    for(int i=0; i<60; ++i) {
        DataDist.push_back(DataDist_Oct2012[i]);
        MCDist.push_back(MCDist_Summer2012_S10[i]);
    }
    lw = reweight::LumiReWeighting(MCDist, DataDist);
}

//___This Function Takes in the nTrueInteractions and Returns the Recalculated Weight Value___
double LumiReWeights(double nTrue){
    float fnTrue = (float)nTrue;
    return lw.weight(fnTrue);
  }
}
