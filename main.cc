#include <iostream>
#include <TCanvas.h>

#include <map>
#include "ECAL.h"
#include "EMShower.h" 
#include "IncGamma.h" 

//#include "EMECALShowerParametrization.h"
#include "GammaFunctionGenerator.h"

#include <fstream>
using namespace std;

int main(){

bool bFixedLength=true;
int nPart=1;
GammaFunctionGenerator* gamma;
std::vector<double> energy_in;
energy_in.push_back(10);
    
ECALProperties *ecalprop;    
EMECALShowerParametrization *myparam = new EMECALShowerParametrization(ecalprop,{100.0,0.1},{1.0,0.1,100.0,1.0},1,1);
ECAL *TheEcal= new ECAL(5,-7.125,7.125,5,-7.125,7.125);   
    
EMShower TheShower(gamma, myparam, TheEcal,bFixedLength,nPart,energy_in);
TheShower.compute();    


    
/*ECALProperties *ecalprop; 
EMECALShowerParametrization *myParam = new EMECALShowerParametrization(ecalprop,{100.0,0.1},{1.0,0.1,100.0,1.0},1,1);
    
    double theMeanT = myParam->meanT(lny);
    double theMeanAlpha = myParam->meanAlpha(lny);
    double theMeanLnT = myParam->meanLnT(lny);
    double meanLnAlpha = myParam->meanLnAlpha(lny);
    double sigmaLnT = myParam->sigmaLnT(lny);
    double sigmaLnAlpha = myParam->sigmaLnAlpha(lny);
    // The correlation matrix
    double theCorrelation = myParam->correlationAlphaT(lny);
    double rhop = std::sqrt((1. + theCorrelation) / 2.);
    double rhom = std::sqrt((1. - theCorrelation) / 2.);
      
      
    TGraph* hist_theMeanLnT = new TGraph(11);
    hist_theMeanLnT->SetPoint(0,10,myParam->meanLnT(std::log(10)));
    hist_theMeanLnT->SetPoint(1,50,myParam->meanLnT(std::log(50)));
    hist_theMeanLnT->SetPoint(2,100,myParam->meanLnT(std::log(100)));
    hist_theMeanLnT->SetPoint(3,250,myParam->meanLnT(std::log(250)));
    hist_theMeanLnT->SetPoint(4,500,myParam->meanLnT(std::log(500)));
    hist_theMeanLnT->SetPoint(5,750,myParam->meanLnT(std::log(750)));
    hist_theMeanLnT->SetPoint(6,1000,myParam->meanLnT(std::log(1000)));
    hist_theMeanLnT->SetPoint(7,2500,myParam->meanLnT(std::log(2500)));
    hist_theMeanLnT->SetPoint(8,5000,myParam->meanLnT(std::log(5000)));
    hist_theMeanLnT->SetPoint(9,7500,myParam->meanLnT(std::log(7500)));
    hist_theMeanLnT->SetPoint(10,10000,myParam->meanLnT(std::log(10000)));
      
    TGraph* hist_meanLnAlpha = new TGraph(11);
    
    hist_meanLnAlpha->SetPoint(0,10,myParam->meanLnAlpha(std::log(10)));
    hist_meanLnAlpha->SetPoint(1,50,myParam->meanLnAlpha(std::log(50)));
    hist_meanLnAlpha->SetPoint(2,100,myParam->meanLnAlpha(std::log(100)));
    hist_meanLnAlpha->SetPoint(3,250,myParam->meanLnAlpha(std::log(250)));
    hist_meanLnAlpha->SetPoint(4,500,myParam->meanLnAlpha(std::log(500)));
    hist_meanLnAlpha->SetPoint(5,750,myParam->meanLnAlpha(std::log(750)));
    hist_meanLnAlpha->SetPoint(6,1000,myParam->meanLnAlpha(std::log(1000)));
    hist_meanLnAlpha->SetPoint(7,2500,myParam->meanLnAlpha(std::log(2500)));
    hist_meanLnAlpha->SetPoint(8,5000,myParam->meanLnAlpha(std::log(5000)));
    hist_meanLnAlpha->SetPoint(9,7500,myParam->meanLnAlpha(std::log(7500)));
    hist_meanLnAlpha->SetPoint(10,10000,myParam->meanLnAlpha(std::log(10000)));
    
    TGraph* hist_sigmaLnT = new TGraph(11);
    
    hist_sigmaLnT->SetPoint(0,10,myParam->sigmaLnT(std::log(10)));
    hist_sigmaLnT->SetPoint(1,50,myParam->sigmaLnT(std::log(50)));
    hist_sigmaLnT->SetPoint(2,100,myParam->sigmaLnT(std::log(100)));
    hist_sigmaLnT->SetPoint(3,250,myParam->sigmaLnT(std::log(250)));
    hist_sigmaLnT->SetPoint(4,500,myParam->sigmaLnT(std::log(500)));
    hist_sigmaLnT->SetPoint(5,750,myParam->sigmaLnT(std::log(750)));
    hist_sigmaLnT->SetPoint(6,1000,myParam->sigmaLnT(std::log(1000)));
    hist_sigmaLnT->SetPoint(7,2500,myParam->sigmaLnT(std::log(2500)));
    hist_sigmaLnT->SetPoint(8,5000,myParam->sigmaLnT(std::log(5000)));
    hist_sigmaLnT->SetPoint(9,7500,myParam->sigmaLnT(std::log(7500)));
    hist_sigmaLnT->SetPoint(10,10000,myParam->sigmaLnT(std::log(10000)));
    
    
    TGraph* hist_sigmaLnAlpha = new TGraph(11);
    
    hist_sigmaLnAlpha->SetPoint(0,10,myParam->sigmaLnAlpha(std::log(10)));
    hist_sigmaLnAlpha->SetPoint(1,50,myParam->sigmaLnAlpha(std::log(50)));
    hist_sigmaLnAlpha->SetPoint(2,100,myParam->sigmaLnAlpha(std::log(100)));
    hist_sigmaLnAlpha->SetPoint(3,250,myParam->sigmaLnAlpha(std::log(250)));
    hist_sigmaLnAlpha->SetPoint(4,500,myParam->sigmaLnAlpha(std::log(500)));
    hist_sigmaLnAlpha->SetPoint(5,750,myParam->sigmaLnAlpha(std::log(750)));
    hist_sigmaLnAlpha->SetPoint(6,1000,myParam->sigmaLnAlpha(std::log(1000)));
    hist_sigmaLnAlpha->SetPoint(7,2500,myParam->sigmaLnAlpha(std::log(2500)));
    hist_sigmaLnAlpha->SetPoint(8,5000,myParam->sigmaLnAlpha(std::log(5000)));
    hist_sigmaLnAlpha->SetPoint(9,7500,myParam->sigmaLnAlpha(std::log(7500)));
    hist_sigmaLnAlpha->SetPoint(10,10000,myParam->sigmaLnAlpha(std::log(10000)));
    

        
    TGraph* hist_correlationAlphaT = new TGraph(11);
    
    hist_correlationAlphaT->SetPoint(0,10,myParam->correlationAlphaT(std::log(10)));
    hist_correlationAlphaT->SetPoint(1,50,myParam->correlationAlphaT(std::log(50)));
    hist_correlationAlphaT->SetPoint(2,100,myParam->correlationAlphaT(std::log(100)));
    hist_correlationAlphaT->SetPoint(3,250,myParam->correlationAlphaT(std::log(250)));
    hist_correlationAlphaT->SetPoint(4,500,myParam->correlationAlphaT(std::log(500)));
    hist_correlationAlphaT->SetPoint(5,750,myParam->correlationAlphaT(std::log(750)));
    hist_correlationAlphaT->SetPoint(6,1000,myParam->correlationAlphaT(std::log(1000)));
    hist_correlationAlphaT->SetPoint(7,2500,myParam->correlationAlphaT(std::log(2500)));
    hist_correlationAlphaT->SetPoint(8,5000,myParam->correlationAlphaT(std::log(5000)));
    hist_correlationAlphaT->SetPoint(9,7500,myParam->correlationAlphaT(std::log(7500)));
    hist_correlationAlphaT->SetPoint(10,10000,myParam->correlationAlphaT(std::log(10000)));
        
    
    
TCanvas * par= new TCanvas("en_lat","en_lat",1000,100,2500,2000); 
par->Divide(3,2);   
par->cd(1);
hist_theMeanLnT->SetLineWidth(6);
hist_theMeanLnT->SetMarkerStyle(117);
hist_theMeanLnT->SetMarkerSize(7);
hist_theMeanLnT->SetTitle("mean LnT");
hist_theMeanLnT->Draw("ACP*");
gPad->SetLogx();
par->cd(2);
hist_meanLnAlpha->SetLineWidth(6);
hist_meanLnAlpha->SetMarkerStyle(117);
hist_meanLnAlpha->SetMarkerSize(7);
hist_meanLnAlpha->SetTitle("mean LnAlpha");
hist_meanLnAlpha->Draw("ACP*");
gPad->SetLogx();
par->cd(3);
hist_sigmaLnT->SetLineWidth(6);
hist_sigmaLnT->SetMarkerStyle(117);
hist_sigmaLnT->SetMarkerSize(7);
hist_sigmaLnT->SetTitle("sigma LnT");
hist_sigmaLnT->Draw("ACP*");
gPad->SetLogx();
par->cd(4);
hist_sigmaLnAlpha->SetLineWidth(6);
hist_sigmaLnAlpha->SetMarkerStyle(117);
hist_sigmaLnAlpha->SetMarkerSize(7);
hist_sigmaLnAlpha->SetTitle("sigma LnAlpha");
hist_sigmaLnAlpha->Draw("ACP*");
gPad->SetLogx();
par->cd(5);
hist_correlationAlphaT->SetLineWidth(6);
hist_correlationAlphaT->SetMarkerStyle(117);
hist_correlationAlphaT->SetMarkerSize(7);
hist_correlationAlphaT->SetTitle("correlation Alpha-T");
hist_correlationAlphaT->Draw("ACP*");
gPad->SetLogx();
 
    
    
    
par->SaveAs("/Users/eugenia/desktop/EMCal/param.png");*/
    

/*ECAL *TheEcal= new ECAL(5,-7.125,7.125,5,-7.125,7.125);   
TH2F* EcalGrid=TheEcal->CreateGrid(5,-7.125,7.125,5,-7.125,7.125);
    
TheEcal->AddHitCoo(1,1,0,0,0.045,EcalGrid);
TheEcal->AddHitCoo(0,1,0,0,0.045,EcalGrid);
TheEcal->AddHitCoo(0,1,0,0,1,EcalGrid);
TheEcal->AddHitCoo(0,1,0,0,95,EcalGrid);
TheEcal->AddHitCoo(0,1,0,0,0.00004,EcalGrid);
    


TheEcal->AddHitCoo(3,-1,0,0,9,EcalGrid);
TheEcal->AddHitCoo(3,-1,0,0,0.34,EcalGrid);

    
TheEcal->Draw_ECAL(EcalGrid);*/
//TheEcal.GiveNcell(1,1,EcalGrid);

    
    
/*TCanvas * a= new TCanvas("a","a",1000,100,2500,2000); 
MyIncompleteGamma->Draw("L");
a->SaveAs("/Users/eugenia/desktop/EMCal/a.png");
    
/IncGamma gamma;
TF1* MyIncompleteGamma=gamma.MyGamma();
gamma.Set_a(0.5,MyIncompleteGamma);
double result=MyIncompleteGamma->Eval(0.4);
    cout << result <<endl;
TCanvas * a= new TCanvas("a","a",1000,100,2500,2000); 
MyIncompleteGamma->Draw("L");
a->SaveAs("/Users/eugenia/desktop/EMCal/a.png");   */

return 0;
}