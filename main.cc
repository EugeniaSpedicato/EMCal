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
int nPart;
double X0depth;
GammaFunctionGenerator* gamma= new GammaFunctionGenerator;
std::vector<double> energy_in;
energy_in.push_back(100);
    
int part_type=1; 
    
ECALProperties *ecalprop= new ECALProperties();    
EMECALShowerParametrization *myparam = new EMECALShowerParametrization(ecalprop,{100.0,0.1},{1.0,0.1,100.0,1.0},1,1);
ECAL *TheEcal= new ECAL(5,-7.125,7.125,5,-7.125,7.125);   
    
    
for (int i=0;i<4000;++i) 
{
    if (part_type==1) {nPart=1; X0depth=0; TheEcal->SetEnergy(energy_in[0]);}
    if (part_type==2) {nPart=2; X0depth=-log(gRandom->Uniform())*(9./7.);
                       double energy=energy_in[0]+energy_in[1];
                       TheEcal->SetEnergy(energy);}
    EMShower TheShower(gamma, myparam, TheEcal,bFixedLength,nPart,X0depth,energy_in);
    TheShower.compute();
} 
    

TheEcal->Print_();


    
/*
ECALProperties *ecalprop; 
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

    TGraph* hist_theT = new TGraph(11);
    hist_theT->SetPoint(0,25,myParam->meanT(std::log(25)));
    hist_theT->SetPoint(1,50,myParam->meanT(std::log(50)));
    hist_theT->SetPoint(2,100,myParam->meanT(std::log(100)));
    hist_theT->SetPoint(3,250,myParam->meanT(std::log(250)));
    hist_theT->SetPoint(4,500,myParam->meanT(std::log(500)));
    hist_theT->SetPoint(5,750,myParam->meanT(std::log(750)));
    hist_theT->SetPoint(6,1000,myParam->meanT(std::log(1000)));
    hist_theT->SetPoint(7,2500,myParam->meanT(std::log(2500)));
    hist_theT->SetPoint(8,5000,myParam->meanT(std::log(5000)));
    hist_theT->SetPoint(9,7500,myParam->meanT(std::log(7500)));
    hist_theT->SetPoint(10,15000,myParam->meanT(std::log(15000)));
    
    TCanvas * mT= new TCanvas("meanT","mean T",1000,100,2500,2000); 
    hist_theT->SetLineWidth(2);
    hist_theT->SetMarkerStyle(22);
    hist_theT->SetMarkerSize(4);
    hist_theT->SetTitle("Mean T hom.");
    hist_theT->GetXaxis()->SetTitle("y");
    hist_theT->GetYaxis()->SetTitle("T_hom [X0]");
    hist_theT->SetMaximum(10.0);
    hist_theT->Draw("ACP");
    gPad->SetLogx();
    mT->SaveAs("/Users/eugenia/desktop/EMCal/parameters/theMeanT.png");  
      
    
    TGraph* hist_theAlpha = new TGraph(11);
    hist_theAlpha->SetPoint(0,25,myParam->meanAlpha(std::log(25)));
    hist_theAlpha->SetPoint(1,50,myParam->meanAlpha(std::log(50)));
    hist_theAlpha->SetPoint(2,100,myParam->meanAlpha(std::log(100)));
    hist_theAlpha->SetPoint(3,250,myParam->meanAlpha(std::log(250)));
    hist_theAlpha->SetPoint(4,500,myParam->meanAlpha(std::log(500)));
    hist_theAlpha->SetPoint(5,750,myParam->meanAlpha(std::log(750)));
    hist_theAlpha->SetPoint(6,1000,myParam->meanAlpha(std::log(1000)));
    hist_theAlpha->SetPoint(7,2500,myParam->meanAlpha(std::log(2500)));
    hist_theAlpha->SetPoint(8,5000,myParam->meanAlpha(std::log(5000)));
    hist_theAlpha->SetPoint(9,7500,myParam->meanAlpha(std::log(7500)));
    hist_theAlpha->SetPoint(10,15000,myParam->meanAlpha(std::log(15000)));
    
    TCanvas * mAl= new TCanvas("meanAl","mean Alpha",1000,100,2500,2000); 
    hist_theAlpha->SetLineWidth(2);
    hist_theAlpha->SetMarkerStyle(22);
    hist_theAlpha->SetMarkerSize(4);
    hist_theAlpha->SetTitle("Mean Alpha hom");
    hist_theAlpha->GetXaxis()->SetTitle("y");
    hist_theAlpha->GetYaxis()->SetTitle("Alpha_hom");
    hist_theAlpha->SetMaximum(6.0);
    hist_theAlpha->SetMinimum(2.0);
    hist_theAlpha->Draw("ACP");
    gPad->SetLogx();
    mAl->SaveAs("/Users/eugenia/desktop/EMCal/parameters/theMeanAlpha.png");  
    
    TGraph* hist_theMeanLnT = new TGraph(11);
    hist_theMeanLnT->SetPoint(0,25,myParam->meanLnT(std::log(25)));
    hist_theMeanLnT->SetPoint(1,50,myParam->meanLnT(std::log(50)));
    hist_theMeanLnT->SetPoint(2,100,myParam->meanLnT(std::log(100)));
    hist_theMeanLnT->SetPoint(3,250,myParam->meanLnT(std::log(250)));
    hist_theMeanLnT->SetPoint(4,500,myParam->meanLnT(std::log(500)));
    hist_theMeanLnT->SetPoint(5,750,myParam->meanLnT(std::log(750)));
    hist_theMeanLnT->SetPoint(6,1000,myParam->meanLnT(std::log(1000)));
    hist_theMeanLnT->SetPoint(7,2500,myParam->meanLnT(std::log(2500)));
    hist_theMeanLnT->SetPoint(8,5000,myParam->meanLnT(std::log(5000)));
    hist_theMeanLnT->SetPoint(9,7500,myParam->meanLnT(std::log(7500)));
    hist_theMeanLnT->SetPoint(10,15000,myParam->meanLnT(std::log(15000)));
      
    TGraph* hist_meanLnAlpha = new TGraph(11);
    hist_meanLnAlpha->SetPoint(0,25,myParam->meanLnAlpha(std::log(25)));
    hist_meanLnAlpha->SetPoint(1,50,myParam->meanLnAlpha(std::log(50)));
    hist_meanLnAlpha->SetPoint(2,100,myParam->meanLnAlpha(std::log(100)));
    hist_meanLnAlpha->SetPoint(3,250,myParam->meanLnAlpha(std::log(250)));
    hist_meanLnAlpha->SetPoint(4,500,myParam->meanLnAlpha(std::log(500)));
    hist_meanLnAlpha->SetPoint(5,750,myParam->meanLnAlpha(std::log(750)));
    hist_meanLnAlpha->SetPoint(6,1000,myParam->meanLnAlpha(std::log(1000)));
    hist_meanLnAlpha->SetPoint(7,2500,myParam->meanLnAlpha(std::log(2500)));
    hist_meanLnAlpha->SetPoint(8,5000,myParam->meanLnAlpha(std::log(5000)));
    hist_meanLnAlpha->SetPoint(9,7500,myParam->meanLnAlpha(std::log(7500)));
    hist_meanLnAlpha->SetPoint(10,15000,myParam->meanLnAlpha(std::log(15000)));
    
    TGraph* hist_sigmaLnT = new TGraph(11);
    
    hist_sigmaLnT->SetPoint(0,25,myParam->sigmaLnT(std::log(25)));
    hist_sigmaLnT->SetPoint(1,50,myParam->sigmaLnT(std::log(50)));
    hist_sigmaLnT->SetPoint(2,100,myParam->sigmaLnT(std::log(100)));
    hist_sigmaLnT->SetPoint(3,250,myParam->sigmaLnT(std::log(250)));
    hist_sigmaLnT->SetPoint(4,500,myParam->sigmaLnT(std::log(500)));
    hist_sigmaLnT->SetPoint(5,750,myParam->sigmaLnT(std::log(750)));
    hist_sigmaLnT->SetPoint(6,1000,myParam->sigmaLnT(std::log(1000)));
    hist_sigmaLnT->SetPoint(7,2500,myParam->sigmaLnT(std::log(2500)));
    hist_sigmaLnT->SetPoint(8,5000,myParam->sigmaLnT(std::log(5000)));
    hist_sigmaLnT->SetPoint(9,7500,myParam->sigmaLnT(std::log(7500)));
    hist_sigmaLnT->SetPoint(10,15000,myParam->sigmaLnT(std::log(15000)));
    
    
    TGraph* hist_sigmaLnAlpha = new TGraph(11);
    
    hist_sigmaLnAlpha->SetPoint(0,25,myParam->sigmaLnAlpha(std::log(25)));
    hist_sigmaLnAlpha->SetPoint(1,50,myParam->sigmaLnAlpha(std::log(50)));
    hist_sigmaLnAlpha->SetPoint(2,100,myParam->sigmaLnAlpha(std::log(100)));
    hist_sigmaLnAlpha->SetPoint(3,250,myParam->sigmaLnAlpha(std::log(250)));
    hist_sigmaLnAlpha->SetPoint(4,500,myParam->sigmaLnAlpha(std::log(500)));
    hist_sigmaLnAlpha->SetPoint(5,750,myParam->sigmaLnAlpha(std::log(750)));
    hist_sigmaLnAlpha->SetPoint(6,1000,myParam->sigmaLnAlpha(std::log(1000)));
    hist_sigmaLnAlpha->SetPoint(7,2500,myParam->sigmaLnAlpha(std::log(2500)));
    hist_sigmaLnAlpha->SetPoint(8,5000,myParam->sigmaLnAlpha(std::log(5000)));
    hist_sigmaLnAlpha->SetPoint(9,7500,myParam->sigmaLnAlpha(std::log(7500)));
    hist_sigmaLnAlpha->SetPoint(10,15000,myParam->sigmaLnAlpha(std::log(15000)));
    

        
    TGraph* hist_correlationAlphaT = new TGraph(11);
    
    hist_correlationAlphaT->SetPoint(0,25,myParam->correlationAlphaT(std::log(25)));
    hist_correlationAlphaT->SetPoint(1,50,myParam->correlationAlphaT(std::log(50)));
    hist_correlationAlphaT->SetPoint(2,100,myParam->correlationAlphaT(std::log(100)));
    hist_correlationAlphaT->SetPoint(3,250,myParam->correlationAlphaT(std::log(250)));
    hist_correlationAlphaT->SetPoint(4,500,myParam->correlationAlphaT(std::log(500)));
    hist_correlationAlphaT->SetPoint(5,750,myParam->correlationAlphaT(std::log(750)));
    hist_correlationAlphaT->SetPoint(6,1000,myParam->correlationAlphaT(std::log(1000)));
    hist_correlationAlphaT->SetPoint(7,2500,myParam->correlationAlphaT(std::log(2500)));
    hist_correlationAlphaT->SetPoint(8,5000,myParam->correlationAlphaT(std::log(5000)));
    hist_correlationAlphaT->SetPoint(9,7500,myParam->correlationAlphaT(std::log(7500)));
    hist_correlationAlphaT->SetPoint(10,15000,myParam->correlationAlphaT(std::log(15000)));
        
    
    
TCanvas * par= new TCanvas("en_lat","en_lat",10000,10000,2500,2500); 
par->Divide(2,3);   
par->cd(1);
hist_theMeanLnT->SetLineWidth(2);
hist_theMeanLnT->SetMarkerStyle(22);
hist_theMeanLnT->SetMarkerSize(4);
hist_theMeanLnT->SetTitle("Mean LnT");
    hist_theMeanLnT->GetXaxis()->SetTitle("y");
    hist_theMeanLnT->GetYaxis()->SetTitle("<ln T_hom>");
    hist_theMeanLnT->SetMaximum(2.5);
    hist_theMeanLnT->SetMinimum(0.5);
hist_theMeanLnT->Draw("ACP");
gPad->SetLogx();
par->cd(2);
hist_meanLnAlpha->SetLineWidth(2);
hist_meanLnAlpha->SetMarkerStyle(22);
hist_meanLnAlpha->SetMarkerSize(4);
hist_meanLnAlpha->SetTitle("mean LnAlpha");
    hist_meanLnAlpha->GetXaxis()->SetTitle("y");
    hist_meanLnAlpha->GetYaxis()->SetTitle("<ln Alpha_hom>");
    hist_meanLnAlpha->SetMaximum(2.0);
    hist_meanLnAlpha->SetMinimum(0.5);
hist_meanLnAlpha->Draw("ACP");
gPad->SetLogx();
par->cd(3);
hist_sigmaLnT->SetLineWidth(2);
hist_sigmaLnT->SetMarkerStyle(22);
hist_sigmaLnT->SetMarkerSize(4);
hist_sigmaLnT->SetTitle("sigma LnT");
    hist_sigmaLnT->GetXaxis()->SetTitle("y");
    hist_sigmaLnT->GetYaxis()->SetTitle("sigma(ln T_hom)");
    hist_sigmaLnT->SetMaximum(0.4);
    hist_sigmaLnT->SetMinimum(0.0);
hist_sigmaLnT->Draw("ACP");
gPad->SetLogx();
par->cd(4);
hist_sigmaLnAlpha->SetLineWidth(2);
hist_sigmaLnAlpha->SetMarkerStyle(22);
hist_sigmaLnAlpha->SetMarkerSize(4);
hist_sigmaLnAlpha->SetTitle("sigma LnAlpha");
    hist_sigmaLnAlpha->GetXaxis()->SetTitle("y");
    hist_sigmaLnAlpha->GetYaxis()->SetTitle("sigma(ln Alpha_hom)");
    hist_sigmaLnAlpha->SetMaximum(0.5);
    hist_sigmaLnAlpha->SetMinimum(0.0);
hist_sigmaLnAlpha->Draw("ACP");
gPad->SetLogx();
par->cd(5);
hist_correlationAlphaT->SetLineWidth(2);
hist_correlationAlphaT->SetMarkerStyle(22);
hist_correlationAlphaT->SetMarkerSize(4);
hist_correlationAlphaT->SetTitle("correlation Alpha-T");
    hist_correlationAlphaT->GetXaxis()->SetTitle("y");
    hist_correlationAlphaT->GetYaxis()->SetTitle("rho(ln T_hom, ln Alpha_hom)");
    hist_correlationAlphaT->SetMaximum(1.0);
    hist_correlationAlphaT->SetMinimum(0.0);
hist_correlationAlphaT->Draw("ACP");
gPad->SetLogx();
 
    
    
    
par->SaveAs("/Users/eugenia/desktop/EMCal/parameters/param.png");*/  
    
    
    
   /* TGraph* RC = new TGraph(16);
    RC->SetPoint(0,0.05,myParam->rC(0.05,75));
    RC->SetPoint(1,0.25,myParam->rC(0.25,75));
    RC->SetPoint(2,0.50,myParam->rC(0.50,75));
    RC->SetPoint(3,0.75,myParam->rC(0.75,75));
    RC->SetPoint(4,1.00,myParam->rC(1.00,75));
    RC->SetPoint(5,1.25,myParam->rC(1.25,75));
    RC->SetPoint(6,1.50,myParam->rC(1.50,75));
    RC->SetPoint(7,1.75,myParam->rC(1.75,75));
    RC->SetPoint(8,2.00,myParam->rC(2.00,75));
    RC->SetPoint(9,2.25,myParam->rC(2.25,75));
    RC->SetPoint(10,2.50,myParam->rC(2.50,75));
    RC->SetPoint(11,2.75,myParam->rC(2.75,75));
    RC->SetPoint(12,3.00,myParam->rC(3.00,75));
    RC->SetPoint(13,3.25,myParam->rC(3.25,75));
    RC->SetPoint(14,3.50,myParam->rC(3.50,75));
    RC->SetPoint(15,3.75,myParam->rC(3.75,75));
    
    TGraph* RT = new TGraph(16);
    RT->SetPoint(0,0.05,myParam->rT(0.05,75));
    RT->SetPoint(1,0.25,myParam->rT(0.25,75));
    RT->SetPoint(2,0.50,myParam->rT(0.50,75));
    RT->SetPoint(3,0.75,myParam->rT(0.75,75));
    RT->SetPoint(4,1.00,myParam->rT(1.00,75));
    RT->SetPoint(5,1.25,myParam->rT(1.25,75));
    RT->SetPoint(6,1.50,myParam->rT(1.50,75));
    RT->SetPoint(7,1.75,myParam->rT(1.75,75));
    RT->SetPoint(8,2.00,myParam->rT(2.00,75));
    RT->SetPoint(9,2.25,myParam->rT(2.25,75));
    RT->SetPoint(10,2.50,myParam->rT(2.50,75));
    RT->SetPoint(11,2.75,myParam->rT(2.75,75));
    RT->SetPoint(12,3.00,myParam->rT(3.00,75));
    RT->SetPoint(13,3.25,myParam->rT(3.25,75));
    RT->SetPoint(14,3.50,myParam->rT(3.50,75));
    RT->SetPoint(15,3.75,myParam->rT(3.75,75));
    
    TGraph* P = new TGraph(16);
    P->SetPoint(0,0.05,myParam->p(0.05,75));
    P->SetPoint(1,0.25,myParam->p(0.25,75));
    P->SetPoint(2,0.50,myParam->p(0.50,75));
    P->SetPoint(3,0.75,myParam->p(0.75,75));
    P->SetPoint(4,1.00,myParam->p(1.00,75));
    P->SetPoint(5,1.25,myParam->p(1.25,75));
    P->SetPoint(6,1.50,myParam->p(1.50,75));
    P->SetPoint(7,1.75,myParam->p(1.75,75));
    P->SetPoint(8,2.00,myParam->p(2.00,75));
    P->SetPoint(9,2.25,myParam->p(2.25,75));
    P->SetPoint(10,2.50,myParam->p(2.50,75));
    P->SetPoint(11,2.75,myParam->p(2.75,75));
    P->SetPoint(12,3.00,myParam->p(3.00,75));
    P->SetPoint(13,3.25,myParam->p(3.25,75));
    P->SetPoint(14,3.50,myParam->p(3.50,75));
    P->SetPoint(15,3.75,myParam->p(3.75,75));
    
TCanvas * pp= new TCanvas("pp","pp",10000,10000,2500,2500); 
RC->SetLineWidth(2);
RC->SetLineColor(kRed);
RC->SetTitle("Rc[RM]");
RC->GetYaxis()->SetTitle("Rc[RM], Rt[RM], p");
RC->GetXaxis()->SetTitle("tau=t/T");
RC->SetMaximum(2.5);
RC->SetMinimum(0.0);
RC->Draw();   
    
RT->SetLineWidth(2);
RT->SetTitle("Rt[RM]");
RT->SetLineColor(kBlue);
RT->Draw("same"); 

P->SetLineWidth(2);
P->SetTitle("p");
P->SetLineColor(kOrange);
P->Draw("same"); 
    
pp->BuildLegend(0.25,0.15,0.25,0.15);
pp->SaveAs("/Users/eugenia/desktop/EMCal/parameters/RcRtP.png");  */


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