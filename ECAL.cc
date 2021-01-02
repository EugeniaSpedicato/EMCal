#include <TH2F.h>
#include <map>
#include <cmath>

#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"


#include "ECAL.h"


using namespace std;
using namespace RooFit;


ECAL::ECAL(double nbinsx, 
    double xlow, 
    double xup, 
    double nbinsy, 
    double ylow, 
    double yup)
    :
    nbinX(nbinsx),nbinY(nbinsy),Xlow(xlow),Xup(xup),Ylow(ylow),Yup(yup) 
    {
        //Queste mappe servono a mappare il numero di bin nel numero vero della cella e viceversa
        //perchè i numeri dei bin sono sballati a causa degli overflow e underflow bins
        number[36]=1; number[37]=2; number[38]=3; number[39]=4; number[40]=5;
        number[29]=6; number[30]=7; number[31]=8; number[32]=9; number[33]=10;
        number[22]=11; number[23]=12; number[24]=13; number[25]=14; number[26]=15;
        number[15]=16; number[16]=17; number[17]=18; number[18]=19; number[19]=20;
        number[8]=21; number[9]=22; number[10]=23; number[11]=24; number[12]=25;
        
        Rev_number[1]=36; Rev_number[2]=37; Rev_number[3]=38; Rev_number[4]=39; Rev_number[5]=40;
        Rev_number[6]=29; Rev_number[7]=30; Rev_number[8]=31; Rev_number[9]=32; Rev_number[10]=33;
        Rev_number[11]=22; Rev_number[12]=23; Rev_number[13]=24; Rev_number[14]=25; Rev_number[15]=26;
        Rev_number[16]=15; Rev_number[17]=16; Rev_number[18]=17; Rev_number[19]=18; Rev_number[20]=19;
        Rev_number[21]=8; Rev_number[22]=9; Rev_number[23]=10; Rev_number[24]=11; Rev_number[25]=12;
      

        
    EnRad_3 = new TProfile("Step3", "Radial Profile 2-3 X0", 20, 0, 4, 0, 5);
    EnRad_3->SetErrorOption("S");
    EnRad_6 = new TProfile("Step6", "Radial Profile 6-7 X0", 20, 0, 4, 0, 5);
    EnRad_6->SetErrorOption("S");
    EnRad_13 = new TProfile("Step13", "Radial Profile 19-20 X0", 20, 0, 4, 0, 5);
    EnRad_13->SetErrorOption("S");
    EnRad_20 = new TProfile("Step20", "Radial Profile 22-23 X0", 20, 0, 4, 0, 5);
    EnRad_20->SetErrorOption("S");
    EnRad_tot = new TProfile("Step20", "Radial Profile Total", 20, 0, 4, 0, 2);
    EnRad_tot->SetErrorOption("S");
        
    EnRad_3ERR = new TProfile("Step3", "Radial Profile 1-2 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_3ERR->SetErrorOption("S");
    EnRad_6ERR = new TProfile("Step6", "Radial Profile 5-6 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_6ERR->SetErrorOption("S");
    EnRad_13ERR = new TProfile("Step13", "Radial Profile 13-14 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_13ERR->SetErrorOption("S");
    EnRad_20ERR = new TProfile("Step20", "Radial Profile 22-23 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_20ERR->SetErrorOption("S");
        
    EnLong = new TProfile("Long", "Longitudinal Profile", 25, 0, 25);
    EnLong->SetErrorOption("S");
    EnLongERR = new TProfile("Long", "Longitudinal Profile with RMS", 25, 0, 25);
    EnLongERR->SetErrorOption("S");
    
    Er= new TProfile("<E(r)>/E", "Mean Energy Fraction in r <E(r)>/E ", 20, 0, 4, 0, 5);
    Er->SetErrorOption("S");
        
    Er2= new TProfile("<E(r)>", "Mean Energy Fraction in r <E(r)> ", 20, 0, 4);
    Er2->SetErrorOption("S");
    
    Energy_dist =new TH1F("Energy", "Energy",100,90,100);
    Energy_dist1 =new TH1F("Energy", "Energy 1 cell",150,90,105);
    Energy_dist3x3 =new TH1F("Energy", "Energy 3x3 cells",150,80,100);
    

    sigma =  new TProfile("Res", "Stochastic term",20, 0, 4, 0, 5);
    sigma->SetErrorOption("S");
    
    Array9=0;
        
    }



// metodo che crea l'istogramma rappresentante il calorimetro
TH2F* ECAL::CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup)
{
    TH2F* EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter with E in GeV",nbinsx,xlow,xup,nbinsy,ylow,yup);
    return EcalGrid;
};


// metodo che assegna il numero della cella che viene colpita dalla particella 
double ECAL::GiveCentralCell(double coox,double cooy,TH2F* a)
{   
    int binx = a->GetXaxis()->FindBin(coox);
    int biny = a->GetYaxis()->FindBin(cooy);
    int nbin = a->GetBin(binx,biny);

    cout <<"Number of the cell:" << number[nbin] << endl;

    return number[nbin];
};

int* ECAL::GiveArray3x3(int n)
{
    if (n==1) {Array9= new int[9]{1,6,7,2,0,0,0,0,0};}
    if (n==2) {Array9= new int[9]{1,2,6,7,8,3,0,0,0};}
    if (n==3) {Array9= new int[9]{2,3,7,8,9,4,0,0,0};}
    if (n==4) {Array9= new int[9]{3,4,8,9,10,0,0,0};}
    if (n==5) {Array9= new int[9]{4,5,9,10,0,0,0};}
    if (n==6) {Array9= new int[9]{1,2,6,7,12,11,0,0,0};}
    if (n==7) {Array9= new int[9]{1,2,3,6,6,8,11,12,13};}
    if (n==8) {Array9= new int[9]{2,3,4,7,8,9,12,13,14};}
    if (n==9) {Array9= new int[9]{3,4,5,8,9,10,13,14,15};}
    if (n==10) {Array9= new int[9]{4,5,9,10,14,15,0,0,0};}
    if (n==11) {Array9= new int[9]{6,7,11,12,16,17,0,0,0};}
    if (n==12) {Array9= new int[9]{6,7,8,11,12,13,16,17,18};}
    if (n==13) {Array9= new int[9]{7,8,9,12,13,14,17,18,19};}
    if (n==14) {Array9= new int[9]{8,9,10,13,14,15,18,19,20};}
    if (n==15) {Array9= new int[9]{9,10,14,15,19,20,0,0,0};}
    if (n==16) {Array9= new int[9]{11,12,16,17,21,22,0,0,0};}
   if (n==17) {Array9= new int[9]{11,12,13,16,17,18,21,22,23};}
   if (n==18) {Array9= new int[9]{12,13,14,17,18,19,22,23,24};}
   if (n==19) {Array9= new int[9]{13,14,15,18,19,20,23,24,25};}
    if (n==20) {Array9= new int[9]{14,15,19,20,24,25,0,0,0};}
    if (n==21) {Array9= new int[9]{16,17,21,22,0,0,0,0,0};}
    if (n==22) {Array9= new int[9]{16,17,18,21,22,23,0,0,0};}
    if (n==23) {Array9= new int[9]{17,18,19,22,23,24,0,0,0};}
    if (n==24) {Array9= new int[9]{18,19,20,23,24,25,0,0,0};}
    if (n==25) {Array9= new int[9]{19,20,24,25,0,0,0,0,0};}
    
    return 0;
}

// metodo che aggiunge il punto di coo(x,y) all'istogramma, quindi al calorimetro e dà numero cella
double ECAL::AddHitCoo(double r, double phi,double xi, double yi, double w, TH2F* a)
{   r *= 2.19;
    double x=r*cos(phi)+xi; // coo x in cm
    double y=r*sin(phi)+yi; // coo y in cm
    a->Fill(x,y,w);   
 
double number=ECAL::GiveCentralCell(x,y,a);
return number;
};

void ECAL::AddHitCooDepth(double r, double phi,double xi, double yi, double w, double depth, double X0depth, TH2F* a)
{   depth += X0depth;
    r *= 2.19;
    double x=r*cos(phi)+xi; // coo x in cm
    double y=r*sin(phi)+yi; // coo y in cm
 if (24.7-X0depth>depth) 
 {a->Fill(x,y,w);   
double number=ECAL::GiveCentralCell(x,y,a); cout <<"è giusto"<< endl;}

};




// metodo che disegna l'evento nel calorimetro e le celle che vengono colpite
void ECAL::Draw_ECAL(TH2F* a){

TCanvas * Ecal_= new TCanvas("Ecal_","Ecal_",1500,100,3500,2000);
Ecal_->Divide(2,1);
Ecal_->cd(1);
gStyle->SetPalette(kAquamarine);
//TColor::InvertPalette();
a->SetXTitle("x (cm)");
a->SetYTitle("y (cm)");
a->Draw("COL");
a->Draw("TEXT SAME");
Ecal_->cd(2);
a->Draw("LEGO");
Ecal_->SaveAs("/Users/eugenia/desktop/EMCal/Ecal.png");


// riempi celle    
int binMax=a->GetMaximumBin();  
int CentralCell=number[binMax];
cout << "cella centrale rev " << Rev_number[CentralCell] <<" and vera " << CentralCell << endl;
Energy_dist1->Fill(a->GetBinContent(binMax));

double energy3x3=0.;    
ECAL::GiveArray3x3(CentralCell);
for (int i=0; i<9; ++i)
{
    if (Array9[i]>0 & Array9[i]<25 & Array9[i]!=0) energy3x3+=a->GetBinContent(Rev_number[Array9[i]]);
    cout << Rev_number[Array9[i]] << " and vera " << Array9[i]<< " c'è energia " << energy3x3 << endl;
}
Energy_dist3x3->Fill(energy3x3);  
};

//inline double getX0back() const { return maxX0_; }

/*void ECAL::Fill_(vector<pair<double,double>> &Spot1,vector<pair<double,double>> &Spot2,vector<pair<double,double>> &Spot3,vector<pair<double,double>> &Spot4,double realTotalEnergy12,double realTotalEnergy56,double realTotalEnergy1314,double realTotalEnergy2223)
{
    for (int i=0; i< Spot1.size(); ++i)
    {EnRad_3->Fill(Spot1[i].first,Spot1[i].second/(1000*realTotalEnergy12));}     
    for (int i=0; i< Spot2.size(); ++i)
    {EnRad_6->Fill(Spot2[i].first,Spot2[i].second/(1000*realTotalEnergy56));}
    for (int i=0; i< Spot3.size(); ++i)
    {EnRad_13->Fill(Spot3[i].first,Spot3[i].second/(1000*realTotalEnergy1314));}
    for (int i=0; i< Spot4.size(); ++i)
    {EnRad_20->Fill(Spot4[i].first,Spot4[i].second/(1000*realTotalEnergy2223)); }         
}*/
    


void ECAL::Fill_(TH1F* &Rad1, TH1F* &Rad2, TH1F* &Rad3, TH1F* &Rad4, TH1F* &RadTot, TH1F* &en_1cell, TH1F* &en_3x3cell)
{


    
double c1=Rad1->Integral();
    
double c2=Rad2->Integral();

double c3=Rad3->Integral();

double c4=Rad4->Integral();
    
double c5=RadTot->Integral();


    // 0.2/2 a metà del bin
    int nx1=Rad1->GetNbinsX();
    for (int i=0; i<nx1+1; ++i) 
    {
       EnRad_3->Fill(Rad1->GetXaxis()->GetBinCenter(i),Rad1->GetBinContent(i)/(0.2*c1));
    }
    int nx2=Rad2->GetNbinsX();
    for (int i=0; i<nx2+1; ++i) 
    {
       EnRad_6->Fill(Rad2->GetXaxis()->GetBinCenter(i),Rad2->GetBinContent(i)/(0.2*c2));
    }
    int nx3=Rad1->GetNbinsX();
    for (int i=0; i<nx3+1; ++i) 
    {
       EnRad_13->Fill(Rad3->GetXaxis()->GetBinCenter(i),Rad3->GetBinContent(i)/(0.2*c3));
    }
    int nx4=Rad1->GetNbinsX();
    for (int i=0; i<nx4+1; ++i) 
    {
       EnRad_20->Fill(Rad4->GetXaxis()->GetBinCenter(i),Rad4->GetBinContent(i)/(0.2*c4));
    }
    int nx5=RadTot->GetNbinsX();
    for (int i=0; i<nx5+1; ++i) 
    {
       EnRad_tot->Fill(RadTot->GetXaxis()->GetBinCenter(i),RadTot->GetBinContent(i)/(0.2*c5));
        //Er->Fill(RadTot->GetXaxis()->GetBinCenter(i),RadTot->Integral(0,i+1)/100);
        //Er2->Fill(RadTot->GetXaxis()->GetBinCenter(i),RadTot->Integral(0,i+1));
        
        //sigma->Fill(RadTot->GetXaxis()->GetBinCenter(i),(1/RadTot->Integral(0,i+1))*(sqrt(10)));
    }
    
    Energy_dist->Fill(RadTot->Integral());  
    //Energy_dist1->Fill(en_1cell->Integral());  
    //Energy_dist3x3->Fill(en_3x3cell->Integral());  
    
    
}



void ECAL::Fill_Lat(TH1F* &Longit)
{
   double c=Longit->Integral();

    // 0.2/2 a metà del bin
    int nx=Longit->GetNbinsX();
    for (int i=0; i<nx+1; ++i) 
    {
       EnLong->Fill(Longit->GetXaxis()->GetBinCenter(i),Longit->GetBinContent(i)/(1*c));
        
    }

}

void ECAL::Print_()
{
    
/*int n=Er->GetNbinsX();
for (int i=0; i<n+1; ++i) 
{   
sigma->Fill(Er2->GetXaxis()->GetBinCenter(i),(1/Er2->GetBinContent(i))*(sqrt(100)));}*/

    int nx1=EnRad_3->GetNbinsX();
    for (int i=0; i<nx1+1; ++i) 
    {
       EnRad_3ERR->Fill(EnRad_3->GetXaxis()->GetBinCenter(i),EnRad_3->GetBinContent(i)+EnRad_3->GetBinError(i));
    } 
    int nx2=EnRad_6->GetNbinsX();
    for (int i=0; i<nx2+1; ++i) 
    {
       EnRad_6ERR->Fill(EnRad_6->GetXaxis()->GetBinCenter(i),EnRad_6->GetBinContent(i)+EnRad_6->GetBinError(i));
    } 

    int nx3=EnRad_13->GetNbinsX();
    for (int i=0; i<nx3+1; ++i) 
    {
       EnRad_13ERR->Fill(EnRad_13->GetXaxis()->GetBinCenter(i),EnRad_13->GetBinContent(i)+EnRad_13->GetBinError(i));
    } 
    
    int nx4=EnRad_20->GetNbinsX();
    for (int i=0; i<nx4+1; ++i) 
    {
       EnRad_20ERR->Fill(EnRad_20->GetXaxis()->GetBinCenter(i),EnRad_20->GetBinContent(i)+EnRad_20->GetBinError(i));
    } 
    
   
TCanvas * en_lat= new TCanvas("en_lat","en_lat",1500,1000,3500,2000); 
en_lat->Divide(2,2);
en_lat->cd(1);
EnRad_3->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_3->SetXTitle("r (RM)");
EnRad_3->SetMarkerColor(kBlue);
EnRad_3->SetMarkerStyle(20);
EnRad_3->SetMarkerSize(2);
EnRad_3->Draw("HIST SAME");
EnRad_3->Draw("HIST SAME P");
EnRad_3ERR->SetMarkerColor(38);
EnRad_3ERR->SetMarkerStyle(20);
EnRad_3ERR->SetMarkerSize(2);
EnRad_3ERR->Draw("P same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);
gPad->SetLogy();
en_lat->cd(2);
EnRad_6->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_6->SetXTitle("r (RM)");
EnRad_6->SetMarkerColor(kRed);
EnRad_6->SetMarkerStyle(20);
EnRad_6->SetMarkerSize(2);
EnRad_6->Draw("HIST SAME");
EnRad_6->Draw("HIST SAME P");
EnRad_6ERR->SetMarkerColor(46);
EnRad_6ERR->SetMarkerStyle(20);
EnRad_6ERR->SetMarkerSize(2);
EnRad_6ERR->Draw("P same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);
gPad->SetLogy();
en_lat->cd(3);
EnRad_13->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_13->SetXTitle("r (RM)");
EnRad_13->SetMarkerColor(kBlack);
EnRad_13->SetMarkerStyle(20);
EnRad_13->SetMarkerSize(2);
EnRad_13->Draw("HIST SAME");
EnRad_13->Draw("HIST SAME P");
EnRad_13ERR->SetMarkerColor(39);
EnRad_13ERR->SetMarkerStyle(20);
EnRad_13ERR->SetMarkerSize(2);
EnRad_13ERR->Draw("P same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);
gPad->SetLogy();
en_lat->cd(4);
EnRad_20->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_20->SetXTitle("r (RM)");
EnRad_20->SetMarkerColor(kOrange);
EnRad_20->SetMarkerStyle(20);
EnRad_20->SetMarkerSize(2);
EnRad_20->Draw("HIST SAME");
EnRad_20->Draw("HIST SAME P");
EnRad_20ERR->SetMarkerColor(41);
EnRad_20ERR->SetMarkerStyle(20);
EnRad_20ERR->SetMarkerSize(2);
EnRad_20ERR->Draw("P same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);
gPad->SetLogy();
en_lat->SaveAs("/Users/eugenia/desktop/EMCal/profRadPerStep.png");
    
    
TCanvas * en_tot= new TCanvas("en_to","Profiles Rad",1000,100,2500,2000); 
    
EnRad_3->SetTitle("Average Radial Profile");
EnRad_3->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_3->SetXTitle("r (RM)");
EnRad_3->SetLineColor(kBlue);
EnRad_3->SetMarkerColor(kBlue);
EnRad_3->SetLineWidth(2);
EnRad_3->SetMarkerStyle(20);
EnRad_3->SetMarkerSize(2);
EnRad_3->Draw("HIST");
EnRad_3->Draw("HIST SAME P");
    
EnRad_20->SetLineColor(kOrange);
EnRad_20->SetMarkerColor(kOrange);
EnRad_20->SetLineWidth(2);
EnRad_20->SetMarkerStyle(20);
EnRad_20->SetMarkerSize(2);
EnRad_20->Draw("HIST SAME");
EnRad_20->Draw("HIST SAME P");
    
EnRad_6->SetLineColor(kRed);
EnRad_6->SetMarkerColor(kRed); 
EnRad_6->SetLineWidth(2);
EnRad_6->SetMarkerStyle(20);
EnRad_6->SetMarkerSize(2);
EnRad_6->Draw("HIST SAME");
EnRad_6->Draw("HIST SAME P");
    
EnRad_13->SetLineColor(kBlack);
EnRad_13->SetMarkerColor(kBlack);
EnRad_13->SetLineWidth(2);
EnRad_13->SetMarkerStyle(20);
EnRad_13->SetMarkerSize(2);
EnRad_13->Draw("HIST SAME");
EnRad_13->Draw("HIST SAME P");
gPad->SetLogy();
en_tot->SaveAs("/Users/eugenia/desktop/EMCal/profRad.png");

    int nx=EnLong->GetNbinsX();
    for (int i=0; i<nx+1; ++i) 
    {
       EnLongERR->Fill(EnLong->GetXaxis()->GetBinCenter(i),EnLong->GetBinContent(i)+EnLong->GetBinError(i));
    } 

TCanvas * en_tot2= new TCanvas("en_to2","Profile Long",1000,100,2500,2000); 
EnLong->SetMaximum(0.15);
EnLong->SetMinimum(0);
EnLong->SetYTitle("E^-1 dE(t)/dt (X0^-1)");
EnLong->SetXTitle("t (X0)");
EnLong->SetTitle("Average Longitudinal Profile");
EnLong->SetMarkerColor(8);
EnLong->SetMarkerStyle(20);
EnLong->SetMarkerSize(2);
EnLong->Draw();
EnLongERR->SetMarkerColor(46);
EnLongERR->SetMarkerStyle(20);
EnLongERR->SetMarkerSize(2);
EnLongERR->Draw("P same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);
    
en_tot2->SaveAs("/Users/eugenia/desktop/EMCal/profiloLong.png");
   
    
TCanvas * en_tot3= new TCanvas("en_to2","Profile Long",1000,100,2500,2000); 
EnRad_tot->SetTitle("Average Radial Profile All X0");
EnRad_tot->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_tot->SetXTitle("r (RM)");
EnRad_tot->SetLineColor(30);
EnRad_tot->SetMarkerColor(30);
EnRad_tot->SetLineWidth(2);
EnRad_tot->SetMarkerStyle(20);
EnRad_tot->SetMarkerSize(2);
EnRad_tot->Draw("HIST");
EnRad_tot->Draw("HIST SAME P");
gPad->SetLogy();
en_tot3->SaveAs("/Users/eugenia/desktop/EMCal/TotRad.png");
    
/*TCanvas * E_r= new TCanvas("E_r","<Er>/E",1000,100,2500,2000);   
Er->SetMaximum(1.1);
Er->SetMinimum(0.5);
Er->SetYTitle("<E(r)>/E");
Er->SetXTitle("r (RM)");
Er->SetLineColor(kPink);
Er->SetMarkerColor(kBlack);
Er->SetLineWidth(2);
Er->SetMarkerStyle(20);
Er->SetMarkerSize(2);
Er->Draw("HIST");
Er->Draw("HIST SAME P");
E_r->SaveAs("/Users/eugenia/desktop/EMCal/Er.png");
    
TCanvas * sig= new TCanvas("Stoc","Stochastic Term",1000,100,2500,2000);   
sigma->SetMaximum(0.5);
sigma->SetMinimum(0.);
sigma->SetYTitle("sigma/<E(r)> sqrt(E)");
sigma->SetXTitle("r (RM)");
sigma->SetLineColor(kOrange);
sigma->SetMarkerColor(kBlue);
sigma->SetLineWidth(2);
sigma->SetMarkerStyle(20);
sigma->SetMarkerSize(2);
sigma->Draw("HIST");
sigma->Draw("HIST SAME P");
sig->SaveAs("/Users/eugenia/desktop/EMCal/StocTerm.png");*/

    
TCanvas * encell= new TCanvas("Energy cells","Energy cells",1000,100,2500,2000);     
encell->Divide(1,3);
encell->cd(1);
/*TF1* f1 = new TF1("f1", "gaus",141.5, 144.5);
TF1* f1pol = new TF1("f1pol", "pol3",140, 141.5);
    
TF1* f2 = new TF1("f2", "gaus",146.5, 148);
TF1* f2pol = new TF1("f2pol", "pol3",144, 146.5);*/
    

    
Energy_dist1->Fit("gaus");
Energy_dist1->SetLineWidth(2);
Energy_dist1->Draw("same");
encell->cd(2);
Energy_dist3x3->Fit("gaus");
/*Energy_dist3x3->Fit("f1","R");
Energy_dist3x3->Fit("f1pol","R+");*/
    
Energy_dist3x3->SetLineWidth(2);
Energy_dist3x3->Draw("same");
encell->cd(3);
Energy_dist->Fit("gaus");
/*Energy_dist->Fit("f2","R");
Energy_dist->Fit("f2pol","R+");*/
    
Energy_dist->SetLineWidth(2);
Energy_dist->Draw("same");    

encell->SaveAs("/Users/eugenia/desktop/EMCal/EnCell.png");

    
// Observable
RooRealVar energy3("energy3","energy3",80,100) ;
RooRealVar mean("mean","mean",94,96.5) ;
RooRealVar sigma("sigma","sigma",0.2,1.6) ;
RooRealVar alpha("alpha","alpha",1,0,20) ;
RooRealVar n("n","n",4,1,8) ;
RooCBShape CrystallBall("CrystallBall", "CrystallBall", energy3, mean, sigma, alpha, n);


RooDataHist en3("en3x3","en3x3",energy3,Import(*Energy_dist3x3));

RooPlot *frame = energy3.frame(Title("energy 3x3 cells"));
en3.plotOn(frame,MarkerStyle(kFullDotMedium));
CrystallBall.fitTo(en3);
CrystallBall.plotOn(frame);
CrystallBall.paramOn(frame,Layout(0.12,0.50));
frame->Draw();
    
TCanvas* cROO= new TCanvas("cROO","cROO",400,10,1100,800);
frame->GetXaxis()->SetTitle("Energy [GeV]");
frame->Draw();
cROO->SaveAs("/Users/eugenia/desktop/EMCal/frame.png");

    
}