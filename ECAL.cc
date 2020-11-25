#include <TH2F.h>
#include <map>
#include <cmath>

#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>


#include "ECAL.h"

using namespace std;

ECAL::ECAL(double nbinsx, 
    double xlow, 
    double xup, 
    double nbinsy, 
    double ylow, 
    double yup)
    :
    nbinX(nbinsx),nbinY(nbinsy),Xlow(xlow),Xup(xup),Ylow(ylow),Yup(yup) 
    {
        //Questa mappa serve a mappare il numero di bin nel numero vero della cella
        //perchè i numeri dei bin sono sballati a causa degli overflow e underflow bins
        number[36]=1; number[37]=2; number[38]=3; number[39]=4; number[40]=5;
        number[29]=6; number[30]=7; number[31]=8; number[32]=9; number[33]=10;
        number[22]=11; number[23]=12; number[24]=13; number[25]=14; number[26]=15;
        number[15]=16; number[16]=17; number[17]=18; number[18]=19; number[19]=20;
        number[8]=21; number[9]=22; number[10]=23; number[11]=24; number[12]=25;
        
    EnRad_3 = new TH1F("Step3", "Radial Profile Step 3", 20, 0, 4);
    EnRad_6 = new TH1F("Step6", "Radial Profile Step 6", 20, 0, 4);
    EnRad_13 = new TH1F("Step13", "Radial Profile Step 13", 20, 0, 4);
    EnRad_20 = new TH1F("Step20", "Radial Profile Step 20", 20, 0, 4);
    EnLong = new TH1F("Long", "Longitudinal Profile", 25, 0, 25);
    
    
    }



// metodo che crea l'istogramma rappresentante il calorimetro
TH2F* ECAL::CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup)
{
    TH2F* EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter with E in GeV",nbinsx,xlow,xup,nbinsy,ylow,yup);
    return EcalGrid;
};


// metodo che aggiunge il punto di coo(x,y) all'istogramma, quindi al calorimetro
void ECAL::AddHitCoo(double r, double phi,double xi, double yi, double w, TH2F* a)
{   r *= 2.19;
    double x=r*cos(phi)+xi; // coo x
    double y=r*sin(phi)+yi; // coo y
    a->Fill(x,y,w);    
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

};

void ECAL::Fill_(vector<pair<double,double>> &Spot1,vector<pair<double,double>> &Spot2,vector<pair<double,double>> &Spot3,vector<pair<double,double>> &Spot4,double realTotalEnergy12,double realTotalEnergy56,double realTotalEnergy1314,double realTotalEnergy2223)
{
    for (int i=0; i< Spot1.size(); ++i)
    {EnRad_3->Fill(Spot1[i].first,Spot1[i].second/realTotalEnergy12);}     
    for (int i=0; i< Spot2.size(); ++i)
    {EnRad_6->Fill(Spot2[i].first,Spot2[i].second/realTotalEnergy56);}
    for (int i=0; i< Spot3.size(); ++i)
    {EnRad_13->Fill(Spot3[i].first,Spot3[i].second/realTotalEnergy1314);}
    for (int i=0; i< Spot4.size(); ++i)
    {EnRad_20->Fill(Spot4[i].first,Spot4[i].second/realTotalEnergy2223); }         
}
    

/*void ECAL::Fill_(double ri, double spote, double iStep,double realTotalEnergy12,double realTotalEnergy56,double realTotalEnergy1314,double realTotalEnergy2223)
{
             if (iStep==0 || iStep==1 )
              {
                EnRad_3->Fill(ri,spote/realTotalEnergy12);
              }
              
             if (iStep==4 || iStep==5)
              {
                EnRad_6->Fill(ri,spote/realTotalEnergy56);
              }
              
              if (iStep==12 || iStep==13)
              {
                EnRad_13->Fill(ri,spote/realTotalEnergy1314);
              }
              
              if (iStep==21 || iStep==22)
              {
                EnRad_20->Fill(ri,spote/realTotalEnergy2223);
              }                
}*/

void ECAL::Fill_Lat(double tt, double stepEn)
{
    EnLong->Fill(tt,stepEn);
}

void ECAL::Print_()
{
/*double c3=EnRad_3->Integral();
EnRad_3->Scale(1/c3);
double c6=EnRad_6->Integral();
EnRad_6->Scale(1/c6);
double c13=EnRad_13->Integral();
EnRad_13->Scale(1/c13);
double c20=EnRad_20->Integral();
EnRad_20->Scale(1/c20);*/

// because of 1/dr in the histo article
EnRad_3->Scale(1/0.2);
EnRad_6->Scale(1/0.2);
EnRad_13->Scale(1/0.2);
EnRad_20->Scale(1/0.2);
EnRad_3->Scale(1/100);
EnRad_6->Scale(1/100);
EnRad_13->Scale(1/100);
EnRad_20->Scale(1/100);
    
double c4=EnLong->Integral();
EnLong->Scale(1/c4);
// because of 1/dt in the histo article
EnLong->Scale(1/1);



    
TCanvas * en_lat= new TCanvas("en_lat","en_lat",1500,1000,3500,2000); 
en_lat->Divide(2,2);
en_lat->cd(1);
EnRad_3->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_3->SetXTitle("r (RM)");
EnRad_3->SetMarkerColor(kBlue);
EnRad_3->SetMarkerStyle(20);
EnRad_3->SetMarkerSize(2);
EnRad_3->Draw();
gPad->SetLogy();
en_lat->cd(2);
EnRad_6->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_6->SetXTitle("r (RM)");
EnRad_6->SetMarkerColor(kRed);
EnRad_6->SetMarkerStyle(20);
EnRad_6->SetMarkerSize(2);
EnRad_6->Draw();
gPad->SetLogy();
en_lat->cd(3);
EnRad_13->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_13->SetXTitle("r (RM)");
EnRad_13->SetMarkerColor(kBlack);
EnRad_13->SetMarkerStyle(20);
EnRad_13->SetMarkerSize(2);
EnRad_13->Draw();
gPad->SetLogy();
en_lat->cd(4);
EnRad_20->SetYTitle("dE(t)^-1 dE(t,r)/dr (RM^-1)");
EnRad_20->SetXTitle("r (RM)");
EnRad_20->SetMarkerColor(kOrange);
EnRad_20->SetMarkerStyle(20);
EnRad_20->SetMarkerSize(2);
EnRad_20->Draw();
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
en_tot2->SaveAs("/Users/eugenia/desktop/EMCal/profiloLong.png");
}




// metodo che assegna il numero della cella che viene colpita dalla particella 
void ECAL::GiveNcell(double coox,double cooy,TH2F* a)
{   
    int binx = a->GetXaxis()->FindBin(coox);
    int biny = a->GetYaxis()->FindBin(cooy);
    int nbin = a->GetBin(binx,biny);

    //return number[nbin];
    cout <<"Number of the bin:" << number[nbin] << endl;
    //return number;
};




