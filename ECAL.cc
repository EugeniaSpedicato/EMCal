#include <TH2F.h>
#include <map>
#include <cmath>

#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TColor.h>
#include <TFile.h>

#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
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
      

Rev_numberX[1]=1; Rev_numberX[2]=2; Rev_numberX[3]=3; Rev_numberX[4]=4; Rev_numberX[5]=5;
Rev_numberX[6]=1; Rev_numberX[7]=2; Rev_numberX[8]=3; Rev_numberX[9]=4; Rev_numberX[10]=5;
Rev_numberX[11]=1; Rev_numberX[12]=2; Rev_numberX[13]=3; Rev_numberX[14]=4; Rev_numberX[15]=5;
Rev_numberX[16]=1; Rev_numberX[17]=2; Rev_numberX[18]=3; Rev_numberX[19]=4; Rev_numberX[20]=5;
Rev_numberX[21]=1; Rev_numberX[22]=2; Rev_numberX[23]=3; Rev_numberX[24]=4; Rev_numberX[25]=5;
    
    
Rev_numberY[1]=5; Rev_numberY[2]=5; Rev_numberY[3]=5; Rev_numberY[4]=5; Rev_numberY[5]=5;
Rev_numberY[6]=4; Rev_numberY[7]=4; Rev_numberY[8]=4; Rev_numberY[9]=4; Rev_numberY[10]=4;
Rev_numberY[11]=3; Rev_numberY[12]=3; Rev_numberY[13]=3; Rev_numberY[14]=3; Rev_numberY[15]=3;
Rev_numberY[16]=2; Rev_numberY[17]=2; Rev_numberY[18]=2; Rev_numberY[19]=2; Rev_numberY[20]=2;
Rev_numberY[21]=1; Rev_numberY[22]=1; Rev_numberY[23]=1; Rev_numberY[24]=1; Rev_numberY[25]=1;
        
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
    
    //Energy_dist =new TH1F("Energy", "Energy",100,90,100);
    Energy_dist1 =new TH1F("Energy", "Energy 1 cell",150,0.67,0.87);//0.77,0.87
    Energy_dist3x3 =new TH1F("Energy", "Energy 3x3 cells",150,0.9,0.98);//0.9,0.98
    E_tot =new TH1F("Energy", "Energy 5x5",150,0.3,1);  
    E_sigma_tot=new TH1F("Energyt", "Energy 5x5",150,0,75);  
    E_sigma_2=new TH1F("Energy2", "Energy 5x5",150,20,40); 
    E_sigma_4=new TH1F("Energy4", "Energy 5x5",150,30,40);  
    

    sigma =  new TProfile("Res", "Stochastic term",20, 0, 4, 0, 5);
    sigma->SetErrorOption("S");
    
    Array9=0;
    En_r1x1 = new TProfile("Enr1x1", "Position impact point VS energy1x1", 30, 0, 1.425);
    En_r1x1->SetErrorOption("S");
    En_r3x3 = new TProfile("En_r3x3", "Position impact point VS energy3x3", 30, 0, 1.425 );
    En_r3x3->SetErrorOption("S");

    
    En_N3x3 = new TProfile("En_r3x3", "Number impact cell VS energy3x3", 25, 0.5, 25.5 );
    En_N3x3->SetErrorOption("S");
     
        X_Y_e  = new TH2F("CooEL" , " X  Vs. Y of the electron",500,-7,7,500,-7,7);
        
        residual = new TH1F("residual", "Residual", 100, -1, 1);

    }



// metodo che crea l'istogramma rappresentante il calorimetro
TH2F* ECAL::CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup)
{
    EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter with E in GeV",nbinsx,xlow,xup,nbinsy,ylow,yup);
    return EcalGrid;
};


TH2F* ECAL::GiveEcalGrid()
{return EcalGrid;};


void ECAL::SetEnergy(double energy)
{
    energy_IN=energy;
}
// metodo che assegna il numero della cella che viene colpita dalla particella 
double ECAL::GiveCentralCell(double coox,double cooy )
{   
    int binx = EcalGrid->GetXaxis()->FindBin(coox);
    int biny = EcalGrid->GetYaxis()->FindBin(cooy);
    int nbin = EcalGrid->GetBin(binx,biny);

    //cout <<"Number of the cell:" << number[nbin] << endl;

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
    if (n==7) {Array9= new int[9]{1,2,3,6,7,8,11,12,13};}
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
double ECAL::AddHitCoo(double r, double phi,double xi, double yi, double w)
{   r *= 2.19; //multiplication moliere radius  * radiuscorrectionfactor  preso da CMS barrel
    double x=r*cos(phi)+xi; // coo x in cm
    double y=r*sin(phi)+yi; // coo y in cm
    EcalGrid->Fill(x,y,w); 
    X_Y_e->Fill(x,y,w);
 
double number=ECAL::GiveCentralCell(x,y);
return number;
};

void ECAL::AddHitCooDepth(double r, double phi,double xi, double yi, double w, double depth, double X0depth)
{   depth += X0depth;
    r *= 2.19;
    double x=r*cos(phi)+xi; // coo x in cm
    double y=r*sin(phi)+yi; // coo y in cm
 cout << "depth: " << depth << ", 24.7-X0depth" << 24.7-X0depth << endl;
 if (24.7-X0depth>depth) 
 {EcalGrid->Fill(x,y,w);   cout << "è aggiunto "<<endl; 
 X_Y_e->Fill(x,y);}

};


// metodo che disegna l'evento nel calorimetro e le celle che vengono colpite
void ECAL::Draw_ECAL(double Xi, double Yi, int i){

TCanvas * Ecal_= new TCanvas("Ecal_","Ecal_",1500,100,3500,2000);
Ecal_->Divide(2,1);
Ecal_->cd(1);
gStyle->SetPalette(kAquamarine);
//TColor::InvertPalette();
EcalGrid->SetXTitle("x (cm)");
EcalGrid->SetYTitle("y (cm)");
EcalGrid->Draw("COL");
EcalGrid->Draw("TEXT SAME");
Ecal_->cd(2);
EcalGrid->Draw("LEGO");
Ecal_->SaveAs("/Users/eugenia/desktop/EMCal/Ecal.png");


// riempi celle    
int binMax=EcalGrid->GetMaximumBin();  
int CentralCell=number[binMax];
//cout << "cella centrale rev " << Rev_number[CentralCell] <<" and vera " << CentralCell << endl;
//Energy_dist1->Fill(EcalGrid->GetBinContent(binMax));
Energy_dist1->Fill(EcalGrid->GetBinContent(binMax)/energy_IN);
double energy3x3=0.;    
ECAL::GiveArray3x3(CentralCell);
for (int i=0; i<9; ++i)
{
    if (Array9[i]>0 && Array9[i]<26) energy3x3+=EcalGrid->GetBinContent(Rev_number[Array9[i]]);
    //cout << Rev_number[Array9[i]] << " and vera " << Array9[i]<< " c'è energia " << energy3x3 << endl;
}
Energy_dist3x3->Fill(energy3x3/energy_IN);
//E_tot->Fill(EcalGrid->GetBinContent(binMax)/energy3x3);
double energy_tot=0.;
for (int i=1;i<26;++i)
{ energy_tot+=EcalGrid->GetBinContent(Rev_number[i]);}
E_tot->Fill(energy_tot/energy_IN);
//E_tot->Fill(EcalGrid->GetBinContent(binMax)/energy3x3);
    
double r=sqrt(Xi*Xi+Yi*Yi);    

 En_r1x1->Fill(r,EcalGrid->GetBinContent(binMax)/energy_IN);
 En_r3x3->Fill(r,energy3x3/energy_IN);

 En_N3x3->Fill(CentralCell,energy3x3/energy_IN);
   

double Ex=0.;
double Ey=0.;
double wtot=0.;

for(int i=0; i<9; ++i)
{
double x = EcalGrid->GetXaxis()->GetBinCenter(Rev_numberX[Array9[i]]);
double y = EcalGrid->GetYaxis()->GetBinCenter(Rev_numberY[Array9[i]]);
    
double w=4.0+log(EcalGrid->GetBinContent(Rev_number[Array9[i]])/energy3x3);
double wi=(w>0)?w:0;
//double wi=EcalGrid->GetBinContent(Rev_number[Array9[i]]);
    
Ex+=wi*x;
Ey+=wi*y;
wtot+=wi;
}
double centroidX=(Ex)/wtot;
double centroidY=(Ey)/wtot;  
    
double ddd=sqrt((centroidX-Xi)*(centroidX-Xi)+(centroidY-Yi)*(centroidY-Yi));    
double r_cal=sqrt((centroidX)*(centroidX)+(centroidY)*(centroidY));  
double r_trak=sqrt((Xi)*(Xi)+(Yi)*(Yi));  
double a=(r_cal>r_trak)?(+1):(-1);
double res=ddd*a;

    cout << " r_cal " << r_cal << ", r_trak " << r_trak << endl;
    cout << " centroidX " << centroidX << ", centroidY " << centroidY << endl;
    
residual->Fill(res);
    
};


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
        Er->Fill(RadTot->GetXaxis()->GetBinCenter(i),RadTot->Integral(0,i+1)/energy_IN);
        Er2->Fill(RadTot->GetXaxis()->GetBinCenter(i),RadTot->Integral(0,i+1));
        //sigma->Fill(RadTot->GetXaxis()->GetBinCenter(i),(1/RadTot->Integral(0,i+1))*(sqrt(10)));
    }    
    
E_sigma_tot->Fill(RadTot->Integral());
E_sigma_2->Fill(RadTot->Integral(0,2));
E_sigma_4->Fill(RadTot->Integral(0,4));

    
    //Energy_dist->Fill(RadTot->Integral());  
    //Energy_dist1->Fill(en_1cell->Integral());  
    //Energy_dist3x3->Fill(en_3x3cell->Integral());  
    
    
}

vector<double> ECAL::EnergyContent()
{
    for (int i=1; i<26 ; ++i)
    {E_cell.push_back(EcalGrid->GetBinContent(Rev_number[i]));}
    return E_cell;
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
E_sigma_tot->Fit("gaus");
TF1 *p1 = E_sigma_tot->GetFunction("gaus");
double si = p1->GetParameter(2);    

/*E_sigma_2->Fit("gaus");
TF1 *p2 = E_sigma_2->GetFunction("gaus");
si[1] = p2->GetParameter(2); 

E_sigma_4->Fit("gaus");
TF1 *p3 = E_sigma_4->GetFunction("gaus");
si[2] = p3->GetParameter(2);    */

    

int ner=Er->GetNbinsX();
for (int i=0; i<ner+1; ++i) 
{sigma->Fill(Er2->GetXaxis()->GetBinCenter(i),(si/Er2->GetBinContent(i))*(sqrt(energy_IN)));} 

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
    
    
TCanvas * sig= new TCanvas("Stoc","Stochastic Term",1000,100,2500,2000);   
sigma->SetMaximum(0.5);
sigma->SetMinimum(0.);
sigma->SetYTitle("sigma/<E(r)> sqrt(E)");
sigma->SetXTitle("r (RM)");
sigma->SetLineColor(46);
sigma->SetMarkerColor(kBlack);
sigma->SetLineWidth(2);
sigma->SetMarkerStyle(20);
sigma->SetMarkerSize(2);
sigma->Draw("HIST L");
sigma->Draw("HIST SAME P");
sig->SaveAs("/Users/eugenia/desktop/EMCal/StocTerm.png");

    
TCanvas * encell= new TCanvas("Energy cells","Energy cells",1000,100,2500,2000);     
encell->Divide(1,3);
encell->cd(1);
/*TF1* f1 = new TF1("f1", "gaus",141.5, 144.5);
TF1* f1pol = new TF1("f1pol", "pol3",140, 141.5);
    
TF1* f2 = new TF1("f2", "gaus",146.5, 148);
TF1* f2pol = new TF1("f2pol", "pol3",144, 146.5);*/
    

Energy_dist1->SetLineWidth(3);
Energy_dist1->Fit("gaus");
Energy_dist1->Draw("same");
Energy_dist1->GetXaxis()->SetTitle("E_rec/E_true");
encell->cd(2);
Energy_dist3x3->Fit("gaus");
Energy_dist3x3->SetLineWidth(3);
Energy_dist3x3->Draw("same");
Energy_dist3x3->GetXaxis()->SetTitle("E_rec/E_true");
encell->cd(3);
E_tot->SetLineWidth(3);
E_tot->GetXaxis()->SetTitle("E_rec/E_true");
E_tot->Draw("same");

encell->SaveAs("/Users/eugenia/desktop/EMCal/EnCell.png");

    
// Observable
RooRealVar energy3("energy3","energy3",0.9,0.98) ;//0.9,0.98
RooRealVar mean("mean","mean",0.94,0.96) ;
RooRealVar sigma0("sigma0","sigma0",0.002,0.015) ;
RooRealVar alpha("alpha","alpha",1,0,20) ;
RooRealVar n("n","n",4,1,8) ;
    
RooCBShape CrystallBall("CrystallBall", "CrystallBall", energy3, mean, sigma0, alpha, n);
    

RooDataHist en3("en3x3","en3x3",energy3,Import(*Energy_dist3x3));

RooPlot *frame = energy3.frame(Title("energy 3x3 cells"));
en3.plotOn(frame,MarkerStyle(kFullDotMedium));
CrystallBall.fitTo(en3);
CrystallBall.plotOn(frame);
cout << "chiiii fram3x3 " << frame->chiSquare() << endl;
CrystallBall.paramOn(frame,Layout(0.12,0.33));
en3.statOn(frame,Layout(0.12,0.33,0.65)) ;
frame->Draw();
    
RooRealVar energy5("energy5","energy5",0.9,1) ;//0.9,0.98
RooRealVar mean5("mean5","mean5",0.95,0.98) ;
RooRealVar sigma5("sigma5","sigma5",0.002,0.015) ;
RooRealVar alpha5("alpha5","alpha5",1,0,20) ;
RooRealVar n5("n5","n5",4,1,8) ;
    
RooCBShape CrystallBall5("CrystallBall5", "CrystallBall5", energy5, mean5, sigma5, alpha5, n5);
    

RooDataHist en5("en5x5","en5x5",energy5,Import(*E_tot));

RooPlot *frame5 = energy5.frame(Title("energy 5x5 cells"));
en5.plotOn(frame5,MarkerStyle(kFullDotMedium));
CrystallBall5.fitTo(en5);
CrystallBall5.plotOn(frame5);
cout << "chiiii fram5x5 " << frame5->chiSquare() << endl;
CrystallBall5.paramOn(frame5,Layout(0.12,0.33));
en5.statOn(frame5,Layout(0.12,0.33,0.65)) ;
frame5->Draw();
        
/*RooRealVar energy1("energy1","energy1",0.30,1) ;
RooRealVar mean1("mean1","mean1",0.75,0.85) ;
RooRealVar sigma1("sigma1","sigma1",0,0.05) ;
RooRealVar alpha1("alpha1","alpha1",0,20) ;
RooRealVar n1("n1","n1",1,10) ;
RooCBShape CrystallBall1("CrystallBall1", "CrystallBall1", energy1, mean1, sigma1, alpha1, n1);
    
RooDataHist en1("en1cell","en1cell",energy1,Import(*Energy_dist1));

RooPlot *frame1 = energy1.frame(Title("energy 1 cell"));
en1.plotOn(frame1,MarkerStyle(kFullDotMedium));
CrystallBall1.fitTo(en1);
CrystallBall1.plotOn(frame1);
cout << "chiiii frame1 " << frame1->chiSquare() << endl;
CrystallBall1.paramOn(frame1,Layout(0.12,0.50));
en1.statOn(frame1,Layout(0.12,0.33,0.65)) ;

frame1->Draw();*/
    
RooRealVar energy1("energy1","energy1",0.77,0.87) ;//0.77,0.87
RooRealVar mean1("mean1","mean1",0.81,0.84);
RooRealVar sigma1("sigma1","sigma1",0,2);
RooGaussian Gauss("Gauss", "Gauss", energy1, mean1, sigma1);
    
RooDataHist en1("en1cell","en1cell",energy1,Import(*Energy_dist1));

RooPlot *frame1 = energy1.frame(Title("energy 1 cell"));
en1.plotOn(frame1,MarkerStyle(kFullDotMedium));
Gauss.fitTo(en1);
Gauss.plotOn(frame1);
cout << "chiiii frame1 " << frame1->chiSquare() << endl;
Gauss.paramOn(frame1,Layout(0.12,0.33));
en1.statOn(frame1,Layout(0.12,0.33,0.65)) ;
frame1->Draw();
    
TCanvas * E_r= new TCanvas("E_r","<Er>/E",1000,100,2500,2000);   
Er->SetMaximum(1.1);
Er->SetMinimum(0.5);
Er->SetYTitle("<E(r)>/E");
Er->SetXTitle("r (RM)");
Er->SetLineColor(9);
Er->SetMarkerColor(kBlack);
Er->SetLineWidth(2);
Er->SetMarkerStyle(20);
Er->SetMarkerSize(2);
gStyle->SetOptStat(0);
TLine *line = new TLine(1.95,0.5,1.95,1.5);
TLine *line2 = new TLine(0.65,0.5,0.65,1.5);

line->SetLineColor(kRed);
line->SetLineStyle(2);
line2->SetLineColor(kRed);
line2->SetLineStyle(2);
Er->Draw("HIST L");
Er->Draw("HIST SAME P");
line->Draw("same");
line2->Draw("same");
E_r->SaveAs("/Users/eugenia/desktop/EMCal/Er.png");
    
TCanvas* cROO= new TCanvas("cROO","cROO",400,10,1100,800);
cROO->Divide(1,3);
cROO->cd(1);
frame1->GetXaxis()->SetTitle("E_rec/E_true");
frame1->Draw();
cROO->cd(2);
frame->GetXaxis()->SetTitle("E_rec/E_true");
frame->Draw();
cROO->cd(3);
frame5->GetXaxis()->SetTitle("E_rec/E_true");
frame5->Draw();
cROO->SaveAs("/Users/eugenia/desktop/EMCal/frame.png");

TCanvas* graph= new TCanvas("gr","gr",400,10,1100,800);
graph->Divide(1,3);
    graph->cd(1);
    En_r1x1->SetTitle("Coordinate impact point-energy 1 cell");
    En_r1x1->GetYaxis()->SetTitle("Erec/Etrue");
    En_r1x1->GetXaxis()->SetTitle("r[cm]");
    En_r1x1->SetLineWidth(2);
    En_r1x1->SetLineColor(kRed);
    En_r1x1->SetMaximum(0.9);
    En_r1x1->SetMinimum(0.45);
    En_r1x1->Draw();
    graph->cd(2);
    En_r3x3->SetMaximum(1);
    En_r3x3->SetMinimum(0.78);
    En_r3x3->SetTitle("Coordinate impact point-energy 3x3 cells");
    En_r3x3->GetYaxis()->SetTitle("Erec/Etrue");
    En_r3x3->GetXaxis()->SetTitle("r[cm]");
    En_r3x3->SetLineWidth(2);
    En_r1x1->SetLineColor(kRed);
    En_r3x3->Draw();
    graph->cd(3);
    residual->SetTitle("Number impact cell-energy 1cell/3x3 cells");
    residual->GetYaxis()->SetTitle("R9");
    residual->GetXaxis()->SetTitle("r[cm]");
    residual->SetLineWidth(2);
    residual->Draw();
    /*En_N3x3->SetMaximum(0.98);
    En_N3x3->SetMinimum(0.);
    En_N3x3->SetTitle("Number impact cell-energy 3x3 cells");
    En_N3x3->GetYaxis()->SetTitle("Erec/Etrue");
    En_N3x3->GetXaxis()->SetTitle("N_cell");
    En_N3x3->SetLineWidth(2);
    En_N3x3->Draw();*/

    graph->SaveAs("/Users/eugenia/desktop/EMCal/En-impact.png");
    
    TCanvas * duedmu= new TCanvas("duedmu","duedmu",1000,100,2500,2000);
    gStyle->SetPalette(kLake);
    TColor::InvertPalette();
    gStyle->SetOptStat();
    X_Y_e->Draw("COLZ");
    X_Y_e->GetXaxis()->SetTitle("x [m]");
    X_Y_e->GetYaxis()->SetTitle("y [m]");
  duedmu->SaveAs("/Users/eugenia/desktop/EMCal/dued.png");  
}