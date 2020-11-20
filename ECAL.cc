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
    }



// metodo che crea l'istogramma rappresentante il calorimetro
TH2F* ECAL::CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup)
{
    TH2F* EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter with E in GeV",nbinsx,xlow,xup,nbinsy,ylow,yup);
    return EcalGrid;
};

// metodo che registra il punto di coo(r,phi)
/*TMatrixD ECAL::addHit(double r, double phi)
{
    TMatrixD a(1,2);
    a[0][0]=r;
    a[0][1]=phi;
        
        return a;};


TMatrixD ECAL::Conversion(TMatrixD c){
    TMatrixD ccart(1,1);
    ccart[0][0]=c[0][0]*cos(c[0][1]); // coo x
    ccart[0][1]=c[0][0]*sin(c[0][1]); // coo y
    
    return ccart;
};
*/

// metodo che aggiunge il punto di coo(x,y) all'istogramma, quindi al calorimetro
void ECAL::AddHitCoo(double r, double phi,double xi, double yi, double w, TH2F* a)
{   r *= 2.19;
    double x=r*cos(phi)+xi; // coo x
    double y=r*sin(phi)+yi; // coo y
    a->Fill(x,y,w);    
};

//void ECAL::GiveCellEn(espote,TH2F* a){}

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




