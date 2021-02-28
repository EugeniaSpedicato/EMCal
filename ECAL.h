#ifndef ECAL_h
#define ECAL_h

#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>


#include <cmath>
#include <map>
#include <TMatrixD.h>
#include <TGraph.h>

#include <iostream>

using namespace std;

class ECAL 
//: public TH2
{
public:

//costruttore
ECAL(double nbinsx, 
    double xlow, 
    double xup, 
    double nbinsy, 
    double ylow, 
    double yup);
//distruttore
~ ECAL(){}


//double radlen;

TH2F* CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup);
TH2F* GiveEcalGrid();

double GiveCentralCell(double coox,double cooy );
void SetEnergy(double energy);
int* GiveArray3x3(int n);
double AddHitCoo(double r,double phi,double xi,double yi,double w);
void AddHitCooDepth(double r, double phi,double xi, double yi, double w, double depth, double deX0depthoffset_pth);
void Draw_ECAL(double Xi, double Yi,int i);
void Fill_(TH1F* &Rad1, TH1F* &Rad2, TH1F* &Rad3, TH1F* &Rad4,TH1F* &RadTot, TH1F* &en_1cell,TH1F* &en_3x3cell);
vector<double> EnergyContent();
void Fill_Lat(TH1F* &Longit);
void Print_();
inline void setSpotEnergy(double e) { spotEnergy = e; }

private:
const double nbinX;
const double nbinY;
const double Xlow;
const double Xup;
const double Ylow;
const double Yup;
double spotEnergy;
double energy_IN;
typedef map<int, int>  n_cell;
n_cell number;
n_cell Rev_number;
n_cell Rev_numberX;
n_cell Rev_numberY;
    
vector<double> E_cell;


TProfile* EnRad_3;
TProfile* EnRad_6;
TProfile* EnRad_13;
TProfile* EnRad_20;
TProfile* EnRad_3ERR;
TProfile* EnRad_6ERR;
TProfile* EnRad_13ERR;
TProfile* EnRad_20ERR;
TProfile* EnRad_tot;
TProfile* EnLong;
TProfile* EnLongERR;
TProfile* Er;
TProfile* Er2;

TProfile* sigma;
//TH1F* Energy_dist;
TH1F* Energy_dist1;
TH1F* Energy_dist3x3;
TH1F* E_tot;
TH2F* EcalGrid;
TH1F* E_sigma_tot;
TH1F* E_sigma_2;
TH1F* E_sigma_4;
TH2F* Emean_out;
TProfile* En_r1x1;
TProfile* En_r3x3;
TProfile* En_N3x3;

TH1F* residual;
TH1F* residuoX;
TH1F* residuoY;

 TH2F* X_Y_e;

int *Array9;
};
#endif