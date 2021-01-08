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
double GiveCentralCell(double coox,double cooy,TH2F* a);
void SetEnergy(double energy);
int* GiveArray3x3(int n);
double AddHitCoo(double r,double phi,double xi,double yi,double w,TH2F* a);
void AddHitCooDepth(double r, double phi,double xi, double yi, double w, double depth, double deX0depthoffset_pth, TH2F* a);
void Draw_ECAL(TH2F* a);
void Fill_(TH1F* &Rad1, TH1F* &Rad2, TH1F* &Rad3, TH1F* &Rad4,TH1F* &RadTot, TH1F* &en_1cell,TH1F* &en_3x3cell);
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
TH1F* Energy_dist;
TH1F* Energy_dist1;
TH1F* Energy_dist3x3;

int *Array9;
};
#endif