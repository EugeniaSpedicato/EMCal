#ifndef ECAL_h
#define ECAL_h

#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>


#include <cmath>
#include <map>
#include <TMatrixD.h>
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
void AddHitCoo(double r,double phi,double xi,double yi,double w,TH2F* a);
void Draw_ECAL(TH2F* a);
void Fill_(TH1F* &Rad1, TH1F* &Rad2, TH1F* &Rad3, TH1F* &Rad4,TH1F* &RadTot);
void Fill_Lat(TH1F* &Longit);
void Print_();
void GiveNcell(double coox,double cooy,TH2F* a);
inline void setSpotEnergy(double e) { spotEnergy = e; }

private:
const double nbinX;
const double nbinY;
const double Xlow;
const double Xup;
const double Ylow;
const double Yup;
double spotEnergy;
typedef map<int, int>  n_cell;
n_cell number;

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


};
#endif