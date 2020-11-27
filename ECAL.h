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
void Fill_(TH1F* &Rad1, TH1F* &Rad2, TH1F* &Rad3, TH1F* &Rad4);
void Fill_Lat(double tt, double stepEn);
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
TH1F* EnLong;

};
#endif