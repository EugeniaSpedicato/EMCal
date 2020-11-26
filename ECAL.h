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
void Fill_(vector<pair<double,double>> &Spot1,vector<pair<double,double>> &Spot2,vector<pair<double,double>> &Spot3,vector<pair<double,double>> &Spot4,double realTotalEnergy12,double realTotalEnergy56,double realTotalEnergy1314,double realTotalEnergy2223);
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

TH1F* EnRad_3;
TH1F* EnRad_6;
TH1F* EnRad_13;
TH1F* EnRad_20;
TH1F* EnLong;

};
#endif