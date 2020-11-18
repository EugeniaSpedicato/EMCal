#ifndef EMSHOWER_h
#define EMSHOWER_h

#include <TH2.h>
#include <cmath>
#include "EMECALShowerParametrization.h"
#include "ECAL.h"
#include "IncGamma.h"


//#include "RandomEngineAndDistribution.h"
#include "GammaFunctionGenerator.h" 
#include "RadialInterval.h" 
#include <iostream>
#include "Math/Vector3D.h"
#include <TRandom3.h>
#include <TMath.h>
#include <TGraph.h>


using namespace std;



class EMShower {
  typedef ROOT::Math::XYZVector XYZPoint;
  typedef pair<XYZPoint, double> Spot;
  typedef pair<unsigned int, double> Step;
  typedef vector<Step> Steps;
  typedef Steps::const_iterator step_iterator;

public:
  EMShower(//TRandom3 const* engine,
           GammaFunctionGenerator* gamma,
           EMECALShowerParametrization* const myParam,
           //std::vector<const RawParticle*>* const myPart,
           ECAL* const myGrid,
           bool bFixedLength,
           int nPart,
           vector<double> energy_in);

  virtual ~EMShower() { ; }

  /// Computes the steps before the real compute
  void prepareSteps();

  /// Compute the shower longitudinal and lateral development
  void compute();

  /// get the depth of the centre of gravity of the shower(s)
  //  inline double getMeanDepth() const {return globalMeanDepth;};

  /// get the depth of the maximum of the shower
  inline double getMaximumOfShower() const { return globalMaximum; }

  /// set the grid address
  void setGrid(ECAL* const myGrid) { theGrid = myGrid; }


private:
  // The longitudinal development ersatzt.
  double gam(double x, double a) const;

  // Energy deposited in the layer t-dt-> t, in units of E0 (initial energy)
  double deposit(double t, double a, double b, double dt);

  // Energy deposited between 0 and t, in units of E0 (initial energy)
  double deposit(double a, double b, double t);

  // Set the intervals for the radial development
  void setIntervals(unsigned icomp, RadialInterval& rad);

  // The parametrization
  EMECALShowerParametrization* const theParam;

  // The Calorimeter properties
  const ECALProperties* theECAL;

  // The incident particle(s)
  //std::vector<const RawParticle*>* const thePart;
  int nPart;
  std::vector<double> energy_in;
    

  // The basic quantities for the shower development.
  std::vector<double> theNumberOfSpots;
  std::vector<double> Etot;
  //std::vector<std::vector<double> > Etot_step;
  std::vector<double> Etot_step;
    
  std::vector<double> E;
  std::vector<double> photos;
  std::vector<double> T;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> Ti;
  std::vector<double> TSpot;
  std::vector<double> aSpot;
  std::vector<double> bSpot;

  // F.B : Use the maximum of the shower rather the center of gravity
  //  std::vector<double> meanDepth;
  //  double globalMeanDepth;
  std::vector<double> maximumOfShower;
  std::vector<std::vector<double> > depositedEnergy;
  std::vector<double> meanDepth;
  double innerDepth, outerDepth;

  double globalMaximum;

  double totalEnergy;

/*std::vector<double> spotE1;
std::vector<double> spotE2;
std::vector<double> spotE3;
std::vector<double> spotE4;*/
std::vector<TH1F*> EnRad_;
TH1F* EnRad_3;
TH1F* EnRad_6;
TH1F* EnRad_13;
TH1F* EnRad_20;
    
    

  // The steps for the longitudinal development
  Steps steps;
  unsigned nSteps;
  bool stepsCalculated;

  // The crystal grid
  ECAL* theGrid;
  TH2F* EcalGrid;

  // Histos

  //  Histos* myHistos;
  IncGamma myIncompleteGamma;

  // Random engine
  //const RandomEngineAndDistribution* random;

  // integer gamma function generator
  GammaFunctionGenerator* myGammaGenerator;

  bool bFixedLength_;
};

#endif