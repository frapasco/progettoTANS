#ifndef TRACK_H
#define TRACK_H

#include "TObject.h"
#include "TString.h"

////////////////////////////////////////////////////////////////////
//class for keeping track and finding intersetion point with matter
////////////////////////////////////////////////////////////////////

class Track: public TObject {

public:
  //constructors and destructor
  Track();
  Track(double theta1, double phi1);
  Track(double x, double y, double z);
  virtual  ~Track();
  
  //getters
  double GetTheta() const;
  double GetPhi() const;
   
  bool Intersection(double x0, double y0, double z0, double &X, double &Y, double &Z, double R, double L); //intersection with a cylinder with radius R
  void Intersection2(double x0, double y0, double z0, double &Z); //intersection of the tracklet with x=0 plane
  void MultipleScattering(bool scattering, double sigma=0.001); //multiple scatter if(bool=true)
 
private:
  double Parameter(double x0, double y0, double R); //t parameter for intersections
  double fC1, fC2, fC3;//directions cosines
  double fTheta,fPhi;//theta and phi
 
  ClassDef(Track,1)
};

#endif
