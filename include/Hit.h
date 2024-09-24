#ifndef HIT_H
#define HIT_H

#include "TObject.h"

////////////////////////////////////////////////////////////////////
//class for the hit point
////////////////////////////////////////////////////////////////////

class Hit: public TObject {

public:
  //constructors and destructor
  Hit();
  Hit(double X, double Y, double Z, int label=-1);
  Hit(double R, double L, int label=-1);
  virtual ~Hit();
   
  //getters
  double GetX() const;
  double GetY() const;
  double GetZ() const;
  double GetPhi() const;
  int GetLabel() const;
  
  void Smearing(double R, double sigz=0.012,double sigt=0.003); //function for the smearing of hit points

private:
  double fX, fY, fZ;//coordinates
  int fLabelTrack;//label for track numbering
  ClassDef(Hit,1)
};

#endif