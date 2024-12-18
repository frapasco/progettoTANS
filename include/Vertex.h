#ifndef VERTEX_H
#define VERTEX_H

#include "TObject.h"
#include "TH1D.h"
#include "TString.h"

//////////////////////////////////////////////////////////////////////////////////
// Class to generate primary vertex and product direction
//////////////////////////////////////////////////////////////////////////////////

class Vertex : public TObject {  

public:  
  // Constructors and destructor
  Vertex();                
  Vertex(TString multV, TH1D *multHist, double sigmaX, double sigmaY, double sigmaZ, int u1=1, int u2=90); // multV=distribution choice for multiplicity
  Vertex(const Vertex &v); // Copy constructor
  virtual ~Vertex();
            
  void InitialDir(double min, double max, TH1D *etaHist); // Direction setter for collision products
  
  // Getters
  double GetX() const;
  double GetY() const;
  double GetZ() const;
  int GetM() const;
  double GetPhi() const; 
  double GetTheta() const; 
  
private:
  double fX, fY, fZ; // Vertex coordinates
  int fMulti; // Multiplicity of particles per vertex
 
  double fPhi, fTheta; // Polar and azimuthal angles
  
  double Var(double s); // Box-Muller extraction for Gaussian
  void MultUnif(int u1, int u2); // Uniform extraction for multiplicity in the interval [u1, u2]
  void MultFunc(TH1D *multHist); // Multiplicity from a given histogram
  
  ClassDef(Vertex, 1);
};  
        
#endif
