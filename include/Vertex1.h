#ifndef VERTEX_H
#define VERTEX_H

#include "TObject.h"
#include "TString.h"

//////////////////////////////////////////////////////////////////////////////////
//Classe per generare vertice primario e direzione prodotti
//////////////////////////////////////////////////////////////////////////////////

class Vertex : public TObject {  

public:  
  //constructors and destructor
  Vertex();                
  Vertex(TString multV, TString fileName, TString multDist, double sigmaX, double sigmaY, double sigmaZ, int u1=1, int u2=90); //multV=distribution choice for multiplicity
  Vertex(const Vertex &v); //copy constructor
  virtual ~Vertex();
            
  void InitialDir(double min, double max, TString fileName, TString etaDist); //direction setter for collision products

  //getters
  double GetX() const;
  double GetY() const;
  double GetZ() const;
  double GetM() const;
  double GetPhi() const; 
  double GetTheta() const; 
  
  
private:
  double fX, fY, fZ; //vertex coordinates
  int fMulti; //multiplicity of particles per vertex
 
  double fPhi,fTheta;//polar and azimuthal angles
  
  double Var(double s); //Box-Muller extraction for Gaussian
  void MultUnif(int u1, int u2); //uniform extraction for multiplicity in the interval [u1,u2]
  void MultFunc(TString fileName, TString multDist); //multiplicity from a given distribution
 
  ClassDef(Vertex,1);
};  
        
#endif
