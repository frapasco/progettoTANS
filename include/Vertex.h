#ifndef VERTEX_H
#define VERTEX_H

#include "TObject.h"
#include "TString.h"

//////////////////////////////////////////////////////////////////////////////////
//Classe per generare vertice primario e direzione prodotti
//////////////////////////////////////////////////////////////////////////////////

class Vertex : public TObject {  

public:  
  //costruttori e distruttore
  Vertex();                
  Vertex(TString multv, TString filename, TString muldis, double sigmax, double sigmay, double sigmaz, int u1=1, int u2=90); //multv=scelta distribuzione molteplicità
  Vertex(const Vertex &v); //copy constructor
  virtual ~Vertex();
            
  void InitialDir(double min, double max, TString filename, TString etadis);//estrazione direzione prodotti di collisione
  
  //get
  double GetX() const;
  double GetY() const;
  double GetZ() const;
  double GetM() const;
  double GetPhi() const; 
  double GetTheta() const; 
  
  
private:
  double fX, fY, fZ;//coordinate vertice
  int fMulti;//numero di particelle cariche da collisione (molteplicità)
  double fPhi,fTheta;//angoli polare e azimutale di prodotti
  
  double Var(double s); //estrazione gaussiana con Box-Muller
  void Multunif(int u1, int u2); //estrazione molteplicità uniforme fra u1 e u2
  void Multfunc(TString filename, TString muldis); //estrazione molteplicità da distribuzione data (istogramma in file)
 
  ClassDef(Vertex,1);
};  
        
#endif
