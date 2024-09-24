#include <iostream>
#include "TRandom3.h"
#include "TMath.h"

#include "Hit.h"

using namespace TMath; 

ClassImp(Hit)

//_____________________________________________________________________________________________
//default constructor
Hit :: Hit() : TObject(){
  fX=0;
  fY=0;
  fZ=0;
  fLabelTrack=0;
} 

//______________________________________________________________________________________________
//standard constructor, for storing hits on the beam pipe and detectors
Hit :: Hit (double X, double Y, double Z, int label) : TObject(){
  fX=X;
  fY=Y;
  fZ=Z;
  fLabelTrack=label;
}

//______________________________________________________________________________________________
//constructor that generates noise uniformly on the two detectors
Hit :: Hit (double R, double L, int label) : TObject(){
  fLabelTrack=label;
  double phi = gRandom->Rndm()*2*Pi();  
  fX=R*Cos(phi);
  fY=R*Sin(phi);
  fZ=-L/2+L*gRandom->Rndm();
}

//______________________________________________________________________________________________
double Hit::GetX() const{ 
  return fX;
}

double Hit::GetY() const{ 
  return fY;
}

double Hit::GetZ() const{ 
  return fZ;
}

//GetPhi: method that gives the value of phi angle of intersection in cylindrical coords
double Hit::GetPhi() const{ 
  if (fX==0&&fY==0) {std::cout<<"Errore: x=0 e y=0 in intersezione"<<std::endl; return -1000;}
  if (fX==0&&fY>0) return Pi()/2;
  if (fX==0&&fY<0) return 3*Pi()/2;
  if (fX>0&&fY>=0) return ATan(fY/fX);
  if (fX>0&&fY<0) return ATan(fY/fX)+Pi()*2;
  if (fX<0&&fY>0) return ATan(fY/fX)+Pi();
  if (fX<0&&fY<=0) return ATan(fY/fX)+Pi();
  std::cout<<"Errore in Hit::GetPhi()"<<std::endl;
  return 0;  
}
 
int Hit::GetLabel() const{ 
  return fLabelTrack;
} 

//______________________________________________________________________________________________
//smearing function
void Hit::Smearing(double R, double sigz, double sigt){
  if (sigz<=0||sigt<=0) {std::cout<<"Error: set the parameters"<<std::endl; return;}
  double sz = sigz; //cm. Std dev ricostruzione in z
  double st = sigt; //cm. Std dev ricostruzione trasversale (R*phi)
   
  double zrec=Sqrt(-2*Log(gRandom->Rndm()))*Cos(2*Pi()*gRandom->Rndm())*sz; //smearing in z
  double arec=Sqrt(-2*Log(gRandom->Rndm()))*Cos(2*Pi()*gRandom->Rndm())*st;// smearing trasversale
  
  double phirec = this->GetPhi()+arec/R;//angolo phi dopo smearing
  
  //coordinate dopo smearing
  fX=R*Cos(phirec);
  fY=R*Sin(phirec);
  fZ=fZ+zrec;
 
}

//______________________________________________________________________________________________
//destructor implementation
Hit:: ~Hit(){ 
}
