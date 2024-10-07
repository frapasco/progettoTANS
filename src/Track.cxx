#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "Track.h"

using namespace TMath; 
using namespace std;

ClassImp(Track)

//_____________________________________________________________________________________________
//default constructor
Track :: Track() : TObject(){
  fC1=0;
  fC2=0;
  fC3=0;
}

//_____________________________________________________________________________________________
//standard constructor 
Track :: Track (double theta1, double phi1) : TObject(){
  fTheta=theta1;
  fPhi=phi1;
  fC1=Sin(fTheta)*Cos(fPhi);
  fC2=Sin(fTheta)*Sin(fPhi);
  fC3=Cos(fTheta);	
}

//_____________________________________________________________________________________________
//constructor used to create a line for reconstruction of tracklets
Track :: Track (double x, double y, double z) : TObject(){
  double norm =1/Sqrt(x*x+y*y+z*z); 
  fC1= x*norm;
  fC2= y*norm;
  fC3= z*norm;	
}  

//_____________________________________________________________________________________
//Getters
double Track :: GetTheta () const{ 
  return fTheta;
}	

double Track :: GetPhi () const{ 
  return fPhi;
}  

//_____________________________________________________________________________________________
//implementation of the t parameter finding method
double Track :: Parameter (double x0, double y0, double R){ 
  double prod = x0*fC1+y0*fC2;
  double C = fC1*fC1+fC2*fC2;
  double Delta = (prod*prod)-C*((x0*x0)+(y0*y0)-(R*R));
  
  if (Delta<0){ //delta<0 if (x0,y0,z0) is outside of the cylinder of radius R
    cout<<"Error: Negative delta."<<endl;	
    return 0;
  }
  double t;
  double t1 = (-prod+Sqrt(Delta))/C;
  double t2 = (-prod-Sqrt(Delta))/C;

  if(t1>0) t=t1; //(t>0) 
  else t=t2;

  return t;
}

//_____________________________________________________________________________________________
//
bool Track :: Intersection(double x0, double y0, double z0, double &X, double &Y, double &Z, double R, double L){  
  double tt;
  tt =  Parameter(x0, y0, R); //tt evaluation (line parameter)
  
  //intersection point
  X=x0+fC1*tt; 
  Y=y0+fC2*tt;
  Z=z0+fC3*tt;
  
  //checking if the intersection's z is in the detector
  if(Abs(Z)>L/2){ 
    return false;
  }
  else return true;
}

//___________________________________________________________________________________________
//intersezione tracklet con piano x=0
void Track::Intersection2(double x0, double y0, double z0, double &Z){
  double tt=-x0/fC1; //tt evaluation (tt=(0-x0)/fC1)	
  Z=z0+fC3*tt; //intersection point of the tracklet with x=0 plane
}

//_____________________________________________________________________________________________
//multiple scattering
void Track:: MultipleScattering(bool scattering, double sigma){
  if(scattering==true){
    if (sigma<=0) {cout<<"Error: trying to set a negative sigma in Gauss"<<endl; return ;} 
    //theta and phi angles w.r.t. the track direction before MS
    double thetaPrimo = gRandom-> Gaus(0,sigma); //Gauss
    double phiPrimo = 2*Pi()*gRandom->Rndm(); //Uniform
    
    //Transform matrix between lab Ref. Frame (beams parallel to z axis) and particle R.F. (z' along particle direction)
    double M[3][3];
    M[0][0]=-Sin(fPhi);
    M[0][1]=-Cos(fTheta)*Cos(fPhi);
    M[0][2]=Sin(fTheta)*Cos(fPhi);
    M[1][0]=Cos(fPhi);
    M[1][1]=-Cos(fTheta)*Sin(fPhi);
    M[1][2]=Sin(fTheta)*Sin(fPhi);
    M[2][0]=0;
    M[2][1]=Sin(fTheta);
    M[2][2]=Cos(fTheta);
    
    //coordinates post scattering in the prime R.F.
    double uSecondo[3];
    uSecondo[0]=Sin(thetaPrimo)*Cos(phiPrimo);
    uSecondo[1]=Sin(thetaPrimo)*Sin(phiPrimo);
    uSecondo[2]=Cos(thetaPrimo);

    //coordinates post scattering in the lab R.F. (director cosines)
    double u[3];
    for(int i=0;i<3;i++){
      u[i]=0; 
      for(int j=0;j<3;j++){
	u[i]+=M[i][j]*uSecondo[j];
      }
    }
  
    //theta and phi post scattering
    if (u[0]==0) {
      if (u[1]==0) {cout<<"Error in director cosines: ATan not defined"<<endl; return;}
      if (u[1]>0) fPhi=Pi()/2;
      if (u[1]<0) fPhi=3*Pi()/2;
    }
    if (u[0]>0){
      if (u[1]>=0)	fPhi=ATan(u[1]/u[0]);
      else fPhi=ATan(u[1]/u[0])+2*Pi();
    }
    if (u[0]<0){
      if(u[1]>0) fPhi=ATan(u[1]/u[0])+Pi(); 
      else fPhi=ATan(u[1]/u[0])+Pi(); 
	  } 
    if(u[2]>0) fTheta=ATan(Sqrt(u[0]*u[0]+u[1]*u[1])/u[2]); 
    else {
      if(u[2]<0) fTheta=ATan(Sqrt(u[0]*u[0]+u[1]*u[1])/u[2])+Pi();
      else fTheta=Pi()/2;
    }
  }
}

//_____________________________________________________________________________________________
//distruttore
Track :: ~Track(){	 
}
