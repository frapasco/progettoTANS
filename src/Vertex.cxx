#include <iostream>

#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"

#include "Vertex.h"

using namespace std;
using namespace TMath ;

ClassImp(Vertex)

//Nota: il file contenente le distribuzioni per eta e molteplicità (e i relativi istogrammi) viene manipolato dalle funzioni interne alla classe, 
//cercando così di ottenere la maggiore protezione possibile del suo contenuto dall'utente (minore rischio che una modifica/un errore da parte di un utente rovini 
//il file e porti a impossibilità (o assurdità) di estrazioni random dagli istogrammi e conseguenti possibili problemi con l'esecuzione del programma

//_____________________________________________________________________________________________
//default constructor (vertex in origin, 0 products)
Vertex :: Vertex() : TObject(){  
  fX=0;     
  fY=0;   
  fZ=0;
  fMulti=0;   
}

//_____________________________________________________________________________________________
//standard constructor 
Vertex :: Vertex(TString multV, TString fileName, TString multDist, double sigmaX, double sigmaY, double sigmaZ, int u1, int u2) : TObject(){
  fMulti=0;
  fX=Var(sigmaX);    //x,y,z extracted with a Gauss centered in the origin
  fY=Var(sigmaY);   
  fZ=Var(sigmaZ);
  if (multV=="No"||multV=="NO"||multV=="no") this->MultUnif(u1, u2); //choice for multiplicity extraction
  else this->MultFunc(fileName, multDist); 
}

//_____________________________________________________________________________________________
//copy constructor 
Vertex :: Vertex(const Vertex &v) : 
  TObject(v),
  fX(v.fX),
  fY(v.fY),
  fZ(v.fZ),
  fMulti(v.fMulti)
  {
  }
  
//_____________________________________________________________________________________________
//Box-Muller for coordinates
double Vertex::Var(double s){ 
  if(s>0){
  	double u1=gRandom->Rndm(); 
  	double u2=gRandom->Rndm();   
  	return Sqrt(-2*Log(u1))*Cos(2*Pi()*u2)*s; //Gauss with average=0, sigma=s
  }
  else{  
    cout<<"Error: sigma <=0"<<endl;
    return 0;
  }
}

//_____________________________________________________________________________________________
double Vertex :: GetX() const{
	return fX;
}

double Vertex :: GetY() const{
	return fY;
}

double Vertex :: GetZ() const{
	return fZ;
}

double Vertex :: GetM() const{
	return fMulti;
}

double Vertex :: GetPhi() const{
	return fPhi;
}

double Vertex :: GetTheta() const{
	return fTheta;
}

//_______________________________________________________________
//uniform multiplicity
void Vertex :: MultUnif(int u1, int u2){
 if (u2>=u1)  fMulti=u1+(u2-u1)*gRandom->Rndm();
 else {
   cout<<"Error in parameters choice"<<endl;
   return;
 }
}

//_______________________________________________________________
//given multiplicity distribution function
void Vertex :: MultFunc(TString fileName, TString multDist){
  TFile *f = new TFile(fileName);
  if(f==NULL){
    cout<<"Error: the file doesn't exists"<<endl; 
    return;
  }	
  TH1D *myh = (TH1D*)f->Get(multDist);
  if(myh==NULL){
    cout<<"Error: the histogram doesn't exists"<<endl; 
    return;
  }	
  fMulti=(int) myh->GetRandom();
  f->Close();
}      
  
//_______________________________________________________________
//initializes the initial direction, min and max are for reducing the range in pseudorapidity and getting more tracks inside the detector
//Remark: reducing too much the range eliminates a relevant fraction of the original distribution; (-2, 2) is a good compromise for keeping the peaks of the original distribution
void Vertex :: InitialDir(double min, double max, TString fileName, TString etaDist){ 
  if(max<min){ //exchanging the limits if they are not in the correct order
    double o=max; 
    max=min; 
    min=o;
  }else if(max==min){
    cout<<"Error: etaMax=etaMin"<<endl; 
    return;
  }
  double eta;
  fPhi = gRandom->Rndm()*2*Pi(); //phi exctraction
  TFile *f = new TFile(fileName);
  if(f==NULL){
    cout<<"Error: the file doesn't exists"<<endl; 
    return;
  }
  TH1D *myh = (TH1D*)f->Get(etaDist); // eta distribution from histogram
  if(myh==NULL){
    cout<<"Error: the histogram doesn't exists"<<endl; 
    return;
  }	
  do{
    eta =(double) myh->GetRandom();
     }
  while (eta>max||eta<min); //it continues until it is in the wanted interval
  fTheta=2*ATan(Exp(-eta));
  f->Close();
} 

//_______________________________________________________
//destructor
Vertex :: ~Vertex(){  
}
