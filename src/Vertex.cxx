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
//costruttore di default (vertice=origine, 0 prodotti)
Vertex :: Vertex() : TObject(){  
  fX=0;     
  fY=0;   
  fZ=0;
  fMulti=0;   
}

//_____________________________________________________________________________________________
//costruttore standard 
Vertex :: Vertex(TString multv, TString filename, TString muldis, double sigmax, double sigmay, double sigmaz, int u1, int u2) : TObject(){
  fMulti=0;
  fX=Var(sigmax);    //x,y,z gaussiani attorno a origine
  fY=Var(sigmay);   
  fZ=Var(sigmaz);
  if (multv=="No"||multv=="NO"||multv=="no") this->Multunif(u1, u2); //con stringa si sceglie come estrarre molteplicità (uniforme o funzione)
  else this->Multfunc(filename, muldis); 
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
//Set delle coordiante (Box-Muller): variazione gaussiana di x, y e z rispetto all'origine
double Vertex::Var(double s){ 
  if(s>0){
  	double u1=gRandom->Rndm(); 
  	double u2=gRandom->Rndm();   
  	return Sqrt(-2*Log(u1))*Cos(2*Pi()*u2)*s; //estrazione random da gaussiana con media=0 (origine), sigma=s
  }
  else{  
    cout<<"Errore: sigma <=0"<<endl; //dovrebbe sempre essere maggiore di 0 per controlli in simulazione, verifica errori di memoria e di battitura e cambiamenti in macro simulazione
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
//metodo per molteplicità da distribuzione uniforme fra u1 e u2
void Vertex :: Multunif(int u1, int u2){
 if (u2>=u1)  fMulti=u1+(u2-u1)*gRandom->Rndm();
 else {
   cout<<"Errore in scelta parametri"<<endl;
   return;
 }
}

//_______________________________________________________________
//metodo per molteplicità da distribuzione assegnata
void Vertex :: Multfunc(TString filename, TString muldis){
  TFile *f = new TFile(filename); //apertura e verifica esistenza file
  if(f==NULL){
    cout<<"Errore, non esiste il file"<<endl; 
    return;
  }	
  TH1D *myh = (TH1D*)f->Get(muldis); //estrazione e verifica esistenza istogramma in file
  if(myh==NULL){
    cout<<"Errore, non esiste l'istogramma"<<endl; 
    return;
  }	
  fMulti=(int) myh->GetRandom(); //estrazione random secondo la distribuzione rappresentata dall'istogramma
  f->Close(); //chiusura file
}      
  
//_______________________________________________________________
//inizializza la direzione iniziale di un prodotto dal vertice primario, min e max per ridurre il range della distribuzione al range di proprio interesse (riducendo il numero di tracce fuori rivelatore a causa del loro theta)
//attenzione in scelta max e min per non perdere una frazione troppo grande della distribuzione (molte estrazioni "a vuoto"); con -2\2 e la distribuzione assegnata, si preservano i significativi picchi centrali e non si perdono troppe estrazioni
void Vertex :: InitialDir(double min, double max, TString filename, TString etadis){ 
  if(max<min){ //scambio estremi in modo che siano in ordine giusto
    double o=max; 
    max=min; 
    min=o;
  }
  if(max==min){
    cout<<"Errore: etamax=etamin"<<endl; 
    return;
  }
  double eta;
  fPhi = gRandom->Rndm()*2*Pi(); //phi uniforme 
  TFile *f = new TFile(filename);
  if(f==NULL){
    cout<<"Errore, non esiste il file"<<endl; 
    return;
  }
  TH1D *myh = (TH1D*)f->Get(etadis);//estrazione istogramma per distribuzione in eta
  if(myh==NULL){
    cout<<"Errore, non esiste l'istogramma"<<endl; 
    return;
  }	
  do{
    eta =(double) myh->GetRandom();
     }
  while (eta>max||eta<min); //estrae finchè non trova un valore nel range max--min voluto
  fTheta=2*ATan(Exp(-eta)); //theta come funzione di eta
  f->Close();
} 

//_______________________________________________________
//distruttore      
Vertex :: ~Vertex(){  
}
