//root includes
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TClonesArray.h"
#include <vector>
#include <fstream>
#include "TDirectory.h" 

//headers include
#include "Vertex.h"
#include "Track.h"
#include "Hit.h"

using namespace std;
const int kVerbosity=1; //verbosity
const int kLimit1=17; //last bin to read for histograms for given multiplicity
const int kLimit2=21; //last bin to read for histograms for uniform multiplicity
const double kNoisefrac=1.2; //ratio for noise/multiplicity points
const double kPhiMax=0.01; //max angle between points for creating tracklets (AAAAAAAA leggitelo meglio)
//sfasamento massimo fra due punti fra cui creare tracklet (valutato per eccesso fra punti di una stessa traccia da simulazione+reconstruction indipendenti)
const double kRange=0.1; //range, attorno al picco dell'istogramma delle intersezioni delle tracklet, entro cui mediare le intersezioni per ricostruire z del vertice

const TString kSim="simulation.root"; //input file
const TString kRic="reconstruction.root"; //output file for reconstruction tree
const TString kHisto="histo.root"; //output file for histograms
const string kData="data.txt";  //config file for the run data

///////////////////////////////////////////////////////////////////
//method for plotting out the results (residues, efficiency vs multiplicity, 
//resolution vs multiplicity, efficiency vs zTrue, resolution vs zTrue)
////////////////////////////////////////////////////////////////////

void plot(const vector <double>* const vzTruep,const  vector <double>* const vzgrecp, const vector <int> * const vzmultip, const vector <double> * const vzrecrmsp, const double * const multibin, const int arraylenghtmulti, TH1D *tot1multi, TH1D *tot3multi, const int limit, const double sigmaZ, const double * const zbin, const int arraylenghtz, TH1D *totz, const int sizetrue);

////////////////////////////////////////////////////////////////////
//method for smearing points, generating noise and reconstruct the primary vertex out of tracklets
////////////////////////////////////////////////////////////////////

void reconstruction()  
  double step=0.05;//step histograms intersections with z
  double R1,R2,L,sigmaZ;
  int events;
  char a;
  //reading data.txt to get L, R1, R2, sigmaZ, events and multiplicity generation choice	
  ifstream In = ifstream(kData);
  if(!In.is_open()){
    cout<<"Error: the file doesn't exist"<<endl; 
    return;
  }
  	
  In>>R1>>R2>>L>>sigmaZ>>events>>a;
  In.close();
  TString mul;
  if (a=='N') mul="No"; else mul="sì";
  
  TFile *file = TFile::Open(kSim, "READ");
  if(file==NULL){
    cout<<"Error: the file doesn't exist"<<endl; 
    return;
  }	
  file->ls();//lettura di contenuto file
  TTree *inputTree = (TTree*)file->Get("tree");
  if(inputTree==NULL){ 
    cout<<"Error: the Tree doesn't exist"<<endl; 
    return;
  }

  //creating output file and Tree
  TFile *fileReco = TFile::Open(kRic,"RECREATE"); 
  TTree *treeReco = new TTree("treeReco", "reconstruction");
  treeReco->SetDirectory(fileReco);
   
  //histograms reconstruction primary vertex in z
  if(step<=0||kVerbosity<=0||kNoisefrac<=0||kPhiMax<=0||kRange<=0){
    cout<<"Error in parameters setting"<<endl; 
    return;
  }
  
  double zExtr=R2*(L/(R2-R1))-(L/2); //maximum z obtainable from tracklets, from the edges of the detectors
  double zMin=-zExtr-step/2;
  double zMaxBin=zExtr+step/2;	
  int nbin=(int)(zMaxBin-zMin)/step; 
  step=(zMaxBin-zMin)/nbin; //new step, in order to have a full covering of the range with nbins
                 
  //creation of histogram and TCloneArray for intersection points
  TH1D *zIntersecHisto = new TH1D("zRec", "Intersections in z; z [cm]; number of events",nbin,zMin,zMaxBin);
  TClonesArray *intersecPointDet1 = new TClonesArray("Hit", 200); //T1
  TClonesArray *intersecPointDet2 = new TClonesArray("Hit", 200); //T2
  TClonesArray &int1 = *intersecPointDet1;  
  TClonesArray &int2 = *intersecPointDet2;
 
  Vertex *vertexRec = NULL; //creation of pointer to a Vertex obj
  vector<double> Zintersection; //vector for intersections of tracklets for a single vertex

  //assignment of branches for the output tree
  treeReco->Branch("Vertexreal",&vertexRec); //primary vertex
  treeReco->Branch("zintersectionvector",&Zintersection); //vector di intersezioni Tracklet 
  treeReco->Branch("int1enoise",&intersecPointDet1);//intersection + noise T1
  treeReco->Branch("int2enoise",&intersecPointDet2);//intersection + noise T2

  //reading from input file TCloneArrays
  TClonesArray *clone1= new TClonesArray("Hit",100);
  TClonesArray *clone2= new TClonesArray("Hit",100);
  for (int i=0;i<10;i++){ 
    new ((*clone1)[i]) Hit (); 
    new ((*clone2)[i]) Hit ();
  } //per non connettere il branch a oggetto di dimensione 0
  
  Vertex *verlec = new Vertex(); //creazione oggetto Vertex 
  
  tree->GetBranch("tclone1")->SetAutoDelete(kFALSE);
  tree->SetBranchAddress("tclone1",&clone1);//indirizzamento per lettura TClone
     
  tree->GetBranch("tclone2")->SetAutoDelete(kFALSE);
  tree->SetBranchAddress("tclone2",&clone2);
 	
  TBranch *b3=tree->GetBranch("Vertext");
  b3->SetAddress(&verlec);//lettura Vertex
	
  clone1->Clear();//pulizia vettori clone1 e clone2
  clone2->Clear();

  //preparazione istogramma per efficenza(molteplicità) e array limiti bin. 
  double multibin[]={-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9.5, 11.5, 14.5, 17.5, 20.5, 25.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5,91.5};
  int arraylenghtmulti= sizeof(multibin)/sizeof(multibin[0])-1;//numero bin isto, deve essere 1 in meno di elementi di vettore per costruttore
  int limit;//dice molteplicità massima isto con multibin[limit]
  if (kLimit1>kLimit2||kLimit2>arraylenghtmulti||kLimit1<=0||kLimit2<=0){
    cout<<"Errore su kLimit1 e kLimit2"<<endl;
    return;
  }
  if(mul=="No") limit=kLimit2; else limit=kLimit1; //ci fermiamo prima in molteplicità se non molteplicità uniforme
  
  TH1D *tot1multi = new TH1D("tot1multi","tot1multi",arraylenghtmulti,multibin);//totale eventi con z simulato entro 1 sigma da 0 per efficenza in funzione di molteplicità
  TH1D *tot3multi= new TH1D("tot3multi","tot3multi",arraylenghtmulti,multibin);//entro 3 sigma
    
  //analogo per efficenza(zTrue=z del vertice primario)	
  double zbin[]={-16, -12, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 12, 16};
  int arraylenghtz= sizeof(zbin)/sizeof(zbin[0])-1;//numero bin isto, deve essere 1 in meno di elementi di vettore per costruttore	
  TH1D *totz = new TH1D("totz","totz",arraylenghtz,zbin);//totale eventi per efficenza con zTrue

  //vector per fare istogrammi finali
  vector <double> *vzTrue = new vector <double>; //z di vertice simulato
  vector <double> *vzrec = new vector <double>;	//z ricostruito
  vector <int> *vzmulti = new vector <int>;	//molteplicità (da simulazione)
  vector <double> *vzrecrms = new vector <double>; //errore (RMS) su z ricostruito	
  vzTrue->reserve(events);//allocazione di memoria adeguata a evitare frammentazioni di memoria
  vzrec->reserve(events);
  vzmulti->reserve(events);
  vzrecrms->reserve(events);
 
  //variabili utili	
  int nlines1, nlines2;
  int out1;//contatori per Hit fuori da T1 o T2  dopo smearing
  int out2; 
  int noTrack=0; //numero di tracklet per cui non si trova intersezione sensata
  int tottraccia=0;//totale tracklet generate
  double x1,y1,z1,x2,y2,z2; 
  double Z; //z di intersezione tracklet - piano x=0
  int element1, element2; //numero elementi in int1 e int 2 (senza noise)
  double zrec, zrecrms;//z ricostruita e suo errore
  double zmax;//z a cui si trova il primo massimo dell'istogramma delle intersezioni delle tracklet
  int sizevec=0;//dimensioni vector allocati su heap
	
  //loop sul numero di eventi (entrate di tree)
  for(int ev=0;ev<tree->GetEntries();ev++){
    if (ev%kVerbosity==0) cout <<"\n \n EVENTO  "<<ev+1<<endl;
    Zintersection.clear();//pulizia vettore intersezioni
    Zintersection.reserve(200);
    out1=0;
    out2=0;
    element1=0;
    element2=0;	  
    delete vertexRec;//eliminazione vertexRec	  
    zIntersecHisto->Reset("ICES");//svuotamento istogramma da evento precedente
    tree->GetEntry(ev);//lettura evento ev 
    nlines1 = clone1->GetEntriesFast();//numero eventi in clone1
    nlines2 = clone2->GetEntriesFast();//numero eventi in clone2
    vertexRec = new Vertex(*verlec);//copia di verlec in vertexRec
 	  
    //lettura punti di intersezione da tree e aggiunta dello smearing
    //loop su tracce simulate 
    for(int i=0;i<nlines1;i++){    
      int t=i-out1; 
      
      //estrazione di un' intersezione per volta da T1 
      Hit *i1 = (Hit*) clone1->At(i); 
      i1->Smearing(R1);
           		
      //riempimento treeReco int 1 con le coordinate con smearing 
      if(TMath::Abs(i1->GetZ())<L/2){ //il punto rimane dentro T1
        x1=i1->GetX();
        y1=i1->GetY() ;
        z1=i1->GetZ();	 
        new (int1[t]) Hit(x1,y1,z1,i1->GetLabel());      
	  element1++;      
      }else out1++;   
    }	
    for(int i=0;i<nlines2;i++){   
      int k=i-out2;  
     
      //estrazione di un'intersezione per volta da T2
      
      Hit *i2 = (Hit*) clone2->At(i); 
      i2->Smearing(R2);
     		
      //riempimento treeReco int 2 con le coordinate con smearing     
      if(TMath::Abs(i2->GetZ())<L/2){//il punto rimane dentro T2
     	  x2=i2->GetX();
        y2=i2->GetY() ;
        z2=i2->GetZ();	 
        new (int2[k]) Hit(x2,y2,z2,i2->GetLabel()); 
	  element2++;	
      }else out2++;      
    }	
    if(ev%kVerbosity==0)cout<<"Applicato smearing"<<endl;
    if(ev%kVerbosity==0)cout<<"Hit usciti da T1: "<<out1<<"  Hit usciti da T2: "<<out2<<endl;
        
    //aggiunta del noise
    int Nnoise = (int) (verlec->GetM())*kNoisefrac;//numero di punti di noise, assunto identico per T1 e T2
    for(int j=element1;j<element1+Nnoise;j++){
      new (int1[j]) Hit(R1,L);
    } 
    for(int k=element2;k<element2+Nnoise; k++){
      new (int2[k]) Hit(R2,L);
    }
    if (ev%kVerbosity==0)cout<<"Applicato noise"<<endl;
       
    //inzio reconstruction PV
    int Nrec1=intersecPointDet1->GetEntriesFast();//numero di elementi (hit+noise) in int1 e int2
    int Nrec2=intersecPointDet2->GetEntriesFast();
      
    int ntraccia=0; //contatore tracklet fatte per un vertice   
    for(int i=0; i<Nrec1; i++){
      for(int j=0; j<Nrec2; j++){
        Hit *i1=(Hit*) intersecPointDet1->At(i); 
        Hit *i2=(Hit*) intersecPointDet2->At(j);
        if(TMath::Abs(i2->GetPhi()-i1->GetPhi())<kPhiMax){//Creazione Tracklet Per Punti Compatibili Entro phimax
          x1=i1->GetX();
          x2=i2->GetX();
          y1=i1->GetY();
          y2=i2->GetY();
          z1=i1->GetZ();
          z2=i2->GetZ();
          Track *t = new Track(x2-x1,y2-y1,z2-z1); 
          ntraccia++; 
          tottraccia++;
	    Z=200000000;  //per controllo di funzionamento funzione intersezione
	    t->Intersection2(x2,y2,z2,Z); // aggiorna Z con intersezione tracklet - piano x=0
          if(Z>-10000000&&Z<10000000){//si riempie vettore e istogramma con intersezioni riuscite
	      zIntersecHisto->Fill(Z);
		Zintersection.push_back(Z);
	    }
	    else noTrack++;   
	    delete t;
        }//fine if su phi1-phi2<kPhiMax
      }
    }//fine for sui due TClone di intersezioni per fare tracklet
      
    if (ev%kVerbosity==0)cout<<"Create tracklet e riempito istogramma"<<endl;
    if (zIntersecHisto->GetEntries()!=0){//if per chiedere istogramma non vuoto
      int binmax= zIntersecHisto->GetMaximumBin(); //numero del bin col primo massimo dell'istogramma
      double zmax=(binmax-1)*step+zMin+step/2;//zMin è limite inferiore di range dell'ascissa di istogramma; centro di bin (in z) dove si ha primo massimo
      
      //calcolo zrec come media dei valori di intersezione (aggiunti precedentemente in vector) in un range attorno a zmax
      double sum=0.;
      int denom=0;
      int size=Zintersection.size();
      for(int i=0;i<size;i++){
	  if(TMath::Abs(Zintersection.at(i)-zmax) <kRange){ //scelta range da osservazione degli istogrammi
	    sum=sum+Zintersection.at(i); 
	    denom++;
	  }
      }
      if(denom>0){
	  zrec=sum/denom; 
	  double scarti=0; //somma di scarti quadratici sui valori mediati
	  for(int i=0;i<size;i++){
	    if(TMath::Abs(Zintersection.at(i)-zmax)<kRange)
	    scarti=scarti+((Zintersection.at(i)-zrec)*(Zintersection.at(i)-zrec));
	  }
	  if (denom>1) zrecrms=TMath::Sqrt(scarti/(denom-1)); //calcolo deviazione standard
	  else {zrecrms=0.05;} //0.05 scelto come errore plausibile nel caso di una media calcolata su un singolo valore
	  if (ev%kVerbosity==0) cout<<"zrec= "<<zrec<<" +- "<<zrecrms<<endl;
	  vzrec->push_back(zrec);//riempimento vector
	  vzTrue->push_back(verlec->GetZ());	
	  vzmulti->push_back(verlec->GetM());
	  vzrecrms->push_back(zrecrms); 
	  sizevec++;//contatore per tenere conto delle dimensioni dei vector
      }//fine if denom>0
    }//fine richiesta istogramma non vuoto
   
    if (ev%kVerbosity==0){new TCanvas; zIntersecHisto->DrawCopy();} 
    if(verlec->GetZ()<1*sigmaZ&&verlec->GetM()<multibin[limit]) tot1multi->Fill(verlec->GetM());  
    if(verlec->GetZ()<3*sigmaZ&&verlec->GetM()<multibin[limit]) tot3multi->Fill(verlec->GetM());    
    totz->Fill(verlec->GetZ());//riempimento istogrammi con eventi totali simulati (efficenza=buoni/totali)
      
    //riempimento tree e pulizia TClone		  
    treeReco->Fill();
    clone1->Clear();
    clone2->Clear();
    intersecPointDet1->Clear();
    intersecPointDet2->Clear();					  	
  
  }//fine loop su eventi
  
  cout<<"\n \n Zintersection non trovata correttamente per "<< noTrack<<" Tracklet su "<<tottraccia<<" Trackelet costruite"<<endl; 
   
  plot(vzTrue, vzrec,  vzmulti, vzrecrms, multibin, arraylenghtmulti, tot1multi, tot3multi, limit, sigmaZ,zbin,  arraylenghtz, totz, sizevec);
  
  delete  vzTrue; //eliminazione vector
  delete  vzrec;
  delete  vzmulti;	
  delete  vzrecrms; 
	
  fileReco->cd(); 
  treeReco->Write(); //scrittura tree su file e chiusura file
	
  delete intersecPointDet1;
  delete intersecPointDet2;
  delete clone1;
  delete clone2;
  delete verlec;
  delete vertexRec;
	
  fileReco->Close();
  file->Close();
  
  cout<<"reconstruction completata"<<endl;	




void plot(const vector <double>* const vzTruep,const  vector <double>* const vzrecp, const vector <int> * const vzmultip, const vector <double> * const vzrecrmsp, const double * const multibin, const int arraylenghtmulti, TH1D *tot1multi, TH1D *tot3multi, const int limit, const double sigmaZ, const double * const zbin, const int arraylenghtz, TH1D *totz, const int sizetrue){

  TFile *histo = TFile::Open(kHisto, "RECREATE");  //file per salvare istogrammi 
  const vector <double> &vzTrue = *vzTruep;
  const vector <double> &vzrec = *vzrecp;
  const vector <int> &vzmulti = *vzmultip;	
  const vector <double> &vzrecrms = *vzrecrmsp;
	
  char nome[200]; //nome e titolo istogrammi
  char titolo[200];
  
  //istogrammi dei residui 
  TH1D *isto_res[3];
  int min_mult1=3; //range molteplicità per istogrammi
  int max_mult1=5;
  int min_mult2=45;
  int max_mult2=55;
  //residui con tutte le molteplicità
  snprintf(nome, sizeof(nome), "isto_restot");
  snprintf(titolo, sizeof(titolo), "Residui con tutte le molteplicità");
  isto_res[0]=new TH1D(nome,titolo,500, -0.1, 0.1);
  //residui con molteplicità tra 3 e 5
  snprintf(nome, sizeof(nome), "histo_res1");
  snprintf(titolo, sizeof(titolo), "%i <= molteplicita' <= %i",min_mult1,max_mult1);
  isto_res[1]=new TH1D(nome,titolo,500, -0.1, 0.1);
  //residui con molteplicità tra 45 e 55
  snprintf(nome, sizeof(nome), "histo_res2");
  snprintf(titolo, sizeof(titolo), "%i <= molteplicita' <= %i",min_mult2,max_mult2);
  isto_res[2]=new TH1D(nome,titolo,500, -0.1, 0.1);
  
  for(int i=0;i<sizetrue;i++){
    int multi = vzmulti[i];
    if (TMath::Abs(vzTrue[i] -vzrec[i])<3*vzrecrms[i]){ //così riempiamo solo con gli eventi ricostruiti discretamente e non con quelli con residui enormi dovuti al prendere il picco sbagliato dell'istogramma
       isto_res[0]->Fill(vzrec[i]-vzTrue[i]);
       if(multi>=min_mult1 && multi<=max_mult1) isto_res[1]->Fill(vzrec[i]-vzTrue[i]);
       if(multi>=min_mult2 && multi<=max_mult2) isto_res[2]->Fill(vzrec[i]-vzTrue[i]);
    }
  }
 
  //istogrammi eventi ben ricostruiti per efficenza vs molteplicità 
  TH1D *isto_ok[2];
  snprintf(nome, sizeof(nome), "isto_ok1sigma");
  snprintf(titolo, sizeof(titolo), "ok 1 sigma");
  isto_ok[0]= new TH1D (nome,titolo,arraylenghtmulti,multibin);
  snprintf(nome, sizeof(nome), "isto_ok3sigma");
  snprintf(titolo, sizeof(titolo), "ok 3 sigma");
  isto_ok[1]= new TH1D (nome,titolo,arraylenghtmulti,multibin);
 
  for(int i=0;i<sizetrue;i++){
    if(TMath::Abs(vzTrue[i]-vzrec[i])<3*vzrecrms[i]&&vzTrue[i]<1*sigmaZ&&vzmulti[i]<multibin[limit]) isto_ok[0]->Fill(vzmulti[i]);
    if(TMath::Abs(vzTrue[i]-vzrec[i])<3*vzrecrms[i]&&vzTrue[i]<3*sigmaZ&&vzmulti[i]<multibin[limit]) isto_ok[1]->Fill(vzmulti[i]);
  }
  	
  //istogrammi efficienza vs molteplicità
  TH1D * isto_eff[2];
  //z fra [-sigma;sigma]
  snprintf(nome, sizeof(nome), "isto_eff1sigma");
  snprintf(titolo, sizeof(titolo), "efficienza 1 sigma");
  isto_eff[0]=new TH1D(nome,titolo,arraylenghtmulti,multibin);
  //z fra [-3 sigma;3 sigma]
  snprintf(nome, sizeof(nome), "isto_eff3sigma");
  snprintf(titolo, sizeof(titolo), "efficienza 3 sigma");
  isto_eff[1]=new TH1D(nome,titolo,arraylenghtmulti,multibin);
  
  for(int i=1;i<arraylenghtmulti+1;i++){
    double e;
    double error;
    if(tot1multi->GetBinContent(i)!=0){
      e=isto_ok[0]->GetBinContent (i)/tot1multi->GetBinContent(i);// efficenza =buoni/totali
      error=TMath::Sqrt(e*(1-e)/tot1multi->GetBinContent(i));
      if (1/tot1multi->GetBinContent(i)>error) error=1/tot1multi->GetBinContent(i);//per evitare errore =0 per efficenza e=1 o e=0
      isto_eff[0]->SetBinContent(i,e);
      isto_eff[0]->SetBinError(i,error); 	  
    }
    if(tot3multi->GetBinContent(i)!=0){
      e=isto_ok[1]->GetBinContent (i)/tot3multi->GetBinContent(i);
      error=TMath::Sqrt(e*(1-e)/tot3multi->GetBinContent(i));
      if (1/tot3multi->GetBinContent(i)>error) error=1/tot3multi->GetBinContent(i);
      isto_eff[1]->SetBinContent(i,e);
      isto_eff[1]->SetBinError(i,error);	  
    }
  }
  
  //istogrammi residui bin per bin per risoluzione vs. molteplicità per zTrue entro 1 sigma da 0
  TH1D * isto_residui1[arraylenghtmulti];
  //analogo per 3 sigma
  TH1D * isto_residui3[arraylenghtmulti];	
  for(int i=0;i<arraylenghtmulti;i++){	
    snprintf(nome, sizeof(nome), "isto_residui1[%i]",i);
    snprintf(titolo, sizeof(titolo), "residui 1 sigma per %f<molteplicita'<%f", multibin[i],multibin[i+1]);
    isto_residui1[i]=new TH1D(nome,titolo,500, -0.1, 0.1);	 
    snprintf(nome, sizeof(nome), "isto_residui3[%i]",i);
    snprintf(titolo, sizeof(titolo), "residui 3 sigma per %f<moteplicita'<%f", multibin[i],multibin[i+1]);
    isto_residui3[i]=new TH1D(nome,titolo,500, -0.1, 0.1);
  
    for(int j=0;j<sizetrue;j++){
      if(TMath::Abs(vzTrue[j]-vzrec[j])<3*vzrecrms[j]&&vzTrue[j]<1*sigmaZ&&vzmulti[j]<multibin[limit]&&vzmulti[j]>multibin[i]&&vzmulti[j]<=multibin[i+1])isto_residui1[i]->Fill(vzrec[j]-vzTrue[j]);
      if(TMath::Abs(vzTrue[j]-vzrec[j])<3*vzrecrms[j]&&vzTrue[j]<3*sigmaZ&&vzmulti[j]<multibin[limit]&&vzmulti[j]>multibin[i]&&vzmulti[j]<=multibin[i+1])isto_residui3[i]->Fill(vzrec[j]-vzTrue[j]);
    }
  }
  
  //istogrammi risoluzione vs molteplicità 1 sigma e 3 sigma
  TH1D * isto_ris[2];
  //z fra [-sigma;sigma]
  snprintf(nome, sizeof(nome), "isto_ris1sigma");
  snprintf(titolo, sizeof(titolo), "risoluzione 1 sigma");
  isto_ris[0]=new TH1D(nome,titolo,arraylenghtmulti,multibin);
  //z fra [-3 sigma;3 sigma]
  snprintf(nome, sizeof(nome), "isto_ris3sigma");
  snprintf(titolo, sizeof(titolo), "risoluzione 3 sigma");
  isto_ris[1]=new TH1D(nome,titolo,arraylenghtmulti,multibin);
  
  for(int i=1;i<arraylenghtmulti+1;i++){
    double RMS;//risoluzione è RMS istogrammi residui bin per bin
    double error;
    if(isto_residui1[i-1]->GetRMS()!=0){
      RMS=isto_residui1[i-1]->GetRMS();
      error=isto_residui1[i-1]->GetRMSError();
     
      isto_ris[0]->SetBinContent(i,RMS);
      isto_ris[0]->SetBinError(i,error); 	  
    }
    if(isto_residui3[i-1]->GetRMS()!=0){
      RMS=isto_residui3[i-1]->GetRMS();
      error=isto_residui3[i-1]->GetRMSError();
     
      isto_ris[1]->SetBinContent(i,RMS);
      isto_ris[1]->SetBinError(i,error); 	  
    }
  }
  
  //analogo per istogrammi in funzione zTrue	
	
  //istogrammi eventi ben ricostruiti per efficenza vs zTrue 
  TH1D *isto_okzTrue[1];
  snprintf(nome, sizeof(nome), "isto_okzTrue");
  snprintf(titolo, sizeof(titolo), "ok zTrue");
  isto_okzTrue[0]= new TH1D (nome,titolo,arraylenghtz,zbin);
 

  for(int i=0;i<sizetrue;i++){
   if(TMath::Abs(vzTrue[i]-vzrec[i])<3*vzrecrms[i]) isto_okzTrue[0]->Fill(vzTrue[i]);
  }
 
  //istogrammi efficienza vs zTrue
  TH1D * isto_effzTrue[1];
  snprintf(nome, sizeof(nome), "isto_effzTrue");
  snprintf(titolo, sizeof(titolo), "efficienza vs zTrue ");
  isto_effzTrue[0]=new TH1D(nome,titolo,arraylenghtz,zbin);
  
  for(int i=1;i<arraylenghtz+1;i++){
    double e;
    double error;
    if(totz->GetBinContent(i)!=0){
      e=isto_okzTrue[0]->GetBinContent(i)/totz->GetBinContent(i);
      error=TMath::Sqrt(e*(1-e)/totz->GetBinContent(i));
      if(1/totz->GetBinContent(i)>error) error=1/totz->GetBinContent(i);
      isto_effzTrue[0]->SetBinContent(i,e);
      isto_effzTrue[0]->SetBinError(i,error); 	  
    }
  }      
  
  //Istogramma Per Verificare Distribuzione residui complessiva su tutti i zTrue
  TH1D * isto_ressomma=new TH1D("isto_ressomma","somma di istogrammi residui bin per bin di zTrue; z_rec-z_true [cm]; numero di eventi",500, -0.1, 0.1);
  
  //istogrammi residui bin per bin per risoluzione vs zTrue 
  TH1D * isto_residuizTrue[arraylenghtz];
  
  for (int i=0;i<arraylenghtz;i++){	
    snprintf(nome, sizeof(nome), "isto_residuizTrue[%i]",i);
    snprintf(titolo, sizeof(titolo), "residui  per %f<zTrue <=%f",zbin[i],zbin[i+1]);
    isto_residuizTrue[i]=new TH1D(nome,titolo,500, -0.1, 0.1);	 
    for(int j=0;j<sizetrue;j++){
      if (TMath::Abs(vzTrue[j]-vzrec[j])<3*vzrecrms[j]&&vzTrue[j]>zbin[i]&&vzTrue[j]<=zbin[i+1]){isto_residuizTrue[i]->Fill(vzrec[j]-vzTrue[j]);}
    }
    isto_ressomma->Add(isto_residuizTrue[i]);
  }
 
  //istogrammi risoluzione vs zTrue
  TH1D * isto_riszTrue[1];
  snprintf(nome, sizeof(nome), "isto_riszTrue");
  snprintf(titolo, sizeof(titolo), "risoluzione vs zTrue ");
  isto_riszTrue[0]=new TH1D(nome,titolo,arraylenghtz,zbin);
 
  for(int i=1;i<arraylenghtz+1;i++){
    double RMS;
    double error;
    if(isto_residuizTrue[i-1]->GetRMS()!=0){
      RMS=isto_residuizTrue[i-1]->GetRMS();
      error=isto_residuizTrue[i-1]->GetRMSError();
      isto_riszTrue[0]->SetBinContent(i,RMS);
      isto_riszTrue[0]->SetBinError(i,error); 	  
    }
  }
  
  histo->Write();//scrittura istogrammi su file
  histo->ls();//controllo contenuto file
  histo->Close();	
   
}

  
