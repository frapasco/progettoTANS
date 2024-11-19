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
#include "TStopwatch.h"

//headers include
#include "Vertex.h"
#include "Track.h"
#include "Hit.h"

using namespace std;
const int kVerbosityReco=10000000; //verbosity
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

void plot(const vector <double>* const vzTruep,const  vector <double>* const vzgrecp, const vector <int> * const vzmultip, const vector <double> * const vzrecrmsp, const double * const multiBin, const int arrayLenghtMulti, TH1D *tot1multi, TH1D *tot3multi, const int limit, const double sigmaZ, const double * const zbin, const int arrayLenghtZ, TH1D *totz, const int sizetrue);

////////////////////////////////////////////////////////////////////
//method for smearing points, generating noise and reconstruct the primary vertex out of tracklets
////////////////////////////////////////////////////////////////////

void reconstruction(){
  TStopwatch clock;
  
  double step=0.05; //step histograms intersections with z
  double R1,R2,L,sigmaZ;
  int events;
  char a;
  //reading data.txt to get L, R1, R2, sigmaZ, events and multiplicity generation choice	
  ifstream In = ifstream(kData);
  if(!In.is_open()){
    cout<<"Error: the file data doesn't exist"<<endl; 
    return;
  }
  
  In>>R1>>R2>>L>>sigmaZ>>events>>a;
  In.close();
  TString mul;
  if (a=='N') mul="No"; else mul="yes";
  TFile *file = TFile::Open(kSim, "READ");
  if(!file || file->IsZombie()){
    cout<<"Error: the file doesn't exist"<<endl; 
    return;
  }
  file->ls(); //reading the content of file
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
   if(step<=0||kVerbosityReco<=0||kNoisefrac<=0||kPhiMax<=0||kRange<=0){
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
   treeReco->Branch("zintersectionvector",&Zintersection); //vector of Tracklet intersections 
   treeReco->Branch("int1enoise",&intersecPointDet1); //intersection + noise T1
   treeReco->Branch("int2enoise",&intersecPointDet2); //intersection + noise T2
   //reading from input file TCloneArrays
   TClonesArray *clone1= new TClonesArray("Hit",100);
   TClonesArray *clone2= new TClonesArray("Hit",100);
   for (int i=0;i<10;i++){ 
     new ((*clone1)[i]) Hit(); 
     new ((*clone2)[i]) Hit();
   } //in order to having a non zero dimension obj

   Vertex *verlec = new Vertex(); //creation of a Vertex obj
   inputTree->GetBranch("tclone1")->SetAutoDelete(kFALSE); //prevents TClonesArray from reusing space allocated by the previous obj, setted for redundancy since it is kFALSE by default
   inputTree->SetBranchAddress("tclone1",&clone1);
   inputTree->GetBranch("tclone2")->SetAutoDelete(kFALSE);
   inputTree->SetBranchAddress("tclone2",&clone2);
   TBranch *b3=inputTree->GetBranch("Vertex");
   b3->SetAddress(&verlec); //reading Vertex
   clone1->Clear(); //cleaning clone1 e clone2
   clone2->Clear();

   //histogram for efficiency (multiplicity)
   double multiBin[]={-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9.5, 11.5, 14.5, 17.5, 20.5, 25.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5, 91.5};
   int arrayLenghtMulti= sizeof(multiBin)/sizeof(multiBin[0])-1; //nbins of histo
   int limit; //maximum multuplicity for histo with multiBin[limit]
   if (kLimit1>kLimit2||kLimit2>arrayLenghtMulti||kLimit1<=0||kLimit2<=0){
     cout<<"Error on kLimit1 and kLimit2"<<endl;
     return;
   }
   if(mul=="No") limit=kLimit2; else limit=kLimit1; //if not uniform multi we stop before

   TH1D *tot1multi = new TH1D("tot1multi", "tot1multi", arrayLenghtMulti, multiBin); //histo for events with z simu between 1 sigma from 0 in function of multiplicity
   TH1D *tot3multi = new TH1D("tot3multi", "tot3multi", arrayLenghtMulti, multiBin); //between 3 sigma

   //same for efficiency (zTrue=z of primary vertex)	
   double zbin[]={-16, -12, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 12, 16};
   int arrayLenghtZ= sizeof(zbin)/sizeof(zbin[0])-1; //nbins
   TH1D *totz = new TH1D("totz","totz",arrayLenghtZ,zbin);

   //final histograms ingredients
   vector <double> *vzTrue = new vector <double>; //z of the simulated vertex
   vector <double> *vzrec = new vector <double>; //z reconstructed
   vector <int> *vzmulti = new vector <int>; //multiplicity (from the simulation)
   vector <double> *vzrecrms = new vector <double>; //error (RMS) on z reconstructed	
   vzTrue->reserve(events); //allocation of memory to avoid memory fragmentation
   vzrec->reserve(events);
   vzmulti->reserve(events);
   vzrecrms->reserve(events);

   //variables
   int nlines1, nlines2;
   int out1; //counters for hits outside T1 and T2 after smeaaring
   int out2; 
   int noTrack=0; //number of tracklets without a good candidate for intersection
   int totTracklets=0;//total numer of tracklets generated
   double x1,y1,z1,x2,y2,z2; 
   double Z; //z of intersection between tracklets (x=0 plane)
   int element1, element2; //number of elements in int1 and int2 (without noise)
   double zrec, zrecrms; //z reco and his uncertainty
   double zmax; //z of the first max between tracklets
   int sizevec=0; //size of the vectors on heap

   //loop on the number of events (inputTree entries)
   for(int ev=0;ev<inputTree->GetEntries();ev++){
     if (ev%kVerbosityReco==0) cout <<"\n \n EVENT  "<<ev+1<<endl;
     
     //preliminary steps
     Zintersection.clear();
     Zintersection.reserve(200);
     out1=0;
     out2=0;
     element1=0;
     element2=0;	  
     delete vertexRec;
     zIntersecHisto->Reset("ICES"); //AAAAAAA why just ICES?
     inputTree->GetEntry(ev); //getting event
     nlines1 = clone1->GetEntriesFast(); //number of events in clone1
     nlines2 = clone2->GetEntriesFast(); //number of events in clone2
     vertexRec = new Vertex(*verlec); //clone of verlec in vertexRec
     
     //reading of intersection points from inputTree and smearing of them
     //loop on simulated tracks
     for(int i=0;i<nlines1;i++){    
       int t=i-out1;
       
       //T1 hit processing
       Hit *i1 = (Hit*) clone1->At(i); 
       i1->Smearing(R1);
       
       if(TMath::Abs(i1->GetZ())<L/2){ //check that after smearing the hit point still lies in T1
	 x1=i1->GetX();
	 y1=i1->GetY() ;
	 z1=i1->GetZ();	 
	 new (int1[t]) Hit(x1,y1,z1,i1->GetLabel());      
	 element1++;      
       }else out1++;   
     }	
     for(int i=0;i<nlines2;i++){   
       int k=i-out2;  
       
       //T2 hit processing      
       Hit *i2 = (Hit*) clone2->At(i); 
       i2->Smearing(R2);
       
       if(TMath::Abs(i2->GetZ())<L/2){ //check that after smearing the hit point still lies in T2
	 x2=i2->GetX();
	 y2=i2->GetY() ;
	 z2=i2->GetZ();	 
	 new (int2[k]) Hit(x2,y2,z2,i2->GetLabel()); 
	 element2++;	
       }else out2++;      
     }	
     if(ev%kVerbosityReco==0)cout<<"smearing of hits done"<<endl;
     if(ev%kVerbosityReco==0)cout<<"Hit outside T1: "<<out1<<"  Hit outside T2: "<<out2<<endl;
     
     //adding noise
     int Nnoise = (int) (verlec->GetM())*kNoisefrac; //assuming the number of hits due to noise is equal on both detectors #T1 = #T2
     for(int j=element1;j<element1+Nnoise;j++){
       new (int1[j]) Hit(R1,L);
     } 
     for(int k=element2;k<element2+Nnoise; k++){
       new (int2[k]) Hit(R2,L);
     }
     if (ev%kVerbosityReco==0)cout<<"noise done"<<endl;
     
     //reconstruction of PV
     int Nrec1=intersecPointDet1->GetEntriesFast(); //is the sum of (hit+noise) in int1 and int2
     int Nrec2=intersecPointDet2->GetEntriesFast();
     
     int ntraccia=0; //tracklets counter for a single vertex
     for(int i=0; i<Nrec1; i++){
       for(int j=0; j<Nrec2; j++){
	 Hit *i1=(Hit*) intersecPointDet1->At(i); 
	 Hit *i2=(Hit*) intersecPointDet2->At(j);
	 if(TMath::Abs(i2->GetPhi()-i1->GetPhi())<kPhiMax){ //Creation of Tracklet for points in phimax range
	   x1=i1->GetX();
	   x2=i2->GetX();
	   y1=i1->GetY();
	   y2=i2->GetY();
	   z1=i1->GetZ();
	   z2=i2->GetZ();
	   Track *t = new Track(x2-x1,y2-y1,z2-z1); 
	   ntraccia++; 
	   totTracklets++;
	   Z=200000000;  //initialized far outside the expected region in order to check the method
	   t->Intersection2(x2,y2,z2,Z); //Z update with intersection of tracklets and x=0 plane
	   if(Z>-10000000 && Z<10000000){ //condition for the intersection method to have worked
	     zIntersecHisto->Fill(Z);
	     Zintersection.push_back(Z);
	   }
	   else noTrack++;   
	   delete t;
	 }
       }
     }
     
     if(ev%kVerbosityReco==0) cout<<"the tracklets are generated and the histo is filled"<<endl;
     if(zIntersecHisto->GetEntries()!=0){ //check that the histo is not empty
       int binmax=zIntersecHisto->GetMaximumBin(); //num bin of the first max of histo
       double zmax=(binmax-1)*step+zMin+step/2; //zmax is the centre of the maximum of the histo (in z), zMin is the lower bound of the histo
       
       //evaluation of zrec as average of the intersection values (in the vector) in a range around zmax
       double sum=0.;
       int denom=0;
       int size=Zintersection.size();
       for(int i=0;i<size;i++){
	 if(TMath::Abs(Zintersection.at(i)-zmax) <kRange){ // range choice is taken from an observation of the histo
	   sum=sum+Zintersection.at(i); 
	   denom++;
	 }
       }
       if(denom>0){
	 zrec=sum/denom; 
	 double deviation=0; //sum of deviation of squared averaged values
	 for(int i=0;i<size;i++){
	   if(TMath::Abs(Zintersection.at(i)-zmax)<kRange)
	     deviation += ((Zintersection.at(i)-zrec)*(Zintersection.at(i)-zrec));
	 }
	 if (denom>1) zrecrms = TMath::Sqrt(deviation/(denom-1)); //standard deviation
	 else {zrecrms=0.05;} //0.05 is chosen as error on an average of a single value
	 if (ev%kVerbosityReco==0) cout<<"zrec= "<<zrec<<" +- "<<zrecrms<<endl;
	 vzrec->push_back(zrec); //filling the vector
	 vzTrue->push_back(verlec->GetZ());	
	 vzmulti->push_back(verlec->GetM());
	 vzrecrms->push_back(zrecrms); 
	 sizevec++; //counter
       }
     }

     if (ev%kVerbosityReco==0){new TCanvas; zIntersecHisto->DrawCopy();} 
     if(verlec->GetZ()<1*sigmaZ && verlec->GetM()<multiBin[limit]) tot1multi->Fill(verlec->GetM());  
     if(verlec->GetZ()<3*sigmaZ && verlec->GetM()<multiBin[limit]) tot3multi->Fill(verlec->GetM());    
     totz->Fill(verlec->GetZ()); //filling the histos with simu events (efficiency=reco/total)

     //writing on TTree and clearing the TClones
     treeReco->Fill();
     clone1->Clear();
     clone2->Clear();
     intersecPointDet1->Clear();
     intersecPointDet2->Clear();					  	

   } //end loop on events

   cout<<"\n \n Zintersection not found correctly for  "<< noTrack<<" Tracklet on "<<totTracklets<<" reconstrucetd Trackelet"<<endl; 

   plot(vzTrue, vzrec,  vzmulti, vzrecrms, multiBin, arrayLenghtMulti, tot1multi, tot3multi, limit, sigmaZ,zbin,  arrayLenghtZ, totz, sizevec);

   delete  vzTrue; //deleting the vectors
   delete  vzrec;
   delete  vzmulti;	
   delete  vzrecrms; 

   fileReco->cd(); 
   treeReco->Write();

   delete intersecPointDet1;
   delete intersecPointDet2;
   delete clone1;
   delete clone2;
   delete verlec;
   delete vertexRec;

   fileReco->Close();
   file->Close();

   cout<<"reconstruction completed"<<endl;
   clock.Stop();
   clock.Print();

 }



 void plot(const vector <double>* const vzTruep,const  vector <double>* const vzrecp, const vector <int> * const vzmultip, const vector <double> * const vzrecrmsp, const double * const multiBin, const int arrayLenghtMulti, TH1D *tot1multi, TH1D *tot3multi, const int limit, const double sigmaZ, const double * const zbin, const int arrayLenghtZ, TH1D *totz, const int sizetrue){

   TFile *histo = TFile::Open(kHisto, "RECREATE");  //histograms file
   const vector <double> &vzTrue = *vzTruep;
   const vector <double> &vzrec = *vzrecp;
   const vector <int> &vzmulti = *vzmultip;	
   const vector <double> &vzrecrms = *vzrecrmsp;

   char name[200]; //names and titles for histo
   char title[200];

   //residual vertexes histo 
   TH1D *histo_res[3];
   int min_mult1=3; //range multiplicity
   int max_mult1=5;
   int min_mult2=45;
   int max_mult2=55;
   
   //Residuals with all multiplicities
   snprintf(name, sizeof(name), "histo_restot");
   snprintf(title, sizeof(title), "Residuals with all multiplicities");
   histo_res[0]=new TH1D(name,title,500, -0.1, 0.1);
   //Residuals with multiplicities between 3 and 5
   snprintf(name, sizeof(name), "histo_res1");
   snprintf(title, sizeof(title), "%i <= multiplicity <= %i",min_mult1,max_mult1);
   histo_res[1]=new TH1D(name,title,500, -0.1, 0.1);
   //Residuals with multiplicities between 45 and 55
   snprintf(name, sizeof(name), "histo_res2");
   snprintf(title, sizeof(title), "%i <= multiplicity <= %i",min_mult2,max_mult2);
   histo_res[2]=new TH1D(name,title,500, -0.1, 0.1);

   for(int i=0;i<sizetrue;i++){
     int multi = vzmulti[i];
     if (TMath::Abs(vzTrue[i] -vzrec[i])<3*vzrecrms[i]){ // so that we fill just with good reconstructions and not with the cases where we took the wrong max
       histo_res[0]->Fill(vzrec[i]-vzTrue[i]);
       if(multi>=min_mult1 && multi<=max_mult1)
	 histo_res[1]->Fill(vzrec[i]-vzTrue[i]);
       if(multi>=min_mult2 && multi<=max_mult2)
	 histo_res[2]->Fill(vzrec[i]-vzTrue[i]);
     }
   }

   //istogrammi eventi ben ricostruiti per efficenza vs molteplicità 
   TH1D *histo_ok[2];
   snprintf(name, sizeof(name), "histo_ok1sigma");
   snprintf(title, sizeof(title), "ok 1 sigma");
   histo_ok[0]= new TH1D (name,title,arrayLenghtMulti,multiBin);
   snprintf(name, sizeof(name), "histo_ok3sigma");
   snprintf(title, sizeof(title), "ok 3 sigma");
   histo_ok[1]= new TH1D (name,title,arrayLenghtMulti,multiBin);

   for(int i=0;i<sizetrue;i++){
     if(TMath::Abs(vzTrue[i]-vzrec[i])<3*vzrecrms[i]&&vzTrue[i]<1*sigmaZ&&vzmulti[i]<multiBin[limit])
       histo_ok[0]->Fill(vzmulti[i]);
     if(TMath::Abs(vzTrue[i]-vzrec[i])<3*vzrecrms[i]&&vzTrue[i]<3*sigmaZ&&vzmulti[i]<multiBin[limit])
       histo_ok[1]->Fill(vzmulti[i]);
   }

   //histo efficiency vs multiplicity
   TH1D * histo_eff[2];
   //z between [-sigma;sigma]
   snprintf(name, sizeof(name), "histo_eff1sigma");
   snprintf(title, sizeof(title), "efficiency 1 sigma");
   histo_eff[0]=new TH1D(name,title,arrayLenghtMulti,multiBin);
   //z between [-3 sigma;3 sigma]
   snprintf(name, sizeof(name), "histo_eff3sigma");
   snprintf(title, sizeof(title), "efficiency 3 sigma");
   histo_eff[1]=new TH1D(name,title,arrayLenghtMulti,multiBin);

   for(int i=1;i<arrayLenghtMulti+1;i++){
     double e;
     double error;
     if(tot1multi->GetBinContent(i)!=0){
       e = histo_ok[0]->GetBinContent(i)/tot1multi->GetBinContent(i); //efficency = good/total
       error=TMath::Sqrt(e*(1-e)/tot1multi->GetBinContent(i));
       if (1/tot1multi->GetBinContent(i)>error) //to avoid error=0 for efficiency=1 || e=0
	 error=1/tot1multi->GetBinContent(i);
       histo_eff[0]->SetBinContent(i,e);
       histo_eff[0]->SetBinError(i,error); 	  
     }
     if(tot3multi->GetBinContent(i)!=0){
       e=histo_ok[1]->GetBinContent (i)/tot3multi->GetBinContent(i);
       error=TMath::Sqrt(e*(1-e)/tot3multi->GetBinContent(i));
       if (1/tot3multi->GetBinContent(i)>error) error=1/tot3multi->GetBinContent(i);
       histo_eff[1]->SetBinContent(i,e);
       histo_eff[1]->SetBinError(i,error);	  
     }
   }
   //histograms with residual 
   //istogrammi residui bin per bin per risoluzione vs. molteplicità per zTrue entro 1 sigma da 0
   TH1D * histo_residui1[arrayLenghtMulti];
   //analogo per 3 sigma
   TH1D * histo_residui3[arrayLenghtMulti];	
   for(int i=0;i<arrayLenghtMulti;i++){	
     snprintf(name, sizeof(name), "histo_residui1[%i]",i);
     snprintf(title, sizeof(title), "residui 1 sigma per %f<molteplicita'<%f", multiBin[i],multiBin[i+1]);
     histo_residui1[i]=new TH1D(name,title,500, -0.1, 0.1);	 
     snprintf(name, sizeof(name), "histo_residui3[%i]",i);
     snprintf(title, sizeof(title), "residui 3 sigma per %f<moteplicita'<%f", multiBin[i],multiBin[i+1]);
     histo_residui3[i]=new TH1D(name,title,500, -0.1, 0.1);

     for(int j=0;j<sizetrue;j++){
       if(TMath::Abs(vzTrue[j]-vzrec[j])<1*vzrecrms[j] && vzTrue[j]<1*sigmaZ && vzmulti[j]<multiBin[limit] && vzmulti[j]>multiBin[i] && vzmulti[j]<=multiBin[i+1])
	 histo_residui1[i]->Fill(vzrec[j]-vzTrue[j]);
       if(TMath::Abs(vzTrue[j]-vzrec[j])<3*vzrecrms[j] && vzTrue[j]<3*sigmaZ && vzmulti[j]<multiBin[limit] && vzmulti[j]>multiBin[i] && vzmulti[j]<=multiBin[i+1])
	 histo_residui3[i]->Fill(vzrec[j]-vzTrue[j]);
     }
   }

   //istogrammi risoluzione vs molteplicità 1 sigma e 3 sigma
   TH1D * histo_ris[2];
   //z fra [-sigma;sigma]
   snprintf(name, sizeof(name), "histo_ris1sigma");
   snprintf(title, sizeof(title), "risoluzione 1 sigma");
   histo_ris[0]=new TH1D(name,title,arrayLenghtMulti,multiBin);
   //z fra [-3 sigma;3 sigma]
   snprintf(name, sizeof(name), "histo_ris3sigma");
   snprintf(title, sizeof(title), "risoluzione 3 sigma");
   histo_ris[1]=new TH1D(name,title,arrayLenghtMulti,multiBin);

   for(int i=1;i<arrayLenghtMulti+1;i++){
     double RMS;//risoluzione è RMS istogrammi residui bin per bin
     double error;
     if(histo_residui1[i-1]->GetRMS()!=0){
       RMS=histo_residui1[i-1]->GetRMS();
       error=histo_residui1[i-1]->GetRMSError();

       histo_ris[0]->SetBinContent(i,RMS);
       histo_ris[0]->SetBinError(i,error); 	  
     }
     if(histo_residui3[i-1]->GetRMS()!=0){
       RMS=histo_residui3[i-1]->GetRMS();
       error=histo_residui3[i-1]->GetRMSError();

       histo_ris[1]->SetBinContent(i,RMS);
       histo_ris[1]->SetBinError(i,error); 	  
     }
   }

   //analogo per istogrammi in funzione zTrue	

   //istogrammi eventi ben ricostruiti per efficenza vs zTrue 
   TH1D *histo_okzTrue[1];
   snprintf(name, sizeof(name), "histo_okzTrue");
   snprintf(title, sizeof(title), "ok zTrue");
   histo_okzTrue[0]= new TH1D (name,title,arrayLenghtZ,zbin);


   for(int i=0;i<sizetrue;i++){
     if(TMath::Abs(vzTrue[i]-vzrec[i])<3*vzrecrms[i]) histo_okzTrue[0]->Fill(vzTrue[i]);
   }

   //istogrammi efficienza vs zTrue
   TH1D * histo_effzTrue[1];
   snprintf(name, sizeof(name), "histo_effzTrue");
   snprintf(title, sizeof(title), "efficienza vs zTrue ");
   histo_effzTrue[0]=new TH1D(name,title,arrayLenghtZ,zbin);

   for(int i=1;i<arrayLenghtZ+1;i++){
     double e;
     double error;
     if(totz->GetBinContent(i)!=0){
       e=histo_okzTrue[0]->GetBinContent(i)/totz->GetBinContent(i);
       error=TMath::Sqrt(e*(1-e)/totz->GetBinContent(i));
       if(1/totz->GetBinContent(i)>error) error=1/totz->GetBinContent(i);
       histo_effzTrue[0]->SetBinContent(i,e);
       histo_effzTrue[0]->SetBinError(i,error); 	  
     }
   }      

   //Istogramma Per Verificare Distribuzione residui complessiva su tutti i zTrue
   TH1D * histo_ressomma=new TH1D("histo_ressomma","somma di istogrammi residui bin per bin di zTrue; z_rec-z_true [cm]; numero di eventi",500, -0.1, 0.1);

   //istogrammi residui bin per bin per risoluzione vs zTrue 
   TH1D * histo_residuizTrue[arrayLenghtZ];

   for (int i=0;i<arrayLenghtZ;i++){	
     snprintf(name, sizeof(name), "histo_residuizTrue[%i]",i);
     snprintf(title, sizeof(title), "residui  per %f<zTrue <=%f",zbin[i],zbin[i+1]);
     histo_residuizTrue[i]=new TH1D(name,title,500, -0.1, 0.1);	 
     for(int j=0;j<sizetrue;j++){
       if (TMath::Abs(vzTrue[j]-vzrec[j])<3*vzrecrms[j]&&vzTrue[j]>zbin[i]&&vzTrue[j]<=zbin[i+1]){histo_residuizTrue[i]->Fill(vzrec[j]-vzTrue[j]);}
     }
     histo_ressomma->Add(histo_residuizTrue[i]);
   }

   //istogrammi risoluzione vs zTrue
   TH1D * histo_riszTrue[1];
   snprintf(name, sizeof(name), "histo_riszTrue");
   snprintf(title, sizeof(title), "risoluzione vs zTrue ");
   histo_riszTrue[0]=new TH1D(name,title,arrayLenghtZ,zbin);

   for(int i=1;i<arrayLenghtZ+1;i++){
     double RMS;
     double error;
     if(histo_residuizTrue[i-1]->GetRMS()!=0){
       RMS=histo_residuizTrue[i-1]->GetRMS();
       error=histo_residuizTrue[i-1]->GetRMSError();
       histo_riszTrue[0]->SetBinContent(i,RMS);
       histo_riszTrue[0]->SetBinError(i,error); 	  
     }
   }

   histo->Write();//scrittura istogrammi su file
   histo->ls();//controllo contenuto file
   histo->Close();	
 }



