//root includes 
#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include <stdio.h>
#include <stdlib.h>

//headers includes
#include "Track.h"
#include "Hit.h"
#include "Vertex.h"

using namespace std;
//measurements are given in cm and radii are taken at half width
//global variables
const double kRbp=3.04; //beam pipe radius (half width)
const double kR1=4.01; //radius of first detector (T1) (half width)
const double kR2=7.01; //radius of second detector (T2) (half width)
const double kL=27; //detector lenght

const double kSigmax=0.01; //Gaussian spread (sigma) of the vertex around (0,0,0)
const double kSigmay=0.01;
const double kSigmaz=5.3;

const int kEvents=20; //number of primary vertexes
const double kEtamin=-2.;//boundaries for the eta distribution histo for generating tracks
const double kEtamax=+2.;

const TString kMul="sì"; //molteplicity options ("sì"=given distrib.; "No","NO","no"=uniform)
const bool kMs=true;//multiple scattering (true=ON)
const int kVerbositySimu=1;//verbosity
const unsigned int kSeed=18; //seed for TRandom3

const TString kFile="simulation.root";//output file

////////////////////////////////////////////////////////////////////

void simulation(){
  if(kR1<=kRbp||kR2<=kR1||kRbp<=0||kEvents<=0||kEtamax<=kEtamin||kL<=0||kSigmaz<=0||kSigmax<=0||kSigmay<=0||kVerbositySimu<=0){
    cout<<"Error in parameters setting"<<endl; 
    return;
  }
  
  gRandom->SetSeed(kSeed);
  
  //file and histograms for initial distributions
  TString fileName="kinem.root";
  TString multDistr="hmul";
  TString etaDistr="heta";
  char a;
  if (kMul=="No"||kMul=="NO"||kMul=="no") a='N'; else a='S';
  
  //storing arrays for intersection points
  TClonesArray *bp = new TClonesArray("Hit",100);//beam pipe
  TClonesArray *clone1 = new TClonesArray("Hit",100);//first detector
  TClonesArray *clone2 = new TClonesArray("Hit",100);//second detector
  TClonesArray &intbp = *bp;
  TClonesArray &int1 = *clone1;  
  TClonesArray &int2 = *clone2;
  
  int nParticle;//number of produced particles
  
  //output file for reconstruction purposes
  FILE *hdata = fopen("data.txt","w");
  fprintf (hdata, "%f %f %f %f %d %c", (float) kR1, (float) kR2, (float) kL, (float) kSigmaz, kEvents, a);	
  fclose (hdata);
  
  //creation of a Vertex
  Vertex *vertex = NULL;
  
  //creation TFile, TTree and Branches
  TFile  *hfile = TFile::Open(kFile,"RECREATE");
  TTree *tree = new TTree("tree", "eventi");
  tree->Branch("Vertext",&vertex); //PV
  tree->Branch("tclonebp",&bp); //intersection BP
  tree->Branch("tclone1",&clone1); //intersection T1
  tree->Branch("tclone2",&clone2); //intersection T2
  
  //hit scoring variables
  double X1,Y1,Z1;
  double X2,Y2,Z2;
  double Xbp,Ybp,Zbp;
  int lost1, lost2, lostAll; //no interactions with: T1, T2, both
  int vertexNoHit=0; //total number of vertexes without interactions
  
  tree->SetDirectory(hfile); //setting for correct file output
  
  //loop on vertexes
  for(int i=0; i<kEvents;i++){   
    lost1=0;
    lost2=0;
    lostAll=0;
    
    if (i%kVerbositySimu==0) cout <<"\n \n  EVENT  "<<i+1<<endl;	   
    vertex=new Vertex(kMul,fileName,multDistr,kSigmax,kSigmay,kSigmaz); //vertex creation
    
    //oooOOOoooOOOoooOOOoooOOOoooOOOooo
    //AAAAAAAAAA write a method that assignes the multiplicity via external loading of the TFILE
    nParticle=vertex->GetM(); 
      if (i%kVerbositySimu==0) cout <<"vertex created"<<endl;
      //loop on charged particles created  
      for (int j=0; j<nParticle; j++){  //loop on created particles
        vertex->InitialDir(kEtamin,kEtamax,fileName,etaDistr);//setting of initial direction
	Track *tbp = new Track(vertex->GetTheta(),vertex->GetPhi()); //propagating towards beam pipe
	
	//intersection with BP
        if (tbp->Intersection(vertex->GetX(),vertex->GetY(),vertex->GetZ(),Xbp,Ybp,Zbp,kRbp,kL)){}; //the lenght of the BP is infinite but the only intersections we keep are the ones in (-L/2,L/2)
	new (intbp[j]) Hit (Xbp,Ybp,Zbp,j) ;
        tbp->MultipleScattering(kMs);//multiple scatter
	Track *t1 = new Track(tbp->GetTheta(),tbp->GetPhi()); //propagating towards T1
     	
	//intersection with T1
	if(t1->Intersection(Xbp,Ybp,Zbp,X1,Y1,Z1,kR1,kL)){  
          new (int1[j-lost1]) Hit (X1,Y1,Z1,j); 
	  t1->MultipleScattering(kMs);
	  Track *t2 = new Track(t1->GetTheta(),t1->GetPhi()); //propating towards T2
	  
	//intersection with T2
	  if(t2->Intersection(X1,Y1,Z1,X2,Y2,Z2,kR2,kL)){   
            new (int2[j-lost2]) Hit (X2,Y2,Z2,j);
	  }else{
	    lost2++;
	  }  
	  
	  delete t2;
	}else{//checking for intersections with T2 but not T1
          if(t1->Intersection(Xbp,Ybp,Zbp,X2,Y2,Z2,kR2,kL)){   
            new (int2[j-lost2]) Hit (X2,Y2,Z2,j);
	  }else {
            lost2++; 
            lostAll++;
          }
	  lost1++;
	}
	
	delete t1;     
	delete tbp;			
      }
      
      if (i%kVerbositySimu==0) cout<<"Hits done"<<endl; 
      
      //scoring hits
      tree->Fill(); 
      if (lostAll==nParticle) vertexNoHit++; 
      clone1->Clear(); //clearing TClones for next iteration
      clone2->Clear(); 
      bp->Clear();
      delete vertex; //deleting the Vertex
  }
  cout<<endl;
  cout<<"Out of "<<kEvents<<" vertexes generated, "<<vertexNoHit<<" didn't interact with any of the detectors"<<endl;	
  hfile->cd();
  tree->Write(); //writing on the output file
  delete bp;
  delete clone1;
  delete clone2;
  hfile->Close();
  
  cout<<"Simulation completed"<<endl;
}
