//root includes 
#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TStopwatch.h"

//headers includes
#include "Track.h"
#include "Hit.h"
#include "Vertex.h"

using namespace std;

// Global constants
const double kRbp = 3.04;       // Beam pipe radius (half width)
const double kR1 = 4.01;        // Radius of first detector (T1) (half width)
const double kR2 = 7.01;        // Radius of second detector (T2) (half width)
const double kL = 27;           // Detector length

const double kSigmax = 0.01;    // Gaussian spread (sigma) of the vertex around (0,0,0)
const double kSigmay = 0.01;
const double kSigmaz = 5.3;

const int kEvents = 1000000;    // Number of primary vertices
const double kEtamin = -2.;     // Boundaries for eta distribution for track generation
const double kEtamax = +2.;

const TString kVar = "no";     //Gaus method for extraction ("yes" = Box-Muller, "no" = ROOT)
const TString kMul = "yes";     // Multiplicity option ("yes" = custom distribution, "no" = uniform)
const bool kMs = true;          // Multiple scattering enabled (true = ON)
const int kVerbositySimu = 1000000;   // Verbosity level - every how many events ti gives a printout
const unsigned int kSeed = 18;  // Random seed
const TString kFile = "simulation.root";  // Output file

////////////////////////////////////////////////////////////////////
TH1D* limit(TH1D *histo, double min, double max){
  //modifying histo limits
  TAxis *x = histo->GetXaxis();
  double step = x->GetBinWidth(1);
  int b1 = x->FindBin(min);
  int b2 = x->FindBin(max); 
  double xLow = x->GetBinLowEdge(b1);
  double xHigh = x->GetBinUpEdge(b2);
  int nBins = b2-b1+1;
  TH1D *modifiedHisto = new TH1D("histo", "histo", nBins, xLow, xHigh);
  int j=1;
  for(int i=b1;i<=b2;i++)
    modifiedHisto->SetBinContent(j++, histo->GetBinContent(i));
  return modifiedHisto;
}

////////////////////////////////////////////////////////////////////

void Simulation(){
  TStopwatch clock;
  if (kR1 <= kRbp || kR2 <= kR1 || kRbp <= 0 || kEvents <= 0 || kEtamax <= kEtamin || kL <= 0 || kSigmaz <= 0 || kSigmax <= 0 || kSigmay <= 0 || kVerbositySimu <= 0) {
    cout << "Error in parameters setting" << endl;
    return;
  }

  gRandom->SetSeed();
#ifdef KINEM
  TFile *inputFile = TFile::Open("kinem.root", "READ");
  if(!inputFile || inputFile->IsZombie()){
    cout << "Error: Could not open kinem.root" << endl;
    return;
  }
  TH1D *etaHist = limit((TH1D*)inputFile->Get("heta"),kEtamin, kEtamax);
  TH1D *multHist = (TH1D*)inputFile->Get("hmul");
  if(!etaHist || !multHist){
    cout << "Error: Could not load histograms" << endl;
    inputFile->Close();
    return;
  }  
#endif

  //output file for reconstruction purposes
  char choice;
  if (kMul=="No"||kMul=="NO"||kMul=="no") choice='N'; else choice='S';
  FILE *hdata = fopen("data.txt","w");
  fprintf (hdata, "%f %f %f %f %d %c", (float) kR1, (float) kR2, (float) kL, (float) kSigmaz, kEvents, choice);	
  fclose (hdata);
  
  cout<<"Simulation Started "<<endl;

    // Output file setup
  TFile *hfile = TFile::Open(kFile, "RECREATE");
  TTree *tree = new TTree("tree", "events");
  Vertex *vertex = nullptr;
  TClonesArray *bp = new TClonesArray("Hit", 100);        // Beam pipe intersections
  TClonesArray *clone1 = new TClonesArray("Hit", 100);    // First detector intersections
  TClonesArray *clone2 = new TClonesArray("Hit", 100);    // Second detector intersections
  
  // Tree branches
  tree->Branch("Vertex", &vertex);
  tree->Branch("tclonebp", &bp);
  tree->Branch("tclone1", &clone1);
  tree->Branch("tclone2", &clone2);

  int vertexNoHit = 0; // Tracks number of vertices without interaction
  tree->SetDirectory(hfile); // Set output file for tree storage
  
  for(int i=0;i<kEvents;i++){
#ifdef VERBOSITY
  cout<<"event: "<<i+1<<" completed"<<endl;
#endif

  //creation of a vertex
  vertex = new Vertex();
  
  
  //storing in TClonesArray

  }

  cout<<"Simulation Completed "<<endl;

  clock.Stop();
  clock.Print();
}
