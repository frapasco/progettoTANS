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

const TString kMul = "yes";     // Multiplicity option ("yes" = custom distribution, "no" = uniform)
const bool kMs = true;          // Multiple scattering enabled (true = ON)
const int kVerbositySimu = 10000000;   // Verbosity level - every how many events ti gives a printout
const unsigned int kSeed = 18;  // Random seed
const TString kFile = "simulation.root";  // Output file

//oooOOOoooOOOoooOOOoooOOOoooOOOoooOOOoooOOOoooOOOoooOOOoooOOOooo

void simulation() {
  TStopwatch clock;
  if (kR1 <= kRbp || kR2 <= kR1 || kRbp <= 0 || kEvents <= 0 || kEtamax <= kEtamin || kL <= 0 || kSigmaz <= 0 || kSigmax <= 0 || kSigmay <= 0 || kVerbositySimu <= 0) {
    cout << "Error in parameters setting" << endl;
    return;
  }
  
  gRandom->SetSeed(kSeed);
  
  // Open the histogram file and retrieve eta and multiplicity histograms
  TFile *inputFile = TFile::Open("kinem.root", "READ");
  if (!inputFile || inputFile->IsZombie()) {
    cout << "Error: Could not open kinem.root" << endl;
    return;
  }
  TH1D *etaHist = (TH1D*)inputFile->Get("heta");
  TH1D *multHist = (TH1D*)inputFile->Get("hmul");
  if (!etaHist || !multHist) {
    cout << "Error: Could not load histograms" << endl;
    inputFile->Close();
    return;
  }
  
  //output file for reconstruction purposes
  char a;
  if (kMul=="No"||kMul=="NO"||kMul=="no") a='N'; else a='S';
  FILE *hdata = fopen("data.txt","w");
  fprintf (hdata, "%f %f %f %f %d %c", (float) kR1, (float) kR2, (float) kL, (float) kSigmaz, kEvents, a);	
  fclose (hdata);
  
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
  
  // Main event loop
  for (int i = 0; i < kEvents; i++) {
    int lost1 = 0, lost2 = 0, lostAll = 0;
    if (i % kVerbositySimu == 0) cout << "\n \n  EVENT  " << i + 1 << endl;
    
    // Initialize vertex with preloaded multiplicity histogram
    vertex = new Vertex(kMul, multHist, kSigmax, kSigmay, kSigmaz);
    
    // Set initial direction using the eta histogram
    vertex->InitialDir(kEtamin, kEtamax, etaHist);
    
    int nParticle = vertex->GetM(); // Retrieve particle multiplicity for the vertex
    if (i % kVerbositySimu == 0) cout << "Vertex created with " << nParticle << " particles" << endl;
    
    // Loop over each generated particle
    for (int j = 0; j < nParticle; j++) {
      // Propagate each track through the beam pipe and detectors
      Track *tbp = new Track(vertex->GetTheta(), vertex->GetPhi());
      double Xbp, Ybp, Zbp, X1, Y1, Z1, X2, Y2, Z2;
      
      // Check for intersection with beam pipe
      if (tbp->Intersection(vertex->GetX(), vertex->GetY(), vertex->GetZ(), Xbp, Ybp, Zbp, kRbp, kL)) {
	new ((*bp)[j]) Hit(Xbp, Ybp, Zbp, j);
	tbp->MultipleScattering(kMs); // Apply multiple scattering
	
	Track *t1 = new Track(tbp->GetTheta(), tbp->GetPhi()); // Propagate to T1
	if (t1->Intersection(Xbp, Ybp, Zbp, X1, Y1, Z1, kR1, kL)) {
	  new ((*clone1)[j - lost1]) Hit(X1, Y1, Z1, j);
	  t1->MultipleScattering(kMs); // Apply scattering for T1
	  
	  Track *t2 = new Track(t1->GetTheta(), t1->GetPhi()); // Propagate to T2
	  if (t2->Intersection(X1, Y1, Z1, X2, Y2, Z2, kR2, kL)) {
	    new ((*clone2)[j - lost2]) Hit(X2, Y2, Z2, j);
	  } else {
	    lost2++;
	  }
	  delete t2;
	} else { // Check for intersection with T2 only
	  if (t1->Intersection(Xbp, Ybp, Zbp, X2, Y2, Z2, kR2, kL)) {
	    new ((*clone2)[j - lost2]) Hit(X2, Y2, Z2, j);
	  } else {
	    lost2++;
	    lostAll++;
	  }
	  lost1++;
	}
	delete t1;
      }
      delete tbp;
    }
    
    if (lostAll == nParticle) vertexNoHit++; // Track vertices without detector interaction
    tree->Fill(); // Fill the tree with the event data
    
    clone1->Clear(); // Clear arrays for the next event
    clone2->Clear();
    bp->Clear();
    delete vertex; // Delete the current vertex
  }
  
  // Summary output
  cout << endl;
  cout << "Out of " << kEvents << " vertices generated, " << vertexNoHit << " didn't interact with any of the detectors" << endl;
  
  // Write the tree to the output file and clean up
  hfile->cd();
  tree->Write();
  hfile->Close();
  inputFile->Close();
  
  // Free memory
  delete bp;
  delete clone1;
  delete clone2;
  delete hfile;
  delete inputFile;
  
  cout << "Simulation completed" << endl;
  clock.Stop();
  clock.Print();
}
