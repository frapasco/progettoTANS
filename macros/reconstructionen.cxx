// Root includes
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

// Headers include
#include "Vertex.h"
#include "Track.h"
#include "Hit.h"

using namespace std;

// Constants
const int kVerbosityReco = 10000000; // Verbosity level
const int kLimit1 = 17;               // Last bin for histograms for given multiplicity
const int kLimit2 = 21;               // Last bin for histograms for uniform multiplicity
const double kNoisefrac = 1.2;        // Ratio for noise/multiplicity points
const double kPhiMax = 0.01;          // Max angle for creating tracklets
const double kRange = 0.1;             // Range for averaging intersections to reconstruct z of the vertex

// File names
const TString kSim = "simulation.root"; // Input file
const TString kRic = "reconstruction.root"; // Output file for reconstruction tree
const TString kHisto = "histo.root"; // Output file for histograms
const string kData = "data.txt";     // Config file for run data

// Method for plotting results (residuals, efficiency vs multiplicity, resolution vs multiplicity)
void plot(const vector<double>* const vzTruep, const vector<double>* const vzgrecp, 
          const vector<int>* const vzmultip, const vector<double>* const vzrecrmsp, 
          const double* const multiBin, const int arrayLengthMulti, 
          TH1D* tot1multi, TH1D* tot3multi, const int limit, const double sigmaZ, 
          const double* const zbin, const int arrayLengthZ, TH1D* totz, const int sizeTrue);

// Method for smearing points, generating noise, and reconstructing the primary vertex from tracklets
void reconstruction() {
    double step = 0.05; // Step for histograms intersections in z
    double R1, R2, L, sigmaZ;
    int events;
    char a;

    // Read parameters from data.txt
    ifstream In(kData);
    if (!In.is_open()) {
        cout << "Error: The file data doesn't exist" << endl; 
        return;
    }
    In >> R1 >> R2 >> L >> sigmaZ >> events >> a;
    In.close();
    TString mul = (a == 'N') ? "No" : "Yes";

    // Open the input simulation file
    TFile* file = TFile::Open(kSim, "READ");
    if (file == NULL) {
        cout << "Error: The file doesn't exist" << endl; 
        return;
    }
    file->ls(); // List contents of the file

    // Get the input tree from the file
    TTree* inputTree = (TTree*)file->Get("tree");
    if (inputTree == NULL) {
        cout << "Error: The Tree doesn't exist" << endl; 
        return;
    }

    // Create output file and tree
    TFile* fileReco = TFile::Open(kRic, "RECREATE"); 
    TTree* treeReco = new TTree("treeReco", "Reconstruction");
    treeReco->SetDirectory(fileReco);

    // Check parameter validity
    if (step <= 0 || kVerbosityReco <= 0 || kNoisefrac <= 0 || kPhiMax <= 0 || kRange <= 0) {
        cout << "Error in parameters setting" << endl; 
        return;
    }

    // Calculate the maximum z obtainable from tracklets
    double zExtr = R2 * (L / (R2 - R1)) - (L / 2);
    double zMin = -zExtr - step / 2;
    double zMaxBin = zExtr + step / 2;	
    int nbin = static_cast<int>((zMaxBin - zMin) / step); 
    step = (zMaxBin - zMin) / nbin; // Adjust step to cover the range

    // Create histogram for intersection points
    TH1D* zIntersecHisto = new TH1D("zRec", "Intersections in z; z [cm]; Number of events", nbin, zMin, zMaxBin);
    TClonesArray* intersecPointDet1 = new TClonesArray("Hit", 200); // T1
    TClonesArray* intersecPointDet2 = new TClonesArray("Hit", 200); // T2
    TClonesArray& int1 = *intersecPointDet1;  
    TClonesArray& int2 = *intersecPointDet2;

    Vertex* vertexRec = nullptr; // Pointer to a Vertex object
    vector<double> Zintersection; // Vector for intersections of tracklets

    // Assign branches for the output tree
    treeReco->Branch("Vertexreal", &vertexRec); // Primary vertex
    treeReco->Branch("zintersectionvector", &Zintersection); // Vector of Tracklet intersections 
    treeReco->Branch("int1enoise", &intersecPointDet1); // Intersection + noise T1
    treeReco->Branch("int2enoise", &intersecPointDet2); // Intersection + noise T2

    // Read from input file TCloneArrays
    TClonesArray* clone1 = new TClonesArray("Hit", 100);
    TClonesArray* clone2 = new TClonesArray("Hit", 100);
    for (int i = 0; i < 10; i++) { 
        new ((*clone1)[i]) Hit(); 
        new ((*clone2)[i]) Hit();
    } // Initialize non-zero dimension objects

    Vertex* verlec = new Vertex(); // Create a Vertex object

    // Set branch addresses
    inputTree->GetBranch("tclone1")->SetAutoDelete(kFALSE);
    inputTree->SetBranchAddress("tclone1", &clone1);
    inputTree->GetBranch("tclone2")->SetAutoDelete(kFALSE);
    inputTree->SetBranchAddress("tclone2", &clone2);

    TBranch* b3 = inputTree->GetBranch("Vertext");
    b3->SetAddress(&verlec); // Reading Vertex

    clone1->Clear(); // Clear clone1 and clone2
    clone2->Clear();

    // Histogram for efficiency (multiplicity)
    double multiBin[] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9.5, 
                         11.5, 14.5, 17.5, 20.5, 25.5, 30.5, 40.5, 
                         50.5, 60.5, 70.5, 80.5, 90.5, 91.5};
    int arrayLengthMulti = sizeof(multiBin) / sizeof(multiBin[0]) - 1; // Number of histogram bins
    int limit; // Maximum multiplicity for histogram with multiBin[limit]
    if (kLimit1 > kLimit2 || kLimit2 > arrayLengthMulti || kLimit1 <= 0 || kLimit2 <= 0) {
        cout << "Error on kLimit1 and kLimit2" << endl;
        return;
    }
    limit = (mul == "No") ? kLimit2 : kLimit1; // Determine limit based on multiplicity choice

    // Create histograms for efficiency
    TH1D* tot1multi = new TH1D("tot1multi", "tot1multi", arrayLengthMulti, multiBin); // Histogram for events with z sim between 1 sigma from 0
    TH1D* tot3multi = new TH1D("tot3multi", "tot3multi", arrayLengthMulti, multiBin); // Between 3 sigma

    // Histogram for efficiency based on zTrue (z of primary vertex)	
    double zbin[] = {-16, -12, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 12, 16};
    int arrayLengthZ = sizeof(zbin) / sizeof(zbin[0]) - 1; // Number of bins
    TH1D* totz = new TH1D("totz", "totz", arrayLengthZ, zbin);

    // Prepare final histogram ingredients
    vector<double>* vzTrue = new vector<double>; // z of the simulated vertex
    vector<double>* vzrec = new vector<double>; // z reconstructed
    vector<int>* vzmulti = new vector<int>; // Multiplicity from simulation
    vector<double>* vzrecrms = new vector<double>; // Error (RMS) on z reconstructed	
    vzTrue->reserve(events); // Preallocate memory to avoid fragmentation
    vzrec->reserve(events);
    vzmulti->reserve(events);
    vzrecrms->reserve(events);
}