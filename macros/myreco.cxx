//root includes
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TEfficiency.h" 
#include "TStopwatch.h"


#include <iostream>
#include <vector>
#include <fstream>

//headers include
#include "Vertex.h"
#include "Track.h"
#include "Hit.h"

using namespace std;

const int kVerbosityReco=10000000; //verbosity

const int kLimit1=17; //last bin to read for histograms for user multiplicity
const int kLimit2=21; //last bin to read for histograms for uniform multiplicity

const double kNoiseFrac=1.2; //ratio for noise/multiplicity points
const double kPhiMax=0.01; //max angle between points for creating tracklets
const double kRange=0.1; //range, attorno al picco dell'istogramma delle intersezioni delle tracklet, entro cui mediare le intersezioni per ricostruire z del vertice

const TString kSim="simulation.root"; //input file
const TString kRic="reconstruction.root"; //output file for reconstruction tree
const TString kHisto="histo.root"; //output file for histograms
const string kConfigFile="data.txt";  //config file for the run data

///////////////////////////////////////////////////////////////////
//method for plotting out the results (residues, efficiency vs multiplicity, 
//resolution vs multiplicity, efficiency vs zTrue, resolution vs zTrue)
////////////////////////////////////////////////////////////////////

void analysis(){} //DEFINE ME

void reconstruction(){
    TStopwatch clock;

    double step = 0.05; //step used by the histograms
    double R1, R2, L, sigmaZ;
    int events;
    char multiChoice;

    //opening config files
    ifstream In = ifstream(kConfigFile);
    if(!In.is_open()){
        cout<<"Error: the file config.txt does not exist"<<endl;
        return;
    }
    //assigning values
    In>>R1>>R2>>L>>sigmaZ>>events>>multiChoice;
    In.close();

    TString mul;
    if (multiChoice=='N') // checking multiplicity generation choice
        mul="No"; else mul="yes";

    //opening simulation result file
    TFile *file = TFile::Open(kSim, "READ");
    if(!file || file->IsZombie()){
        cout<<"Error: the file doesn't exist"<<endl; 
        return;
    }
    file->ls(); //reading the content of file
    TTree *simuTree = (TTree*)file->Get("tree");
    if(simuTree==NULL){ 
        cout<<"Error: the Tree doesn't exist"<<endl; 
        return;
    }
    //creating output file and Tree
    TFile *fileReco = TFile::Open(kRic,"RECREATE");
    TTree *treeReco = new TTree("treeReco", "reconstruction");
    treeReco->SetDirectory(fileReco);
    
    double zMax=R2*(L/(R2-R1))-(L/2); //maximum z obtainable from tracklets, from the edges of the detectors
    double zMin= -zMax; //symmetrical w.r.t. the origin
    double zMinBin= -zMax -step/2; //centre of the bin of the minimum
    double zMaxBin= zMax+step/2;
    int nBin=(int)(zMaxBin-zMinBin)/step; 
    step=(zMaxBin-zMinBin)/nBin; //correct step, in order to have a full covering of the range with nbins

    //creation of histogram and TCloneArray for intersection points (preallocation of memory)
    TH1D *zIntersectionHisto = new TH1D("zRec", "Intersections in z; z [cm]; number of events",nbin,zMin,zMaxBin);
    TClonesArray *intersectionPointDet1 = new TClonesArray("Hit", 200); //T1
    TClonesArray *intersectionPointDet2 = new TClonesArray("Hit", 200); //T2
    TClonesArray &int1 = *intersectionPointDet1;  
    TClonesArray &int2 = *intersectionPointDet2;

    //working variables
    Vertex *vertexReco = NULL; //creation of pointer to a Vertex obj
    vector<double> zIntersection; //vector for intersections of tracklets for a single vertex
    
    //assignment of branches for the output tree
    treeReco->Branch("vertexReco",&vertexReco); //primary vertex
    treeReco->Branch("zIntersectionVector",&zIntersection); //vector of Tracklet intersections
    treeReco->Branch("int1+noise",&intersecPointDet1); //intersection + noise T1
    treeReco->Branch("int2+noise",&intersecPointDet2); //intersection + noise T2

    //reading from the input file
    TClonesArray *clone1= new TClonesArray("Hit",100);
    TClonesArray *clone2= new TClonesArray("Hit",100);
    for (int i=0;i<10;i++){ 
     new ((*clone1)[i]) Hit(); 
     new ((*clone2)[i]) Hit();
    } //filling of the Hit arrays
    
    Vertex *vertexLec = new Vertex(); //vertex obj used for reading the ttree
    inputTree->GetBranch("tclone1")->SetAutoDelete(kFALSE); //prevents TClonesArray from reusing space allocated by the previous obj, setted for redundancy since it is kFALSE by default
    inputTree->SetBranchAddress("tclone1",&clone1);
    inputTree->GetBranch("tclone2")->SetAutoDelete(kFALSE);
    inputTree->SetBranchAddress("tclone2",&clone2);
    TBranch *b3=inputTree->GetBranch("Vertex");
    b3->SetAddress(&verlec); //reading Vertex
    clone1->Clear(); //cleaning clone1 e clone2
    clone2->Clear();
    

}
