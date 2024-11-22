#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include "TH1D.h"
#include "Vertex.h"
#include "TAxis.h"

using namespace std;
using namespace TMath;

ClassImp(Vertex)

// Default constructor (vertex in origin, 0 products)
Vertex::Vertex() : TObject() {  
    fX = 0;     
    fY = 0;   
    fZ = 0;
    fMulti = 0;   
}

// Standard constructor with multiplicity options
Vertex::Vertex(TString varV, TString multV, TH1D *multHist, double sigmaX, double sigmaY, double sigmaZ, int u1, int u2) : TObject() {
    fMulti = 0;
    // Choose gaus var extraction based on user input
    if (varV == "No" || varV == "NO" || varV == "no"){
        fX = VarROOT(sigmaX);  // x, y, z extracted with a Gaussian centered in the origin
        fY = VarROOT(sigmaY);   
        fZ = VarROOT(sigmaZ);
    } else {
        fX = VarBM(sigmaX);  // x, y, z extracted with a Gaussian centered in the origin
        fY = VarBM(sigmaY);   
        fZ = VarBM(sigmaZ);
    }
    // Choose multiplicity extraction based on user input
    if (multV == "No" || multV == "NO" || multV == "no") {
        this->MultUnif(u1, u2); // Uniform multiplicity extraction
    } else if (multV == "fixed"){
        fMulti = 12;
    } else {
        this->MultFunc(multHist); // Use distribution from histogram
    }
}

// Copy constructor
Vertex::Vertex(const Vertex &v) :
    TObject(v),
    fX(v.fX),
    fY(v.fY),
    fZ(v.fZ),
    fMulti(v.fMulti)
{
}

// Destructor
Vertex::~Vertex() {}

// Box-Muller method for Gaussian extraction
double Vertex::VarBM(double s) { 
    if (s > 0) {
        double u1 = gRandom->Rndm(); 
        double u2 = gRandom->Rndm();   
        return Sqrt(-2 * Log(u1)) * Cos(2 * Pi() * u2) * s; // Gaussian with mean 0, sigma = s
    } else {  
        cout << "Error: sigma <= 0" << endl;
        return 0;
    }
}
// ROOT Gaus method, using W. Hoermann and G. Derflinger ACR
double Vertex::VarROOT(double s) { 
    if (s > 0) {
        return gRandom->Gaus(0,s);
    } else {  
        cout << "Error: sigma <= 0" << endl;
        return 0;
    }
}


// Getters for vertex properties
double Vertex::GetX() const { return fX; }
double Vertex::GetY() const { return fY; }
double Vertex::GetZ() const { return fZ; }
int Vertex::GetM() const { return fMulti; }
double Vertex::GetPhi() const { return fPhi; }
double Vertex::GetTheta() const { return fTheta; }

// Uniform multiplicity distribution
void Vertex::MultUnif(int u1, int u2) {
    if (u2 >= u1) {
        fMulti = u1 + (u2 - u1) * gRandom->Rndm();
    } else {
        cout << "Error in multiplicity parameters choice" << endl;
    }
}

// Given multiplicity distribution function using histogram
void Vertex::MultFunc(TH1D *multHist) {
    if (multHist == nullptr) {
        cout << "Error: multiplicity histogram is null" << endl;
        return;
    }	
    fMulti = static_cast<int>(multHist->GetRandom());
}

// Initializes the initial direction; min and max are for reducing the range in pseudorapidity
void Vertex::InitialDir(double min, double max, TH1D *etaHist) { 
    if (max < min) { // Swap limits if out of order
        swap(max,min);
    } else if (max == min) {
        cout << "Error: etaMax == etaMin" << endl; 
        return;
    }
    
    fPhi = gRandom->Rndm() * 2 * Pi(); // Phi extraction
    if (etaHist == nullptr) {
        cout << "Error: eta histogram is null" << endl;
        return;
    }
    
    double eta;
    do {
        eta = etaHist->GetRandom();
    } while (eta > max || eta < min); // Retry until within desired range
    
    fTheta = 2 * ATan(Exp(-eta));
}

void Vertex::InitialDir1(TH1D *etaHist) { 
    fPhi = gRandom->Rndm() * 2 * Pi(); // Phi extraction
    if (etaHist == nullptr) {
        cout << "Error: eta histogram is null" << endl;
        return;
    }
    
    double eta;
    eta = etaHist->GetRandom();
    fTheta = 2 * ATan(Exp(-eta));
}
