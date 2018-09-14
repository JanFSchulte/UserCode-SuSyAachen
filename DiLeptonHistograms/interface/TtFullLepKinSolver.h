#ifndef TtFullLepKinSolver_h
#define TtFullLepKinSolver_h

#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"

class TF2;

/*
  \class   TtFullLepKinSolver TtFullLepKinSolver.h "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
  
  \brief   Class to calculate solutions for neutrino momenta in dileptonic ttbar events and related probability weights

  Class to calculate solutions for neutrino momenta in dileptonic ttbar events and related probability weights.
  A fourth-order polynomial in p_x(nu) is used with coefficients that are functions of the top-quark mass.
  This class is based on a code by Jan Valenta.
  
**/

class TtFullLepKinSolver {

 public:

  ///
  struct NeutrinoSolution {
    double weight;
    int smearingsWithSols;
    reco::LeafCandidate neutrino;
    reco::LeafCandidate neutrinoBar; 
    TLorentzVector LV_t;
    TLorentzVector LV_t_;
  };

  /// default constructor
  TtFullLepKinSolver();
  /// constructor with parameters to configure the top-mass scan and the neutrino spectrum
  TtFullLepKinSolver(const std::string="", const bool=true, const double=80.4, const double=4.8);
  /// destructor
  ~TtFullLepKinSolver();
  ///
  void SetConstraints(const double xx=0, const double yy=0);
  ///
  NeutrinoSolution getNuSolution(const TLorentzVector& LV_l, 
                                 const TLorentzVector& LV_l_, 
               const TLorentzVector& LV_b, 
               const TLorentzVector& LV_b_,
               const TLorentzVector& met);        
  
           
 private:

  ///
  void FindCoeff(const TLorentzVector& al, 
     const TLorentzVector& l,
     const TLorentzVector& b_al,
     const TLorentzVector& b_l,
     const double mt, const double mat, const double pxboost, const double pyboost,
     double* q_coeff);
  ///
  void TopRec(const TLorentzVector& al, 
        const TLorentzVector& l,
        const TLorentzVector& b_al,
        const TLorentzVector& b_l, const double sol);
  ///
  double WeightSolfromMass() const;
  ///
  double WeightSolfromMLB(double mlb1, double mlb2);
  ///
  int quartic(double* q_coeff, double* q_sol) const;
  ///
  int cubic(const double* c_coeff, double* c_sol) const;
  ///
  double sqr(const double x) const {return (x*x);}
  ///
  void SWAP(double& realone, double& realtwo) const;
    
 private:
  ///
  double mw;
  ///
  const double mb;
  ///
  const std::string smearingFilename;
  ///
  double pxmiss_, pymiss_;
  ///
  bool doSmearing;
  
  double C;
  double D;
  double F;
  double pom;
  double k16;
  double k26;
  double k36;
  double k46;
  double k56;
  double k51;
  double k61;
  double m1;
  double m2;
  double m3;
  double n1;
  double n2;
  double n3;
  
  ///
  TLorentzVector LV_n, LV_n_, LV_t, LV_t_, LV_tt_t, LV_tt_t_;  
  TLorentzVector LV_L, LV_L_, LV_B, LV_B_;
  /// provisional
  TLorentzVector genLV_n, genLV_n_;  
    
  TH1F leptonEnergyHist, leptonEtaHist, leptonPhiHist;
  TH1F jetEnergyHist, jetEtaHist, jetPhiHist;
  TH1F wMassHist, lbMassHist;
};


#endif
