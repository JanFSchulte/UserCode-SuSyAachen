// -*- C++ -*-
//
// Package:    Histograms
// Class:      pdfUncertaintyNominatorTrees
// 
/**\class pdfUncertaintyNominatorTrees pdfUncertaintyNominatorTrees.cc brot/pdfUncertaintyNominatorTrees/src/pdfUncertaintyNominatorTrees.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  matthias edelhoff
//         Created:  Tue Oct 27 13:50:40 CET 2009
// $Id: pdfUncertaintyNominatorTrees.cc,v 1.31 2012/09/17 17:38:58 sprenger Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/Lepton.h>

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>

#include <DataFormats/Provenance/interface/EventID.h>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctorAbsEta.h>
#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctorFastSim.h>
#include <SuSyAachen/DiLeptonHistograms/interface/WeightFunctorFastSimDifferentLeptons.h>
#include <SuSyAachen/DiLeptonHistograms/interface/PdgIdFunctor.h>
#include <SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h>
#include <SuSyAachen/TagAndProbeTreeWriter/interface/IsolationFunctor.h>

//ROOT
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

//
// class decleration
//

class pdfUncertaintyNominatorTrees : public edm::EDAnalyzer {
public:
  explicit pdfUncertaintyNominatorTrees(const edm::ParameterSet&);
  ~pdfUncertaintyNominatorTrees();

private:
  //  typedef reco::Candidate candidate;
  typedef pat::Lepton<reco::Candidate> candidate;
  typedef edm::View<candidate> collection;

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void initFloatBranch( const std::string &name);
  void initIntBranch( const std::string &name);
 
  template<class aT, class bT> void fillTree( const std::string &treeName, const aT &a, const bT &b, const std::vector<reco::PFCandidate>&pfCands, const pat::MET &patMet, const TLorentzVector &MHT, const int &nGenLeptons, const TLorentzVector &genLept1, const TLorentzVector &genLept2, const TLorentzVector &genLept3, const TLorentzVector &genLept4);

 

  edm::InputTag pdfInfo_;
  std::vector<edm::InputTag> pdfs_;


  //data
  std::map<std::string, TTree*> trees_;  
  std::map<std::string, std::map< std::string, float*> > floatBranches_; 
  std::map<std::string, std::map< std::string, unsigned int*> > intBranches_; 




  bool debug;
};

// constructors and destructor
pdfUncertaintyNominatorTrees::pdfUncertaintyNominatorTrees(const edm::ParameterSet& iConfig){
  debug = false;

  
  // read config

  pdfInfo_ = iConfig.getParameter<edm::InputTag>("pdfInfo");


  // init trees
  edm::Service<TFileService> file;
  trees_["Tree"] = file->make<TTree>("Tree", "Tree");
 
  initFloatBranch( "pdf_x1" );
  initFloatBranch( "pdf_x2" );
  initFloatBranch( "pdf_scale" );
  initFloatBranch( "pdf_id1" );
  initFloatBranch( "pdf_id2" );
  initFloatBranch( "pdf_weight" );

}



void 
pdfUncertaintyNominatorTrees::initFloatBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    floatBranches_[(*it).first][name] = new float;
    (*it).second->Branch(name.c_str(), floatBranches_[(*it).first][name], (name+"/F").c_str());
  }
}

void 
pdfUncertaintyNominatorTrees::initIntBranch(const std::string &name)
{
  for( std::map<std::string, TTree*>::const_iterator it = trees_.begin();
       it != trees_.end(); ++it){
    if(debug) std::cout << (*it).first <<" - "<< name << std::endl;
    intBranches_[(*it).first][name] = new unsigned int;
    (*it).second->Branch(name.c_str(), intBranches_[(*it).first][name], (name+"/I").c_str());
  }
}

pdfUncertaintyNominatorTrees::~pdfUncertaintyNominatorTrees()
{ 
  for( std::map<std::string, std::map< std::string, float*> >::const_iterator it = floatBranches_.begin();
       it != floatBranches_.end(); ++it){
    for( std::map< std::string, float*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug)std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }
  for( std::map<std::string, std::map< std::string, unsigned int*> >::const_iterator it = intBranches_.begin();
       it != intBranches_.end(); ++it){
    for( std::map< std::string, unsigned int*>::const_iterator it2 = (*it).second.begin();
	 it2 != (*it).second.end(); ++it2){
      if(debug) std::cout << "deleting: " << (*it).first << " - "<< (*it2).first << std::endl;
      delete (*it2).second;
    }
  }


}

// member functions
// ------------ method called to for each event  ------------
void
pdfUncertaintyNominatorTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  std::map<std::string, int> intEventProperties;
  std::map<std::string, float> floatEventProperties;
  std::map<std::string, TLorentzVector> tLorentzVectorEventProperties;

  floatEventProperties["pdf_x1"] = -10.;
  floatEventProperties["pdf_x2"] = -10.;
  floatEventProperties["pdf_scale"] = -1.;
  floatEventProperties["pdf_id1"] = -1.;
  floatEventProperties["pdf_id2"] = -1.;
  floatEventProperties["pdf_weight"] = -1.;

  edm::Handle<GenEventInfoProduct> pdfstuff;
  if (iEvent.getByLabel(pdfInfo_, pdfstuff)) {
  	floatEventProperties["pdf_x1"] = pdfstuff->pdf()->x.first;
  	floatEventProperties["pdf_x2"] = pdfstuff->pdf()->x.second;
  	floatEventProperties["pdf_scale"] = pdfstuff->pdf()->scalePDF;
  	floatEventProperties["pdf_id1"] = pdfstuff->pdf()->id.first;
  	floatEventProperties["pdf_id2"] = pdfstuff->pdf()->id.second;
  	floatEventProperties["pdf_weight"] = -1.;


      }
  for(std::map<std::string, int>::const_iterator it = intEventProperties.begin(); it != intEventProperties.end(); ++it){
    assert(intBranches_["Tree"].find((*it).first) != intBranches_["Tree"].end());
    *(intBranches_["Tree"][(*it).first]) = (*it).second;
  }
  for(std::map<std::string, float>::const_iterator it = floatEventProperties.begin(); it != floatEventProperties.end(); ++it){
    assert(floatBranches_["Tree"].find((*it).first) != floatBranches_["Tree"].end());
    *(floatBranches_["Tree"][(*it).first]) = (*it).second;
  }
  trees_["Tree"]->Fill();
}





// ------------ Method called once each job just before starting event loop  ------------
void 
pdfUncertaintyNominatorTrees::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
pdfUncertaintyNominatorTrees::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(pdfUncertaintyNominatorTrees);
