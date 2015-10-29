// -*- C++ -*-
//
// Package:    HBHENoiseFilterMiniAOD
// Class:      HBHENoiseFilterMiniAOD
// 
/**\class HBHENoiseFilterMiniAOD HBHENoiseFilterMiniAOD.cc SuSyAachen/Skimming/src/HBHENoiseFilterMiniAOD.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/



// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaPhi.h"
//
// class declaration
//

class HBHENoiseFilterMiniAOD : public edm::EDFilter {
public:
  explicit HBHENoiseFilterMiniAOD(const edm::ParameterSet&);
  ~HBHENoiseFilterMiniAOD();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::InputTag filterTag_;



  bool debug;
};

// constructors and destructor
HBHENoiseFilterMiniAOD::HBHENoiseFilterMiniAOD(const edm::ParameterSet& iConfig)
{
  filterTag_    = iConfig.getParameter < edm::InputTag > ("src");
 
  debug = false;
}

HBHENoiseFilterMiniAOD::~HBHENoiseFilterMiniAOD(){}


// member functions
// ------------ method called on each new Event  ------------
bool
HBHENoiseFilterMiniAOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<bool> HBHEFilterResult;
  iEvent.getByLabel(filterTag_, HBHEFilterResult);
  bool filter = true;
  
  if (!*HBHEFilterResult){
	  filter = false;
  }

  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HBHENoiseFilterMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HBHENoiseFilterMiniAOD::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HBHENoiseFilterMiniAOD);
