// -*- C++ -*-
//
// Package:    InverseMETFilter
// Class:      InverseMETFilter
// 
/**\class InverseMETFilter InverseMETFilter.cc SuSyAachen/Skimming/src/InverseMETFilter.cc

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

class InverseMETFilter : public edm::EDFilter {
public:
  explicit InverseMETFilter(const edm::ParameterSet&);
  ~InverseMETFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::EDGetTokenT<bool> filterToken_;



  bool debug;
};

// constructors and destructor
InverseMETFilter::InverseMETFilter(const edm::ParameterSet& iConfig):
  filterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("src")))
{
  debug = false;
}

InverseMETFilter::~InverseMETFilter(){}


// member functions
// ------------ method called on each new Event  ------------
bool
InverseMETFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< bool > filterHandle;
  iEvent.getByToken(filterToken_, filterHandle);
  Bool_t filter = (Bool_t)(*filterHandle);

 

  if (filter == true){ 
     return false;
  }
  else{
     return true;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
InverseMETFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
InverseMETFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(InverseMETFilter);
