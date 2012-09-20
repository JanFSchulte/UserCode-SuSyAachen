// -*- C++ -*-
//
// Package:    METFilter
// Class:      METFilter
// 
/**\class METFilter METFilter.cc SuSyAachen/Skimming/src/METFilter.cc

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

class METFilter : public edm::EDFilter {
public:
  explicit METFilter(const edm::ParameterSet&);
  ~METFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::InputTag filterTag_;



  bool debug;
};

// constructors and destructor
METFilter::METFilter(const edm::ParameterSet& iConfig)
{
  filterTag_    = iConfig.getParameter < edm::InputTag > ("src");

  debug = false;
}

METFilter::~METFilter(){}


// member functions
// ------------ method called on each new Event  ------------
bool
METFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< bool > filterHandle;
  iEvent.getByLabel(filterTag_, filterHandle);
  Bool_t filter = (Bool_t)(*filterHandle);

 

 
  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
METFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METFilter);
