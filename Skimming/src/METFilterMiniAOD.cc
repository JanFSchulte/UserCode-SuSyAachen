// -*- C++ -*-
//
// Package:    METFilterMiniAOD
// Class:      METFilterMiniAOD
// 
/**\class METFilterMiniAOD METFilterMiniAOD.cc SuSyAachen/Skimming/src/METFilterMiniAOD.cc

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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
//
// class declaration
//

class METFilterMiniAOD : public edm::EDFilter {
public:
  explicit METFilterMiniAOD(const edm::ParameterSet&);
  ~METFilterMiniAOD();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::InputTag filterTag_;



  bool debug;
};

// constructors and destructor
METFilterMiniAOD::METFilterMiniAOD(const edm::ParameterSet& iConfig)
{
  filterTag_    = iConfig.getParameter < edm::InputTag > ("src");

  debug = false;
}

METFilterMiniAOD::~METFilterMiniAOD(){}


// member functions
// ------------ method called on each new Event  ------------
bool
METFilterMiniAOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<edm::TriggerResults> filterBits;
  iEvent.getByLabel(filterTag_, filterBits);
	bool filter = true;
    for (unsigned int i = 0, n = filterBits->size(); i < n; ++i) {
    
    	if (!filterBits->accept(i)){
    		filter = false;
    	}

    }
 

  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
METFilterMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METFilterMiniAOD::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METFilterMiniAOD);
