// -*- C++ -*-
//
// Package:    MhtFilter
// Class:      MhtFilter
// 
/**\class MhtFilter MhtFilter.cc SuSyAachen/Skimming/src/MhtFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Niklas Mohr
//
//


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

//
// class declaration
//

class MhtFilter : public edm::EDFilter {
public:
  explicit MhtFilter(const edm::ParameterSet&);
  ~MhtFilter();
  
private:
  typedef reco::Candidate cand;
  typedef edm::View<cand> collection;
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  //~ edm::InputTag inputTag_;
  edm::EDGetTokenT< collection > inputToken_;

  double minMHT_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MhtFilter::MhtFilter(const edm::ParameterSet& iConfig):
  inputToken_(consumes< collection >(iConfig.getParameter<edm::InputTag>("src")))
{
  //~ inputTag_ = iConfig.getParameter<edm::InputTag> ("src");
  minMHT_ = iConfig.getParameter<double> ("minMHT");
}


MhtFilter::~MhtFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MhtFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle< collection > candidates;
  //~ iEvent.getByLabel(inputTag_, candidates);
  iEvent.getByToken(inputToken_, candidates);

  reco::Candidate::LorentzVector cand(0.,0.,0.,0.);
  for(collection::const_iterator it = candidates->begin(); it != candidates->end() ; ++it){
     cand += (*it).p4();
  }
  double mht = cand.pt();
  bool result = mht >=minMHT_;

  return result;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MhtFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MhtFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MhtFilter);
