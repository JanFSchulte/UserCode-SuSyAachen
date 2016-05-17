// -*- C++ -*-
//
// Package:    MetSqrtHtFilter
// Class:      MetSqrtHtFilter
// 
/**\class MetSqrtHtFilter MetSqrtHtFilter.cc SuSyAachen/Skimming/src/MetSqrtHtFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Matthias Edelhoff
//         Created:  Mon Nov 16 11:26:19 CET 2009
// $Id: MetSqrtHtFilter.cc,v 1.1 2011/01/17 11:01:11 nmohr Exp $
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

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//

class MetSqrtHtFilter : public edm::EDFilter {
public:
  explicit MetSqrtHtFilter(const edm::ParameterSet&);
  ~MetSqrtHtFilter();
  
private:
  typedef reco::Candidate cand;
  typedef edm::View<cand> collection;
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::EDGetTokenT< collection > inputToken_;
  edm::EDGetTokenT< std::vector<pat::MET> > inputTokenMet_;

  double minCut_;
  double maxCut_;

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
MetSqrtHtFilter::MetSqrtHtFilter(const edm::ParameterSet& iConfig):
  inputToken_(consumes< collection >(iConfig.getParameter<edm::InputTag>("src"))),
  inputTokenMet_(consumes< std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("metSrc")))
{
  minCut_ = iConfig.getParameter<double> ("minCut");
  maxCut_ = iConfig.getParameter<double> ("maxCut");
}


MetSqrtHtFilter::~MetSqrtHtFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MetSqrtHtFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle< collection > candidates;
  iEvent.getByToken(inputToken_, candidates);
  
  edm::Handle< std::vector<pat::MET> > mets;
  iEvent.getByToken(inputTokenMet_, mets);

  bool result = false;
  double ht = 0.0;
  double met = mets->front().pt();
  for(collection::const_iterator it = candidates->begin(); it != candidates->end() ; ++it){
    ht += (*it).pt();
  }
  double val = 0.;
  if (ht>0.) val = met/sqrt(ht);
  if (maxCut_ > 0. ) result = (val >=minCut_ && val <= maxCut_);
  else result = val >=minCut_;

  return result;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MetSqrtHtFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MetSqrtHtFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MetSqrtHtFilter);
