// -*- C++ -*-
//
// Package:    DiLeptonFilter
// Class:      DiLeptonFilter
// 
/**\class DiLeptonFilter DiLeptonFilter.cc SuSyAachen/Skimming/src/DiLeptonFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Matthias Edelhoff
//         Created:  Mon Nov 16 11:26:19 CET 2009
// $Id: DiLeptonFilter.cc,v 1.4 2010/02/19 15:38:22 edelhoff Exp $
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
#include "DataFormats/Math/interface/deltaR.h"

//include "CommonTools/UtilAlgos/interface/DeltaR.h"


//
// class declaration
//

class DiLeptonFilter : public edm::EDFilter {
public:
  explicit DiLeptonFilter(const edm::ParameterSet&);
  ~DiLeptonFilter();
  
private:
  typedef reco::Candidate cand;
  typedef edm::View<cand> collection;
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  bool filterChannel(const collection &primary, const collection &secondary);  
  bool inclusiveFilter( const collection &primary, const collection &secondary);
  bool different(const cand &a, const cand &b);

  // ----------member data ---------------------------
  std::vector<std::string> combinations_;
  
  edm::InputTag primaryTag_;
  edm::InputTag secondaryTag_;
  edm::InputTag tertiaryTag_;

  bool sameSign_;
  bool matching_;
  std::string method_;

  double minDR_;
  double minDpt_;

  bool useSecondary_;
  bool useTertiary_;
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
DiLeptonFilter::DiLeptonFilter(const edm::ParameterSet& iConfig)
{
  combinations_ = iConfig.getParameter< std::vector<std::string> > ("combinations");
  primaryTag_ = iConfig.getParameter<edm::InputTag> ("primarySrc");
  secondaryTag_ = iConfig.getParameter<edm::InputTag> ("secondarySrc");
  tertiaryTag_ = iConfig.getParameter<edm::InputTag> ("tertiarySrc");
  sameSign_ = iConfig.getParameter<bool> ("sameSign");
  matching_ = iConfig.getParameter<bool> ("matching");
  method_ = iConfig.getParameter<std::string> ("method");
  minDR_ = iConfig.getParameter<double> ("minDR");
  minDpt_ = iConfig.getParameter<double> ("minDpt");
  useSecondary_ = false;
  useTertiary_ = false;


  for(std::vector< std::string >::iterator it = combinations_.begin(); it != combinations_.end(); ++it){
    if( (*it).size() != 2
	||( (*it)[0] != 'p' && (*it)[0] != 's' && (*it)[0] != 't') 
	|| ((*it)[1] != 'p' && (*it)[1] != 's' && (*it)[1] != 't') )
      throw new cms::Exception("malformed combination: "+(*it));    
    useSecondary_ |= (*it)[0] == 's' || (*it)[1] == 's';
    useTertiary_ |= (*it)[0] == 't' || (*it)[1] == 't';
  }
}


DiLeptonFilter::~DiLeptonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DiLeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool result = false;

  std::map< char, collection > collectionMap;
  
  edm::Handle< collection > primary;
  iEvent.getByLabel(primaryTag_, primary);
  collectionMap['p'] = *primary;

  if(useSecondary_){
    edm::Handle< collection > secondary;
    iEvent.getByLabel(secondaryTag_, secondary);
    collectionMap['s'] = *secondary;
  }
  

  if(useTertiary_){
    edm::Handle< collection > tertiary;
    iEvent.getByLabel(tertiaryTag_, tertiary);
    collectionMap['t'] = *tertiary;
   }

  for(std::vector< std::string >::iterator it = combinations_.begin(); it != combinations_.end(); ++it){
    bool fRes = filterChannel( collectionMap[(*it)[0]], collectionMap[(*it)[1]] );
    //    std::cout<<" -> " << fRes << " " << *it<<std::endl;
    result |= fRes;
  }

  return result;
}

bool
DiLeptonFilter::filterChannel(const collection &primary, const collection &secondary)
{
  bool result = false;
  if(method_ == "inclusive")
    result = inclusiveFilter( primary, secondary);
  else if (method_ == "exclusive"){
    int overlap = 0;
    for(collection::const_iterator iPrime = primary.begin(); iPrime != primary.end() ; ++iPrime)
      for(collection::const_iterator iSec = secondary.begin(); iSec != secondary.end() ; ++iSec)
	if(!different( (*iPrime), (*iSec)))
	  overlap++;
    result = inclusiveFilter( primary, secondary);
    result = result && (( primary.size() + secondary.size() - overlap) == 2);
  }else
    throw new cms::Exception("unknown method: "+method_);
  return result;
}

bool DiLeptonFilter::inclusiveFilter( const collection &primary, const collection &secondary)
{
   bool result = false;
  for(collection::const_iterator iPrime = primary.begin(); iPrime != primary.end() ; ++iPrime){
    for(collection::const_iterator iSec = secondary.begin(); iSec != secondary.end() ; ++iSec){
      if( ( different( (*iPrime), (*iSec) ) ) 
	  &&( fabs((*iPrime).charge()) == 1 && fabs((*iSec).charge()) == 1 ) ){
	if( sameSign_ ){
	  result = result || ( (*iPrime).charge() == (*iSec).charge() );
	}
	else
	  result = result || ( (*iPrime).charge() != (*iSec).charge() );
	//	std::cout << (*iPrime).charge() << (*iSec).charge() <<"| ";
      }else{
	//		std::cout << "O| ";
      }
    }
    // std::cout << "   ";
  }
  //  std::cout <<" ("<< primary.size()<<", "<< secondary.size()<<") "<<" result: "<< result;
  return result;
}

bool DiLeptonFilter::different(const cand &a, const cand &b)
{
  if(!matching_)
    return &a != &b;
  else
    return reco::deltaR<const cand, const cand>( a, b) > minDR_ 
      && fabs( a.pt() - b.pt() ) > minDpt_ ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DiLeptonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonFilter);
