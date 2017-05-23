#ifndef isolationSelector_h
#define isolationSelector_h

//DataFormats
//include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SuSyAachen/DiLeptonHistograms/interface/IsolationFunctor.h"

//STL
#include <vector>

template<typename rhoType, typename collectionType, typename containerType>
struct isolationSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  isolationSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector iC ):
    fctIsolation_ (cfg.getParameter<edm::ParameterSet>("isolationDefinitions"), (edm::ConsumesCollector&&)iC),
    isoMin_( cfg.getParameter<double>( "isoMin") ),
    isoMax_( cfg.getParameter<double>( "isoMax" ) ),
	method_( cfg.getParameter<string>( "method" ) ),    
    isoMaxEE_(0.1) { }

  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    fctIsolation_.init(ev);


    selected_.clear();
    double iso = 0;
    for(typename collection::const_iterator it = col.product()->begin(); 
	it != col.product()->end(); ++it ){
      iso = 0;
      iso = fctIsolation_(*(it),method_) / (*it).pt();
		
      //std::cout << iso << std::endl;
      if (iso > isoMin_ && iso < isoMax_){
	

	  	selected_.push_back( & (*it) );
      }
   }
  }

  size_t size() const { return selected_.size(); }
  
  private:
    
    container selected_;
    IsolationFunctor fctIsolation_;
    double isoMin_;
    double isoMax_;
    string method_;    
    double isoMaxEE_;

};

#endif
