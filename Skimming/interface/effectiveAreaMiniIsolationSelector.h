#ifndef isolationSelector_h
#define isolationSelector_h

//DataFormats
//include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SuSyAachen/DiLeptonHistograms/interface/IsolationFunctor.h"

//STL
#include <vector>

template<typename rhoType, typename collectionType, typename containerType>
struct isolationSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  isolationSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
    fctIsolation_ (cfg.getParameter<edm::ParameterSet>("isolationDefinitions")),
    isoMin_( cfg.getParameter<double>( "isoMin") ),
    isoMax_( cfg.getParameter<double>( "isoMax" ) ),
	method_( cfg.getParameter<string>( "method" ) ),    
    isoMaxEE_(0.1),

  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    fctIsolation_.init(ev);


    selected_.clear();
    int i = 0;
    double iso = 0;
    for(typename collection::const_iterator it = col.product()->begin(); 
	it != col.product()->end(); ++it ){
      iso = 0;
      iso = fctIsolation_(e,method_)
		
      
      if( (iso > isoMin_ && iso < isoMax_){
	

	  	selected_.push_back( & (*it) );
      }
    }
  }

  size_t size() const { return selected_.size(); }
  
    
    container selected_;
    IsolationFunctor fctIsolation_;
    double isoMin_;
    string method_;
    double isoMax_;
    double isoMaxEE_;

};

#endif
