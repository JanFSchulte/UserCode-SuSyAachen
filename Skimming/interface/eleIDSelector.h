#ifndef eleIDSelector_h
#define eleIDSelector_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"


//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

template<typename T, typename collectionType, typename containerType>
struct eleIDSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  eleIDSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
    idMapSrc_( cfg.getParameter<edm::InputTag>( "idMapSource" ) )  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle<edm::ValueMap<bool> > id_decisions;
    ev.getByLabel(idMapSrc_, id_decisions);


    selected_.clear();

    //~ double etaValue;
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
		const edm::Ptr<pat::Electron> elPtr(col, it - col->begin() );
		bool isPass  = (*id_decisions)[ elPtr ];
		if (isPass){
			selected_.push_back( & (*it) );
    	}
    }
  }

    

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  edm::InputTag idMapSrc_;
};

#endif
