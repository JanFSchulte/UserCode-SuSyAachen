#ifndef tightMuonSelector_h
#define tightMuonSelector_h

//DataFormats
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

template<typename T, typename collectionType, typename containerType>
struct tightMuonSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  tightMuonSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
    vertexSrc_( cfg.getParameter<edm::InputTag>( "vertexSource" ) )  { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    edm::Handle<reco::VertexCollection> vertexHandle;
    ev.getByLabel(vertexSrc_, vertexHandle);
   
    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); 
	 it != col.product()->end(); ++it ){		 
      if ( (*it).isTightMuon(vertexHandle->at(0)))
	selected_.push_back( & (*it) );
    }
  }

  size_t size() const { return selected_.size(); }
private:
  container selected_;
  edm::InputTag vertexSrc_;
};

#endif
