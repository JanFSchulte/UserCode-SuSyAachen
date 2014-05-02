#ifndef TauToTauMatchSelector_h
#define TauToTauMatchSelector_h

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/Math/interface/deltaR.h"

//STL
#include <vector>

template<typename collectionType, typename containerType>
  struct TauToTauMatchSelector {
    typedef collectionType collection;
    typedef containerType container;
    typedef typename container::const_iterator const_iterator;
    TauToTauMatchSelector ( const edm::ParameterSet & cfg, edm::ConsumesCollector ):
      dRMin_( cfg.getParameter<double>( "dRMin") ),
      otherTauId_(cfg.getParameter<std::string>( "otherTauId")),
      otherTauSrc_( cfg.getParameter<edm::InputTag>( "otherTauSource" ) ) { }
  
    const_iterator begin() const { return selected_.begin(); }
    const_iterator end() const { return selected_.end(); }
    void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {

      edm::Handle<std::vector<pat::Tau> > otherTaus;
      ev.getByLabel(otherTauSrc_, otherTaus);

      selected_.clear();
      for(typename collection::const_iterator it = col.product()->begin(); 
	  it != col.product()->end(); ++it ){
	double  minDeltaR = dRMin_;
	for(pat::TauCollection::const_iterator other = otherTaus->begin(); other != otherTaus->end(); ++other){
	  if((*other).tauID(otherTauId_) > 0.5 && reco::deltaR<const pat::Tau, const pat::Tau>( *other, *it) < minDeltaR){
	    //std::cout << "d0 = " << dS.first << "  dZ = " << dS.second << std::endl;
	    minDeltaR = reco::deltaR<const pat::Tau, const pat::Tau>( *other, *it);
	    //std::cout << "Selected it d0 = " << dS.first << "  dZ = " << dS.second << std::endl;
	  }
	}
	if(minDeltaR < dRMin_)
	  selected_.push_back( & (*it) );
      }
    }

    size_t size() const { return selected_.size(); }
    private:
    container selected_;
    double dRMin_;
    std::string otherTauId_;
    edm::InputTag otherTauSrc_;
  };

#endif
