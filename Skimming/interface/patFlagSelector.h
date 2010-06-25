#ifndef patFlagSelector_h
#define patFlagSelector_h

//DataFormats
#include "DataFormats/PatCandidates/interface/Flags.h"

//Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"

//STL
#include <vector>

template<typename collectionType, typename containerType>
struct patFlagSelector {
  typedef collectionType collection;
  typedef containerType container;
  typedef typename container::const_iterator const_iterator;
  patFlagSelector ( const edm::ParameterSet & cfg ):
   flag_( cfg.getParameter<std::string>( "cut" ) ) { }
  
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle< collection > &col , const edm::Event &ev , const edm::EventSetup &setup ) {
    
    selected_.clear();
    for(typename collection::const_iterator it = col.product()->begin(); it != col.product()->end(); ++it ){
        if(pat::Flags::test( *it, pat::Flags::get(flag_))) selected_.push_back( & (*it) );
    }
  }
  
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  std::string flag_;
};

#endif
