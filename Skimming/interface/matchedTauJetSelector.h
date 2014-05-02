#ifndef SuSyAachen_Skimming_MatchedTauJetSelector_h
#define SuSyAachen_Skimming_MatchedTauJetSelector_h

#include "DataFormats/PatCandidates/interface/Tau.h"
struct MatchedTauJetSelector {
MatchedTauJetSelector() { }
MatchedTauJetSelector( const edm::ParameterSet & cfg, edm::ConsumesCollector ) { }

  bool operator()( const pat::Tau & p ) const 
  { 
    return p.genJet() != 0;
  }
};

#endif
