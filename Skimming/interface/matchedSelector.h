#ifndef SuSyAachen_Skimming_MatchedSelector_h
#define SuSyAachen_Skimming_MatchedSelector_h

#include <vector>

template<typename particleType>
struct MatchedSelector {
MatchedSelector( int pdgId, unsigned int status=0, bool autoCharge=false) : pdgId_(pdgId), status_(status), autoCharge_(autoCharge) { }

MatchedSelector( const edm::ParameterSet & cfg ) :
    pdgId_( cfg.getParameter<int>( "pdgId" ) ),
    status_( cfg.getParameter<unsigned  int>( "status" ) ),
    autoCharge_( cfg.getParameter<bool>( "autoCharge" ) ){ }

  bool operator()( const particleType & p ) const 
  { 
    return p.genParticleById(pdgId_, status_, autoCharge_).isNonnull(); 
  }

private:
  int pdgId_;
  unsigned int status_;
  bool autoCharge_;

};

#endif


