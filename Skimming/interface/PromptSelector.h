#ifndef SuSyAachen_Skimming_PromptSelector_h
#define SuSyAachen_Skimming_PromptSelector_h

#include <vector>

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
//include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>

struct PromptSelector {
  typedef reco::CompositeRefCandidateT<reco::GenParticleRefVector>::daughters daughters;

PromptSelector( std::vector<int> promptMothers, bool bsm) : motherIds_( promptMothers ), bsm_(bsm) { }
  PromptSelector( const edm::ParameterSet & cfg ) :
  motherIds_( cfg.getParameter< std::vector<int> >( "promptMotherIds" ) ),
    bsm_( cfg.getParameter< bool >( "bsm" )){ }

  bool operator()( const reco::GenParticle & p ) const 
  {   
    bool result = false;
    if(p.status() == 1 && p.numberOfMothers() == 1){
      //      std::cout << "p "<< p.pdgId()<<" | "<<p.mother()->pdgId()<<": ";
      for( std::vector<int>::const_iterator i = motherIds_.begin();
	   i != motherIds_.end(); ++i ){
	if( p.mother()->pdgId() == (*i)) result |= true;
	std::cout << result<<" ";
      }
      
      result |= bsm_ && p.mother()->pdgId() > 1000000;
      //std::cout << result<<" ";
      result |= p.mother()->pdgId() == p.pdgId();
      //std::cout << result<<std::endl;
    }
    return result; 
  }

private:
  std::vector<int> motherIds_;
  bool bsm_;
};

#endif


