#ifndef SuSyAachen_Skimming_PromptSelector_h
#define SuSyAachen_Skimming_PromptSelector_h

#include <vector>
#include <math.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>

struct PromptSelector {
  typedef reco::CompositeRefCandidateT<reco::GenParticleRefVector>::daughters daughters;

PromptSelector( std::vector<int> promptMothers, bool bsm) : motherIds_( promptMothers ), bsm_(bsm) { }
  PromptSelector( const edm::ParameterSet & cfg ) :
  motherIds_( cfg.getParameter< std::vector<int> >( "promptMotherIds" ) ),
  leptonicDecayIds_( cfg.getParameter< std::vector<int> >( "leptonicDecayIds" ) ),
    bsm_( cfg.getParameter< bool >( "bsm" )){ }

  bool operator()( const reco::GenParticle & p ) const 
  {   
    bool result = false;
    if(p.status() == 1 && p.numberOfMothers() == 1){
      //      std::cout << "p "<< p.pdgId()<<" | "<<p.mother()->pdgId()<<": ";
      for( std::vector<int>::const_iterator i = motherIds_.begin();
	   i != motherIds_.end(); ++i ){
	if( p.mother()->pdgId() == (*i)) result |= true;
	//	std::cout << result<<" ";
      }
      
      result |= bsm_ && p.mother()->pdgId() > 1000000;
      //std::cout << result<<" ";
      result |= p.mother()->pdgId() == p.pdgId();
      //std::cout << result<<std::endl;
    }
    
    if( abs(p.pdgId()) == 15 && p.status() ==2 
	&& abs(p.mother()->pdgId()) == 15 && p.mother()->status()==3){
      bool leptonicDecay = false;
      const reco::GenParticleRefVector& daughters = p.daughterRefVector();
      for( reco::GenParticleRefVector::const_iterator it = daughters.begin(); it != daughters.end(); ++it){
	for( std::vector<int>::const_iterator i = leptonicDecayIds_.begin();
	     i != leptonicDecayIds_.end(); ++i ){
	  leptonicDecay |= (*it)->pdgId() == (*i) && (*it)->status() == 1;
	}
      }
      result |= !leptonicDecay;
    }
    //if(result)     std::cout << abs(p.pdgId()) << " "<< (abs(p.pdgId()) == 15) << std::endl;
    return result; 
  }

private:
  std::vector<int> motherIds_;
  std::vector<int> leptonicDecayIds_;
  bool bsm_;
};

#endif


