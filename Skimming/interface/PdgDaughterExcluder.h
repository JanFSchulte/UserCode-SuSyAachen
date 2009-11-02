#ifndef SuSyAachen_Skimming_PdgDaughterExcluder_h
#define SuSyAachen_Skimming_PdgDaughterExcluder_h

#include <vector>

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
//include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>

struct PdgDaughterExcluder {
  typedef reco::CompositeRefCandidateT<reco::GenParticleRefVector>::daughters daughters;

  PdgDaughterExcluder( std::vector<int> daughters) : daughterIds_( daughters ) { }
  PdgDaughterExcluder( const edm::ParameterSet & cfg ) :
  daughterIds_( cfg.getParameter< std::vector<int> >( "daughterIds" ) ) { }

  bool operator()( const reco::GenParticle & p ) const 
  { 
    bool result = true;
    //    std::cout << "id: "<< p.pdgId() << std::endl;
    if(p.status() == 1){
      for( std::vector<int>::const_iterator i = daughterIds_.begin();
	   i != daughterIds_.end(); ++i ){
	if( p.pdgId() == (*i)) result = false;
	//	std::cout << "---->" << (*i) <<" = "<<result<<std::endl;
      }
    }else{
      for( daughters::const_iterator i = p.daughterRefVector().begin();
	   i != p.daughterRefVector().end(); ++i ){
	result = result & operator()(*(*i));
      }
    }
    //    std::cout << "result: "<< result << std::endl;
    return result; 
  }

private:
  std::vector<int> daughterIds_;

};

#endif


