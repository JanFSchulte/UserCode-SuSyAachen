#include "SuSyAachen/TagAndProbeTreeWriter/interface/LeptonKindFunctor.h"

#include <iostream>

LeptonKindFunctor::LeptonKindFunctor()
{
}

const int LeptonKindFunctor::PromptCategory(const reco::Candidate * genParticle){
    int value = 0;
    if(genParticle->numberOfMothers()==1){
      //Check if lepton is promt (itself,Z,W,SUSY)
      const reco::Candidate * mom = genParticle->mother();
      if( (mom->pdgId()==genParticle->pdgId() && ( abs(genParticle->pdgId()) != 15 || PromptCategory(mom) == 3 ))  // itself (taus need to come from a prompt mom)
	  ||abs(mom->pdgId())==23 ||abs(mom->pdgId())==24 // Z or W 
	  ||abs(mom->pdgId())>1000000){ // SUSY
	if(genParticle->status()==1 || abs(genParticle->pdgId())==15 )//only status 1 particles are prompt (excpet for taus)
	  value = 3;
      } else if( abs(mom->pdgId())==15 ){ //if the mother is a tau the category is that of the tau
	int motherCat = PromptCategory(mom);
	//	std::cout << "tau mother: "<<motherCat<<" | ";
	value = motherCat;
      }
      else if( mom->pdgId()==443 || mom->pdgId()==553 ||mom->pdgId()==100553) // J/psi (s1) upsilon(s1) upsilon(s2)
	value = 5;
      else
	value = 4;
      //      std::cout <<"mom: " << mom->pdgId();
    }
    else value = 4;
    //    std::cout << " result: "<< value <<std::endl;
    return value;
}



