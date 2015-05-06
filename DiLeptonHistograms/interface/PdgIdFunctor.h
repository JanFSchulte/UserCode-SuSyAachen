/*
 * PdgIdFunctor.h
 *
 *  Created on: 21.12.2010
 *      Author: heron
 */

#ifndef PDGIDFUNCTOR_H_
#define PDGIDFUNCTOR_H_

#include <vector>
#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//#include <DataFormats/PatCandidates/interface/Electron.h>
//#include <DataFormats/PatCandidates/interface/Muon.h>
//#include <DataFormats/PatCandidates/interface/Tau.h>

using namespace std;

class PdgIdFunctor
{
public:
 PdgIdFunctor(edm::ParameterSet const & params):genParticles_(){
    hasGen_ = false;
    deltaR_ = params.getParameter<double>("deltaR");
    genTag_ = params.getParameter<edm::InputTag> ("genSrc");
  }
  
  void loadGenParticles( const edm::Event& ev){
    edm::Handle<reco::GenParticleCollection> genParticles;
    hasGen_ = ev.getByLabel(genTag_, genParticles);
    if(hasGen_) genParticles_ = *genParticles;
  }
  
  bool hasGen(){return hasGen_;}
  template<class T> std::vector<int> operator()(const T& lepton);
  
 private:
  bool hasGen_;
  reco::GenParticleCollection genParticles_; 
  
  double deltaR_;
  edm::InputTag genTag_;  
};

template<class T>
std::vector<int>
PdgIdFunctor::operator()(const T& lepton)
{

	double minDR = 999.;
	double deltaR = 0.;
	const reco::GenParticle * matchedGenPart = new reco::GenParticle();
	const reco::GenParticle * matchedGenPartMother = new reco::GenParticle();
	const reco::GenParticle * matchedGenPartGrandMother = new reco::GenParticle();	
	bool matched = false;

	int pdgId = -9999;
	int motherPdgId = -9999;
	int grandMotherPdgId = -9999;
	
	 std::vector<int> res;
	 
	for (std::vector<reco::GenParticle>::const_iterator itGenParticle = genParticles_.begin(); itGenParticle != genParticles_.end(); itGenParticle++) {
	    if (itGenParticle->status() != 1) continue;
	  
	    deltaR = reco::deltaR(lepton.eta(),lepton.phi(), itGenParticle->eta(), itGenParticle->phi());
	    if (deltaR > 0.1) continue;
	    
	    double ndpt = fabs(lepton.pt() - itGenParticle->pt())/itGenParticle->pt();
	    if(ndpt > 2.) continue;
	    
	    if(deltaR > minDR) continue;
	    
	    minDR = deltaR;
	    
	    matched = true;
	    matchedGenPart = &(*itGenParticle);

	}
	if (matched == true){
	  pdgId = matchedGenPart->pdgId();
	  if (!matchedGenPart->mother()){
	      motherPdgId = -9999;
	      grandMotherPdgId = -9999;
	  }
	  else{
	      matchedGenPartMother = static_cast<const reco::GenParticle*> (matchedGenPart->mother());
	      if (matchedGenPartMother->pdgId() != pdgId){
		    motherPdgId = matchedGenPartMother->pdgId();
	      }
	      else{
		    int tempid = matchedGenPartMother->pdgId();
		    int loop_counter = 0;
		    while(tempid == pdgId){
		       loop_counter++;
		       if(loop_counter>=10){
		         matchedGenPartMother = static_cast<const reco::GenParticle*>(matchedGenPart->mother());
		         break;
		       }
		  	   if (!matchedGenPartMother->mother()){
		  	   	  break;
		  	   }
		  	   else{
				   matchedGenPart = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());
				   matchedGenPartMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());		  
				   tempid = matchedGenPartMother->pdgId();
			   }   
		     }
		
	      }
	      motherPdgId = matchedGenPartMother->pdgId();
	  
	      if (!matchedGenPartMother->mother()){
		  grandMotherPdgId = -9999;
	      }	  
	      else{
		  matchedGenPartGrandMother = static_cast<const reco::GenParticle*> (matchedGenPartMother->mother());
		  if (matchedGenPartGrandMother->pdgId() != motherPdgId){
		    grandMotherPdgId = matchedGenPartGrandMother->pdgId();
		  }
		  else{
		    
		    int tempid = matchedGenPartGrandMother->pdgId();
		    int loop_counter = 0;
		    while(tempid == motherPdgId){
		      loop_counter++;
		      if(loop_counter>=10){
			matchedGenPartGrandMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());
			break;
		      }
		      if (!matchedGenPartMother->mother()){
		         break;
		      }
		      else{
				  matchedGenPartMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());
				  matchedGenPartGrandMother = static_cast<const reco::GenParticle*>(matchedGenPartMother->mother());		  
				  tempid = matchedGenPartGrandMother->pdgId();
			  }
		    }
		    
		  }
		  grandMotherPdgId = matchedGenPartGrandMother->pdgId();
	      }
	  }      
	}
	res.push_back(pdgId);
	res.push_back(motherPdgId);
	res.push_back(grandMotherPdgId);

	return res;
}


#endif /* PDGIDFUNCTOR_H_ */
