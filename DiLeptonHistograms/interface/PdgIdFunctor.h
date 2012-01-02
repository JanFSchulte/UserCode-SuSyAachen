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

class PdgIdFunctor
{
public:
        PdgIdFunctor(edm::ParameterSet const & params):genParticles_(){
	  deltaR_ = params.getParameter<double>("deltaR");
	  genTag_ = params.getParameter<edm::InputTag> ("genSrc");
	}
	
	void loadGenParticles( const edm::Event& ev){
	  edm::Handle<reco::GenParticleCollection> genParticles;
	  ev.getByLabel(genTag_, genParticles);
	  genParticles_ = *genParticles;
	}

	bool isUseable(){return initialized_;}
	template<class T> int operator()(const T& lepton);

private:
	bool initialized_;
	reco::GenParticleCollection genParticles_; 

	double deltaR_;
	edm::InputTag genTag_;

};

template<class T>
int
PdgIdFunctor::operator()(const T& lepton)
{
   int result = -9999;
   if(lepton.genLepton() == NULL){
     double  minDeltaR = deltaR_;
     for(reco::GenParticleCollection::const_iterator it = genParticles_.begin(); it != genParticles_.end(); ++it){
       if((*it).status() == 3 && reco::deltaR<const T, const reco::GenParticle>( lepton, *it) < minDeltaR){
	 minDeltaR = reco::deltaR<const T, const reco::GenParticle>( lepton, *it);
	 result = (*it).pdgId();
       }
     }
   }else{
     reco::GenParticle p(*lepton.genLepton());
     if(p.status() == 3)
       result = p.pdgId();
     else if(p.mother() != NULL){
       if(abs(p.mother()->pdgId()) == 11 || abs(p.mother()->pdgId()) == 13 || abs(p.mother()->pdgId()) == 15)
	 result = p.mother()->pdgId();
       else
	 result = p.pdgId();
     }
   }   
   return result;
}


#endif /* PDGIDFUNCTOR_H_ */
