// -*- C++ -*-
//
// Package:    GenParticlePatProducer
// Class:      GenParticlePatProducer
// 
/**\class GenParticlePatProducer bJetLeptonProducer.cc SuSyAachen/Skimming/src/bJetLeptonProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Wed Aug 18 15:37:34 CEST 2010
// $Id: GenParticleProducer.cc,v 1.1 2011/05/18 10:55:37 nmohr Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonTools/Utils/interface/PtComparator.h"

//
// class declaration
//
template< typename T >
class GenParticlePatProducer : public edm::EDProducer {
   public:
      explicit GenParticlePatProducer(const edm::ParameterSet&);
      ~GenParticlePatProducer();

   private:
      edm::InputTag leptonSrc;
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      GreaterByPt<reco::GenParticle>    pTComparator_;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
template< typename T >
GenParticlePatProducer<T>::GenParticlePatProducer(const edm::ParameterSet& iConfig)
{ 
   produces< std::vector< reco::GenParticle > > ();
   
   leptonSrc      = iConfig.getParameter<edm::InputTag> ("src");
}


template< typename T >
GenParticlePatProducer<T>::~GenParticlePatProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template< typename T >
void
GenParticlePatProducer<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle< std::vector< T > > leptons;
   iEvent.getByLabel(leptonSrc, leptons);

   std::auto_ptr<std::vector< reco::GenParticle > > theGens ( new std::vector< reco::GenParticle >() );
   
   for (typename std::vector< T >::const_iterator l_i = leptons->begin(); l_i != leptons->end(); ++l_i){
     if (l_i->genLepton()){
        reco::GenParticle genLepton = *l_i->genLepton();
        theGens->push_back(genLepton);
     }
   }
   std::sort(theGens->begin(), theGens->end(), pTComparator_);
   iEvent.put( theGens );
 
}

// ------------ method called once each job just before starting event loop  ------------
template< typename T >
void 
GenParticlePatProducer<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T >
void 
GenParticlePatProducer<T>::endJob() {
}

typedef GenParticlePatProducer< pat::Muon > GenParticleMuonProducer;
typedef GenParticlePatProducer< pat::Electron > GenParticleElectronProducer;

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleMuonProducer);
DEFINE_FWK_MODULE(GenParticleElectronProducer);
