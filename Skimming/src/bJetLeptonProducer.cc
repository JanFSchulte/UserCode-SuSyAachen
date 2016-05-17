// -*- C++ -*-
//
// Package:    bJetLeptonProducer
// Class:      bJetLeptonProducer
// 
/**\class bJetLeptonProducer bJetLeptonProducer.cc SuSyAachen/Skimming/src/bJetLeptonProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Wed Aug 18 15:37:34 CEST 2010
// $Id: bJetLeptonProducer.cc,v 1.1 2010/08/27 14:48:04 nmohr Exp $
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
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


//
// class declaration
//
class PtGreater {
    public:
        template<typename T> bool operator ()(const T& i, const T& j) {
            return (i.pt() > j.pt());
    }
};

template< typename T >
class bJetLeptonProducer : public edm::EDProducer {
   public:
      explicit bJetLeptonProducer(const edm::ParameterSet&);
      ~bJetLeptonProducer();

   private:
      edm::EDGetTokenT< std::vector<pat::Jet> > jetToken_;
      edm::EDGetTokenT< std::vector< T > > leptonToken_;
      double dRCut;
      double dPhiCut;
      double bTagCut;
      std::string bJetAlgo;
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
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
bJetLeptonProducer<T>::bJetLeptonProducer(const edm::ParameterSet& iConfig):
   jetToken_(consumes< std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src"))),
   leptonToken_(consumes< std::vector< T > >(iConfig.getParameter<edm::InputTag>("leptSrc")))
{ 
   produces< std::vector< T > > ();
   
   dRCut          = iConfig.getParameter<double> ("dRJetLepton");
   dPhiCut        = iConfig.getParameter<double> ("dPhiOppositeJetLepton");
   bJetAlgo       = iConfig.getParameter<std::string> ("user_bJetAlgo");
   bTagCut        = iConfig.getParameter<double> ("user_bTagDiscriminator");
}


template< typename T >
bJetLeptonProducer<T>::~bJetLeptonProducer()
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
bJetLeptonProducer<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByToken(jetToken_, jets);
   
   edm::Handle< std::vector< T > > leptons;
   iEvent.getByToken(leptonToken_, leptons);

   std::auto_ptr<std::vector< T > > theLeptons ( new std::vector< T >() );
   double dR_ = 999999999.;
   double dRLJ = 999999999.;
   double dPhi_ = 0.;
   double dPhiLJ = 0.;
   
   for (typename std::vector< T >::const_iterator l_i = leptons->begin(); l_i != leptons->end(); ++l_i){
     dR_ = 999999999.;
     dPhi_ = 0.;
     for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){ 
        dRLJ = reco::deltaR(jet_i->eta(),jet_i->phi(),l_i->eta(),l_i->phi());
        dPhiLJ = reco::deltaPhi(jet_i->phi(), l_i->phi()); 
        if ( dRLJ < dR_ ) dR_ = dRLJ;
        if ( jet_i->bDiscriminator(bJetAlgo) > bTagCut && dPhiLJ > dPhi_  ) dPhi_ = dPhiLJ;
     }
     if (dR_ <= dRCut && dPhi_ >= dPhiCut){
        T bJetLepton = *l_i;
        theLeptons->push_back(bJetLepton);
     }
   }
   std::sort(theLeptons->begin(), theLeptons->end(), PtGreater());
   iEvent.put( theLeptons );
 
}

// ------------ method called once each job just before starting event loop  ------------
template< typename T >
void 
bJetLeptonProducer<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T >
void 
bJetLeptonProducer<T>::endJob() {
}

typedef bJetLeptonProducer< pat::Muon > bJetMuonProducer;
typedef bJetLeptonProducer< pat::Electron > bJetElectronProducer;

//define this as a plug-in
DEFINE_FWK_MODULE(bJetMuonProducer);
DEFINE_FWK_MODULE(bJetElectronProducer);
