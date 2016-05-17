// -*- C++ -*-
//
// Package:    UnCorrJetsProducer
// Class:      UnCorrJetsProducer
// 
/**\class UnCorrJetsProducer UnCorrJetsProducer.cc SuSyAachen/Skimming/src/UnCorrJetsProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Wed Aug 18 15:37:34 CEST 2010
// $Id: UnCorrJetsProducer.cc,v 1.2 2011/06/20 09:57:39 sprenger Exp $
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

class UnCorrJetsProducer : public edm::EDProducer {
   public:
      explicit UnCorrJetsProducer(const edm::ParameterSet&);
      ~UnCorrJetsProducer();

   private:
      edm::EDGetTokenT< std::vector<pat::Jet> > jetToken_;
      std::string jetCorrections;
      FactorizedJetCorrector* JEC;
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
UnCorrJetsProducer::UnCorrJetsProducer(const edm::ParameterSet& iConfig):
  jetToken_(consumes< std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src")))
{ 
   produces< std::vector<pat::Jet> > ();
}


UnCorrJetsProducer::~UnCorrJetsProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
UnCorrJetsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByToken(jetToken_, jets);

   std::auto_ptr<std::vector<pat::Jet> > theJets ( new std::vector<pat::Jet>() );
   for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){
        pat::Jet rescaledJet = jet_i->correctedJet("Uncorrected");
        theJets->push_back(rescaledJet);
    }
   std::sort(theJets->begin(), theJets->end(), PtGreater());
   iEvent.put( theJets );
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
UnCorrJetsProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
UnCorrJetsProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(UnCorrJetsProducer);
