// -*- C++ -*-
//
// Package:    FastJetUnCorrJetsProducer
// Class:      FastJetUnCorrJetsProducer
// 
/**\class FastJetUnCorrJetsProducer FastJetUnCorrJetsProducer.cc SuSyAachen/Skimming/src/FastJetUnCorrJetsProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Daniel Sprenger
//         Created:  Sat Jun 18 17:53:34 CEST 2010
// $Id: FastJetUnCorrJetsProducer.cc,v 1.1 2011/06/18 17:53:59 sprenger Exp $
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

class FastJetUnCorrJetsProducer : public edm::EDProducer {
   public:
      explicit FastJetUnCorrJetsProducer(const edm::ParameterSet&);
      ~FastJetUnCorrJetsProducer();

   private:
      edm::InputTag jetSrc;
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
FastJetUnCorrJetsProducer::FastJetUnCorrJetsProducer(const edm::ParameterSet& iConfig)
{ 
   produces< std::vector<pat::Jet> > ();
   
   jetSrc            = iConfig.getParameter<edm::InputTag> ("src");
}


FastJetUnCorrJetsProducer::~FastJetUnCorrJetsProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
FastJetUnCorrJetsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   std::auto_ptr<std::vector<pat::Jet> > theJets ( new std::vector<pat::Jet>() );
   for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){
        float factor = 1.0/jet_i->jecFactor("L1FastJet");
	//pat::Jet rescaledJet = jet_i->correctedJet("L1FastJet");
        pat::Jet uncorrectedJet = jet_i->correctedJet("Uncorrected");
	uncorrectedJet.scaleEnergy(factor);
        theJets->push_back(uncorrectedJet);
    }
   std::sort(theJets->begin(), theJets->end(), PtGreater());
   iEvent.put( theJets );
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
FastJetUnCorrJetsProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FastJetUnCorrJetsProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(FastJetUnCorrJetsProducer);