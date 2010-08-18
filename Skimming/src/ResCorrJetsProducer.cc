// -*- C++ -*-
//
// Package:    ResCorrJetsProducer
// Class:      ResCorrJetsProducer
// 
/**\class ResCorrJetsProducer ResCorrJetsProducer.cc SuSyAachen/Skimming/src/ResCorrJetsProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Wed Aug 18 15:37:34 CEST 2010
// $Id: ResCorrJetsProducer.cc,v 1.1 2010/08/18 19:23:06 nmohr Exp $
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

class ResCorrJetsProducer : public edm::EDProducer {
   public:
      explicit ResCorrJetsProducer(const edm::ParameterSet&);
      ~ResCorrJetsProducer();

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
ResCorrJetsProducer::ResCorrJetsProducer(const edm::ParameterSet& iConfig)
{ 
   produces< std::vector<pat::Jet> > ();
   
   jetSrc            = iConfig.getParameter<edm::InputTag> ("src");
   jetCorrections    = iConfig.getParameter<std::string> ("jetCorrections");

   edm::FileInPath fipRes(jetCorrections);
   JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
   std::vector<JetCorrectorParameters> vParam;
   vParam.push_back(*ResJetCorPar);
   JEC = new FactorizedJetCorrector(vParam);
}


ResCorrJetsProducer::~ResCorrJetsProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ResCorrJetsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);
   bool isRealData = iEvent.isRealData();

   std::auto_ptr<std::vector<pat::Jet> > theJets ( new std::vector<pat::Jet>() );
   for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){
        JEC->setJetEta(jet_i->eta());
        JEC->setJetPt(jet_i->pt());
        pat::Jet rescaledJet = *jet_i;
        if (isRealData) rescaledJet.scaleEnergy(JEC ? JEC->getCorrection() : 1.);
        theJets->push_back(rescaledJet);
    }
   std::sort(theJets->begin(), theJets->end(), PtGreater());
   iEvent.put( theJets );
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
ResCorrJetsProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ResCorrJetsProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ResCorrJetsProducer);
