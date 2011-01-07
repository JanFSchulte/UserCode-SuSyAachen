// -*- C++ -*-
//
// Package:    PfElectronProducer
// Class:      PfElectronProducer
// 
/**\class PfElectronProducer PfElectronProducer.cc SuSyAachen/Skimming/src/PfElectronProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Wed Aug 18 15:37:34 CEST 2010
// $Id: PfElectronProducer.cc,v 1.2 2010/08/18 20:55:13 nmohr Exp $
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
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "CommonTools/Utils/interface/PtComparator.h"


//
// class declaration
//
class PfElectronProducer : public edm::EDProducer {
   public:
      explicit PfElectronProducer(const edm::ParameterSet&);
      ~PfElectronProducer();

   private:
      edm::InputTag eleSrc;
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      GreaterByPt<pat::Electron>       pTComparator_;
      
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
PfElectronProducer::PfElectronProducer(const edm::ParameterSet& iConfig)
{ 
   produces< std::vector<pat::Electron> > ();
   
   eleSrc            = iConfig.getParameter<edm::InputTag> ("src");
}


PfElectronProducer::~PfElectronProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PfElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle< std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(eleSrc, electrons);

   std::auto_ptr<std::vector<pat::Electron> > theElectrons ( new std::vector<pat::Electron>() );
   for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ++ele_i){
        pat::Electron pfEle = *ele_i;
        pfEle.setP4( ele_i->pfCandidateRef()->p4() );
        theElectrons->push_back(pfEle);
    }
   std::sort(theElectrons->begin(), theElectrons->end(), pTComparator_);
   iEvent.put( theElectrons );
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
PfElectronProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PfElectronProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PfElectronProducer);
