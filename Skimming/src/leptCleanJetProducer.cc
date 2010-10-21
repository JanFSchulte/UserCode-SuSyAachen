// -*- C++ -*-
//
// Package:    jetLeptonCleaner
// Class:      jetLeptonCleaner
// 
/**\class jetLeptonCleaner jetLeptonCleaner.cc SuSyAachen/Skimming/src/jetLeptonCleaner.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Wed Aug 18 15:37:34 CEST 2010
// $Id: jetLeptonCleaner.cc,v 1.1 2010/08/27 14:48:04 nmohr Exp $
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

#include "CommonTools/Utils/interface/PtComparator.h"

//
// class declaration
//
template< typename T >
class jetLeptonCleaner : public edm::EDProducer {
   public:
      explicit jetLeptonCleaner(const edm::ParameterSet&);
      ~jetLeptonCleaner();

   private:
      edm::InputTag jetSrc;
      edm::InputTag leptonSrc;
      double dRCut;
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      GreaterByPt<pat::Jet>       pTComparator_;
      
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
jetLeptonCleaner<T>::jetLeptonCleaner(const edm::ParameterSet& iConfig)
{ 
   produces< std::vector< pat::Jet > > ();
   
   leptonSrc      = iConfig.getParameter<edm::InputTag> ("leptSrc");
   jetSrc         = iConfig.getParameter<edm::InputTag> ("src");
   dRCut          = iConfig.getParameter<double> ("dRJetLepton");
}


template< typename T >
jetLeptonCleaner<T>::~jetLeptonCleaner()
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
jetLeptonCleaner<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle< std::vector< T > > leptons;
   iEvent.getByLabel(leptonSrc, leptons);

   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   std::auto_ptr<std::vector< pat::Jet > > theJets ( new std::vector< pat::Jet >() );
   double dR_ = 999999999.;
   double dRLJ = 999999999.;
   
   for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){ 
     dR_ = 999999999.;
     for (typename std::vector< T >::const_iterator l_i = leptons->begin(); l_i != leptons->end(); ++l_i){
        dRLJ = reco::deltaR(jet_i->eta(),jet_i->phi(),l_i->eta(),l_i->phi());
        if ( dRLJ < dR_ ) dR_ = dRLJ;
     }
     if (dR_ >= dRCut){
        pat::Jet goodJet = *jet_i;
        theJets->push_back(goodJet);
     }
   }
   std::sort(theJets->begin(), theJets->end(), pTComparator_);
   iEvent.put( theJets );
 
}

// ------------ method called once each job just before starting event loop  ------------
template< typename T >
void 
jetLeptonCleaner<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template< typename T >
void 
jetLeptonCleaner<T>::endJob() {
}

typedef jetLeptonCleaner< pat::Muon > jetMuonCleaner;
typedef jetLeptonCleaner< pat::Electron > jetElectronCleaner;

//define this as a plug-in
DEFINE_FWK_MODULE(jetMuonCleaner);
DEFINE_FWK_MODULE(jetElectronCleaner);
