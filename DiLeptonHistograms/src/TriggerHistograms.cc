// -*- C++ -*-
//
// Package:    TriggerHistograms
// Class:      TriggerHistograms
// 
/**\class TriggerHistograms TriggerHistograms.cc SuSyAachen/TriggerHistograms/src/TriggerHistograms.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Niklas Mohr
//         Created:  Mon Jun 28 13:06:46 CEST 2010
// $Id$
//
//


// system include files
#include <memory>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"
//
// class declaration
//

class TriggerHistograms : public edm::EDAnalyzer {
   public:
      explicit TriggerHistograms(const edm::ParameterSet&);
      ~TriggerHistograms();

   private:
      edm::InputTag trgSrc_;
      double weight_;
      TFile *theFile;
      TH1F *hTrigger;

      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
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
TriggerHistograms::TriggerHistograms(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   
   trgSrc_            = iConfig.getParameter<edm::InputTag> ("triggerSource");
   weight_   = iConfig.getUntrackedParameter<double>   ("external_Weight",1.);
   edm::Service<TFileService> theFile;
   hTrigger = theFile->make<TH1F>( "Trigger paths", "Trigger paths", 160, 0, 160);


}


TriggerHistograms::~TriggerHistograms()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TriggerHistograms::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle< edm::TriggerResults > trigger;
    iEvent.getByLabel(trgSrc_, trigger);
    const edm::TriggerNames & names = iEvent.triggerNames(*trigger);

    std::vector< std::string > hlNames;
    for (edm::TriggerNames::Strings::const_iterator j = names.triggerNames().begin(); j !=names.triggerNames().end(); ++j ) {
        hlNames.push_back(*j);
    }
    //LogPrint("Trigger") << "Event" << "\n";
    for (unsigned int i=0; i<trigger->size()-1; ++i){
        hTrigger->Fill(hlNames[i].c_str(),0);
        if((*trigger)[i].accept()){
                hTrigger->Fill(hlNames[i].c_str(),weight_);
        }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerHistograms::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerHistograms::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerHistograms);
