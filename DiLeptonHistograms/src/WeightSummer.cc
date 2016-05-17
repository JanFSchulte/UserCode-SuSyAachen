/** \class WeightSummer
 *
 *  
 *  This class is an EDAnalyzer for PAT
 *  Layer 0 and Layer 1 output
 *
 *  $Date: 2011/12/05 13:47:01 $
 *  $Revision: 1.1 $
 *  for CMSSW_4_2_8
 *  \author Daniel Sprenger  --  daniel.sprenger@cern.ch
 *
 */


#include "SuSyAachen/DiLeptonHistograms/interface/WeightSummer.h"


//Constructor
WeightSummer::WeightSummer(const edm::ParameterSet &iConfig):
fctVtxWeight_    (iConfig.getParameter<edm::ParameterSet>("vertexWeights"),consumesCollector()  )
{
    // Create the root file
    edm::Service<TFileService> theFile;

    hWeights = theFile->make<TH1F>( "Weights", "Weights", 1, 0.5, 1.5);
}

//Destructor
WeightSummer::~WeightSummer()
{
} 


//Event loop
void WeightSummer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  float weight = fctVtxWeight_(iEvent);
  hWeights->Fill(1.0, weight);
}


DEFINE_FWK_MODULE(WeightSummer);
