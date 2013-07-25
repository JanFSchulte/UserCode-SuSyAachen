#ifndef WeightSummer_h
#define WeightSummer_h

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

// system include files
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <vector>
#include <map>
#include <utility> 
#include <memory>
#include <iostream>
#include <fstream>
#include <string>

// user include files
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SuSyAachen/DiLeptonHistograms/interface/Combinations.h"
#include "SuSyAachen/DiLeptonHistograms/interface/WeightFunctor.h"
#include "SuSyAachen/DiLeptonHistograms/interface/VertexWeightFunctor.h"


class WeightSummer : public edm::EDAnalyzer {
    public:
    
    /// Constructor
    WeightSummer(const edm::ParameterSet &iConfig);

    /// Destructor
    virtual ~WeightSummer();

    /// Perform the real analysis
    void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

    private:

    //Histograms
    TH1F * hWeights;

    // The file which will store the histos
    TFile * theFile;

    //weight definition
    VertexWeightFunctor fctVtxWeight_;
};

#endif
