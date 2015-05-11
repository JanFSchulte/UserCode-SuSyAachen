#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "SuSyAachen/Skimming/interface/effectiveAreaIsolationSelector.h"
#include "SuSyAachen/Skimming/interface/matchedConversionSelector.h"
#include "SuSyAachen/Skimming/interface/eleIDSelector.h"
// define your producer name
typedef ObjectSelector< effectiveAreaIsolationSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronEffectiveAreaSelector;

typedef ObjectSelector< matchedConversionSelector<reco::ConversionCollection, reco::BeamSpot, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronMatchedConversionSelector;

typedef ObjectSelector< eleIDSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronIDSelector;

// declare the module as plugin
DEFINE_FWK_MODULE( PATElectronEffectiveAreaSelector );
DEFINE_FWK_MODULE( PATElectronMatchedConversionSelector );
DEFINE_FWK_MODULE( PATElectronIDSelector );

