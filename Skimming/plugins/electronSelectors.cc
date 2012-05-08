#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "SuSyAachen/Skimming/interface/effectiveAreaIsolationSelector.h"

// define your producer name
typedef ObjectSelector< effectiveAreaIsolationSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronEffectiveAreaSelector;


// declare the module as plugin
DEFINE_FWK_MODULE( PATElectronEffectiveAreaSelector );

