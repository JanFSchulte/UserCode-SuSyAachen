#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "SuSyAachen/Skimming/interface/isolationSelector.h"
#include "SuSyAachen/Skimming/interface/eleMVAIDSelector.h"
// define your producer name

typedef ObjectSelector< isolationSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronIsolationSelector;

typedef ObjectSelector< eleMVAIDSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronMVAIDSelector;


// declare the module as plugin
DEFINE_FWK_MODULE( PATElectronIsolationSelector );
DEFINE_FWK_MODULE( PATElectronMVAIDSelector );

