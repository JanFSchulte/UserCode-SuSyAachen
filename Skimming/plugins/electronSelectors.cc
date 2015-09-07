#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "SuSyAachen/Skimming/interface/effectiveAreaIsolationSelector.h"
#include "SuSyAachen/Skimming/interface/isolationSelector.h"
#include "SuSyAachen/Skimming/interface/matchedConversionSelector.h"
#include "SuSyAachen/Skimming/interface/eleIDSelector.h"
#include "SuSyAachen/Skimming/interface/eleMVAIDSelector.h"
#include "SuSyAachen/Skimming/interface/eleLooseMVAIDSelector.h"
// define your producer name
typedef ObjectSelector< effectiveAreaIsolationSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronEffectiveAreaSelector;

typedef ObjectSelector< isolationSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronIsolationSelector;

typedef ObjectSelector< matchedConversionSelector<reco::ConversionCollection, reco::BeamSpot, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronMatchedConversionSelector;

typedef ObjectSelector< eleIDSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronIDSelector;

typedef ObjectSelector< eleMVAIDSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronMVAIDSelector;

typedef ObjectSelector< eleLooseMVAIDSelector<double, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronLooseMVAIDSelector;

// declare the module as plugin
DEFINE_FWK_MODULE( PATElectronEffectiveAreaSelector );
DEFINE_FWK_MODULE( PATElectronIsolationSelector );
DEFINE_FWK_MODULE( PATElectronMatchedConversionSelector );
DEFINE_FWK_MODULE( PATElectronIDSelector );
DEFINE_FWK_MODULE( PATElectronMVAIDSelector );
DEFINE_FWK_MODULE( PATElectronLooseMVAIDSelector );

