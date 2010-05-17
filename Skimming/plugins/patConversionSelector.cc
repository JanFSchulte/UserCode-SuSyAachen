#include "FWCore/Framework/interface/MakerMacros.h"

// include the definition of custom cloning procedure 
// for track collection, that clones together with Tracks, 
// also TrackExtras and RecHits
//#include "PhysicsTools/RecoAlgos/interface/TrackSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "SuSyAachen/Skimming/interface/patConversionSelector.h"

// define your producer name
typedef ObjectSelector< patConversionSelector<pat::ElectronCollection, std::vector<const pat::Electron* > > > PATElectronConversionSelector;

// declare the module as plugin
DEFINE_FWK_MODULE(PATElectronConversionSelector);
