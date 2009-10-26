#include "FWCore/Framework/interface/MakerMacros.h"

// include the definition of custom cloning procedure 
// for track collection, that clones together with Tracks, 
// also TrackExtras and RecHits
//#include "PhysicsTools/RecoAlgos/interface/TrackSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "SuSyAachen/Skimming/interface/d0Selector.h"

// define your producer name
typedef ObjectSelector< d0Selector<pat::MuonCollection, std::vector<const pat::Muon *> > > PATMuonD0Selector;
typedef ObjectSelector< d0Selector<pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronD0Selector;

// declare the module as plugin
DEFINE_FWK_MODULE(PATMuonD0Selector);
DEFINE_FWK_MODULE(PATElectronD0Selector);
