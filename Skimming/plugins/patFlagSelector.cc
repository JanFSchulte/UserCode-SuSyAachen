#include "FWCore/Framework/interface/MakerMacros.h"

// include the definition of custom cloning procedure 
// for track collection, that clones together with Tracks, 
// also TrackExtras and RecHits
//#include "PhysicsTools/RecoAlgos/interface/TrackSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "SuSyAachen/Skimming/interface/patFlagSelector.h"

// define your producer name
typedef ObjectSelector< patFlagSelector<pat::JetCollection, std::vector<const pat::Jet* > > > PATJetFlagSelector;

// declare the module as plugin
DEFINE_FWK_MODULE(PATJetFlagSelector);
