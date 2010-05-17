#include "FWCore/Framework/interface/MakerMacros.h"

// include the definition of custom cloning procedure 
// for track collection, that clones together with Tracks, 
// also TrackExtras and RecHits
//#include "PhysicsTools/RecoAlgos/interface/TrackSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SuSyAachen/Skimming/interface/d0Selector.h"

// define your producer name
typedef ObjectSelector< d0Selector<reco::BeamSpot, pat::MuonCollection, std::vector<const pat::Muon *> > > PATMuonD0BSSelector;
typedef ObjectSelector< d0Selector<reco::BeamSpot, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronD0BSSelector;
typedef ObjectSelector< d0Selector<reco::VertexCollection, pat::MuonCollection, std::vector<const pat::Muon *> > > PATMuonD0PVSelector;
typedef ObjectSelector< d0Selector<reco::VertexCollection, pat::ElectronCollection, std::vector<const pat::Electron *> > > PATElectronD0PVSelector;

// declare the module as plugin
DEFINE_FWK_MODULE(PATMuonD0BSSelector);
DEFINE_FWK_MODULE(PATElectronD0BSSelector);
DEFINE_FWK_MODULE(PATMuonD0PVSelector);
DEFINE_FWK_MODULE(PATElectronD0PVSelector);
