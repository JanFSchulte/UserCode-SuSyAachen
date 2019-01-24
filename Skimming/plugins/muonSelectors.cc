#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SuSyAachen/Skimming/interface/mediumMuonSelector.h"
#include "SuSyAachen/Skimming/interface/isolationSelector.h"

// define your producer name
typedef ObjectSelector< mediumMuonSelector<reco::BeamSpot, pat::MuonCollection, std::vector<const pat::Muon *> > > PATMuonMediumIDSelector;
typedef ObjectSelector< isolationSelector<double, pat::MuonCollection, std::vector<const pat::Muon *> > > PATMuonIsolationSelector;

// declare the module as plugin
DEFINE_FWK_MODULE( PATMuonMediumIDSelector );
DEFINE_FWK_MODULE( PATMuonIsolationSelector );


