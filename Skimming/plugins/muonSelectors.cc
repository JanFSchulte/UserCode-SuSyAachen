#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SuSyAachen/Skimming/interface/tightMuonSelector.h"

// define your producer name
typedef ObjectSelector< tightMuonSelector<reco::BeamSpot, pat::MuonCollection, std::vector<const pat::Muon *> > > PATMuonTightIDSelector;


// declare the module as plugin
DEFINE_FWK_MODULE( PATMuonTightIDSelector );


