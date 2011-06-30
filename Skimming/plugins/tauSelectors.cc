#include "FWCore/Framework/interface/MakerMacros.h"


#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "SuSyAachen/Skimming/interface/TauToTauMatchSelector.h"

// define your producer name
typedef ObjectSelector< TauToTauMatchSelector<pat::TauCollection, std::vector<const pat::Tau *> > > PATTauToTauMatchSelector;


// declare the module as plugin
DEFINE_FWK_MODULE(PATTauToTauMatchSelector);
