#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SuSyAachen/Skimming/interface/PdgDaughterExcluder.h"

typedef SingleObjectSelector< reco::GenParticleCollection, PdgDaughterExcluder > GenDaughterExcluder;

DEFINE_FWK_MODULE(GenDaughterExcluder);
