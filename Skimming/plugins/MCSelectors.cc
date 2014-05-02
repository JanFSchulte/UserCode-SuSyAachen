#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SuSyAachen/Skimming/interface/PdgDaughterExcluder.h"
#include "SuSyAachen/Skimming/interface/PromptSelector.h"
#include "SuSyAachen/Skimming/interface/matchedSelector.h"
#include "SuSyAachen/Skimming/interface/matchedTauJetSelector.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

// typedef SingleObjectSelector< reco::GenParticleCollection, PdgDaughterExcluder > GenDaughterExcluder;
// typedef SingleObjectSelector< reco::GenParticleCollection, PromptSelector > GenPromptSelector;
// 
// typedef SingleObjectSelector< pat::ElectronCollection, MatchedSelector<pat::Electron> > PATElectronMatchedSelector;
// typedef SingleObjectSelector< pat::MuonCollection, MatchedSelector<pat::Muon> > PATMuonMatchedSelector;
// typedef SingleObjectSelector< pat::TauCollection, MatchedSelector<pat::Tau> > PATTauMatchedSelector;
// 
// typedef SingleObjectSelector< pat::TauCollection, MatchedTauJetSelector > PATTauJetMatchedSelector;
// 
// 
// DEFINE_FWK_MODULE(GenDaughterExcluder);
// DEFINE_FWK_MODULE(GenPromptSelector);
// //--
// DEFINE_FWK_MODULE(PATElectronMatchedSelector);
// DEFINE_FWK_MODULE(PATMuonMatchedSelector);
// DEFINE_FWK_MODULE(PATTauMatchedSelector);
// //
// DEFINE_FWK_MODULE(PATTauJetMatchedSelector);
