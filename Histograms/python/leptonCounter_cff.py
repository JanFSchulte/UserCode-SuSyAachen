import FWCore.ParameterSet.Config as cms
import SuSyAachen.Histograms.leptonCounter_cfi 

############# reco slevel #####################
#---- All reco particles
anyRecoLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "cleanLayer1Electrons",
  muonSource = "cleanLayer1Muons",
  tauSource = "cleanLayer1Taus",
) 

#---- basic reco particles
basicRecoLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "basicElectrons",
  muonSource = "basicMuons",
  tauSource = "basicTaus",
) 

#---- isolated reco particles
isoRecoLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "isoElectrons",
  muonSource = "isoMuons",
  tauSource = "isoTaus",
) 

################ gen Level #####################
#---- All gen particles
anyGenLeptonExclusiveCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  method = 'exclusive',
  electronSource = "electronGenParticles",
  muonSource = "muonGenParticles",
  tauSource = "tauGenParticles",
)

anyGenLeptonInclusiveCounter = anyGenLeptonExclusiveCounter.clone(
  method = "inclusive"
)

anyMatchedLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "anyElectronMatchedParticles",
  muonSource = "anyMuonMatchedParticles",
  tauSource = "anyTauJetMatchedParticles"#anyTauMatchedParticles",
)

basicMatchedLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "basicElectronMatchedParticles",
  muonSource = "basicMuonMatchedParticles",
  tauSource = "basicTauJetMatchedParticles"
)

isoMatchedLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "isoElectronMatchedParticles",
  muonSource = "isoMuonMatchedParticles",
  tauSource = "isoTauJetMatchedParticles",
)


#---- Cut Gen Particles
basicGenLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "electronBasicGenParticles",
  muonSource = "muonBasicGenParticles",
  tauSource = "tauBasicGenParticles",
)

seqLeptonCounter = cms.Sequence( 
  anyGenLeptonExclusiveCounter 
  + anyGenLeptonInclusiveCounter 
  + basicGenLeptonCounter
  + anyMatchedLeptonCounter
  + basicMatchedLeptonCounter 
  + isoMatchedLeptonCounter

  + anyRecoLeptonCounter 
  + isoRecoLeptonCounter 
  + basicRecoLeptonCounter
)
