import FWCore.ParameterSet.Config as cms

from SuSyAachen.Skimming.diLeptonFilter_cfi import diLeptonFilter

SSEE = diLeptonFilter.clone(
    primarySrc = "cleanLayer1Electrons",
    secondarySrc = "cleanLayer1Electrons",
)

SSMuMu = diLeptonFilter.clone(
    primarySrc = "cleanLayer1Muons",
    secondarySrc = "cleanLayer1Muons",
)

SSTauTau = diLeptonFilter.clone(
    primarySrc = "cleanLayer1Taus",
    secondarySrc = "cleanLayer1Taus",
)

SSEMu = diLeptonFilter.clone(
    primarySrc = "cleanLayer1Muons",
    secondarySrc = "cleanLayer1Electrons",
)

SSETau = diLeptonFilter.clone(
    primarySrc = "cleanLayer1Electrons",
    secondarySrc = "cleanLayer1Taus",
)

SSMuTau = diLeptonFilter.clone(
    primarySrc = "cleanLayer1Muons",
    secondarySrc = "cleanLayer1Taus",
)

seqSSDiLeptons = cms.Sequence( 
    cms.ignore( SSEE ) + cms.ignore( SSMuMu ) + cms.ignore( SSTauTau )
    + cms.ignore( SSEMu ) + cms.ignore( SSETau ) + cms.ignore( SSMuTau )
    )
