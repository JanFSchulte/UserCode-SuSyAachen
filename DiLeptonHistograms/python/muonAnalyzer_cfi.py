import FWCore.ParameterSet.Config as cms
from MuonAnalysis.Examples.inclusiveMuonPlots_cfi import *

muonPlots = cms.EDAnalyzer("InclusiveMuonPlots",
                           makeInclusiveMuonPlots(1),  # change "1" to your rebin factor; it can be a fractionary number, like 0.5
                           muons     = cms.InputTag('muons'),
                           selection = cms.string("isGlobalMuon && muonID('GlobalMuonPromptTight')"), # or anything else
                           primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                           )
