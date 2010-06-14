import FWCore.ParameterSet.Config as cms

from SuSyAachen.Skimming.leptonSelectors_cfi import patLeptonCountFilter as patLeptonCountFilterOrig
patLeptonCountFilter = patLeptonCountFilterOrig.clone(
           filter = cms.bool(True),
           minNumber = 1
)
