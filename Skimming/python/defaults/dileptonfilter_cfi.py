import FWCore.ParameterSet.Config as cms
from SuSyAachen.Skimming.diLeptonFilter_cfi import diLeptonFilter as diLeptonFilterFilterOrig

defaultSelector = diLeptonFilterFilterOrig.clone(
           primarySrc = "cleanLayer1Electrons",
           secondarySrc = "cleanLayer1Muons",
)
