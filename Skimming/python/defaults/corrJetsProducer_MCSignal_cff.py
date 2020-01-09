# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.Skimming.defaults.corrJetsProducer_cff import corrJetsProd
from SuSyAachen.DiLeptonHistograms.jecToUse_cfi import *

def corrJetsProducer_MCSignal16(process):
        process.seqcorrJetsProducer_MCSignal16 = corrJetsProd(process, "2016", False, fastSimUseDB["2016"], fastSimDB["2016"], special="AK8PFchs")
        process.seqcorrJetsPath_MCSignal = cms.Path(process.seqcorrJetsProducer_MCSignal16)

def corrJetsProducer_MCSignal17(process):
        process.seqcorrJetsProducer_MCSignal17 = corrJetsProd(process, "2017", False, fastSimUseDB["2017"], fastSimDB["2017"])
        process.seqcorrJetsPath_MCSignal = cms.Path(process.seqcorrJetsProducer_MCSignal17)

def corrJetsProducer_MCSignal18(process):
        process.seqcorrJetsProducer_MCSignal18 = corrJetsProd(process, "2018", False, fastSimUseDB["2018"], fastSimDB["2018"])
        process.seqcorrJetsPath_MCSignal = cms.Path(process.seqcorrJetsProducer_MCSignal18)
