# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.Skimming.corrJetsProducer_cff import corrJetsProd
from SuSyAachen.DiLeptonHistograms.jecToUse_cfi import *

def corrJetsProducer_MC16(process):
        process.seqcorrJetsProducer_MC16 = corrJetsProd(process, "2016", False, mcUseDB["2016"], mcDB["2016"])
        process.seqcorrJetsPath_MC = cms.Path(process.seqcorrJetsProducer_MC16)

def corrJetsProducer_MC17(process):
        process.seqcorrJetsProducer_MC17 = corrJetsProd(process, "2017", False, mcUseDB["2017"], mcDB["2017"])
        process.seqcorrJetsPath_MC = cms.Path(process.seqcorrJetsProducer_MC17)

def corrJetsProducer_MC18(process):
        process.seqcorrJetsProducer_MC18 = corrJetsProd(process, "2018", False, mcUseDB["2018"], mcDB["2018"])
        process.seqcorrJetsPath_MC = cms.Path(process.seqcorrJetsProducer_MC18)
