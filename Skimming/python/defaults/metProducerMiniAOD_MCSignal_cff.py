# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.Skimming.defaults.metProducerMiniAOD_cff import metProducerMiniAOD
from SuSyAachen.DiLeptonHistograms.jecToUse_cfi import *


def metProducerMiniAOD_MCSignal16(process):
        process.seqmetProducerMiniAOD_MCSignal16 = metProducerMiniAOD(process, "2016", runOnData=False, usePrivateSQlite=fastSimUseDB["2016"], era=fastSimDB["2016"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MCSignal16)

def metProducerMiniAOD_MCSignal17(process):
        process.seqmetProducerMiniAOD_MCSignal17 = metProducerMiniAOD(process, "2017", runOnData=False, usePrivateSQlite=fastSimUseDB["2017"], era=fastSimDB["2017"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MCSignal17)

def metProducerMiniAOD_MCSignal18(process):
        process.seqmetProducerMiniAOD_MCSignal18 = metProducerMiniAOD(process, "2018", runOnData=False, usePrivateSQlite=fastSimUseDB["2018"], era=fastSimDB["2018"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MCSignal18)

