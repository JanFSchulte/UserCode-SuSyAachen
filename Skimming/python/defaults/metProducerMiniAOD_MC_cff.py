# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.Skimming.metProducerMiniAOD_cff import metProducerMiniAOD
from SuSyAachen.DileptonHistograms.jecToUse_cfi import *


def metProducerMiniAOD_MC16(process):
        process.seqmetProducerMiniAOD_MC16 = metProducerMiniAOD(process, "2016", runOnData=False, usePrivateSQlite=mcUseDB["2016"], era=mcDB["2016"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MC16)

def metProducerMiniAOD_MC17(process):
        process.seqmetProducerMiniAOD_MC17 = metProducerMiniAOD(process, "2017", runOnData=False, usePrivateSQlite=mcUseDB["2017"], era=mcDB["2017"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MC17)

def metProducerMiniAOD_MC18(process):
        process.seqmetProducerMiniAOD_MC18 = metProducerMiniAOD(process, "2018", runOnData=False, usePrivateSQlite=mcUseDB["2018"], era=mcDB["2018"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MC18)
