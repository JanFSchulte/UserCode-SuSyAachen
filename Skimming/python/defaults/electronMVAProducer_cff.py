import FWCore.ParameterSet.Config as cms

def electronMVAProducer16(process):
        from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
        setupEgammaPostRecoSeq( process,
                                isMiniAOD=True,
                                applyEnergyCorrections=False,
                                applyVIDOnCorrectedEgamma=False,
                                era="2016-Legacy"
        )
        process.seqelectronMVAProducer16 = cms.Sequence(process.egammaPostRecoSeq)
        process.electronMVAPath = cms.Path(process.seqelectronMVAProducer16)
        
def electronMVAProducer17(process):
        from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
        setupEgammaPostRecoSeq( process,
                                era="2017-Nov17ReReco"
        )
        process.seqelectronMVAProducer17 = cms.Sequence(process.egammaPostRecoSeq)
        process.electronMVAPath = cms.Path(process.seqelectronMVAProducer17)

def electronMVAProducer18(process):
        from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
        setupEgammaPostRecoSeq( process,
                                isMiniAOD=True,
                                runEnergyCorrections=False,
                                era="2018-Prompt")
        process.seqelectronMVAProducer18 = cms.Sequence(process.egammaPostRecoSeq)
        process.electronMVAPath = cms.Path(process.seqelectronMVAProducer18)

