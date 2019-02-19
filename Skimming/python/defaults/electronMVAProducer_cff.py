import FWCore.ParameterSet.Config as cms

def electronMVAProducer17(process):
        #process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
        #process.seqelectronMVAProducer = cms.Sequence(process.electronMVAValueMapProducer)
        #process.electronMVAPath = cms.Path(process.seqelectronMVAProducer)
        from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
        setupEgammaPostRecoSeq( process,
                                isMiniAOD=True,
                                applyEnergyCorrections=False,
                                applyVIDOnCorrectedEgamma=False,
                                era="2017-Nov17ReReco"
        )
        process.seqelectronMVAProducer17 = cms.Sequence(process.egammaScaleSmearSeq*process.egammaPostRecoSeq)
        process.electronMVAPath = cms.Path(process.seqelectronMVAProducer17)



def electronMVAProducer18(process):
        #process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
        #process.seqelectronMVAProducer = cms.Sequence(process.electronMVAValueMapProducer)
        #process.electronMVAPath = cms.Path(process.seqelectronMVAProducer)
        from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
        setupEgammaPostRecoSeq( process,
                                isMiniAOD=True,
                                runEnergyCorrections=False,
                                applyEnergyCorrections=False,
                                applyVIDOnCorrectedEgamma=False,
                                era="2018-Prompt")
        process.seqelectronMVAProducer18 = cms.Sequence(process.egammaScaleSmearSeq*process.egammaPostRecoSeq)
        process.electronMVAPath = cms.Path(process.seqelectronMVAProducer18)

