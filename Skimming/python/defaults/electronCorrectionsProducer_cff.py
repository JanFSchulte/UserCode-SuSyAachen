import FWCore.ParameterSet.Config as cms

def electronCorrectionsProducer(process):
        from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
        setupEgammaPostRecoSeq(process,
                               runVID=False, #saves CPU time by not needlessly re-running VID
                               era='2017-Nov17ReReco')  

        
        process.seqelectronCorrectionsPath = cms.Path(process.seqseqelectronCorrections)
