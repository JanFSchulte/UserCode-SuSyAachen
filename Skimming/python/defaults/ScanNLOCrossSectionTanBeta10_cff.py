import FWCore.ParameterSet.Config as cms
def ScanNLOCrossSectionTanBeta10(process):  
    from SuSyAachen.Skimming.genSelection_cff import susyCrossSectionProducer
    from SuSyAachen.Skimming.scale_xsection_nlo10_tandat10_cfi import scale_xsection_nlo10_tandat10
    from SuSyAachen.Skimming.scale_xsection_nlo20_tandat10_cfi import scale_xsection_nlo20_tandat10
    from SuSyAachen.Skimming.scale_xsection_nlo05_tandat10_cfi import scale_xsection_nlo05_tandat10
    from SuSyAachen.Skimming.scale_kfactordat10_cfi import scale_kfactordat10
    process.susyScanNLOCrossSection = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo10_tandat10
                )
    process.susyScanNLOCrossSectionScale2 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo20_tandat10
                )
    process.susyScanNLOCrossSectionScale05 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo05_tandat10
                )
    process.susyScankFactor = susyCrossSectionProducer.clone(
                crossSections = scale_kfactordat10
                )

    process.seqScanNLOCrossSectionTanBeta10 = cms.Sequence(process.susyScanNLOCrossSection+process.susyScanNLOCrossSectionScale2+process.susyScanNLOCrossSectionScale05+process.susyScankFactor)
                         
