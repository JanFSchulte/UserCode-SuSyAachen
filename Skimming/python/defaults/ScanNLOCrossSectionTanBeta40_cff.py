import FWCore.ParameterSet.Config as cms
def ScanNLOCrossSectionTanBeta40(process):  
    from SuSyAachen.Skimming.genSelection_cff import susyCrossSectionProducer
    from SuSyAachen.Skimming.scale_xsection_nlo10_tandat40_cfi import scale_xsection_nlo10_tandat40
    from SuSyAachen.Skimming.scale_xsection_nlo20_tandat40_cfi import scale_xsection_nlo20_tandat40
    from SuSyAachen.Skimming.scale_xsection_nlo05_tandat40_cfi import scale_xsection_nlo05_tandat40
    from SuSyAachen.Skimming.scale_kfactordat40_cfi import scale_kfactordat40
    process.susyScanNLOCrossSection = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo10_tandat40
                )
    process.susyScanNLOCrossSectionScale2 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo20_tandat40
                )
    process.susyScanNLOCrossSectionScale05 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo05_tandat40
                )
    process.susyScankFactor = susyCrossSectionProducer.clone(
                crossSections = scale_kfactordat40
                )

    process.seqScanNLOCrossSectionTanBeta40 = cms.Sequence(process.susyScanNLOCrossSection+process.susyScanNLOCrossSectionScale2+process.susyScanNLOCrossSectionScale05+process.susyScankFactor)
                         
