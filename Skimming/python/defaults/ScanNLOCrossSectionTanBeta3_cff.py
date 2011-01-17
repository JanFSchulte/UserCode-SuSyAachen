import FWCore.ParameterSet.Config as cms
def ScanNLOCrossSectionTanBeta3(process):  
    from SuSyAachen.Skimming.genSelection_cff import susyCrossSectionProducer
    from SuSyAachen.Skimming.scale_xsection_nlo10_tandat3_cfi import scale_xsection_nlo10_tandat3
    from SuSyAachen.Skimming.scale_xsection_nlo20_tandat3_cfi import scale_xsection_nlo20_tandat3
    from SuSyAachen.Skimming.scale_xsection_nlo05_tandat3_cfi import scale_xsection_nlo05_tandat3
    from SuSyAachen.Skimming.scale_kfactordat3_cfi import scale_kfactordat3
    process.susyScanNLOCrossSection = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo10_tandat3
                )
    process.susyScanNLOCrossSectionScale2 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo20_tandat3
                )
    process.susyScanNLOCrossSectionScale05 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo05_tandat3
                )
    process.susyScankFactor = susyCrossSectionProducer.clone(
                crossSections = scale_kfactordat3
                )

    process.seqScanNLOCrossSectionTanBeta3 = cms.Sequence(process.susyScanNLOCrossSection+process.susyScanNLOCrossSectionScale2+process.susyScanNLOCrossSectionScale05+process.susyScankFactor)
                         
