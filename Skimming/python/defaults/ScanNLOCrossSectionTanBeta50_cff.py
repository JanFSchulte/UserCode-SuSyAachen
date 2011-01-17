import FWCore.ParameterSet.Config as cms
def ScanNLOCrossSectionTanBeta50(process):  
    from SuSyAachen.Skimming.genSelection_cff import susyCrossSectionProducer
    from SuSyAachen.Skimming.scale_xsection_nlo10_tandat50_cfi import scale_xsection_nlo10_tandat50
    from SuSyAachen.Skimming.scale_xsection_nlo20_tandat50_cfi import scale_xsection_nlo20_tandat50
    from SuSyAachen.Skimming.scale_xsection_nlo05_tandat50_cfi import scale_xsection_nlo05_tandat50
    from SuSyAachen.Skimming.scale_kfactordat50_cfi import scale_kfactordat50
    process.susyScanNLOCrossSection = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo10_tandat50
                )
    process.susyScanNLOCrossSectionScale2 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo20_tandat50
                )
    process.susyScanNLOCrossSectionScale05 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo05_tandat50
                )
    process.susyScankFactor = susyCrossSectionProducer.clone(
                crossSections = scale_kfactordat50
                )

    process.seqScanNLOCrossSectionTanBeta50 = cms.Sequence(process.susyScanNLOCrossSection+process.susyScanNLOCrossSectionScale2+process.susyScanNLOCrossSectionScale05+process.susyScankFactor)
                         
