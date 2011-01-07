import FWCore.ParameterSet.Config as cms
def ScanNLOCrossSectionTanBeta10(process):  
    from SuSyAachen.Skimming.genSelection_cff import susyCrossSectionProducer
    from SuSyAachen.Skimming.scale_xsection_nlo10_tanssdat10_cfi import scale_xsection_nlo10_tanssdat10
    from SuSyAachen.Skimming.scale_xsection_nlo20_tanssdat10_cfi import scale_xsection_nlo20_tanssdat10
    from SuSyAachen.Skimming.scale_xsection_nlo05_tanssdat10_cfi import scale_xsection_nlo05_tanssdat10
    process.susyScanNLOCrossSection = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo10_tanssdat10
                )
    process.susyScanNLOCrossSectionScale2 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo20_tanssdat10
                )
    process.susyScanNLOCrossSectionScale05 = susyCrossSectionProducer.clone(
                crossSections = scale_xsection_nlo05_tanssdat10
                )

    process.seqScanNLOCrossSectionTanBeta10 = cms.Sequence(process.susyScanNLOCrossSection+process.susyScanNLOCrossSectionScale2+process.susyScanNLOCrossSectionScale05)
                         
