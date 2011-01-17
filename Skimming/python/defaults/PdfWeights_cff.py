import FWCore.ParameterSet.Config as cms
def PdfWeights(process):
    from SuSyAachen.Skimming.genSelection_cff import pdfWeights
    process.susyScanPdfWeights = pdfWeights.clone(
    )

    process.seqPdfWeights = cms.Sequence(process.susyScanPdfWeights)




