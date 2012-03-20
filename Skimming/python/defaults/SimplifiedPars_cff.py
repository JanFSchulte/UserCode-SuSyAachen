import FWCore.ParameterSet.Config as cms
def SimplifiedPars(process):  
    process.load('SUSYBSMAnalysis.SusyCAF.SusyCAF_Simplified_cfi')

    from SuSyAachen.Histograms.scanCounter_cfi import scanCounter
    from SuSyAachen.DiLeptonHistograms.DiLeptonTrees_cfi import simplifedVars
    process.scanCounter = scanCounter.clone(
        susyVars = simplifedVars,
        )
    process.seqSimplifiedPars = cms.Sequence(process.susycafscan * process.scanCounter)


                         
