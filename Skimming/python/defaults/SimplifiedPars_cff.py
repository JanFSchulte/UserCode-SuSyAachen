import FWCore.ParameterSet.Config as cms
def SimplifiedPars(process):  
    process.load('SUSYBSMAnalysis.SusyCAF.SusyCAF_Simplified_cfi')
    process.seqSimplifiedPars = process.susycafscan
                         
