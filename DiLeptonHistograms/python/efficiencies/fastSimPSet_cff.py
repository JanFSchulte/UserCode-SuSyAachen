#created at Thu Sep 22 15:59:28 2011
#created from  '['/Users/niklas/SUSY/Histos/pfTnP/MC.pfTnP.TTJets_madgraph_Fall10.root']'  

if not "cms" in globals(): import FWCore.ParameterSet.Config as cms


eeScaleFactors =  cms.VPSet(
                 cms.PSet(
                        weight = cms.double(0.970959),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.006370),
                ),

 
                cms.PSet(
                        weight = cms.double(0.924484),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.014489),
                ),

 
                cms.PSet(
                        weight = cms.double(0.893397),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.013062),
                ),

 
                cms.PSet(
                        weight = cms.double(0.799996),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.025873),
                ),

 
                cms.PSet(
                        weight = cms.double(0.980203),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.009002),
                ),

 
                cms.PSet(
                        weight = cms.double(0.910504),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.021174),
                ),

 
                cms.PSet(
                        weight = cms.double(0.890994),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.018817),
                ),

 
                cms.PSet(
                        weight = cms.double(0.771308),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.039578),
                ),

 
                cms.PSet(
                        weight = cms.double(0.998869),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.013118),
                ),

 
                cms.PSet(
                        weight = cms.double(0.893160),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.031124),
                ),

 
                cms.PSet(
                        weight = cms.double(0.935214),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.028894),
                ),

 
                cms.PSet(
                        weight = cms.double(0.930187),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.072956),
                ),

)

emScaleFactors = cms.VPSet(

                 cms.PSet(
                        weight = cms.double(0.986340),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.004384),
                ),

 
                cms.PSet(
                        weight = cms.double(1.006434),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.009087),
                ),

 
                cms.PSet(
                        weight = cms.double(0.911217),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.009332),
                ),

 
                cms.PSet(
                        weight = cms.double(0.942950),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.018109),
                ),

 
                cms.PSet(
                        weight = cms.double(1.001630),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.006233),
                ),

 
                cms.PSet(
                        weight = cms.double(0.997904),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.013088),
                ),

 
                cms.PSet(
                        weight = cms.double(0.936552),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.013889),
                ),

 
                cms.PSet(
                        weight = cms.double(0.944355),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.028444),
                ),

 
                cms.PSet(
                        weight = cms.double(1.024684),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.009153),
                ),

 
                cms.PSet(
                        weight = cms.double(1.025384),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.019904),
                ),

 
                cms.PSet(
                        weight = cms.double(0.922190),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.020659),
                ),

 
                cms.PSet(
                        weight = cms.double(0.944144),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.043338),
                ),




)


mmScaleFactors = cms.VPSet(
                cms.PSet(
                        weight = cms.double(0.998287),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.005981),
                ),

 
                cms.PSet(
                        weight = cms.double(1.014286),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.012968),
                ),

 
                cms.PSet(
                        weight = cms.double(1.004746),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.011817),
                ),

 
                cms.PSet(
                        weight = cms.double(1.042044),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(2.000000),
                        nJetsMax = cms.double(3.000000),
                        uncert = cms.double(0.024128),
                ),

 
                cms.PSet(
                        weight = cms.double(1.010975),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.008435),
                ),

 
                cms.PSet(
                        weight = cms.double(1.042208),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.019318),
                ),

 
                cms.PSet(
                        weight = cms.double(1.017324),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.016922),
                ),

 
                cms.PSet(
                        weight = cms.double(0.943732),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(3.000000),
                        nJetsMax = cms.double(4.000000),
                        uncert = cms.double(0.034000),
                ),

 
                cms.PSet(
                        weight = cms.double(1.036016),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.012460),
                ),

 
                cms.PSet(
                        weight = cms.double(1.076845),
                        eta1Min = cms.double(0.000000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(1.400000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.030477),
                ),

 
                cms.PSet(
                        weight = cms.double(1.055968),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(0.000000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(1.400000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.026010),
                ),

 
                cms.PSet(
                        weight = cms.double(1.019435),
                        eta1Min = cms.double(1.400000),
                        pt1Min = cms.double(20.000000),
                        pt1Max = cms.double(1000.000000),
                        eta1Max = cms.double(2.500000),
                        eta2Min = cms.double(1.400000),
                        pt2Min = cms.double(20.000000),
                        pt2Max = cms.double(1000.000000),
                        eta2Max = cms.double(2.500000),
                        nJetsMin = cms.double(4.000000),
                        nJetsMax = cms.double(100.000000),
                        uncert = cms.double(0.056865),
                ),



)