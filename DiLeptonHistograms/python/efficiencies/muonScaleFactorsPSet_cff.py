#created at Thu Sep 22 16:05:38 2011
#created from  '['/Users/niklas/SUSY/Histos/pfTnP/MC.pfTnP.TTJets_madgraph_Fall10.root']'  

if not "cms" in globals(): import FWCore.ParameterSet.Config as cms


muonTrackingScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(0.9869),
			etaMin = cms.double(-2.4),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-2.1),
		),
		cms.PSet(
			weight = cms.double(0.9948),
			etaMin = cms.double(-2.1),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-1.6),
		),
		cms.PSet(
			weight = cms.double(0.9967),
			etaMin = cms.double(-1.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-1.2),
		),
		cms.PSet(
			weight = cms.double(0.9974),
			etaMin = cms.double(-1.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.9),
		),	
		cms.PSet(
			weight = cms.double(0.9980),
			etaMin = cms.double(-0.9),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.6),
		),	
		cms.PSet(
			weight = cms.double(0.9980),
			etaMin = cms.double(-0.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.3),
		),
		cms.PSet(
			weight = cms.double(0.9972),
			etaMin = cms.double(-0.3),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.2),
		),	
		cms.PSet(
			weight = cms.double(0.9963),
			etaMin = cms.double(-0.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
				etaMax = cms.double(0.2),
			),	
		cms.PSet(
			weight = cms.double(0.9978),
			etaMin = cms.double(0.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.3),
		),    
		cms.PSet(
			weight = cms.double(0.9977),
			etaMin = cms.double(0.3),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.6),
		),	
		cms.PSet(
			weight = cms.double(0.9976),
			etaMin = cms.double(0.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.9),
		),
		cms.PSet(
			weight = cms.double(0.9968),
			etaMin = cms.double(0.9),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.2),
		),
		cms.PSet(
			weight = cms.double(0.9959),
			etaMin = cms.double(1.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.6),
		),
		cms.PSet(
			weight = cms.double(0.9970),
			etaMin = cms.double(1.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.1),
		),
		cms.PSet(
			weight = cms.double(0.9836),
			etaMin = cms.double(2.1),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.4),
		),				    
)	

muonLowerTrackingScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(0.9862),
			etaMin = cms.double(-2.4),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-2.1),
		),
		cms.PSet(
			weight = cms.double(0.9946),
			etaMin = cms.double(-2.1),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-1.6),
		),
		cms.PSet(
			weight = cms.double(0.9965),
			etaMin = cms.double(-1.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-1.2),
		),					
		 cms.PSet(
			weight = cms.double(0.9972),
			etaMin = cms.double(-1.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.9),
		),	
		cms.PSet(
			weight = cms.double(0.9979),
			etaMin = cms.double(-0.9),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.6),
		),	
		cms.PSet(
			weight = cms.double(0.9979),
			etaMin = cms.double(-0.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.3),
		),
		cms.PSet(
			weight = cms.double(0.9970),
			etaMin = cms.double(-0.3),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.2),
		),	
		cms.PSet(
			weight = cms.double(0.9962),
			etaMin = cms.double(-0.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.2),
		),	
		cms.PSet(
			weight = cms.double(0.9976),
			etaMin = cms.double(0.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.3),
		),    
		cms.PSet(
			weight = cms.double(0.9976),
			etaMin = cms.double(0.3),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.6),
		),	
		cms.PSet(
			weight = cms.double(0.9975),
			etaMin = cms.double(0.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.9),
		),
		cms.PSet(
			weight = cms.double(0.9966),
			etaMin = cms.double(0.9),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.2),
		),
		cms.PSet(
			weight = cms.double(0.9956),
			etaMin = cms.double(1.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.6),
		),
		cms.PSet(
			weight = cms.double(0.9968),
			etaMin = cms.double(1.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.1),
		),
		cms.PSet(
			weight = cms.double(0.9828),
			etaMin = cms.double(2.1),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.4),
		),												    
							
)	


muonUpperTrackingScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(0.9876),
			etaMin = cms.double(-2.4),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-2.1),
		),
		cms.PSet(
			weight = cms.double(0.9950),
			etaMin = cms.double(-2.1),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-1.6),
		),
		cms.PSet(
			weight = cms.double(0.9969),
			etaMin = cms.double(-1.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-1.2),
		),					
		cms.PSet(
			weight = cms.double(0.9976),
			etaMin = cms.double(-1.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.9),
		),	
		cms.PSet(
			weight = cms.double(0.9981),
			etaMin = cms.double(-0.9),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.6),
		),	
		cms.PSet(
			weight = cms.double(0.9981),
			etaMin = cms.double(-0.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.3),
		),
		cms.PSet(
			weight = cms.double(0.9974),
			etaMin = cms.double(-0.3),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(-0.2),
		),	
		cms.PSet(
			weight = cms.double(0.9964),
			etaMin = cms.double(-0.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.2),
		),	
		cms.PSet(
			weight = cms.double(0.9980),
			etaMin = cms.double(0.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.3),
		),    
		cms.PSet(
			weight = cms.double(0.9978),
			etaMin = cms.double(0.3),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.6),
		),	
		cms.PSet(
			weight = cms.double(0.9977),
			etaMin = cms.double(0.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.9),
		),
		cms.PSet(
			weight = cms.double(0.9970),
			etaMin = cms.double(0.9),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.2),
		),
		cms.PSet(
			weight = cms.double(0.9962),
			etaMin = cms.double(1.2),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.6),
		),
		cms.PSet(
			weight = cms.double(0.9972),
			etaMin = cms.double(1.6),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.1),
		),
		cms.PSet(
			weight = cms.double(0.9844),
			etaMin = cms.double(2.1),
			ptMin = cms.double(0.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.4),
		),												    
							
)

muonIDScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(0.9939),
			etaMin = cms.double(0.000000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.900000),
		),
		cms.PSet(
			weight = cms.double(0.9902),
			etaMin = cms.double(0.900000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.200000),
		),
		cms.PSet(
			weight = cms.double(0.9970),
			etaMin = cms.double(1.200000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.400000),
		),
			
		
 
)

muonLowerIDScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(0.9889),
			etaMin = cms.double(0.000000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.900000),
		),
		cms.PSet(
			weight = cms.double(0.9852),
			etaMin = cms.double(0.900000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.200000),
		),
		cms.PSet(
			weight = cms.double(0.9920),
			etaMin = cms.double(1.200000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.400000),
		),		
		
 
)
muonUpperIDScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(0.9989),
			etaMin = cms.double(0.000000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.900000),
		),
		cms.PSet(
			weight = cms.double(0.9952),
			etaMin = cms.double(0.900000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.200000),
		),
		cms.PSet(
			weight = cms.double(1.0020),
			etaMin = cms.double(1.200000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.400000),
		),		
		
 
)

muonIsolationScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(1.0004),
			etaMin = cms.double(0.000000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.900000),
		),
		cms.PSet(
			weight = cms.double(1.0031),
			etaMin = cms.double(0.900000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.200000),
		),
		cms.PSet(
			weight = cms.double(1.0050),
			etaMin = cms.double(1.200000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.400000),
		),		
		
 
)

muonLowerIsolationScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(0.9984),
			etaMin = cms.double(0.000000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.900000),
		),
		cms.PSet(
			weight = cms.double(1.0011),
			etaMin = cms.double(0.900000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.200000),
		),
		cms.PSet(
			weight = cms.double(1.0030),
			etaMin = cms.double(1.200000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.400000),
		),		

)

muonUpperIsolationScaleFactors =  cms.VPSet(
		cms.PSet(
			weight = cms.double(1.0024),
			etaMin = cms.double(0.000000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(0.900000),
		),
		cms.PSet(
			weight = cms.double(1.0051),
			etaMin = cms.double(0.900000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(1.200000),
		),
		cms.PSet(
			weight = cms.double(1.0070),
			etaMin = cms.double(1.200000),
			ptMin = cms.double(20.000000),
			ptMax = cms.double(1000.000000),
			etaMax = cms.double(2.400000),
		),		
	
 
)
