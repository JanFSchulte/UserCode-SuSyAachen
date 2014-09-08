#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
ANA_JSON=./FullJSON2012.txt
BlockA_JSON = ./BlockAJSON.txt
BlockB_JSON = ./BlockBJSON.txt
#/home/home4/institut_1b/edelhoff/superSymmetry/SubmitScripts/Input/Cert_190456-195016_8TeV_PromptReco_Collisions12_JSON.txt
#/home/home4/institut_1b/edelhoff/superSymmetry/SubmitScripts/Input/Cert_190456-193557_8TeV_PromptReco_Collisions12_JSON.txt
#/home/home4/institut_1b/edelhoff/superSymmetry/SubmitScripts/Input/Cert_190456-191276_8TeV_PromptReco_Collisions12_JSON.txt
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr.txt
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-195165.txt
#5fbi? /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-196531_patch2.txt
#1p6fbi/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-193575.txt
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-191276.txt
PU_CALC=pileupCalc.py --calcMode true  --maxPileupBin 60 --numPileupBins 60 
MAKE_HIST=python makeMChisto.py

mc_PU_3D_dist.root:
	ls

mc_PU_dist_S10.root:
	$(MAKE_HIST) $@

PU_dist.root:
	$(PU_CALC) --minBiasXsec 69400 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) PU_dist_Full2012.root
	$(PU_CALC) --minBiasXsec 72870 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) PU_dist_Full2012_Up.root
	$(PU_CALC) --minBiasXsec 65930 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) PU_dist_Full2012_Down.root

PU_dist_BlockA.root:
	$(PU_CALC) --minBiasXsec 69400 -i $(BlockA_JSON) --inputLumiJSON $(PU_JSON) PU_dist_BlockA.root
	$(PU_CALC) --minBiasXsec 72870 -i $(BlockA_JSON) --inputLumiJSON $(PU_JSON) PU_dist_BlockA_Up.root
	$(PU_CALC) --minBiasXsec 65930 -i $(BlockA_JSON) --inputLumiJSON $(PU_JSON) PU_dist_BlockA_Down.root

PU_dist_BlockB.root:
	$(PU_CALC) --minBiasXsec 69400 -i $(BlockB_JSON) --inputLumiJSON $(PU_JSON) PU_dist_BlockB.root
	$(PU_CALC) --minBiasXsec 72870 -i $(BlockB_JSON) --inputLumiJSON $(PU_JSON) PU_dist_BlockB_Up.root
	$(PU_CALC) --minBiasXsec 65930 -i $(BlockB_JSON) --inputLumiJSON $(PU_JSON) PU_dist_BlockB_Down.root


PU_3D_dist.root:
	ls
