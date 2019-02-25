#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
# JSON: 36.459/fb https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2766.html
#ANA_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt 
ANA_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt 
#PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt 
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt 
NORMTAG_JSON=/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json



PU_CALC=pileupCalc.py --calcMode true --maxPileupBin 100 --numPileupBins 100 -i ana_json.txt --inputLumiJSON pu_json.txt
BRIL_CALC=brilcalc lumi --normtag $(NORMTAG_JSON) -u /pb -i $(ANA_JSON)

all: data/puWeights.root

data:
	mkdir -p $@

mcHist: makeMChisto.py
	python2.7 $^ mc_PU_dist_Fall2017.root

# For the cross section and its uncertainty, see https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
dataHist: Makefile
#	scp mteroerd@lxplus.cern.ch:$(ANA_JSON) ana_json.txt
#	scp mteroerd@lxplus.cern.ch:$(PU_JSON) pu_json.txt
	$(PU_CALC) --minBiasXsec 69200 --pileupHistName pileup PU_dist_2018.root
	$(PU_CALC) --minBiasXsec 72383 --pileupHistName pileupUp PU_dist_2018_Up.root
	$(PU_CALC) --minBiasXsec 66017 --pileupHistName pileupDown PU_dist_2018_Down.root

# need to do this where brilws is installed, for now it only works on lxplus: https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html
# alias brilcalc in .bashrc, source .bashrc in .bash_login and add 'shopt -s expand_aliases' to .bashrc so that the brilcalc command is usable when it is executed in ssh
lumi:
	ssh mteroerd@lxplus.cern.ch $(BRIL_CALC) > lumi_summary.txt
