#!/usr/bin/python
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Submit some dang jobs')
parser.add_argument('-v','--verbose', dest='verbose', action='store_true',
		    help='verbosity')
parser.add_argument('--submit', dest='submit', action='store_true',
		    help='submit')
parser.add_argument('--metRate', dest='metRate', action='store_true',
                    help='metRate')
parser.add_argument('--sample', dest='sample', default='dy_140',
                    help='data sample to run on')
parser.add_argument('--tkptcut', dest='tkptcut', default='2.0',
                    help='data sample to run on')

args = parser.parse_args()

def main():

	print args.sample;
	if args.verbose: print "verbose!"

	tplFileIn  = 'runNtupleProducer_cfg_tmp.py'
	tplFileOut = 'runNtupleProducer_cfg_tmpO.py'

	cursample = args.sample

	samples = {};
	samples['dy_140'] = '/store/relval/CMSSW_9_1_0_pre1/RelValZMM_14/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v9_D12PU140-v1/00000';
	samples['tt_140'] = '/store/relval/CMSSW_9_1_0_pre1/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v9_D11PU140-v1/00000';
	samples['qcdmu_140'] = '/store/relval/CMSSW_9_1_0_pre1/RelValQCD_Pt-20toInf_MuEnrichedPt15_14TeV/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v9_D12PU140-v1/00000'
	samples['mugun_140'] = '/store/relval/CMSSW_9_1_0_pre1/RelValSingleMuPt1Extended/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v9_D12PU140-v1/00000'
	oput = os.listdir('/eos/cms'+samples[cursample]);
	# print oput

	sMetRate = "False"
	if args.metRate: sMetRate = "True"
	sTkPtCut = args.tkptcut;

	fIn = open(tplFileIn, "r") ;
	fOut = open(tplFileOut, "wr");
	for line in fIn: 
		newline  = line.replace('YYYY',sTkPtCut)
		newline2 = newline.replace('ZZZZ',sMetRate)
		fOut.write(newline2[:-1] + "\n");
		print newline2[:-1]

	for i,n in enumerate(oput):
		label = cursample+"_"+sTkPtCut
		if args.metRate: label += "_metrate"
		command = "bsub -q cmscaf1nh -o out.%%J run.sh %s %s $PWD %s %s" % (samples[cursample],n,label,tplFileOut)
		print args.submit, args.sample, label, command
		if args.submit: os.system(command);

##---------------------------------------
if __name__ == '__main__':
		main();
##---------------------------------------  
