#!/usr/bin/env python3.9
import os, math, optparse, ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from array import array
from root_numpy import hist2array
import numpy as np

import joblib
#from sklearn.externals import joblib
from xgboost import XGBClassifier
import glob
import xgboost as xgb 
from shutil import copy2

ROOT.gStyle.SetOptStat(111111)
ROOT.gROOT.SetBatch(True)

def computeBdtBatch(isPFPF, ijob, therootfile, theeosdir):

    # -------------------------------------------------
    # XGBoost model in .pkl
    if (isPFPF==1): 
        fmodel = '/afs/cern.ch/user/c/crovelli/public/4RK/xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_0.pkl'
        print('dumping BDT output for PFPF reading from...'.format(fmodel))

    if (isPFPF==0): 
        fmodel = '/afs/cern.ch/user/c/crovelli/public/4RK/xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.2_0.pkl'
        print('dumping BDT output for PFLPT reading from...'.format(fmodel))

    # --------------------------------------------------
    # Load the BDT model 'weights'
    bdt = joblib.load(fmodel)
    print(bdt)

    booster = bdt._Booster
    print(booster)

    # --------------------------------------------------
    # Copy the input file from EOS to here
    startfile = theeosdir+therootfile
    endfile = therootfile
    print('startfile = ', startfile)
    print('endfile = ', endfile)
    copy2(startfile, endfile)

    # Read the original root tree (local copy)
    tf = ROOT.TFile.Open(endfile, "update")
    tftree = tf.Get("tnpAna/fitter_tree")
    entries = tftree.GetEntries();
    print("entries = ", entries)


    # ---------------------------------------------------
    # Prepare an extra branch with the analysis BDT output - computed event by event
    xgb_arr = array("f",[0])
    branch = tftree.Branch("xgb",xgb_arr,"xgb/F")    


    # ---------------------------------------------------  
    ## now loop over the original tree entries and fill the extra branch
    for ie, event in enumerate(tftree):

        if ie%500==1: print ("Processing event ie = ",ie)

        # inputs to BDT
        f0  = tftree.B_svprob
        f1  = tftree.B_xysig
        f2  = tftree.B_cos2d
        f3  = tftree.tagPt/tftree.B_mass
        f4  = tftree.probePt/tftree.B_mass
        f5  = tftree.kPt/tftree.B_mass
        f6  = tftree.theLKdz
        f7  = tftree.theL1L2dr_raw
        f8  = tftree.theLKdr_raw
        f9  = tftree.tagId 
        f10 = tftree.probeId
        f11 = tftree.probeIso04Rel
        f12 = tftree.kIso04Rel
        f13 = tftree.theBBDPhi
        f14 = tftree.theBTrkdxy2
        f15 = tftree.thePtAsym
        f16 = tftree.B_i3dsig

        # format the BDT inputs
        if (isPFPF==1): thefile = os.environ['PWD']+'/dumpThisEntryPFPF_'+ijob+'.txt'
        if (isPFPF==0): thefile = os.environ['PWD']+'/dumpThisEntryPFLPT_'+ijob+'.txt'
        with open(thefile, 'w') as f:
            line = "1 0:{} 1:{} 2:{} 3:{} 4:{} 5:{} 6:{} 7:{} 8:{} 9:{} 10:{} 11:{} 12:{} 13:{} 14:{} 15:{} 16:{}".format(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16) 
            f.write(line)
            f.write('\n')
        f.close()    
        dtest = xgb.DMatrix(thefile, silent =1)

        # compute the BDT score   
        ypred = booster.predict(dtest)

        # Store the BDT score in the tree as an extra branch
        xgb_arr[0] = float(ypred)
        branch.Fill()
        

    # --------------------------------------------------
    # Store the updated tree    
    tf.Write("", ROOT.TFile.kOverwrite)
    tf.Close()


if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')    
    parser.add_option('', '--isPFPF' , type='string' , default='1' , help='isPFPFP (options = 1 (default) or 0)')
    parser.add_option('', '--ijob' , type='string' , default='0' , help='which job')

    parser.add_option('', '--rootfile' , type='string' , default='pippo.root' , help='root file to be analyzed (wo dir)')
    parser.add_option('', '--eosdir' , type='string' , default='/eos/cms/store/user/crovelli/LowPtEle/TnpDataTest/')
    (options, args) = parser.parse_args()

    import os, sys

    if options.isPFPF in ['0']: computeBdtBatch(0, options.ijob, options.rootfile, options.eosdir)
    if options.isPFPF in ['1']: computeBdtBatch(1, options.ijob, options.rootfile, options.eosdir)
