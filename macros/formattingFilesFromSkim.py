#!/usr/bin/env python3.9 
import os, math, optparse, ROOT
from sklearn.externals import joblib
import glob
from shutil import copy2

# compile code for formatting
if "/ tnpForBTreeFormat_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/RKAnalysis/macros/tnpForBTreeFormat.cc++" % os.environ['CMSSW_BASE']);

if __name__ == "__main__":
  directory = '/eos/cms/store/user/crovelli/LowPtEle/TnpDataTest/'

  # Loop over merged files
  filename = directory+'Merged*.root'

  # copy files from EOS to this dir and act
  for file in glob.glob(filename):
    #print(file)
    endfile = file.replace(directory,"")
    copy2(file, endfile)
    print('file ',file,' copied')

    # format: PFPF
    print('Formatting for PFPF')
    ROOT.tnpForBTreeFormat(endfile, 1)

    # format: PF-LPT
    print('Formatting for PFLPT')
    ROOT.tnpForBTreeFormat(endfile, 0)

    # remove files before formatting
    os.system("rm {}".format(endfile))

    # copy formatted file to EOS
    tocp1 = "cp FormattedTnPForB_PFPF_{fullname} {directory}".format(fullname=endfile, directory=directory)
    os.system(tocp1)
    tocp2 = "cp FormattedTnPForB_PFLPT_{fullname} {directory}".format(fullname=endfile, directory=directory)
    os.system(tocp2)



