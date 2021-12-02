### source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
import os, math, optparse, ROOT
from sklearn.externals import joblib
import glob
from shutil import copy2

if __name__ == "__main__":
  directory = '/eos/cms/store/user/crovelli/LowPtEle/TnpDataTest/'

  # Loop over files ending with 0.root, 1.root... 9.root 
  fullnamestart = ''
  namestart = ''
  for theend in range(10):

    filenameA = '*{}.root'.format(theend)
    filename = directory+filenameA
    #print(filename) 

    # copy files from EOS to this dir
    for file in glob.glob(filename):
      #print(file)
      endfile = file.replace(directory,"")
      copy2(file, endfile)

      # save the beginning of the file names, to be used later
      if namestart=='':
        namestart = endfile[0:5]

      if fullnamestart=='':
        extfileend = '_{}.root'.format(theend)    
        if endfile[-7:]==extfileend:
          fullnamestartA = endfile
          fullnamestart = fullnamestartA.replace(extfileend,'')

    # print ('fullnamestart=',fullnamestart)

    # prepare the string to add all the files which have been copied
    filesA = [x.strip() for x in glob.glob('*.root')]
    #print('files before cleaning = ', filesA)

    # clean up the list and keep only the wanted files (i.e. copied and ending with 0.root, 1.root...)
    fileend = '{}.root'.format(theend)   
    files = [x for x in filesA if x[0:5]==namestart and x[-6:]==fileend ]
    #print('files after cleaning = ', files)   

    toadd1 = "hadd Merged_{fullname}_{number}.root {myfiles}".format(fullname=fullnamestart, number=theend, myfiles=files)
    toadd2 = toadd1.replace("[", "")
    toadd3 = toadd2.replace("]", "")
    toadd4 = toadd3.replace("'", "")
    toadd = toadd4.replace(",", "")
    print(toadd)

    # now add
    os.system(toadd)

    # copy back to EOS
    tocp = "cp Merged_{fullname}_{number}.root {directory}".format(fullname=fullnamestart, number=theend, directory=directory)
    os.system(tocp)

    # finally remove the files before merging
    # print('now cancel')
    for filetc in glob.glob('{fullname}*{number}.root'.format(fullname=fullnamestart,  number=theend)):
      # print(filetc)
      os.system("rm {}".format(filetc))
