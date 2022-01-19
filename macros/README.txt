# For more recent Otto's models:

source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.csh

# pkl->txt conversion:
python pkl2txt_py3.py

# offline bdt values dump:
python checkFromPkl_py3.py

# to add a branch with the bdt value to Formatted trees:
python computeBdtOutput.py --isPFPF=0
or
python computeBdtOutput.py --isPFPF=1


# ---------------------------------------------------------------

## to merge and format (before adding the BDT)
## No need for python3
## From e.g. /tmp

cmsenv

python mergingFilesFromSkim.py
python formattingFilesFromSkim.py

# ---------------------------------------------------------------
# 
## to add the BDT output in batch
# 

cmsenv

voms-proxy-init --voms cms --valid 168:00
mv /tmp/x509up_u7291 ~/X509_USER_PROXY
setenv X509_USER_PROXY ~/X509_USER_PROXY
echo $X509_USER_PROXY 

python submit_batch_BDTout.py -c -N -1 -n 1 -p testtnp --eos=ok --ispfpf=1 --inputfile=myfiles.txt
