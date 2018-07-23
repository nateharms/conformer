#!/usr/bin/env python

# This example use python to count number of confomers each system has 
# below a threshold.
#

import os

aalist = ["ala","arg","argH","asn","asp","aspH","cys","gln","glu","gluH",
          "gly","hisD","hisE","hisH","ile","leu","lys","lysH","met","phe",
          "pro","ser","thr","trp","tyr","val"]

caps = ["uncapped","dipeptide"]

thrs = 0.020 # 20 meV threshold

for aa in aalist:
    for cap in caps:
        # get list of ions
        ipath = os.path.join(aa,cap)
        ions = []
        for x in os.listdir(ipath):
            if os.path.isdir(os.path.join(ipath,x)):
                ions.append(x)

        # Go though each hierarchy file and find number of conformers below threshold
        for ion in ions:
            hfname = os.path.join(ipath,ion,"hier_rel_PDB+vdW.dat")
            hfile = open(hfname)
            
            count = 0
            for line in hfile:
                if not '#' == line[0]:
                    n,e = line.split()
                    n = int(n)
                    e = float(e)
                    if e > thrs:
                        break
            hfile.close()
            # print system information and number of conformers under threshold
            print("%s %s %s %i" %(aa,cap,ion, n-1))

