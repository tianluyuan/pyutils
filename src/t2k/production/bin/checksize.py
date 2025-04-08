#!/usr/bin/env python

import sys
import ROOT
import subprocess

def getMean(sizes):
    return sum(sizes)/len(sizes)

filelist = sys.argv[1:]
cmd = 'ls -l '
sizes = []
for file in filelist:
    fileLS = subprocess.Popen(cmd+file, shell=True, stderr=subprocess.PIPE,
    stdout=subprocess.PIPE).communicate()[0]
    sizes.append(int(fileLS.split()[4]))

print getMean(sizes)
print max(sizes)
print min(sizes)

h1 = ROOT.TH1D('h1', '', 30, min(sizes), max(sizes))
for size in sizes:
    h1.Fill(size)

h1.Draw()
#fPois = ROOT.TF1('fPois', '[0]*TMath::Poisson(x, [1])')
#fPois.SetParameter(0, 1)
#fPois.SetParameter(1, getMean(h1))
h1.Fit('gaus',"L")
raw_input()
