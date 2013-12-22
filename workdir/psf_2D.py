#!/ScienceTools/i386-apple-darwin10.8.0/bin/python

import math
import ROOT
import numpy as np
import pylab as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

def pass8_mom(Ene,ZDir,Sigma):

  c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )

  hist = ROOT.TH1F("hist","hist",500,0,90)

  cut = "DirPsf < 45 && CovPsf < 45 && NumCoreXtals_mom > 10  && ChiSq_mom > 0 && ChiSq_mom < 10 && log10(McEnergy) > %f && log10(McEnergy) < (%f + %f) && abs(McZDir) > %f && abs(McZDir) < (%f + %f)" % (float(Ene), float(Ene), float(dEne), float(ZDir),float(ZDir), float(dZDir))

  tree.Draw("DirPsf>>hist",cut)

  fig_name = "Dir_%f_%f.png" % (float(Ene),float(ZDir))

  if hist.GetEntries() > nent:
   
#    c1.SaveAs(fig_name)

    bin = 0
    sum = 0

    while sum < (Sigma/100.)*(hist.GetEntries()):
      sum = sum + hist.GetBinContent(bin)
      bin = bin + 1
    return hist.GetBinCenter(bin)

  else:
    return -1

def pass8_momerr(Ene,ZDir,Sigma):

  hist = ROOT.TH1F("hist","hist",500,0,90)

  c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )

  cut = "DirPsf < 45 && CovPsf < 45 && NumCoreXtals_mom > 10 && ChiSq_mom > 0 && ChiSq_mom < 10 && log10(McEnergy) > %f && log10(McEnergy) < (%f + %f) && abs(McZDir) > %f && abs(McZDir) < (%f + %f)" % (float(Ene), float(Ene), float(dEne), float(ZDir),float(ZDir), float(dZDir))

  fig_name = "Cov_%f_%f.png" % (float(Ene),float(ZDir))

  tree.Draw("CovPsf>>hist",cut)

  if hist.GetEntries() > nent:

#    c1.SaveAs(fig_name)

    bin = 0
    sum = 0

    while sum < (Sigma/100.)*(hist.GetEntries()):
      sum = sum + hist.GetBinContent(bin)
      bin = bin + 1
    return hist.GetBinCenter(bin)

  else:
    return -1

def pass8_old(Ene,ZDir,Sigma):

  hist = ROOT.TH1F("hist","hist",500,0,90)

  c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )

  cut = "DirPsf < 45 && CovPsf < 45 && NumCoreXtals_mom > 10 && ChiSq_mom > 0 && ChiSq_mom < 10 && log10(McEnergy) > %f && log10(McEnergy) < (%f + %f) && abs(McZDir) > %f && abs(McZDir) < (%f + %f)" % (float(Ene), float(Ene), float(dEne), float(ZDir),float(ZDir), float(dZDir))

  fig_name = "Old_%f_%f.png" % (float(Ene),float(ZDir))

  tree.Draw("OldPsf>>hist",cut)

  if hist.GetEntries() > nent:

#    c1.SaveAs(fig_name)

    bin = 0
    sum = 0

    while sum < (Sigma/100.)*(hist.GetEntries()):
      sum = sum + hist.GetBinContent(bin)
      bin = bin + 1
    return hist.GetBinCenter(bin)

  else:
    return -1

def pass8Plots_momerr(Sigma):

  global tree

  filename = "CalMomP8.root"
  file = ROOT.TFile(filename)
  tree = file.Get("newtree")

  clf()

  X = []
  Y = []
  Z = []

  Ene = np.arange(4,7+dEne,dEne)
  ZDir = np.arange(0,1+dZDir,dZDir)

  for i in Ene:
    for j in ZDir:
      X = np.append(X,i )
      Y = np.append(Y,j )
      Z = np.append(Z,pass8_momerr(i,j,Sigma))
      print i,j

  zi = griddata(X,Y,Z,Ene,ZDir)
  xim, yim = meshgrid(Ene,ZDir)

  cs = pcolormesh(xim,yim,zi)
  cs.cmap.set_under('w')
  cs.set_clim(0,45)
  cbar = colorbar(cs)

  cbar.ax.set_ylabel("|Pass8Dir_momerr| (Degrees)")

  xlim(4.25,7)
  ylim(0,1)
  xlabel("Log10(McEnergy)")
  ylabel("McZDir")
  titlename = "%i Pass8_mom Error for McEnergy and McZDir" % (int(Sigma))
  title(titlename)
  figname = "psf%i_2D_pass8_momerr.pdf" % (int(Sigma))
  print figname
  plt.savefig(figname)

def pass8Plots_old(Sigma):

  global tree

  filename = "CalMomP8.root"
  file = ROOT.TFile(filename)
  tree = file.Get("newtree")

  clf()

  X = []
  Y = []
  Z = []

  Ene = np.arange(4,7+dEne,dEne)
  ZDir = np.arange(0,1+dZDir,dZDir)

  for i in Ene:
    for j in ZDir:
      X = np.append(X,i )
      Y = np.append(Y,j )
      Z = np.append(Z,pass8_old(i,j,Sigma))
      print i,j

  zi = griddata(X,Y,Z,Ene,ZDir)
  xim, yim = meshgrid(Ene,ZDir)

  cs = pcolormesh(xim,yim,zi)
  cs.cmap.set_under('w')
  cs.set_clim(0,45)
  cbar = colorbar(cs)

  cbar.ax.set_ylabel("|CalDir - McDir| (Degrees)")

  xlim(4.25,7)
  ylim(0,1)
  xlabel("Log10(Incident Energy) (MeV)")
  ylabel("Cos(Incident Angle)")
  titlename = "%i%% Pass8 |CalDir - McDir| for McEnergy and McZDir" % (int(Sigma))
  title(titlename)
  figname = "psf%i_2D_pass8_old.pdf" % (int(Sigma))
  print figname
  plt.savefig(figname)

def pass8Plots_mom(Sigma):

  global tree

  filename = "CalMomP8.root"
  file = ROOT.TFile(filename)
  tree = file.Get("newtree")

  clf()

  X = []
  Y = []
  Z = []

  Ene = np.arange(4,7+dEne,dEne)
  ZDir = np.arange(0,1+dZDir,dZDir)

  for i in Ene:
    for j in ZDir:
      X = np.append(X,i )
      Y = np.append(Y,j )
      Z = np.append(Z,pass8_mom(i,j,Sigma))
      print i,j

  zi = griddata(X,Y,Z,Ene,ZDir)
  xim, yim = meshgrid(Ene,ZDir)

  cs = pcolormesh(xim,yim,zi)
  cs.cmap.set_under('w')
  cs.set_clim(0,45)
  cbar = colorbar(cs)

  cbar.ax.set_ylabel("|CalDir - McDir| (Degrees)")

  xlim(4.25,7)
  ylim(0,1)
  xlabel("Log10(Incident Energy) (MeV)")
  ylabel("Cos(Incident Angle)")
  titlename = "%i%% |CalDir-McDir| for McEnergy and McZDir" % (int(Sigma))
  title(titlename)
  figname = "psf%i_2D_pass8_mom.pdf" % (int(Sigma))
  print figname
  plt.savefig(figname)

global dEne
global dZDir
global nent 

dEne = (7-4.25)/10.
dZDir = 1/10.
nent = 10

con = (68,95,99)
for i in con:
  pass8Plots_mom(i)
  pass8Plots_momerr(i)
  pass8Plots_old(i)

