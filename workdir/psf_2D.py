#!/usr/bin/python

import sys
import math
import ROOT
import numpy as np
import pylab as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

def getBins(Ene,ZDir,Sigma):

  c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )

  histOld = ROOT.TH1F("histOld","histOld",500,0,90)
  histCov = ROOT.TH1F("histCov","histCov",500,0,90)
  histDir = ROOT.TH1F("histDir","histDir",500,0,90)

  cut = "DirPsf < 45 && CovPsf < 45 && NumCoreXtals_mom > 10  && ChiSq_mom > 0 && ChiSq_mom < 10 && log10(McEnergy) > %f && log10(McEnergy) < (%f + %f) && abs(McZDir) > %f && abs(McZDir) < (%f + %f)" % (float(Ene), float(Ene), float(dEne), float(ZDir),float(ZDir), float(dZDir))

  tree.Draw("DirPsf>>histDir",cut)
  tree.Draw("CovPsf>>histCov",cut)
  tree.Draw("OldPsf>>histOld",cut)

  psfDir = -1
  psfOld = -1
  psfCov = -1

  if histDir.GetEntries() > nent:
    binDir = 0
    sumDir = 0
    while sumDir < (Sigma/100.)*(histDir.GetEntries()):
      sumDir = sumDir + histDir.GetBinContent(binDir)
      binDir = binDir + 1
    psfDir =  histDir.GetBinCenter(binDir)

  if histCov.GetEntries() > nent:
    binCov = 0
    sumCov = 0
    while sumCov < (Sigma/100.)*(histCov.GetEntries()):
      sumCov = sumCov + histCov.GetBinContent(binCov)
      binCov = binCov + 1
    psfCov =  histCov.GetBinCenter(binCov)

  if histOld.GetEntries() > nent:
    binOld = 0
    sumOld = 0
    while sumOld < (Sigma/100.)*(histOld.GetEntries()):
      sumOld = sumOld + histOld.GetBinContent(binOld)
      binOld = binOld + 1
    psfOld =  histOld.GetBinCenter(binOld)

  return psfDir, psfCov, psfOld

def Plots(Sigma):

  global tree

  file = ROOT.TFile(inputFile)
  tree = file.Get("newtree")

  X = []
  Y = []
  ZDir = []
  ZCov = []
  ZOld = []

  Ene = np.arange(4,7+dEne,dEne)
  zDir = np.arange(0,1+dZDir,dZDir)

  for i in Ene:
    for j in zDir:
      psf = getBins(i,j,Sigma)
      X = np.append(X,i )
      Y = np.append(Y,j )
      ZDir = np.append(ZDir,psf[0])
      ZCov = np.append(ZCov,psf[1])
      ZOld = np.append(ZOld,psf[2])

  clf()

  ##
  #Make the Covariance plots
  ##
  zi = griddata(X,Y,ZCov,Ene,zDir)
  xim, yim = meshgrid(Ene,zDir)

  cs = pcolormesh(xim,yim,zi)
  cs.cmap.set_under('w')
  cs.set_clim(0,45)
  cbar = colorbar(cs)

  cbar.ax.set_ylabel("|Covariance Error| (Degrees)")

  xlim(4.25,7)
  ylim(0,1)
  xlabel("Log10(Incident Energy) (MeV)")
  ylabel("Cos(Incident Angle)")
  titlename = "%i |Covariance Error| for McEnergy and McZDir" % (int(Sigma))
  title(titlename)
  pdfname = "Cov%i_2D_" % (int(Sigma)) + input + ".pdf"
  pngname = "Cov%i_2D_" % (int(Sigma)) + input + ".png"
  print pdfname
  plt.savefig(pngname)
  plt.savefig(pdfname)

  clf()
  ##
  #Make the modified direction plots
  ##
  zi = griddata(X,Y,ZDir,Ene,zDir)
  xim, yim = meshgrid(Ene,zDir)

  cs = pcolormesh(xim,yim,zi)
  cs.cmap.set_under('w')
  cs.set_clim(0,45)
  cbar = colorbar(cs)

  cbar.ax.set_ylabel("|CalDir - McDir| (Degrees)")

  xlim(4.25,7)
  ylim(0,1)
  xlabel("Log10(Incident Energy) (MeV)")
  ylabel("Cos(Incident Angle)")
  titlename = "%i |CalDir-McDir| for McEnergy and McZDir" % (int(Sigma))
  title(titlename)
  pdfname = "PSF%i_2D_" % (int(Sigma)) + input + ".pdf"
  pngname = "PSF%i_2D_" % (int(Sigma)) + input + ".png"
  print pdfname
  plt.savefig(pngname)
  plt.savefig(pdfname)

  clf()
  ##
  #Make the old direction plots
  ##
  zi = griddata(X,Y,ZOld,Ene,zDir)
  xim, yim = meshgrid(Ene,zDir)

  cs = pcolormesh(xim,yim,zi)
  cs.cmap.set_under('w')
  cs.set_clim(0,45)
  cbar = colorbar(cs)

  cbar.ax.set_ylabel("|CalDir - McDir| (Degrees)")

  xlim(4.25,7)
  ylim(0,1)
  xlabel("Log10(Incident Energy) (MeV)")
  ylabel("Cos(Incident Angle)")
  titlename = "%i |CalDir - McDir| for McEnergy and McZDir" % (int(Sigma))
  title(titlename)
  pdfname = "Old%i_2D_" % (int(Sigma)) + input + ".pdf"
  pngname = "Old%i_2D_" % (int(Sigma)) + input + ".png"
  print pdfname
  plt.savefig(pngname)
  plt.savefig(pdfname)

global dEne
global dZDir
global nent 
global inputFile

dEne = (7-4.25)/10.
dZDir = 1/10.
nent = 10

input = sys.argv[1]
inputFile = input + ".root"

con = (68,95,99)
for i in con:
  Plots(i)

