#!/usr/bin/python

import math,sys
import numpy as np
import ROOT

def getLongitudinalPositionErr(longPos,hit,ene):
  axis = (m_axis[1,0], m_axis[1,1], m_axis[1,2])
  x = math.fabs(longPos)
  p0 = 0.1 + math.exp(0.045*x-5.)
  diffVec = np.subtract(hit,m_centroid)
  crossProd = np.cross(axis,diffVec)

  c = np.dot(diffVec,axis)
  d = math.sqrt(np.dot(crossProd,crossProd))

  p0 = min(8., p0*(1. + 0.0073*c + 0.11*d - 0.0003*c*d))
  p1 = max(4.2 - math.exp(0.09*x-11.), 0.)
  err = p0 + p1/math.sqrt(ene/10)
  return err

def getWeightErr(longPos,ene):
  x = math.fabs(longPos)
  p0 = 0.000557 + math.exp(0.0483*x -11.6)
  p1 = max(0.00821 - math.exp(0.0358*x -10.5), 0.)
  err = p0 + p1/math.sqrt(ene/10.)
  return err

def point2Pos(dataVec):
 
  idXtal = (dataVec[0],dataVec[1],dataVec[2])
  pXtal  = (dataVec[3],dataVec[4],dataVec[5])
 
  towerIdy = int(idXtal[0]/4.)
  towerIdx = int(idXtal[0] - 4*towerIdy)
 
  if idXtal[2] %2 == 0:
    minCoord = -1.5*towerPitch + towerPitch*towerIdx - CsILength/2
    maxCoord = -1.5*towerPitch + towerPitch*towerIdx + CsILength/2
    xLong = pXtal[0] - (minCoord + maxCoord)/2
  else:
    minCoord = -1.5*towerPitch + towerPitch*towerIdy - CsILength/2
    maxCoord = -1.5*towerPitch + towerPitch*towerIdy + CsILength/2
    xLong = pXtal[1] - (minCoord + maxCoord)/2

  return xLong

def calcCoordAlongAxis(point,centroid,axis):
  diffVec = np.subtract(point,centroid)
  dotProd = np.dot(diffVec,axis)
  return dotProd

def calcDistToAxis(point,centroid,axis):
  diffVec = np.subtract(centroid,point)
  crossProd = np.cross(axis,diffVec)
  return math.sqrt(np.dot(crossProd,crossProd))

def getEntry(initPos,initDir):

  calZTop = -48.12
  calZBot = -217.47
  calXLo = -728.22
  calXHi = 728.22
  calYLo = -728.22
  calYHi = 728.22

  x0 = initPos[0]
  y0 = initPos[1]
  z0 = initPos[2]
 
  ux0 = initDir[0]
  uy0 = initDir[1]
  uz0 = initDir[2]

  if uz0 != 0:
    calEntry = (initPos + np.multiply((calZTop - z0)/uz0,initDir))

    if calEntry[0] > calXHi:
      calEntry = (initPos + np.multiply((calXHi - x0)/ux0,initDir))
    if calEntry[0] < calXLo:
      calEntry = (initPos + np.multiply((calXLo - x0)/ux0,initDir))
    if calEntry[1] > calYHi:
      calEntry = (initPos + np.multiply((calYHi - y0)/uy0,initDir))
    if calEntry[1] < calYLo:
      calEntry = (initPos + np.multiply((calYLo - y0)/uy0,initDir))
 
    return calEntry
  else:
    return np.zeros(3)

def getExit(initPos,initDir):

  calZTop = -48.12
  calZBot = -217.47
  calXLo = -728.22
  calXHi = 728.22
  calYLo = -728.22
  calYHi = 728.22

  x0 = initPos[0]
  y0 = initPos[1]
  z0 = initPos[2]
 
  ux0 = initDir[0]
  uy0 = initDir[1]
  uz0 = initDir[2]

  if uz0 != 0:
    calExit  = (initPos + np.multiply((calZBot - z0)/uz0,initDir))

    if calExit[0] > calXHi:
      calExit = (initPos + np.multiply((calXHi - x0)/ux0,initDir))
    if calExit[0] < calXLo:
      calExit = (initPos + np.multiply((calXLo - x0)/ux0,initDir))
    if calExit[1] > calYHi:
      calExit = (initPos + np.multiply((calYHi - y0)/uy0,initDir))
    if calExit[1] < calYLo:
      calExit = (initPos + np.multiply((calYLo - y0)/uy0,initDir))
    
    return calExit
  else:
    return np.zeros(3)

def doMomentsAnalysis(dataVec,iniCentroid,coreRadius):

  global m_centroid
  global m_axis
  global m_centroid
  global m_transRms
  global m_weightSum
  global m_moment
  global m_longSkewness
  global m_fullLength
  global m_mcfullLength
  global m_longRms
  global m_longRmsAsym
  global m_coreEnergyFrac

  chisq = -1

  if len(dataVec) < 2:
    return -1
 
  m_moment = np.zeros(3)
  m_axis = np.zeros(9).reshape(3,3)
 
  Ixx = 0
  Iyy = 0
  Izz = 0
  Ixy = 0
  Ixz = 0
  Iyz = 0
  weightSum = 0
  centroid = [0,0,0]
  
  coreRadius *= 35.7
  
  for nXtal in range(len(dataVec)):

    #Get all of the relavent data from the dataVec
    point = (dataVec[nXtal,3],dataVec[nXtal,4],dataVec[nXtal,5])
    weight = dataVec[nXtal,6]
    hit = np.subtract(point,iniCentroid)
    Rsq = np.dot(hit,hit)
    xprm = hit[0]
    yprm = hit[1]
    zprm = hit[2]

    #Fill the elements of the inertia tensor
    Ixx += (Rsq - xprm*xprm)*weight
    Iyy += (Rsq - yprm*yprm)*weight
    Izz += (Rsq - zprm*zprm)*weight
    Ixy -= xprm*yprm*weight
    Ixz -= xprm*zprm*weight
    Iyz -= yprm*zprm*weight

    #Don't forget about the centroid
    weightSum += weight
    centroid[0] += weight*point[0]
    centroid[1] += weight*point[1]
    centroid[2] += weight*point[2]

  #This is some cray cray method for solving the identiy equation and finding the eigenvectors and values
  p = - (Ixx + Iyy + Izz)
  q = Iyy*Izz + Iyy*Ixx + Izz*Ixx - (Ixy*Ixy + Iyz*Iyz + Ixz*Ixz)
  r = - Ixx*Iyy*Izz + Ixx*Iyz*Iyz + Iyy*Ixz*Ixz + Izz*Ixy*Ixy - 2.*Ixy*Iyz*Ixz

  a = (3.*q - p*p)/3.
  b = (2.*p*p*p - 9.*p*q + 27.*r)/27.
  
  rad_test = b*b/4. + a*a*a/27.
  if rad_test < 0 and Ixy != 0 and Ixz != 0 and Iyz !=0:
    m_weightSum = weightSum
    m_centroid = centroid
    m_centroid /= m_weightSum

    m = 2. * math.sqrt(-a/3.)
    psi = ROOT.acos(3.*b/(a*m))/3.
  
    m_moment[0] = m*math.cos(psi) - p/3.
    m_moment[1] = m*math.cos(psi + 2.*math.pi/3.) - p/3.
    m_moment[2] = m*math.cos(psi + 4.*math.pi/3.) - p/3. 

    #Now we find the direction cosines. The dircos for the middle root is parallel to the longest principal axis
    for iroot in range(Pos):
      A = Iyz * (Ixx - m_moment[iroot]) - Ixy*Ixz
      B = Ixz * (Iyy - m_moment[iroot]) - Ixy*Iyz
      C = Ixy * (Izz - m_moment[iroot]) - Ixz*Iyz
      D = math.sqrt(1./(1./(A*A) + 1./(B*B) + 1/(C*C)))/C
    
      m_axis[iroot,0] = D*C/A
      m_axis[iroot,1] = D*C/B
      m_axis[iroot,2] = D

      if m_axis[iroot,2] < 0:
        m_axis[iroot] = np.multiply(-1,m_axis[iroot])

    #Calculate the covariance matrix on the primary axis
    covCalcStatus = calcCovariance(dataVec,iniCentroid)

    mcaxis  = (oldtree.McXDir, oldtree.McYDir, oldtree.McZDir)
    mcpoint = (oldtree.McX0,   oldtree.McY0,   oldtree.McZ0)
    mcxdir = mcaxis[0]
    mcydir = mcaxis[1]
    absmcxdir = math.fabs(mcxdir)
    absmcydir = math.fabs(mcydir)
    mctmin = 9999.
    mctmax = -9999.

    xdir = m_axis[1,0]
    ydir = m_axis[1,1]
    absxdir = math.fabs(xdir)
    absydir = math.fabs(ydir)
    chisq = 0
    tmin = 9999.
    tmax = -9999.
    tave = 0.
    tvar = 0.
    tskew = 0.
    coreEnergy = 0.
  
    for nXtal in range(len(dataVec)):
      point = (dataVec[nXtal,3],dataVec[nXtal,4],dataVec[nXtal,5])
      distToAxis = calcDistToAxis(point,m_centroid,m_axis[1])
      dataVec[nXtal,8] = distToAxis

      t = calcCoordAlongAxis(point,m_centroid,m_axis[1])
      tmc = calcCoordAlongAxis(point,mcpoint,mcaxis)

      weight = dataVec[nXtal,6]

      #Rough attempt to correct for the gaps in the CAL modules
      tower = dataVec[nXtal,0]

      if absxdir > 0.05:
        t -= towerGap*(tower % 4)/xdir
      if absydir > 0.05:
        t -= towerGap*(tower/4)/ydir
      if t < tmin:
        tmin = t
      if t > tmax:
        tmax = t

      if absmcxdir > 0.05:
        tmc -= towerGap*(tower % 4)/mcxdir
      if absmcydir > 0.05:
        tmc -= towerGap*(tower/4)/mcydir
      if tmc < mctmin:
        mctmin = tmc
      if tmc > mctmax:
        mctmax = tmc

      chisq += weight*distToAxis*distToAxis
      tave += weight*t
      tvar += weight*t*t
      tskew += weight*t*t*t
   
      if distToAxis < coreRadius:
        coreEnergy += weight

    #Normalize everything
    tave /= m_weightSum
    tvar /= m_weightSum
    tvar -= tave*tave
    tskew /= m_weightSum
    tskew -= (3*tvar*tave + tave*tave*tave)
    if tvar > 0:
      tskew /= math.pow(tvar,1.5)
 
    chisq /= weightSum*len(dataVec)

    m_fullLength = (tmax - tmin)/nucInterLen
    m_mcfullLength = (mctmax - mctmin)/nucInterLen
    longMag1 = math.fabs(m_moment[0])
    longMag2 = math.fabs(m_moment[2])
    transMag = math.fabs(m_moment[1])
    m_longRms = math.sqrt((longMag1 + longMag2)/(2.*m_weightSum))
    m_transRms = math.sqrt(transMag/m_weightSum)
    m_longRmsAsym = (longMag1 - longMag2)/(longMag1 + longMag2)
    m_longSkewness = tskew
    m_coreEnergyFrac = coreEnergy/m_weightSum

  else:
    chisq = -1

  return chisq

def calcCovariance(dataVec,iniCentroid):
  global cov
  global centroidCovMatrix

  L0 = m_moment[0]
  L1 = m_moment[1]
  L2 = m_moment[2]

  S = m_axis
  Sp = np.zeros(27).reshape(9,3)

  Sp[0,1] = S[2,0]
  Sp[0,2] = -1*S[1,0]
  Sp[1,0] = -1*S[2,0]
  Sp[1,2] = S[0,0]
  Sp[2,0] = S[1,0]
  Sp[2,1] = -1*S[0,0]
  Sp[3,1] = S[2,1]
  Sp[3,2] = -1*S[1,1]
  Sp[4,0] = -1*S[2,1]
  Sp[4,2] = S[0,1]
  Sp[5,0] = S[1,1]
  Sp[5,1] = -1*S[0,1]
  Sp[6,1] = S[2,2]
  Sp[6,2] = -1*S[1,2]
  Sp[7,0] = -1*S[2,2]
  Sp[7,2] = S[0,2]
  Sp[8,0] = S[1,2]
  Sp[8,1] = -1*S[0,2]

  D = np.zeros(54).reshape(6,9)

  D[0,0] = 1.
  D[1,4] = 1.
  D[2,8] = 1.
  D[3,1] = 0.5
  D[3,3] = 0.5
  D[4,2] = 0.5
  D[4,6] = 0.5
  D[5,5] = 0.5
  D[5,7] = 0.5

  Dplus = np.zeros(54).reshape(9,6)
  Dplus[0,0] = 1
  Dplus[1,3] = 1
  Dplus[2,4] = 1
  Dplus[3,3] = 1
  Dplus[4,1] = 1
  Dplus[5,5] = 1
  Dplus[6,4] = 1
  Dplus[7,5] = 1
  Dplus[8,2] = 1

  g1 = 0.5/(L2 - L1)
  g2 = 0.5/(L0 - L2)
  g3 = 0.5/(L1 - L0)
  Gplus = np.zeros(54).reshape(6,9)
  Gplus[0,0] = 1
  Gplus[1,4] = 1
  Gplus[2,8] = 1
  Gplus[3,5] = g1
  Gplus[3,7] = g1
  Gplus[4,2] = g2
  Gplus[4,6] = g2
  Gplus[5,1] = g3
  Gplus[5,3] = g3

  Skron = np.kron(S,S)

  Finverse = np.asmatrix(Gplus)*np.asmatrix(Skron)*np.asmatrix(Dplus)

  K_left = np.zeros(72).reshape(12,6)

  for i in range(9):
    for j in range(3):
       K_left[i,j+3] = Sp[i,j]
  K_left[9,0] = 1
  K_left[10,1] = 1
  K_left[11,2] = 1

  K = np.asmatrix(K_left)*np.asmatrix(Finverse)

  cIxx_xx = 0.0;
  cIxx_yy = 0.0;
  cIxx_zz = 0.0;
  cIxx_xy = 0.0;
  cIxx_xz = 0.0;
  cIxx_yz = 0.0;

  cIyy_yy = 0.0;
  cIyy_zz = 0.0;
  cIyy_xy = 0.0;
  cIyy_xz = 0.0;
  cIyy_yz = 0.0;

  cIzz_zz = 0.0;
  cIzz_xy = 0.0;
  cIzz_xz = 0.0;
  cIzz_yz = 0.0;

  cIxy_xy = 0.0;
  cIxy_xz = 0.0;
  cIxy_yz = 0.0;

  cIxz_xz = 0.0;
  cIxz_yz = 0.0;

  cIyz_yz = 0.0;

  cCxx = 0.0;
  cCyy = 0.0;
  cCzz = 0.0;
  wtot = 0.0;

  tlen = 27.
  llen = 326.
  zlen = 20.
  sqrt12 = math.sqrt(12)

  wtot = m_weightSum

  for nXtal in range(len(dataVec)):
    pos = (dataVec[nXtal,3],dataVec[nXtal,4],dataVec[nXtal,5])
    hit = np.subtract(pos,m_centroid)
    layer = dataVec[nXtal,3]
    xLong = dataVec[nXtal,7]
    x = hit[0]
    y = hit[1]
    z = hit[2]
    w = dataVec[nXtal,6]
    x2 = x*x
    y2 = y*y
    z2 = z*z
    w2 = w*w

    if layer %2 == 0:
      dx = getLongitudinalPositionErr(xLong,pos,w)
      dy = tlen/sqrt12
      dw = getWeightErr(xLong,w)
    else:
      dx = tlen/sqrt12
      dy = getLongitudinalPositionErr(xLong,pos,w)
      dw = getWeightErr(xLong,w)
    dz = zlen/sqrt12

    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz
    dw2 = dw*dw
    ww  = (1. - w/wtot)*(1. - w/wtot)

    cIxx_xx +=  4*w2 * (y2 * dy2 + z2 * dz2) + pow((y2 + z2),2)  * dw2;
    cIxx_yy +=  4*w2 * z2 * dz2         + (y2 + z2)*(x2 + z2) * dw2;
    cIxx_zz +=  4*w2 * y2 * dy2         + (y2 + z2)*(x2 + y2) * dw2;
    cIxx_xy += -2*w2 * x*y  * dy2         - (y2 + z2)*(x*y)        * dw2;
    cIxx_xz += -2*w2 * x*z  * dz2         - (y2 + z2)*(x*z)        * dw2;
    cIxx_yz += -2*w2 * y*z  * (dy2+dz2) - (y2 + z2)*(y*z)        * dw2;

    cIyy_yy +=  4*w2 * (x2 * dx2 + z2 * dz2) + pow((x2 + z2),2)   * dw2;
    cIyy_zz +=  4*w2 * x2 * dx2             + (x2 + z2)*(x2 + y2) * dw2;
    cIyy_xy += -2*w2 * x*y  * dx2           - (x2 + z2)*(x*y)     * dw2;
    cIyy_xz += -2*w2 * x*z  * (dx2+dz2) - (x2 + z2)*(x*z)         * dw2;
    cIyy_yz += -2*w2 * y*z  * dz2           - (x2 + z2)*(y*z)     * dw2;

    cIzz_zz += 4*w2 * (x2 * dx2 + y2 * dy2) + pow((x2 + y2),2)    * dw2;
    cIzz_xy += -2*w2 * x*y  * (dx2+dy2) - (x2 + y2)*(x*y)         * dw2;
    cIzz_xz += -2*w2 * x*z  * dx2           - (x2 + y2)*(x*z)     * dw2;
    cIzz_yz += -2*w2 * y*z  * dy2           - (x2 + y2)*(y*z)     * dw2;

    cIxy_xy += w2 * (y2 * dx2 + x2 * dy2)   + x2 * y2         * dw2;
    cIxy_xz += w2 * y*z * dx2                     + x2 * y*z          * dw2;
    cIxy_yz += w2 * x*z * dy2                     + x * y2 * z        * dw2;

    cIxz_xz += w2 * (z2 * dx2 + x2 * dz2)   + x2 * z2         * dw2;
    cIxz_yz += w2 * x*y * dz2                     + x * y * z2        * dw2;

    cIyz_yz += w2 * (y2 * dz2 + z2 * dy2)   + y2 * z2         * dw2;

    cCxx += w2 * dx2 + x2 * dw2 * ww;
    cCyy += w2 * dy2 + y2 * dw2 * ww;
    cCzz += w2 * dz2 + z2 * dw2 * ww;

  centroidCovMatrix = np.zeros(9).reshape(3,3)
  centroidCovMatrix[0,0] = cCxx
  centroidCovMatrix[1,1] = cCyy
  centroidCovMatrix[2,2] = cCzz

  VdICovMatrix = np.zeros(36).reshape(6,6)
  VdICovMatrix[0,0] = cIxx_xx;
  VdICovMatrix[0,1] = cIxx_yy;
  VdICovMatrix[0,2] = cIxx_zz;
  VdICovMatrix[0,3] = cIxx_xy;
  VdICovMatrix[0,4] = cIxx_xz;
  VdICovMatrix[0,5] = cIxx_yz;
  VdICovMatrix[1,0] = cIxx_yy;
  VdICovMatrix[1,1] = cIyy_yy;
  VdICovMatrix[1,2] = cIyy_zz;
  VdICovMatrix[1,3] = cIyy_xy;
  VdICovMatrix[1,4] = cIyy_xz;
  VdICovMatrix[1,5] = cIyy_yz;
  VdICovMatrix[2,0] = cIxx_zz;
  VdICovMatrix[2,1] = cIyy_zz;
  VdICovMatrix[2,2] = cIzz_zz;
  VdICovMatrix[2,3] = cIzz_xy;
  VdICovMatrix[2,4] = cIzz_xz;
  VdICovMatrix[2,5] = cIzz_yz;
  VdICovMatrix[3,0] = cIxx_xy;
  VdICovMatrix[3,1] = cIyy_xy;
  VdICovMatrix[3,2] = cIzz_xy;
  VdICovMatrix[3,3] = cIxy_xy;
  VdICovMatrix[3,4] = cIxy_xz;
  VdICovMatrix[3,5] = cIxy_yz;
  VdICovMatrix[4,0] = cIxx_xz;
  VdICovMatrix[4,1] = cIyy_xz;
  VdICovMatrix[4,2] = cIzz_xz;
  VdICovMatrix[4,3] = cIxy_xz;
  VdICovMatrix[4,4] = cIxz_xz;
  VdICovMatrix[4,5] = cIxz_yz;
  VdICovMatrix[5,0] = cIxx_yz;
  VdICovMatrix[5,1] = cIyy_yz;
  VdICovMatrix[5,2] = cIzz_yz;
  VdICovMatrix[5,3] = cIxy_yz;
  VdICovMatrix[5,4] = cIxz_yz;
  VdICovMatrix[5,5] = cIyz_yz;

  errorCovarianceMatrix = np.asmatrix(K)*np.asmatrix(VdICovMatrix)*np.asmatrix(K.transpose())

  cov = np.zeros(9).reshape(3,3)
  cov[0,0] = errorCovarianceMatrix[1,1]
  cov[0,1] = errorCovarianceMatrix[1,4]
  cov[0,2] = errorCovarianceMatrix[1,7]
  cov[1,0] = errorCovarianceMatrix[4,1]
  cov[1,1] = errorCovarianceMatrix[4,4]
  cov[1,2] = errorCovarianceMatrix[4,7]
  cov[2,0] = errorCovarianceMatrix[7,1]
  cov[2,1] = errorCovarianceMatrix[7,4]
  cov[2,2] = errorCovarianceMatrix[7,7]

  return 0

def doIterativeMomentAnalysis(dataVec,inputCentroid,transScaleFactor,transScaleFactorBoost,coreRadius):

  global data
  global m_numIterations
  global m_numDroppedPoints
  m_numIterations = 0
  m_numDroppedPoints = 0
 
  data = dataVec
 
  chisq = -1
  iterate = True

  centroid = inputCentroid
 
  while iterate:
  
    localChisq = doMomentsAnalysis(data,centroid,coreRadius)
    if localChisq < 0:
      break
    chisq = localChisq

    centroid = m_centroid
    transRms = m_transRms
  
    data = data[data[:,8].argsort()] 
    
    iterate = False
    m_numIterations += 1

    while len(data) != 0:
      if data[len(data)-1,8] > transScaleFactor * transRms:
        data = np.delete(data,len(data)-1,0)
        iterate = True
        m_numDroppedPoints += 1
      else:
        break
    transScaleFactor *= transScaleFactorBoost

    if len(data) < 3:
      break

  return chisq

def fillMomentsData(dataVec):

  momCentroid = np.zeros(Pos) 
  momAxis = np.zeros(Pos) 

  ene = 0
  pCluster = (0,0,0)

  for nXtal in range(len(dataVec)):
    eneXtal = dataVec[nXtal,6]
    pXtal = (dataVec[nXtal,3],dataVec[nXtal,4],dataVec[nXtal,5])
    ptmp = np.multiply(eneXtal,pXtal)
    ene += eneXtal
    pCluster = np.add(pCluster,ptmp)

  if ene > 0:
    pCluster = np.divide(pCluster,ene)
  else:
    pCluster = (-1000,-1000,-1000)

  inputCentroid = pCluster 

  numXtals = len(dataVec)
  chiSq = doIterativeMomentAnalysis(dataVec,inputCentroid,transScaleFactor,transScaleFactorBoost,coreRadius)
  dataVec = data
  if chiSq >=0 and m_axis[1,2] > 0.001 and m_numIterations > 1:
    momCentroid[0] = m_centroid[0]
    momCentroid[1] = m_centroid[1]
    momCentroid[2] = m_centroid[2]
    momAxis[0] = m_axis[1,0]
    momAxis[1] = m_axis[1,1]
    momAxis[2] = m_axis[1,2]
    NumIter_mom[0] = m_numIterations
    NumCoreXtals_mom[0] = numXtals - m_numDroppedPoints
    CalLSkew_mom[0] = m_longSkewness 
    CalFullLen_mom[0] = m_fullLength 
    McFullLen_mom[0] = m_mcfullLength 
    CalEne_mom[0] = m_weightSum
 
  momCentroidErr = centroidCovMatrix 
  momAxisErr =  cov
  
  if m_numIterations > 1:
    ChiSq_mom[0] = doMomentsAnalysis(dataVec,momCentroid,coreRadius)
    CalLongRms_mom[0] = m_longRms 
    CalTransRms_mom[0] = m_transRms
    CalLRmsAsym_mom[0] = m_longRmsAsym 
    CoreEneFrac_mom[0] = m_coreEnergyFrac 

    calEntry = getEntry(m_centroid,m_axis[1])
    calExit  = getExit(m_centroid,m_axis[1])

    CalDist[0] = np.linalg.norm(calEntry - calExit)/nucInterLen
 
    CalXDir_mom[0] = m_axis[1,0] 
    CalYDir_mom[0] = m_axis[1,1] 
    CalZDir_mom[0] = m_axis[1,2] 

    CalXEcntr_mom[0] = m_centroid[0]
    CalYEcntr_mom[0] = m_centroid[1]
    CalZEcntr_mom[0] = m_centroid[2]

    CentCovXX[0] = momCentroidErr[0,0]
    CentCovYY[0] = momCentroidErr[1,1]
    CentCovZZ[0] = momCentroidErr[2,2]

    DirCovXX[0] = momAxisErr[0,0]
    DirCovYY[0] = momAxisErr[1,1]
    DirCovZZ[0] = momAxisErr[2,2]
    DirCovXY[0] = momAxisErr[0,1]
    DirCovXZ[0] = momAxisErr[0,2]
    DirCovYZ[0] = momAxisErr[1,2]

  else:
    ChiSq_mom[0] = -1

  return 

def buildIdVec():

  XtalId = []
  for iTow in range(Tow):
    for iLay in range(Lay):
      for iLog in range(Log):
        tmp = (int(iTow),int(iLog),int(iLay))
        XtalId = np.append(XtalId,tmp)
  XtalId = XtalId.reshape(Tow*Lay*Log,Pos)

  return XtalId

def buildEneVec(tree):

  XtalEne = np.ndarray(Tow*Log*Lay)
  for nXtal in range(Tow*Log*Lay):
    XtalEne[nXtal] = tree.CalXtalEnePos[nXtal]

  return XtalEne

def buildPosVec(tree):

  XtalPos = np.ndarray(Tow*Log*Lay*Pos)
  for nXtal in range(Tow*Log*Lay*Pos):
    XtalPos[nXtal] = tree.CalXtalPos[nXtal]
  XtalPos = XtalPos.reshape(Tow*Lay*Log,Pos)

  return XtalPos

def buildDataVec(tree,XtalId):

  XtalEne = buildEneVec(tree)
  XtalPos = buildPosVec(tree)
  
  dataVec = np.concatenate((XtalId,XtalPos),axis=1)
  dataVec = np.insert(dataVec,6,XtalEne,axis=1)
  dataVec = np.insert(dataVec,7,0,axis=1)
  for n in range(len(dataVec)):
    dataVec[n,7] = math.fabs(point2Pos(dataVec[n]))
  dataVec = np.insert(dataVec,8,0,axis=1)

  return dataVec

def clearTreeVar():

  CalEne_mom = np.zeros(1,dtype=float)

  DirPsf = np.zeros(1,dtype=float)
  CovPsf = np.zeros(1,dtype=float)
  OldPsf = np.zeros(1,dtype=float)

  CalTransRms_mom = np.zeros(1,dtype=float)
  CalLongRms_mom = np.zeros(1,dtype=float)
  CalLRmsAsym_mom = np.zeros(1,dtype=float)
  CalLSkew_mom = np.zeros(1,dtype=float)

  DirCovXX = np.zeros(1,dtype=float)
  DirCovXY = np.zeros(1,dtype=float)
  DirCovXZ = np.zeros(1,dtype=float)
  DirCovYY = np.zeros(1,dtype=float)
  DirCovYZ = np.zeros(1,dtype=float)
  DirCovZZ = np.zeros(1,dtype=float)

  CentCovXX = np.zeros(1,dtype=float)
  CentCovYY = np.zeros(1,dtype=float)
  CentCovZZ = np.zeros(1,dtype=float)

  CalXDir = np.zeros(1,dtype=float)
  CalYDir = np.zeros(1,dtype=float)
  CalZDir = np.zeros(1,dtype=float)

  CalXDir_mom = np.zeros(1,dtype=float)
  CalYDir_mom = np.zeros(1,dtype=float)
  CalZDir_mom = np.zeros(1,dtype=float)

  CalXEcntr_mom = np.zeros(1,dtype=float)
  CalYEcntr_mom = np.zeros(1,dtype=float)
  CalZEcntr_mom = np.zeros(1,dtype=float)

  VarTheta = np.zeros(1,dtype=float)
  VarPhi = np.zeros(1,dtype=float)

  ChiSq_mom = np.zeros(1,dtype=float)

  CalEnergyRaw = np.zeros(1,dtype=float)

  McCharge = np.zeros(1,dtype=float)
  McEnergy = np.zeros(1,dtype=float)

  McXDir = np.zeros(1,dtype=float)
  McYDir = np.zeros(1,dtype=float)
  McZDir = np.zeros(1,dtype=float)

  McX0 = np.zeros(1,dtype=float)
  McY0 = np.zeros(1,dtype=float)
  McZ0 = np.zeros(1,dtype=float)

  CoreEneFrac_mom = np.zeros(1,dtype=float)
  NumIter_mom = np.zeros(1,dtype=float)
  NumCoreXtals_mom = np.zeros(1,dtype=float)
  CalFullLen_mom = np.zeros(1,dtype=float)
  McFullLen_mom = np.zeros(1,dtype=float)


def buildNewTree():

  global newfile,newtree,CalTransRms_mom,CalLongRms_mom,CalLRmsAsym_mom,CalLSkew_mom,DirCovXX,DirCovXY,DirCovXZ,DirCovYY,DirCovYZ,DirCovZZ,CentCovXX,CentCovYY,CentCovZZ,CalXDir_mom,CalYDir_mom,CalZDir_mom,CalXEcntr_mom,CalYEcntr_mom,CalZEcntr_mom,ChiSq_mom,CalEnergyRaw,McCharge,McEnergy,McXDir,McYDir,McZDir,McX0,McY0,McZ0,CoreEneFrac_mom,NumIter_mom,NumCoreXtals_mom,CalFullLen_mom, DirPsf, CovPsf, CalEne_mom, VarPhi, VarTheta, OldPsf, CalXDir, CalYDir, CalZDir, McFullLen_mom, CalDist, McDist;

  outputName = "CalMom" + inputName 
  newfile_name = outputName
  newfile = ROOT.TFile(newfile_name,"RECREATE")
  newtree = ROOT.TTree("newtree","newtree")

  CalEne_mom = np.zeros(1,dtype=float)

  CalDist = np.zeros(1,dtype=float)
  McDist = np.zeros(1,dtype=float)

  DirPsf = np.zeros(1,dtype=float)
  CovPsf = np.zeros(1,dtype=float)
  OldPsf = np.zeros(1,dtype=float)

  CalTransRms_mom = np.zeros(1,dtype=float)
  CalLongRms_mom = np.zeros(1,dtype=float)
  CalLRmsAsym_mom = np.zeros(1,dtype=float)
  CalLSkew_mom = np.zeros(1,dtype=float)

  DirCovXX = np.zeros(1,dtype=float)
  DirCovXY = np.zeros(1,dtype=float)
  DirCovXZ = np.zeros(1,dtype=float)
  DirCovYY = np.zeros(1,dtype=float)
  DirCovYZ = np.zeros(1,dtype=float)
  DirCovZZ = np.zeros(1,dtype=float)

  CentCovXX = np.zeros(1,dtype=float)
  CentCovYY = np.zeros(1,dtype=float)
  CentCovZZ = np.zeros(1,dtype=float)

  CalXDir_mom = np.zeros(1,dtype=float)
  CalYDir_mom = np.zeros(1,dtype=float)
  CalZDir_mom = np.zeros(1,dtype=float)

  CalXDir = np.zeros(1,dtype=float)
  CalYDir = np.zeros(1,dtype=float)
  CalZDir = np.zeros(1,dtype=float)

  CalXEcntr_mom = np.zeros(1,dtype=float)
  CalYEcntr_mom = np.zeros(1,dtype=float)
  CalZEcntr_mom = np.zeros(1,dtype=float)

  VarTheta = np.zeros(1,dtype=float)
  VarPhi = np.zeros(1,dtype=float)

  ChiSq_mom = np.zeros(1,dtype=float)

  CalEnergyRaw = np.zeros(1,dtype=float)

  McCharge = np.zeros(1,dtype=float)
  McEnergy = np.zeros(1,dtype=float)

  McXDir = np.zeros(1,dtype=float)
  McYDir = np.zeros(1,dtype=float)
  McZDir = np.zeros(1,dtype=float)

  McX0 = np.zeros(1,dtype=float)
  McY0 = np.zeros(1,dtype=float)
  McZ0 = np.zeros(1,dtype=float)

  CoreEneFrac_mom = np.zeros(1,dtype=float)
  NumIter_mom = np.zeros(1,dtype=float)
  NumCoreXtals_mom = np.zeros(1,dtype=float)
  CalFullLen_mom = np.zeros(1,dtype=float)
  McFullLen_mom = np.zeros(1,dtype=float)

  newtree.Branch('CalEne_mom',CalEne_mom,'CalEne_mom/D') 
  newtree.Branch('CalDist',CalDist,'CalDist/D') 
  newtree.Branch('McDist',McDist,'McDist/D') 

  newtree.Branch('DirPsf',DirPsf,'DirPsf/D') 
  newtree.Branch('CovPsf',CovPsf,'CovPsf/D') 
  newtree.Branch('OldPsf',OldPsf,'OldPsf/D') 

  newtree.Branch('CalTransRms_mom',CalTransRms_mom,'CalTransRms_mom/D') 
  newtree.Branch('CalLongRms_mom',CalLongRms_mom,'CalLongRms_mom/D') 
  newtree.Branch('CalLRmsAsym_mom',CalLRmsAsym_mom,'CalLRmsAsym_mom/D') 
  newtree.Branch('CalLSkew_mom',CalLSkew_mom,'CalLSkew_mom/D') 

  newtree.Branch('DirCovXX',DirCovXX,'DirCovXX/D') 
  newtree.Branch('DirCovXY',DirCovXY,'DirCovXY/D') 
  newtree.Branch('DirCovXZ',DirCovXZ,'DirCovXZ/D') 
  newtree.Branch('DirCovYY',DirCovYY,'DirCovYY/D') 
  newtree.Branch('DirCovYZ',DirCovYZ,'DirCovYZ/D') 
  newtree.Branch('DirCovZZ',DirCovZZ,'DirCovZZ/D') 

  newtree.Branch('CentCovXX',CentCovXX,'CentCovXX/D') 
  newtree.Branch('CentCovYY',CentCovYY,'CentCovYY/D') 
  newtree.Branch('CentCovZZ',CentCovZZ,'CentCovZZ/D') 

  newtree.Branch('CalXDir_mom',CalXDir_mom,'CalXDir_mom/D') 
  newtree.Branch('CalYDir_mom',CalYDir_mom,'CalYDir_mom/D') 
  newtree.Branch('CalZDir_mom',CalZDir_mom,'CalZDir_mom/D') 

  newtree.Branch('CalXDir',CalXDir,'CalXDir/D') 
  newtree.Branch('CalYDir',CalYDir,'CalYDir/D') 
  newtree.Branch('CalZDir',CalZDir,'CalZDir/D') 

  newtree.Branch('CalXEcntr_mom',CalXEcntr_mom,'CalXEcntr_mom/D') 
  newtree.Branch('CalYEcntr_mom',CalYEcntr_mom,'CalYEcntr_mom/D') 
  newtree.Branch('CalZEcntr_mom',CalZEcntr_mom,'CalZEcntr_mom/D') 

  newtree.Branch('VarTheta',VarTheta,'VarTheta/D') 
  newtree.Branch('VarPhi',VarPhi,'VarPhi/D') 

  newtree.Branch('ChiSq_mom',ChiSq_mom,'ChiSq_mom/D') 

  newtree.Branch('CalEnergyRaw',CalEnergyRaw,'CalEnergyRaw/D') 

  newtree.Branch('McCharge',McCharge,'McCharge/D') 
  newtree.Branch('McEnergy',McEnergy,'McEnergy/D') 
  newtree.Branch('McXDir',McXDir,'McXDir/D') 
  newtree.Branch('McYDir',McYDir,'McYDir/D') 
  newtree.Branch('McZDir',McZDir,'McZDir/D') 
  newtree.Branch('McX0',McX0,'McX0/D') 
  newtree.Branch('McY0',McY0,'McY0/D') 
  newtree.Branch('McZ0',McZ0,'McZ0/D') 

  newtree.Branch('CoreEneFrac_mom',CoreEneFrac_mom,'CoreEneFrac_mom/D') 
  newtree.Branch('NumIter_mom',NumIter_mom,'NumIter_mom/D') 
  newtree.Branch('NumCoreXtals_mom',NumCoreXtals_mom,'NumCoreXtals_mom/D') 
  newtree.Branch('CalFullLen_mom',CalFullLen_mom,'CalFullLen_mom/D') 
  newtree.Branch('McFullLen_mom',McFullLen_mom,'McFullLen_mom/D') 

  return


#Set up a couple of constants that used in the moment analysis
Tow = 16
Log = 12
Lay = 8
Pos = 3
towerGap = 45.6
CsILength = 326.0
nucInterLen = 368.653
towerPitch = 374.5
transScaleFactor = 1
transScaleFactorBoost = 2
coreRadius = 0.75

inputName = sys.argv[1]
oldfile = ROOT.TFile(inputName)
oldtree = oldfile.Get("tuple")

XtalCut = 400

nent = oldtree.GetEntries()
#nent = 1000

XtalId = buildIdVec()
buildNewTree()


for nEvent in range(nent):
  oldtree.GetEntry(nEvent)


  clearTreeVar()

  #Build the data Vector with all of the needed information
  dataVec = buildDataVec(oldtree,XtalId)

  #Trim the data vector using the Crystal energy cut
  dataVec = dataVec[~(dataVec[:,6]<XtalCut)]

  #Preform the moment analysis
  fillMomentsData(dataVec)

  CalXDir[0] = oldtree.CalXDir 
  CalYDir[0] = oldtree.CalYDir 
  CalZDir[0] = oldtree.CalZDir 

  McCharge[0] = oldtree.McCharge
  McEnergy[0] = oldtree.McEnergy
  CalEnergyRaw[0] = oldtree.CalEnergyRaw

  McDir = (oldtree.McXDir, oldtree.McYDir, oldtree.McZDir)
  McPoint = (oldtree.McX0, oldtree.McY0, oldtree.McZ0)

  McXDir[0] = McDir[0]
  McYDir[0] = McDir[1]
  McZDir[0] = McDir[2]
  McX0[0] = McPoint[0]
  McY0[0] = McPoint[1]
  McZ0[0] = McPoint[2]

  mcEntry = getEntry(McPoint,McDir)
  mcExit  = getExit(McPoint,McDir)
  McDist[0] = np.linalg.norm(mcEntry - mcExit)/nucInterLen

  DirPsf[0] = 180/3.14159*ROOT.sqrt((math.acos(abs(CalZDir_mom[0])) - math.acos(abs(McZDir[0])))**2 + (math.atan(abs(CalYDir_mom[0])/abs(CalXDir_mom[0])) - math.atan(abs(McYDir[0])/abs(McXDir[0])))**2)

  VarPhi[0] = (CalYDir_mom[0]*CalYDir_mom[0]*DirCovXX[0] + CalXDir_mom[0]*CalXDir_mom[0]*DirCovYY[0] -  2*CalXDir_mom[0]*CalYDir_mom[0]*DirCovXY[0])/(CalXDir_mom[0]*CalXDir_mom[0] + CalYDir_mom[0]*CalYDir_mom[0])**2 
 
  VarTheta[0] = DirCovZZ[0]/(1-CalZDir_mom[0]*CalZDir_mom[0])

  CovPsf[0] = ROOT.sqrt(VarPhi[0] + VarTheta[0])*180/3.14159

  OldPsf[0] = 180/3.14159*ROOT.sqrt((math.acos(abs(CalZDir[0])) - math.acos(abs(McZDir[0])))**2 + (math.atan(abs(CalYDir[0])/abs(CalXDir[0])) - math.atan(abs(McYDir[0])/abs(McXDir[0])))**2)

  newtree.Fill()
  
  if nEvent % 1000 == 0:
    print "Event ", nEvent

newtree.Write()
