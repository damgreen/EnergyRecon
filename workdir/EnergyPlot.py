#!/usr/bin/python

import sys,math,ROOT

inputFile = argv[1]
file = ROOT.TFile(inputFile)
tree = file.Get("newtree")

cut = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10"
cut01 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom < 1"
cut12 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom > 1 && CalFullLen_mom < 2"
cut23 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom > 2 && CalFullLen_mom < 3"
cut34 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom > 3"
