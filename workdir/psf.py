import ROOT
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)

file = ROOT.TFile("CalMomP8.root")
tree = file.Get("newtree")

c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )

cut = "DirPsf < 45 && CovPsf < 45 && NumCoreXtals_mom > 10  && ChiSq_mom > 0 && ChiSq_mom < 10"

hist = ROOT.TH1F("hist","hist",90,0,25)
tree.Draw("DirPsf>>hist",cut)

bin = 0
sum = 0
while sum < (0.68)*(hist.GetEntries()):
  sum = sum + hist.GetBinContent(bin)
  bin = bin + 1
psf68 =  hist.GetBinCenter(bin)

bin = 0
sum = 0
while sum < (0.95)*(hist.GetEntries()):
  sum = sum + hist.GetBinContent(bin)
  bin = bin + 1
psf95 =  hist.GetBinCenter(bin)

bin = 0
sum = 0
while sum < (0.99)*(hist.GetEntries()):
  sum = sum + hist.GetBinContent(bin)
  bin = bin + 1
psf99 =  hist.GetBinCenter(bin)

max = 1.2*(hist.GetMaximum())
hist.SetMaximum(max)

line68 = ROOT.TLine(psf68,0,psf68,max)
line68.Draw("same")

line95 = ROOT.TLine(psf95,0,psf95,max)
line95.Draw("same")

line99 = ROOT.TLine(psf99,0,psf99,max)
line99.Draw("same")


hist.SetXTitle("|CalDir - McDir| (degrees)")
hist.SetYTitle("Events")
hist.SetTitle("Histogram for |CalDir - McDir|")

text68 = ROOT.TText(psf68+0.25,2000,"68%")
text95 = ROOT.TText(psf95+0.25,800,"95%")
text99 = ROOT.TText(psf99+0.25,200,"99%")

text68.Draw()
text95.Draw()
text99.Draw()
ROOT.gPad.SetLogy()
c1.SaveAs("psf.pdf")
