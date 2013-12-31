{
gStyle->SetOptStat(0);
gStyle->SetPalette(1);
gStyle->SetPadRightMargin(0.13);

TFile f("CalMomPass8.root");

TCut cut = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10";
TCut cut01 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom < 1";
TCut cut12 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom > 1 && CalFullLen_mom < 2";
TCut cut23 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom > 2 && CalFullLen_mom < 3";
TCut cut34 = "DirPsf < 45 && CovPsf < 45 && ChiSq_mom > 0 && ChiSq_mom < 10 && NumCoreXtals_mom > 10 && CalFullLen_mom > 3";

Int_t bin = 25;

TH2F *hist = new TH2F("hist","hist",bin,4.5,7,bin,3.75,6.5);
newtree->Draw("log10(CalEnergyRaw/CalFullLen_mom):log10(McEnergy)>>hist",cut01);
hist->SetTitle("Incident Energy vs Energy Deposited/Path Length");
hist->SetXTitle("Log10(Incident Energy) (MeV)");
hist->SetYTitle("Log10(Deposited Energy/(Path Length)) (MeV/#\lambda_{i})");
hist->Draw("colz");
text = new TText(4.75,6,"Preliminary");
text.SetTextColor(2);
c1.Update();
pal = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette"); 
pal->GetAxis()->SetTitle("Events"); 
hist->SetMinimum(0);
gPad->SetLogz(0);
text->Draw();
c1->SaveAs("CalEnergyRaw-CalFullLen_mom01.pdf");
c1->SaveAs("CalEnergyRaw-CalFullLen_mom01.png");
gPad->SetLogz(0);

hist.Delete();
TH2F *hist = new TH2F("hist","hist",bin,4.5,7,bin,3.75,6.5);
newtree->Draw("log10(CalEnergyRaw/CalFullLen_mom):log10(McEnergy)>>hist",cut12);
hist->SetTitle("Incident Energy vs Energy Deposited/Path Length");
hist->SetXTitle("Log10(Incident Energy) (MeV)");
hist->SetYTitle("Log10(Deposited Energy/(Path Length)) (MeV/#\lambda_{i})");
hist->Draw("colz");
c1.Update();
pal = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
pal->GetAxis()->SetTitle("Events");       
hist->SetMinimum(0);
text->Draw();
gPad->SetLogz(0);
c1->SaveAs("CalEnergyRaw-CalFullLen_mom12.pdf");
c1->SaveAs("CalEnergyRaw-CalFullLen_mom12.png");
gPad->SetLogz(0);


hist.Delete();
TH2F *hist = new TH2F("hist","hist",bin,4.5,7,bin,3.75,6.5);
newtree->Draw("log10(CalEnergyRaw/CalFullLen_mom):log10(McEnergy)>>hist",cut23);
hist->SetTitle("Incident Energy vs Energy Deposited/Path Length");
hist->SetXTitle("Log10(Incident Energy) (MeV)");
hist->SetYTitle("Log10(Deposited Energy/(Path Length)) (MeV/#\lambda_{i})");
hist->Draw("colz");
c1.Update();
pal = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
pal->GetAxis()->SetTitle("Events");       
hist->SetMinimum(0);
text->Draw();
gPad->SetLogz(0);
c1->SaveAs("CalEnergyRaw-CalFullLen_mom23.pdf");
c1->SaveAs("CalEnergyRaw-CalFullLen_mom23.png");
gPad->SetLogz(0);

hist.Delete();
TH2F *hist = new TH2F("hist","hist",bin,4.5,7,bin,3.75,6.5);
newtree->Draw("log10(CalEnergyRaw/CalFullLen_mom):log10(McEnergy)>>hist",cut34);
hist->SetTitle("Incident Energy vs Energy Deposited/Path Length");
hist->SetXTitle("Log10(Incident Energy) (MeV)");
hist->SetYTitle("Log10(Deposited Energy/(Path Length)) (MeV/#\lambda_{i})");
hist->Draw("colz");
c1.Update();
pal = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
pal->GetAxis()->SetTitle("Events");       
hist->SetMinimum(0);
text->Draw();
gPad->SetLogz(0);
c1->SaveAs("CalEnergyRaw-CalFullLen_mom34.pdf");
c1->SaveAs("CalEnergyRaw-CalFullLen_mom34.png");
gPad->SetLogz(0);

hist.Delete();
TH2F *hist = new TH2F("hist","hist",bin,4.5,7,bin,3.75,6.5);
newtree->Draw("log10(CalEnergyRaw/CalFullLen_mom):log10(McEnergy)>>hist",cut);
hist->SetTitle("Incident Energy vs Energy Deposited/Path Length");
hist->SetXTitle("Log10(Incident Energy) (MeV)");
hist->SetYTitle("Log10(Deposited Energy/(Path Length)) (MeV/#\lambda_{i})");
hist->Draw("colz");
hist->Draw("colz");
c1.Update();
pal = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
pal->GetAxis()->SetTitle("Events");       
hist->SetMinimum(0);
text->Draw();
gPad->SetLogz(0);
c1->SaveAs("CalEnergyRaw-CalFullLen_mom.pdf");
c1->SaveAs("CalEnergyRaw-CalFullLen_mom.png");
gPad->SetLogz(0);

}
