#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

void overlay() {

    TFile *f1 = TFile::Open("plot_NH3_1D_rgc_sum22_kaonp.root");
    TFile *f2 = TFile::Open("plot_NH3_1D_rgc_fall22_kaonp.root");
    TFile *f3 = TFile::Open("plot_NH3_1D_rgc_spring23_kaonp.root");

    TH1D *h1 = (TH1D*)f1->Get("_el_vz");
    TH1D *h2 = (TH1D*)f2->Get("_el_vz");
    TH1D *h3 = (TH1D*)f3->Get("_el_vz");

    // Normalizza le aree a 1
    if (h1->Integral() != 0)h1->Scale(1.0 / h1->Integral());
    if (h2->Integral() != 0) h2->Scale(1.0 / h2->Integral());
    if (h3->Integral() != 0) h3->Scale(1.0 / h3->Integral());
    // Stile
    h1->SetLineColor(kBlue+1);
    h2->SetLineColor(kRed+1);
    h3->SetLineColor(kGreen+2);

    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h3->SetLineWidth(2);

    // Canvas
    TCanvas *c = new TCanvas("c", "Overlay", 800, 600);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");

    // Legenda
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(h1, "Summer 2022", "l");
    leg->AddEntry(h2, "Fall 2022", "l");
    leg->AddEntry(h3, "Spring 2023", "l");
    leg->Draw();

    // Crea file ROOT di output
    TFile *fout = new TFile("overlay_normalized.root", "RECREATE");

    // Salva gli istogrammi con nomi diversi
    h1->Write("el_vz_sum22");
    h2->Write("el_vz_fall22");
    h3->Write("el_vz_spring23");
    // Salva anche il canvas completo
    c->Write("c_el_vz");

    fout->Close();

    // Salva anche PDF, opzionale
    //c->SaveAs("el_vz_overlay_normalized.pdf");
}