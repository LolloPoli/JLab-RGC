// ROOT: root -l -b -q 'DF_plot_Kplus_overlay.C()'
// Legge i tre CSV e salva tre canvas + nove grafici in un file ROOT.
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TColor.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

namespace DFOverlay {
struct Point { int bin; double x, df, err; };
struct Data { std::vector<Point> xb, xbq2, zpt; };

Data readCSV(const char* filename) {
    Data d;
    const double pfScale = 1.0; //5.86 / 5.0;
    std::ifstream in(filename);
    if (!in) { std::cerr << "ERRORE: non posso aprire " << filename << '\n'; return d; }
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<std::string> c;
        std::string s;
        while (std::getline(ss, s, ',')) c.push_back(s);
        if (c.size() < 9) continue;
        try {
            const int bin = std::stoi(c[2]);
            const double xmin = std::stod(c[3]);
            const double xmax = std::stod(c[4]);
            const double df = std::stod(c[7]) * pfScale;
            const double err = std::stod(c[8]) * pfScale;
            if (!std::isfinite(df) || !std::isfinite(err) || err < 0) continue;
            if (c[1] == "xB") {
                if (bin == 1 || bin == 2 || bin == 20) continue;
                d.xb.push_back({bin, (xmin+xmax)/2., df, err});
            } else if (c[1] == "xB_Q2") {
                if (bin >= 1 && bin <= 14) d.xbq2.push_back({bin, double(bin), df, err});
            } else if (c[1] == "z_Pt") {
                if (bin >= 1 && bin <= 25) d.zpt.push_back({bin, double(bin), df, err});
            }
        } catch (...) { std::cerr << "Riga ignorata in " << filename << ": " << line << '\n'; }
    }
    std::cout << filename << ": xB=" << d.xb.size()
              << " xB_Q2=" << d.xbq2.size() << " z_Pt=" << d.zpt.size() << '\n';
    return d;
}

TGraphErrors* makeGraph(const std::vector<Point>& pts, const char* name,
                        Color_t color, Style_t marker, double markerSize) {
    auto* g = new TGraphErrors((int)pts.size());
    g->SetName(name);
    for (int i=0; i<(int)pts.size(); ++i) {
        g->SetPoint(i, pts[i].x, pts[i].df);
        g->SetPointError(i, 0., pts[i].err); // solo errori Y
    }
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(markerSize);
    g->GetXaxis()->SetRangeUser(0.0, 0.9);
    //g->SetLineWidth(2);
    return g;
}
}

void DF_plot_Kplus_overlay(const char* summerCSV = "DF_Kplus_summer22.csv",
                           const char* fallCSV   = "DF_Kplus_fall22.csv",
                           const char* springCSV = "DF_Kplus_spring23.csv",
                           const char* rootFile  = "DF_Kplus_overlay.root") {
    // =============== STILE (modifica qui) ===============
    const Color_t colors[3] = {kBlue-3, kRed-3, kGreen-3};
    const Style_t markers[3] = {20, 21, 22};
    const double markerSize = 1.0;
    const double yMin = 0.0;
    const double yMax = 0.4; 
    const double legendX1 = 0.13, legendY1 = 0.70;
    const double legendX2 = 0.35, legendY2 = 0.88;
    const bool gridY = false;
    // ===================================================

    const char* files[3] = {summerCSV, fallCSV, springCSV};
    const char* names[3] = {"Summer 2022", "Fall 2022", "Spring Inb. 2023"};
    const char* keys[3] = {"summer22", "fall22", "spring23"};
    DFOverlay::Data data[3];
    for (int j=0; j<3; ++j) data[j] = DFOverlay::readCSV(files[j]);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(4);
    TFile out(rootFile, "RECREATE");
    if (out.IsZombie()) { std::cerr << "Impossibile creare " << rootFile << '\n'; return; }

    const char* canvasNames[3] = {"c_DF_xB", "c_DF_xB_Q2", "c_DF_zPt"};
    const char* plotNames[3] = {"xB", "xB_Q2", "z_Pt"};
    const char* titles[3] = {
        "Dilution factor in the e^{-}K^{+}X channel; x_{B}; D_{F}^{K^{+}}",
        "e^{-}K^{+} dilution factor; Bin (x_{B}, Q^{2}); DF_{K^{+}}",
        "e^{-}K^{+} dilution factor; Bin (z, P_{T}); DF_{K^{+}}"
    };
    for (int k=0; k<3; ++k) {
        auto* c = new TCanvas(canvasNames[k], titles[k], k==0 ? 950 : 1150, 700);
        if (gridY) c->SetGridy();
        //c->SetLeftMargin(0.12);
        //c->SetBottomMargin(0.12);
        auto* legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        bool first = true;
        for (int j=0; j<3; ++j) {
            const std::vector<DFOverlay::Point>* pts =
                k==0 ? &data[j].xb : (k==1 ? &data[j].xbq2 : &data[j].zpt);
            if (pts->empty()) continue;
            auto* g = DFOverlay::makeGraph(*pts, Form("g_DF_%s_%s", plotNames[k], keys[j]),
                                           colors[j], markers[j], markerSize);
            g->SetTitle(titles[k]);
            if (first) {
                g->Draw("AP");
                //g->GetXaxis()->SetLimits(k==0 ? 0.0 : 0.5, k==0 ? 0.8 : (k==1 ? 14.5 : 25.5));
                g->GetYaxis()->SetRangeUser(yMin, yMax);
                g->GetXaxis()->SetRangeUser(0.0, 0.9);
                g->GetXaxis()->SetTitleSize(0.040);
                g->GetYaxis()->SetTitleSize(0.040);
                first = false;
            } else {
                g->Draw("P SAME");
            }
            legend->AddEntry(g, names[j], "p");
            g->Write();
        }
        if (first) {
            std::cerr << "Nessun dato valido per " << plotNames[k] << '\n';
        } else {
            legend->Draw();
            c->Modified();
            c->Update();
        }
        c->Write();
    }
    out.Close();
    std::cout << "Creato " << rootFile << " (3 canvas: xB, xB_Q2, z_Pt)\n";
}
