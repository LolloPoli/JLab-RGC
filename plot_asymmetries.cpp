#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

struct Row {
    int    bin_zPt;
    double mean_z;
    double AUL_sinPhi,    AUL_sinPhi_err;
    double AUL_sin2Phi,   AUL_sin2Phi_err;
    double ALL,           ALL_err;
    double ALL_cosPhi,    ALL_cosPhi_err;
    double ALU_sinPhi,    ALU_sinPhi_err;
};

std::vector<Row> readCSV(const std::string& filename) {
    std::vector<Row> data;
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Cannot open " << filename << std::endl;
        return data;
    }
    std::string line;
    std::getline(f, line);
    std::map<std::string, int> col;
    std::stringstream ss(line);
    std::string token;
    int idx = 0;
    while (std::getline(ss, token, ',')) {
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        col[token] = idx++;
    }
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<double> vals;
        std::stringstream ss2(line);
        std::string tok;
        while (std::getline(ss2, tok, ',')) vals.push_back(std::stod(tok));
        Row r;
        r.bin_zPt          = (int)vals[col["bin_zPt"]];
        r.mean_z           = vals[col["mean_z"]];
        r.AUL_sinPhi       = vals[col["AUL_sinPhi"]];
        r.AUL_sinPhi_err   = vals[col["AUL_sinPhi_err"]];
        r.AUL_sin2Phi      = vals[col["AUL_sin2Phi"]];
        r.AUL_sin2Phi_err  = vals[col["AUL_sin2Phi_err"]];
        r.ALL              = vals[col["ALL"]];
        r.ALL_err          = vals[col["ALL_err"]];
        r.ALL_cosPhi       = vals[col["ALL_cosPhi"]];
        r.ALL_cosPhi_err   = vals[col["ALL_cosPhi_err"]];
        r.ALU_sinPhi       = vals[col["ALU_sinPhi"]];
        r.ALU_sinPhi_err   = vals[col["ALU_sinPhi_err"]];
        data.push_back(r);
    }
    return data;
}

std::vector<Row> filterBins(const std::vector<Row>& data, int bin_min, int bin_max) {
    std::vector<Row> out;
    for (auto& r : data)
        if (r.bin_zPt >= bin_min && r.bin_zPt < bin_max)
            out.push_back(r);
    std::sort(out.begin(), out.end(),
              [](const Row& a, const Row& b){ return a.mean_z < b.mean_z; });
    return out;
}

TGraphErrors* makeGraph(const std::vector<Row>& rows,
                        double Row::* val, double Row::* err) {
    auto* g = new TGraphErrors(rows.size());
    for (int i = 0; i < (int)rows.size(); i++) {
        g->SetPoint(i, rows[i].mean_z, rows[i].*val);
        g->SetPointError(i, 0, rows[i].*err);
    }
    return g;
}

void styleGraph(TGraphErrors* g, int color, int marker) {
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(0.9);
    g->SetLineWidth(2);
}

// ── disegna un singolo pad con sovrapposizione ───────────────────────────────
void drawOverlay(TVirtualPad* pad,
                 TGraphErrors* g1, TGraphErrors* g2,TGraphErrors* g3,
                 const std::string& asym_label,
                 const std::string& pt_title,
                 const std::string& leg1, const std::string& leg2, const std::string& leg3,
                 double ymin, double ymax,
                 bool draw_ylabel) {
    pad->cd();
    gPad->SetTickx(); gPad->SetTicky();
    gPad->SetLeftMargin(draw_ylabel ? 0.15 : 0.08);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.15);

    auto* frame = gPad->DrawFrame(0.2, ymin, 1.0, ymax);
    frame->GetXaxis()->SetTitle("z");
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    if (draw_ylabel) {
        frame->GetYaxis()->SetTitle(asym_label.c_str());
        frame->SetTitle((asym_label + " | " + pt_title + " | K^{+}").c_str());
        frame->GetYaxis()->SetTitleSize(0.04);
        frame->GetYaxis()->SetTitleOffset(1.1);
    }

    TLine* zero = new TLine(0.2, 0, 1.0, 0);
    zero->SetLineStyle(2);
    zero->SetLineColor(kGray+1);
    zero->Draw();

    g1->Draw("P same");
    g2->Draw("P same");
    g3->Draw("P same");
    /*
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(22);
    latex.DrawLatex(0.55, 0.92, pt_title.c_str());
    */

    TLegend* leg = new TLegend(0.18, 0.77, 0.38, 0.87);
    //leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(g1, leg1.c_str(), "p");
    leg->AddEntry(g2, leg2.c_str(), "p");
    leg->AddEntry(g3, leg3.c_str(), "p");
    leg->Draw();
}

// ────────────────────────────────────────────────────────────────────────────
void plot_asymmetries() {

    gStyle->SetOptStat(0);
    gStyle->SetFrameLineWidth(1);

    auto sum22 = readCSV("output_RGC_NH3_asymmetries_sum22.csv");
    auto fall22 = readCSV("output_RGC_NH3_asymmetries_fall22.csv");
    auto spring23 = readCSV("output_RGC_NH3_asymmetries_spring23.csv");

    if (sum22.empty() || fall22.empty() || spring23.empty()) {
        std::cerr << "Error loading CSV files." << std::endl;
        return;
    }

    TFile* outfile = TFile::Open("rgc_asymmetries_comparison_NH3.root", "RECREATE");

    struct AsymDef {
        std::string name;
        std::string label;
        double Row::* val;
        double Row::* err;
        double ymin, ymax;
    };

    std::vector<AsymDef> asyms = {
        {"AUL_sinPhi",  "A_{UL}^{sin#Phi}",  &Row::AUL_sinPhi,  &Row::AUL_sinPhi_err,  -0.3,  0.3},
        {"AUL_sin2Phi", "A_{UL}^{sin2#Phi}", &Row::AUL_sin2Phi, &Row::AUL_sin2Phi_err, -0.3,  0.3},
        {"ALL",         "A_{LL}",            &Row::ALL,         &Row::ALL_err,         -0.2,  0.8},
        {"ALL_cosPhi",  "A_{LL}^{cos#Phi}",  &Row::ALL_cosPhi,  &Row::ALL_cosPhi_err,  -0.4,  0.4},
        {"ALU_sinPhi",  "A_{LU}^{sin#Phi}",  &Row::ALU_sinPhi,  &Row::ALU_sinPhi_err,  -0.15, 0.15}
    };

    struct PtBin {
        std::string name;
        std::string title;
        int bmin, bmax;
    };

    std::vector<PtBin> ptbins = {
        {"PhT_bin1", "0.0 < P_{hT} < 0.25 GeV",  1,  8},
        {"PhT_bin2", "0.25 < P_{hT} < 0.5 GeV",  8,  15},
        {"PhT_bin3", "0.5 < P_{hT} < 0.8 GeV",   15, 22},
        {"PhT_bin4", "0.8 < P_{hT} < 1.4 GeV",   22, 26}
    };

    for (auto& asym : asyms) {

        // ── 4 canvas singoli (uno per bin PhT) ──────────────────────────────
        for (int ipt = 0; ipt < 4; ipt++) {

            auto sum22_sel = filterBins(sum22, ptbins[ipt].bmin, ptbins[ipt].bmax);
            auto fall22_sel = filterBins(fall22, ptbins[ipt].bmin, ptbins[ipt].bmax);
            auto spring23_sel = filterBins(spring23, ptbins[ipt].bmin, ptbins[ipt].bmax);
            auto* g1 = makeGraph(sum22_sel, asym.val, asym.err);
            auto* g2 = makeGraph(fall22_sel, asym.val, asym.err);
            auto* g3 = makeGraph(spring23_sel, asym.val, asym.err);
            styleGraph(g1, kAzure-5,  20);
            styleGraph(g2, kOrange+1, 21);
            styleGraph(g3, kGreen+2, 22);

            std::string cname_single = asym.name + "_" + ptbins[ipt].name + "_overlay";
            TCanvas* c_single = new TCanvas(
                cname_single.c_str(),
                (asym.label + " | " + ptbins[ipt].title + " | K^{+}").c_str(),
                600, 500
            );

            drawOverlay(c_single,
                        g1, g2, g3,
                        asym.label,
                        ptbins[ipt].title,
                        "NH3 sum22", "NH3 fall22", "NH3 spring23",
                        asym.ymin, asym.ymax,
                        true);

            outfile->cd();
            c_single->Write();
            delete c_single;
        }

        // ── canvas 4-pad riepilogativo ───────────────────────────────────────
        TCanvas* c4 = new TCanvas(
            asym.name.c_str(),
            (asym.label + "  |  K^{+}").c_str(),
            1800, 400
        );
        c4->Divide(4, 1, 0.0001, 0.2);

        for (int ipt = 0; ipt < 4; ipt++) {

            auto sum22_sel = filterBins(sum22, ptbins[ipt].bmin, ptbins[ipt].bmax);
            auto fall22_sel = filterBins(fall22, ptbins[ipt].bmin, ptbins[ipt].bmax);
            auto spring23_sel = filterBins(spring23, ptbins[ipt].bmin, ptbins[ipt].bmax);

            auto* g1 = makeGraph(sum22_sel, asym.val, asym.err);
            auto* g2 = makeGraph(fall22_sel, asym.val, asym.err);
            auto* g3 = makeGraph(spring23_sel, asym.val, asym.err);
            styleGraph(g1, kAzure-5,  20);
            styleGraph(g2, kOrange+1, 21);
            styleGraph(g3, kGreen+2, 22);

            g1->SetName((asym.name + "_sum22_" + ptbins[ipt].name).c_str());
            g2->SetName((asym.name + "_fall22_" + ptbins[ipt].name).c_str());
            g3->SetName((asym.name + "_spring23_" + ptbins[ipt].name).c_str());

            drawOverlay(c4->GetPad(ipt+1),
                        g1, g2, g3,
                        asym.label,
                        ptbins[ipt].title,
                        "NH3 sum22", "NH3 fall22", "NH3 spring23",
                        asym.ymin, asym.ymax,
                        ipt == 0);

            outfile->cd();
            g1->Write();
            g2->Write();
            g3->Write();
        }

        outfile->cd();
        c4->Write();
        delete c4;
    }

    outfile->Close();
    std::cout << "Done! Output: rgc_asymmetries_comparison_NH3.root" << std::endl;
}

int main() {
    plot_asymmetries();
    return 0;
}