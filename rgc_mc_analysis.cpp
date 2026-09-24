#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TAxis.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <TLine.h>
#include <cmath>
#include <TMultiGraph.h>

using namespace std;
namespace fs = std::filesystem;
gROOT->SetBatch(kTRUE);
// to download the data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rgc/fall22/ND3/kaon_plus/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/fall22_ND3_data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/w/hallb-scshelf2102/clas12/lpolizzi/rgc/output_prova.root /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/fall22_ND3_data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rgc/mc_output/fall22_neg/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/MC_fall22_neg
// 
double MeanVect(const vector<double>& v) {
    if (v.empty()) return 0.0;
    double sum = 0.0;
    for (double x : v) sum += x;
    return sum / v.size();
}


struct PhiHSysResult {
    double sys_sin, sys_sin2, sys_cos, sys_cos2;
    double sys_sin_err, sys_sin2_err, sys_cos_err, sys_cos2_err;
    double auu_cos_mc,   auu_cos_mc_err,   auu_cos2_mc,   auu_cos2_mc_err;
    double auu_cos_reco, auu_cos_reco_err, auu_cos2_reco, auu_cos2_reco_err;
    int validBins_mc, validBins_reco;
    bool ok; // false se non ci sono abbastanza bin validi per il fit
};

// ---------------------------------------------------------------------
// Helper: fit Auu = p0*(1 + p1*cos(x) + p2*cos(2x)) su un istogramma
// COMBINATO (plus + minus, indipendente da helicity per definizione di Auu)
// ---------------------------------------------------------------------
void FitAuu(TH1D* h_plus, TH1D* h_minus, const char* tag,
            double& auu_cos, double& auu_cos_err,
            double& auu_cos2, double& auu_cos2_err) {
 
    TH1D* h_sum = (TH1D*)h_plus->Clone(Form("h_sum_auu_%s", tag));
    h_sum->Add(h_minus);
 
    if (h_sum->Integral() > 0) h_sum->Scale(1.0 / h_sum->Integral());
 
    TF1* f_auu = new TF1(Form("f_auu_%s", tag),
                          "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))",
                          -TMath::Pi(), TMath::Pi());
    f_auu->SetParLimits(0, 0.0, 1000.0);
    f_auu->SetParLimits(1, -1.0, 1.0);
    f_auu->SetParLimits(2, -1.0, 1.0);
    f_auu->SetParameters(1.0, 0.0, 0.0);
 
    h_sum->Fit(f_auu, "Q");
 
    auu_cos      = f_auu->GetParameter(1);
    auu_cos_err  = f_auu->GetParError(1);
    auu_cos2     = f_auu->GetParameter(2);
    auu_cos2_err = f_auu->GetParError(2);
}

PhiHSysResult ComputePhiHSystematic(TH1D* h_plus_mc, TH1D* h_minus_mc,
                                     TH1D* h_plus_reco, TH1D* h_minus_reco,
                                     double r_inv, const char* tag) {
 
    PhiHSysResult res{};
    res.ok = false;
 
    int nBins = h_plus_mc->GetNbinsX();
 
    TH1D* h_A_mc   = (TH1D*)h_plus_mc->Clone(Form("h_A_mc_%s", tag));
    TH1D* h_A_reco = (TH1D*)h_plus_reco->Clone(Form("h_A_reco_%s", tag));
    h_A_mc->Reset();
    h_A_reco->Reset();
 
    int validBins_mc = 0, validBins_reco = 0;
 
    for (int ib = 1; ib <= nBins; ib++) {
        // MC
        double Npos_mc = h_plus_mc->GetBinContent(ib);
        double Nneg_mc = h_minus_mc->GetBinContent(ib);
        double denom_mc = (Npos_mc + r_inv * Nneg_mc);
        if (denom_mc != 0) {
            validBins_mc++;
            double A_mc = (Npos_mc - r_inv * Nneg_mc) / denom_mc;
            double err_mc = sqrt((1 - A_mc * A_mc) / denom_mc);
            h_A_mc->SetBinContent(ib, A_mc);
            h_A_mc->SetBinError(ib, err_mc);
        }
 
        // RECO
        double Npos_reco = h_plus_reco->GetBinContent(ib);
        double Nneg_reco = h_minus_reco->GetBinContent(ib);
        double denom_reco = (Npos_reco + r_inv * Nneg_reco);
        if (denom_reco != 0) {
            validBins_reco++;
            double A_reco = (Npos_reco - r_inv * Nneg_reco) / denom_reco;
            double err_reco = sqrt((1 - A_reco * A_reco) / denom_reco);
            h_A_reco->SetBinContent(ib, A_reco);
            h_A_reco->SetBinError(ib, err_reco);
        }
    }
 
    res.validBins_mc = validBins_mc;
    res.validBins_reco = validBins_reco;

    // --- AUU ---
    FitAuu(h_plus_mc, h_minus_mc, Form("mc_%s", tag),
           res.auu_cos_mc, res.auu_cos_mc_err, res.auu_cos2_mc, res.auu_cos2_mc_err);
    FitAuu(h_plus_reco, h_minus_reco, Form("reco_%s", tag),
           res.auu_cos_reco, res.auu_cos_reco_err, res.auu_cos2_reco, res.auu_cos2_reco_err);
 
    if (validBins_mc <= 6 || validBins_reco <= 6) {
        return res; // ok resta false, non fittiamo
    }
 
    TF1* f_mc = new TF1(Form("f_mc_%s", tag),
                         "[4] + [0]*sin(x) + [1]*sin(2*x) + [2]*cos(x) + [3]*cos(2*x)",
                         -TMath::Pi(), TMath::Pi());
    f_mc->SetParLimits(0, -1.0, 1.0);
    f_mc->SetParLimits(1, -1.0, 1.0);
    f_mc->SetParLimits(2, -1.0, 1.0);
    f_mc->SetParLimits(3, -1.0, 1.0);
    h_A_mc->Fit(f_mc, "Q");
 
    TF1* f_reco = new TF1(Form("f_reco_%s", tag),
                          "[4] + [0]*sin(x) + [1]*sin(2*x) + [2]*cos(x) + [3]*cos(2*x)",
                          -TMath::Pi(), TMath::Pi());
    f_reco->SetParLimits(0, -1.0, 1.0);
    f_reco->SetParLimits(1, -1.0, 1.0);
    f_reco->SetParLimits(2, -1.0, 1.0);
    f_reco->SetParLimits(3, -1.0, 1.0);
    h_A_reco->Fit(f_reco, "Q");
 
    double A0_mc = f_mc->GetParameter(0), A1_mc = f_mc->GetParameter(1);
    double A2_mc = f_mc->GetParameter(2), A3_mc = f_mc->GetParameter(3);
    double A0_reco = f_reco->GetParameter(0), A1_reco = f_reco->GetParameter(1);
    double A2_reco = f_reco->GetParameter(2), A3_reco = f_reco->GetParameter(3);
 
    res.sys_sin  = std::abs(A0_mc - A0_reco);
    res.sys_sin2 = std::abs(A1_mc - A1_reco);
    res.sys_cos  = std::abs(A2_mc - A2_reco);
    res.sys_cos2 = std::abs(A3_mc - A3_reco);
 
    res.sys_sin_err  = sqrt(pow(f_mc->GetParError(0), 2) + pow(f_reco->GetParError(0), 2));
    res.sys_sin2_err = sqrt(pow(f_mc->GetParError(1), 2) + pow(f_reco->GetParError(1), 2));
    res.sys_cos_err  = sqrt(pow(f_mc->GetParError(2), 2) + pow(f_reco->GetParError(2), 2));
    res.sys_cos2_err = sqrt(pow(f_mc->GetParError(3), 2) + pow(f_reco->GetParError(3), 2));
 
    res.ok = true;
    return res;
}


double AUL_loglike_withLL(const double* Aul,
                          const vector<double>& phi_h,
                          const vector<double>& spinstate,   // +1 / -1 (target helicity)
                          const vector<double>& depol,       // eps
                          const vector<double>& Ptarget,     // target polarization
                          const vector<double>& helicity,    // beam helicity (+1 / -1)
                          const vector<double>& sintheta,
                          bool use_spinstate_for_helicity = true)
{
    double logLike = 0.0;

    // -------- Fit parameters --------
    const double A_sin      = Aul[0];
    const double A_sin2     = Aul[1];
    const double A_LL_0     = Aul[2];
    const double A_LL_cos1  = Aul[3];
    const double A_UU_cos1  = Aul[4];
    const double A_UU_cos2  = Aul[5];
    const double A_LU_sin   = Aul[6];
    //const double A_sin3     = Aul[7];
    //const double A_LL_cos2  = Aul[8];

    const double dilution = 0.24; // 0.24 for NH3 and 0.32 for ND3

    // -------- Relative luminosity weights --------
    double Nplus = 0.0, Nminus = 0.0;
    for (size_t i = 0; i < spinstate.size(); ++i) {
        if (spinstate[i] > 0) Nplus++;
        else                  Nminus++;
    }

    // protect against division by zero
    if (Nplus <= 0 || Nminus <= 0) return 1e12;

    const double w_plus  = 0.5 * (Nplus + Nminus) / Nplus;
    const double w_minus = 0.5 * (Nplus + Nminus) / Nminus;

    // -------- Event loop --------
    for (size_t i = 0; i < phi_h.size(); ++i) {

        const double eps = depol[i];

        // --- Depolarization factors ---
        const double D_sin   = sqrt(2.0 * eps * (1.0 + eps));
        const double D_2sin  = eps;

        const double D_LL_0    = sqrt(1.0 - eps * eps);
        const double D_LL_cos  = sqrt(2.0 * eps * (1.0 - eps));
        const double D_LL_cos2 = 1.0;

        const double D_UU_cos1 = sqrt(2.0 * eps * (1.0 + eps));
        const double D_UU_cos2 = eps;

        const double D_LU_sin = sqrt(2.0 * eps * (1.0 - eps));

        // --- UL modulation ---
        const double UL_mod =
              D_sin  * A_sin  * sin(phi_h[i])
            + D_2sin * A_sin2 * sin(2.0 * phi_h[i]);
            //- A_sin3 * sin(3.0 * phi_h[i]);

        // --- LL modulation ---
        const int lam = helicity[i];
        const double LL_mod =
              D_LL_0    * A_LL_0
            + D_LL_cos  * A_LL_cos1 * cos(phi_h[i]);
            //- D_LL_cos2 * A_LL_cos2 * cos(2.0 * phi_h[i]);

        // --- UU modulation ---
        const double UU_mod =
            D_UU_cos1 * A_UU_cos1 * cos(phi_h[i])
            + D_UU_cos2 * A_UU_cos2 * cos(2.0 * phi_h[i]);
        // --- LU modulation ---
        const double LU_mod =
            D_LU_sin * A_LU_sin * sin(phi_h[i]);


        // --- Effective polarization ---
        const double pol = dilution * fabs(Ptarget[i]);

        // Full spin-dependent term
        const double Seff = pol * spinstate[i] * (UL_mod + lam * LL_mod);

        // --- Likelihood argument ---
        const double arg = 1.0 + UU_mod + lam*LU_mod + Seff;

        // Physical protection
        if (arg <= 0.0) return 1e12;

        // --- Luminosity weight ---
        const double w = (spinstate[i] > 0 ? w_plus : w_minus);

        logLike += w * log(arg);
    }

    return -logLike;
}



bool KinematicPID_Vertex(double vz, int torus, int pid){
  bool signal = false;
  if(torus == -1){  // INBENDING
    if(pid > 20){
      if(-10 < vz && vz < 2.5) signal = true;
    } 
    if(pid < 20){
      if(-8 < vz && vz < 3) signal = true;
    }
  }
  if(torus == 1){   // OUTBENDING
    if(pid > 20){
      if(-8 < vz && vz < 3) signal = true;
    } 
    if(pid < 20){
      if(-10 < vz && vz < 2.5) signal = true;
    }
  }
  return signal;
}

std::vector<double> CreateLogBinning(int nbins, double xmin, double xmax) {
    std::vector<double> bin_edges(nbins + 1);
    double logxmin = std::log10(xmin);
    double logxmax = std::log10(xmax);
    double bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = std::pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}

void SetStatsBox(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.75);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.89);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.68);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.88);  // Posizione pannello (alto)
    }
}
void SetStatsBox2(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.12);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.26);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.68);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.88);  // Posizione pannello (alto)
    }
}


bool passRichSelection(double RL, double RQ, int nTot, int rich_pid, int track_pid, double rich_chi2, double Ph, double phi){
    double mA = 0.04, qA = -0.02;
    double mB = 0.04, qB = -0.02;
    double mC = 0.06, qC = -0.08;
    double mD = 0.06, qD = -0.08;
    if(Ph <= 3.0){
        if(track_pid == 321) return true;
        else return false;
    }
    if(nTot > 0){
        if(std::fabs(phi) >= 2.5 || std::fabs(phi) <= 0.65){
            if(Ph >= 1.5){
                // RICH: kaon and CLAS: hadron
                if (rich_pid == 321 && track_pid != 321) {
                    if (RL > 6 || RQ < 0.1 || nTot <= 2) return false;

                    double m = 0.0, q = 0.0;
                    if (rich_chi2 < 50 && Ph <= 3.0) {
                        m = mA; q = qA;
                    } else if (rich_chi2 < 50 && Ph > 3.0) {
                        m = mB; q = qB;
                    } else if (rich_chi2 >= 50 && Ph <= 3.0) {
                        m = mC; q = qC;
                    } else { 
                        m = mD; q = qD;
                    }
                    double y_limit = m * RL + q;
                    if (RQ < y_limit) return false;

                    return true;
                }
                // both RICH and CLAS see a kaon
                if (rich_pid == 321 && track_pid == 321) {
                    if (RL > 6 || RQ < 0.1 || nTot <= 2) return false;
                    return true;
                }
                // RICH: hadron and CLAS: kaon
                if (rich_pid != 321 && Ph <= 3.0) { // since we cannot trust CLAS above 3.0 GeV
                    if (RL > 6.5 && nTot <= 5) return true; // the oppisite, if the RICH is bad, we want these data, we trust CLAS
                    if (RQ < 0.1) return true;

                    double m = 0.0, q = 0.0;
                    if (rich_chi2 < 50 && Ph <= 3.0) {
                        m = mA; q = qA;
                    } else if (rich_chi2 < 50 && Ph > 3.0) {
                        m = mB; q = qB;
                    } else if (rich_chi2 >= 50 && Ph <= 3.0) {
                        m = mC; q = qC;
                    } else { // rich_chi2 >= 50 && Ph > 3.0
                        m = mD; q = qD;
                    }
                    double y_limit = m * RL + q;
                    double x_limit = (RQ - q)/m;
                    if (RL > x_limit) return true; 

                    return false;
                }
            }
        }
        else if(Ph < 3.0){
            if(track_pid == 321) return true;
            else return false;
        }
    }

    return false; // other cases
}

int getBinIndex_xQ2(double xB, double Q2){
    double bin_Q2[][2] = {{1,3}, {3, 5}, {5, 7}, {7,11}};
    vector<vector<array<double,2>>> binning_xB_for_Q2 = {
        // in sequenza i bin di Q2 per i rispettivi bin di xB
        {{0.0, 0.12}, {0.12, 0.18}, {0.18, 0.24}, {0.24, 0.3}, {0.3, 0.8}},
        {{0.0, 0.24}, {0.24, 0.3}, {0.3, 0.36}, {0.36, 0.44}, {0.44, 0.8}},
        {{0.0, 0.44}, {0.44, 0.8}},
        {{0.0, 0.55}, {0.55, 0.8}}
    };
    int binIndex = 1;
    for (int ix = 0; ix < 4; ix++){
        if (Q2 >= bin_Q2[ix][0] && Q2 <= bin_Q2[ix][1]){
            for (size_t iq = 0; iq < binning_xB_for_Q2[ix].size(); iq++){
                double xBmin = binning_xB_for_Q2[ix][iq][0];
                double xBmax = binning_xB_for_Q2[ix][iq][1];
                if (xB >= xBmin && xB < xBmax){
                    return binIndex;
                }
                binIndex++;
            }
            // nel caso non si trovi un valore di Q2
            return -1;
        } else {
            binIndex += binning_xB_for_Q2[ix].size();
        }
    }
    return -1;
}

int getBinIndex_zPt(double z, double Pt){
    double bin_Pt[][2] = {{0, 0.25}, {0.25, 0.5}, {0.5, 0.8}, {0.8, 1.4}};
    vector<vector<array<double,2>>> binning_z_for_Pt = {
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 1.0}}
    };
    int binIndex = 1;
    for (int ipt = 0; ipt < 4; ipt++){
        if (Pt >= bin_Pt[ipt][0] && Pt <= bin_Pt[ipt][1]){
            for (size_t iz = 0; iz < binning_z_for_Pt[ipt].size(); iz++){
                double z_min = binning_z_for_Pt[ipt][iz][0];
                double z_max = binning_z_for_Pt[ipt][iz][1];
                if (z >= z_min && z < z_max){
                    return binIndex;
                }
                binIndex++;
            }
            // nel caso non si trovi un valore di Pt
            return -1;
        } else {
            binIndex += binning_z_for_Pt[ipt].size();
        }
    }

    return -1;
}

int getBinIndex_z(double z){
    double bin_z[][2] = {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}};
    int binIndex = 1;
    for (int iz = 0; iz < 7; iz++){
        if (z >= bin_z[iz][0] && z <= bin_z[iz][1]){
            return binIndex;
            binIndex++;
            // nel caso non si trovi un valore di Pt
            return -1;
        } 
    }

    return -1;
}
int getBinIndex_Pt(double Pt){
    double bin_Pt[][2] = {{0, 0.25}, {0.25, 0.5}, {0.5, 0.8}, {0.8, 1.4}};
    int binIndex = 1;
    for (int iPt = 0; iPt < 7; iPt++){
        if (Pt >= bin_Pt[iPt][0] && Pt <= bin_Pt[iPt][1]){
            return binIndex;
            binIndex++;
            // nel caso non si trovi un valore di Pt
            return -1;
        } 
    }
    return -1;
}

struct PhiRatioFitResult {
    double p0, p1, p2, p3, p4, p5;
    double e0, e1, e2, e3, e4, e5;
    double chi2ndf;
};

PhiRatioFitResult FitPhiRatio(
    TH1D* h_plus,
    TH1D* h_minus,
    const TString& name) {
    PhiRatioFitResult res{};

    if (!h_plus || !h_minus) return res;
    if (h_plus->Integral() == 0 || h_minus->Integral() == 0) return res;

    // --- Clone & normalize
    TH1D* h_ratio = (TH1D*)h_plus->Clone(name);
    h_ratio->Scale(1.0 / h_plus->Integral());

    TH1D* h_minus_norm = (TH1D*)h_minus->Clone(name + "_minus");
    h_minus_norm->Scale(1.0 / h_minus->Integral());

    h_ratio->Divide(h_minus_norm);

    // --- Fit function
    TF1* f = new TF1(
        name + "_fit",
        "[0] + [1]*sin(x) + [2]*sin(2*x) + [3]*sin(3*x) + [4]*cos(x) + [5]*cos(2*x)",
        -TMath::Pi(), TMath::Pi()
    );
    f->SetParameters(1, 0, 0, 0, 0, 0);

    h_ratio->Fit(f, "QER");

    // --- Store results
    res.p0 = f->GetParameter(0); res.e0 = f->GetParError(0);
    res.p1 = f->GetParameter(1); res.e1 = f->GetParError(1);
    res.p2 = f->GetParameter(2); res.e2 = f->GetParError(2);
    res.p3 = f->GetParameter(3); res.e3 = f->GetParError(3);
    res.p4 = f->GetParameter(4); res.e4 = f->GetParError(4);
    res.p5 = f->GetParameter(5); res.e5 = f->GetParError(5);
    res.chi2ndf = f->GetChisquare() / f->GetNDF();

    

    return res;
}


// FUNZIONI UTILI ==================================================================================================

    struct SysWithCovResult {
        double A_A, A_B;
        double delta;          // |A_A - A_B|, valore grezzo
        double sigma_delta_naive;     // sqrt(sigma_A^2 + sigma_B^2), quello che avevi ora
        double sigma_delta_corrected; // con -2*Cov(A_A,A_B) incluso
        double delta_corrected;       // quadCorrect(delta, sigma_delta_corrected)
    };
    
    SysWithCovResult ComputeSysWithCov(double N_A_pos, double N_A_neg,
                                        double N_B_pos, double N_B_neg,
                                        double N_Omega_pos, double N_Omega_neg,
                                        double r_inv) {
        SysWithCovResult res{};
    
        double D_A = N_A_pos + r_inv * N_A_neg;
        double D_B = N_B_pos + r_inv * N_B_neg;
    
        res.A_A = (N_A_pos - r_inv * N_A_neg) / D_A;
        res.A_B = (N_B_pos - r_inv * N_B_neg) / D_B;
        res.delta = std::abs(res.A_A - res.A_B);
    
        double sigma_A = sqrt((1 - res.A_A * res.A_A) / D_A);
        double sigma_B = sqrt((1 - res.A_B * res.A_B) / D_B);
        res.sigma_delta_naive = sqrt(sigma_A * sigma_A + sigma_B * sigma_B);
    
        // --- covarianza annidata (B sottoinsieme di A, entrambi sottoinsiemi di Omega) ---
        // Cov(n_A, n_B) = n_B * (1 - n_A/N_Omega), calcolata separatamente per +/-
        double cov_pos = (N_Omega_pos > 0) ? N_B_pos * (1.0 - N_A_pos / N_Omega_pos) : 0.0;
        double cov_neg = (N_Omega_neg > 0) ? N_B_neg * (1.0 - N_A_neg / N_Omega_neg) : 0.0;
    
        // derivate di A_A e A_B rispetto a N_pos, N_neg (delta method)
        double dAA_dNpos =  2 * r_inv * N_A_neg / (D_A * D_A);
        double dAA_dNneg = -2 * r_inv * N_A_pos / (D_A * D_A);
        double dAB_dNpos =  2 * r_inv * N_B_neg / (D_B * D_B);
        double dAB_dNneg = -2 * r_inv * N_B_pos / (D_B * D_B);
    
        // Nota: N_A e N_B condividono gli stessi eventi "pos" e "neg" separatamente
        // (B è sottoinsieme di A per costruzione), quindi la covarianza tra le
        // derivate va presa componente per componente (pos con pos, neg con neg) -
        // non c'e' covarianza incrociata pos-neg perche' sono eventi fisicamente
        // distinti (elicita' diverse, mai lo stesso evento).
        double Cov_AB = dAA_dNpos * dAB_dNpos * cov_pos + dAA_dNneg * dAB_dNneg * cov_neg;
    
        double sigma_delta2_corrected = sigma_A * sigma_A + sigma_B * sigma_B - 2 * Cov_AB;
        res.sigma_delta_corrected = sqrt(std::max(0.0, sigma_delta2_corrected));
    
        double diff2 = res.delta * res.delta - res.sigma_delta_corrected * res.sigma_delta_corrected;
        res.delta_corrected = (diff2 > 0.0) ? sqrt(diff2) : 0.0;
    
        return res;
    }



void rgc_mc_analysis(const char* period) {
    // Istanza di TDatabasePDG per accedere alle proprietà delle particelle, non va in conflitto con clas12database, sono indipendenti
    auto db2 = TDatabasePDG::Instance();
    // Variables
    double helicity;
    double beta_CMS_Double;
    int _torus;
    TVector3 beta_CMS;
    int electron_Nphe, electron_status, electron_sector;
    double N_up, N_down;
    double electron_PCAL, electron_ECAL, electron_ECIN, electron_CAL_Tot, electron_vz;
    double electron_edge1, electron_edge2, electron_edge3;
    double Ptarget, Pbeam;
    // elettrone
    double electron_Phi;
    double electron_px, electron_py, electron_pz, electron_mom, electron_Theta, electron_ThetaDeg, electron_E, electron_W, electron_Q2, electron_ass;
    // gamma
    double gamma_px, gamma_py, gamma_pz, gamma_nu;
    // kaon + 
    //double kaonp_px, kaonp_py, kaonp_pz;
    // kaone +
    double kaonp_px, kaonp_py, kaonp_pz;
    double kaonp_xF, kaonp_xB, kaonp_Q2, kaonp_z, kaonp_Ph, kaonp_Pt;
    double kaonp_Phi_h, kaonp_Phi, kaonp_Theta, kaonp_y, kaonp_s, kaonp_E, kaonp_W, kaonp_Phi_hDeg;
    double kaonp_Ph_x, kaonp_Ph_y, kaonp_PhT, kaonp_eta, kaonp_Mx, kaonp_chi2pid;
    double kaonp_edge1, kaonp_edge2, kaonp_edge3, kaonp_vz, kaonp_SinPhi;
    double kaonp_rich_Id, kaonp_rich_PID = 0;
    double kaonp_phi_cambio, kaonp_phi_cambio2;
    double kaonp_helicity, kaonp_Phi_Hup, kaonp_Phi_Hdw;
    double kaonp_rich_id, kaonp_rich_pid, kaonp_rich_RQ, kaonp_rich_nTot, kaonp_rich_ch;
    double kaonp_PID_event, kaonp_tracker_chi2, kaonp_Pol, kaonp_epsilon;
    double kaonp_rich_tr1_x, kaonp_rich_tr1_y, kaonp_rich_tr1_z, kaonp_rich_tr1_path;
    double kaonp_rich_tr2_x, kaonp_rich_tr2_y, kaonp_rich_tr2_z, kaonp_rich_tr2_path;
    double kaonp_rich_tr3_x, kaonp_rich_tr3_y, kaonp_rich_tr3_z, kaonp_rich_tr3_path;
    double kaonp_rich_tr4_x, kaonp_rich_tr4_y, kaonp_rich_tr4_z, kaonp_rich_tr4_path;
    double kaonp_rich_chi2, kaonp_rich_Mchi2, kaonp_rich_RL, kaonp_track_pid;
    double kaonp_beta, kaonp_m, kaonp_bestMass, kaonp_Pt_zQ, kaonp_gamma;
    double kaonp_rich_tr1_edge, kaonp_rich_tr2_edge, kaonp_rich_tr3_edge, kaonp_rich_tr4_edge;
    double proton_spin, kaonp_sintheta, kaonp_sintheta2;
    double kaonp_t;
    double mc_PID;
    double C_LSA_TSA;
    double bin_Q2[][2] = {{1,3}, {3, 5}, {5, 7}, {7,11}};
    vector<vector<array<double,2>>> binning_xB_for_Q2 = {
        // in sequenza i bin di Q2 per i rispettivi bin di xB
        {{0.0, 0.12}, {0.12, 0.18}, {0.18, 0.24}, {0.24, 0.3}, {0.3, 0.8}},
        {{0.0, 0.24}, {0.24, 0.3}, {0.3, 0.36}, {0.36, 0.44}, {0.44, 0.8}},
        {{0.0, 0.44}, {0.44, 0.8}},
        {{0.0, 0.55}, {0.55, 0.8}}
    };
    double bin_Pt[][2] = {{0, 0.25}, {0.25, 0.5}, {0.5, 0.8}, {0.8, 1.4}};
    vector<vector<array<double,2>>> binning_z_for_Pt = {
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 1.0}}
    };

    //string inputDir = "fall22_NH3_data";
    //string inputDir = "fall22_NH3_data";
    TString ifile_pos = Form("MC_%s_pos", period);
    TString ifile_neg = Form("MC_%s_neg", period);
    vector<string> rootFiles = {ifile_pos.Data(), ifile_neg.Data()};

    TChain chainKaonP("Kaon+");
    TChain MC_chainKaonP("MC_Kaon+");
    int fileCount = 0;
    /*
    for (const auto &entry : fs::directory_iterator(inputDir)) {
        if (entry.path().extension() == ".root") {
            string filePath = entry.path().string();
            chainKaonP.Add(Form("%s/Kaon+", filePath.c_str()));

            fileCount++;
        }
    }
    */
    
    for (const auto& dir : rootFiles) {

        if (!fs::exists(dir)) {
            cerr << "Directory non trovata: " << dir << endl;
            continue;
        }

        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.path().extension() == ".root") {

                string filePath = entry.path().string();
                chainKaonP.Add(Form("%s/Kaon+", filePath.c_str()));
                MC_chainKaonP.Add(Form("%s/MC_Kaon+", filePath.c_str()));
                fileCount++;
            }
        }
    }
    
    if (fileCount == 0) {
        cerr << "Nessun file .root trovato nelle directory di input" << endl;
        return;
    }
    // creo un output root 
    //const char* outputFile = "studies_mc_summer22_kaonp.root"; 
    TString outputFile = Form("studies_mc_%s_kaonp.root", period);
    double torus = -1;

    //const char* csv_filename = "table_RGC_MC_summer22.csv";
    TString csv_filename = Form("table_RGC_MC_%s.csv", period);
    TString csv_filename_xQ2 = Form("table_RGC_MC_%s_xQ2.csv", period);
    TString csv_filename_zPt = Form("table_RGC_MC_%s_zPt.csv", period);
    std::ofstream csvFile(csv_filename.Data());
    std::ofstream csvFile_xQ2(csv_filename_xQ2.Data());
    std::ofstream csvFile_zPt(csv_filename_zPt.Data());
    csvFile << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err,Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err, Phih_sys_sinx, Phih_sys_sinx_err, Phih_sys_sin2x, Phih_sys_sin2x_err, Phih_sys_cosx, Phih_sys_cosx_err, Phih_sys_cos2x, Phih_sys_cos2x_err, Auu_cos_mc, Auu_cos_err_mc, Auu_cos2_mc, Auu_cos2_err_mc, Auu_cos_reco, Auu_cos_err_reco, Auu_cos2_reco, Auu_cos2_err_reco\n";
    csvFile_xQ2 << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err, Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err, Phih_sys_sinx, Phih_sys_sinx_err, Phih_sys_sin2x, Phih_sys_sin2x_err, Phih_sys_cosx, Phih_sys_cosx_err, Phih_sys_cos2x, Phih_sys_cos2x_err, Auu_cos_mc, Auu_cos_err_mc, Auu_cos2_mc, Auu_cos2_err_mc, Auu_cos_reco, Auu_cos_err_reco, Auu_cos2_reco, Auu_cos2_err_reco\n";
    csvFile_zPt << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err, Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err, Phih_sys_sinx, Phih_sys_sinx_err, Phih_sys_sin2x, Phih_sys_sin2x_err, Phih_sys_cosx, Phih_sys_cosx_err, Phih_sys_cos2x, Phih_sys_cos2x_err, Auu_cos_mc, Auu_cos_err_mc, Auu_cos2_mc, Auu_cos2_err_mc, Auu_cos_preID, Auu_cos_err_preID, Auu_cos2_preID, Auu_cos2_err_preID, Auu_cos_true, Auu_cos_err_true, Auu_cos2_true, Auu_cos2_err_true, Auu_cos_reco, Auu_cos_err_reco, Auu_cos2_reco, Auu_cos2_err_reco\n";

    TFile outFile(outputFile.Data(), "RECREATE");  // File di output ROOT
    TTree treeKaonP("Kaon+", "");
    TTree MC_treeKaonP("MC_Kaon+", "");

    TDirectory* dir_xQ2 = outFile.mkdir("Binning xB-Q2");
    TDirectory* dir_zPt = outFile.mkdir("Binning z-Pt");
    TDirectory* dir_eff = outFile.mkdir("Efficiency");
    TDirectory* dir_bin_mig = outFile.mkdir("Bin migration");
    TDirectory* dir_Phih_sin = outFile.mkdir("Phi_{h} systematics of sin(x)");
    TDirectory* dir_Phih_sin2 = outFile.mkdir("Phi_{h} systematics of sin(2x)");
    TDirectory* dir_Phih_cos = outFile.mkdir("Phi_{h} systematics of cos(x)");
    TDirectory* dir_Phih_cos2 = outFile.mkdir("Phi_{h} systematics of cos(2x)");

    double kaonp_xB_recoMC, kaonp_Q2_recoMC, kaonp_z_recoMC, kaonp_Pt_recoMC;
    double kaonp_Phi_h_recoMC, kaonp_Theta_recoMC, kaonp_Phi_recoMC;

    // To save all the variables
    // Kaon +
    //chainKaonP.SetBranchAddress("E", &kaonp_E);
    chainKaonP.SetBranchAddress("el_px", &electron_px);
    chainKaonP.SetBranchAddress("el_py", &electron_py);
    chainKaonP.SetBranchAddress("el_pz", &electron_pz);
    chainKaonP.SetBranchAddress("el_mom", &electron_mom);
    chainKaonP.SetBranchAddress("el_theta", &electron_Theta);
    chainKaonP.SetBranchAddress("el_phi", &electron_Phi);
    chainKaonP.SetBranchAddress("el_W", &electron_W);
    chainKaonP.SetBranchAddress("kaon_px", &kaonp_px);
    chainKaonP.SetBranchAddress("kaon_py", &kaonp_py);
    chainKaonP.SetBranchAddress("kaon_pz", &kaonp_pz);
    chainKaonP.SetBranchAddress("kaon_mom", &kaonp_Ph);
    chainKaonP.SetBranchAddress("beta", &kaonp_beta);
    chainKaonP.SetBranchAddress("Ptarget", &Ptarget);
    chainKaonP.SetBranchAddress("Pbeam", &Pbeam);
    chainKaonP.SetBranchAddress("Polarization", &kaonp_Pol);
    chainKaonP.SetBranchAddress("Helicity", &helicity);
    chainKaonP.SetBranchAddress("gamma", &kaonp_gamma);
    chainKaonP.SetBranchAddress("epsilon", &kaonp_epsilon);
    chainKaonP.SetBranchAddress("W", &kaonp_W);
    chainKaonP.SetBranchAddress("Q2", &kaonp_Q2);
    chainKaonP.SetBranchAddress("xF", &kaonp_xF);
    chainKaonP.SetBranchAddress("xB", &kaonp_xB);
    chainKaonP.SetBranchAddress("y", &kaonp_y);
    chainKaonP.SetBranchAddress("z", &kaonp_z);
    chainKaonP.SetBranchAddress("xB_mc", &kaonp_xB_recoMC);
    chainKaonP.SetBranchAddress("Q2_mc", &kaonp_Q2_recoMC);
    chainKaonP.SetBranchAddress("z_mc", &kaonp_z_recoMC);
    chainKaonP.SetBranchAddress("Pt_mc", &kaonp_Pt_recoMC);
    chainKaonP.SetBranchAddress("Pt", &kaonp_PhT);
    chainKaonP.SetBranchAddress("Pt_over_zQ", &kaonp_Pt_zQ);
    chainKaonP.SetBranchAddress("kaon_phi_lab", &kaonp_Phi);
    chainKaonP.SetBranchAddress("kaon_theta", &kaonp_Theta);
    chainKaonP.SetBranchAddress("eta", &kaonp_eta);
    chainKaonP.SetBranchAddress("kaon_phi_h", &kaonp_Phi_h);
    chainKaonP.SetBranchAddress("kaon_phi_lab_mc", &kaonp_Phi_recoMC);
    chainKaonP.SetBranchAddress("kaon_theta_mc", &kaonp_Theta_recoMC);
    chainKaonP.SetBranchAddress("kaon_phi_h_mc", &kaonp_Phi_h_recoMC);
    chainKaonP.SetBranchAddress("Mx", &kaonp_Mx);
    chainKaonP.SetBranchAddress("clas_chi2", &kaonp_chi2pid);
    chainKaonP.SetBranchAddress("mc_PID", &mc_PID);
    chainKaonP.SetBranchAddress("clas_pid", &kaonp_track_pid);
    chainKaonP.SetBranchAddress("Rich_Id", &kaonp_rich_id);
    chainKaonP.SetBranchAddress("Rich_PID", &kaonp_rich_pid); // PID_rich == 321 e RQ > 0.1 con 1.2<P<8 GeV sono i tagli da richiedere per il rich
    chainKaonP.SetBranchAddress("Rich_RQ", &kaonp_rich_RQ);   // PID_rich != 321 con 1.2<P<3 GeV + HadronPID_Chi2Pid + KinematicPID_Vertex per event builder
    chainKaonP.SetBranchAddress("Rich_mass", &kaonp_bestMass);
    chainKaonP.SetBranchAddress("Rich_RL", &kaonp_rich_RL);
    chainKaonP.SetBranchAddress("Rich_nTot", &kaonp_rich_nTot);
    chainKaonP.SetBranchAddress("Rich_ch", &kaonp_rich_ch);
    chainKaonP.SetBranchAddress("Rich_chi2", &kaonp_rich_chi2);
    chainKaonP.SetBranchAddress("Rich_Mchi2", &kaonp_rich_Mchi2);
    chainKaonP.SetBranchAddress("Rich_PMT_edge", &kaonp_rich_tr1_edge);
    chainKaonP.SetBranchAddress("Rich_PMT_x", &kaonp_rich_tr1_x);
    chainKaonP.SetBranchAddress("Rich_PMT_y", &kaonp_rich_tr1_y);
    chainKaonP.SetBranchAddress("Rich_aerogel1_edge", &kaonp_rich_tr2_edge);
    chainKaonP.SetBranchAddress("Rich_aerogel1_x", &kaonp_rich_tr2_x);
    chainKaonP.SetBranchAddress("Rich_aerogel1_y", &kaonp_rich_tr2_y);
    chainKaonP.SetBranchAddress("Rich_aerogel2_edge", &kaonp_rich_tr3_edge);
    chainKaonP.SetBranchAddress("Rich_aerogel2_x", &kaonp_rich_tr3_x);
    chainKaonP.SetBranchAddress("Rich_aerogel2_y", &kaonp_rich_tr3_y);
    chainKaonP.SetBranchAddress("Rich_aerogel3_edge", &kaonp_rich_tr4_edge);
    chainKaonP.SetBranchAddress("Rich_aerogel3_x", &kaonp_rich_tr4_x);
    chainKaonP.SetBranchAddress("Rich_aerogel3_y", &kaonp_rich_tr4_y);

    double electron_px_mc, electron_py_mc, electron_pz_mc, electron_mom_mc, electron_Theta_mc, electron_Phi_mc, electron_W_mc;
    double kaonp_px_mc, kaonp_py_mc, kaonp_pz_mc, kaonp_Ph_mc, kaonp_Theta_mc, kaonp_Phi_mc, kaonp_xF_mc, kaonp_xB_mc, kaonp_Q2_mc;
    double kaonp_y_mc, kaonp_z_mc, kaonp_PhT_mc, kaonp_Pt_zQ_mc, kaonp_eta_mc, kaonp_Mx_mc, kaonp_gamma_mc, kaonp_epsilon_mc;
    double helicity_mc, Ptarget_mc, Pbeam_mc, kaonp_Phi_h_mc;

    // MC Kaon +
    MC_chainKaonP.SetBranchAddress("el_px_mc", &electron_px_mc);
    MC_chainKaonP.SetBranchAddress("el_py_mc", &electron_py_mc);
    MC_chainKaonP.SetBranchAddress("el_pz_mc", &electron_pz_mc);
    MC_chainKaonP.SetBranchAddress("el_theta_mc", &electron_Theta_mc);
    MC_chainKaonP.SetBranchAddress("el_phi_mc", &electron_Phi_mc);
    MC_chainKaonP.SetBranchAddress("ep_W_mc", &electron_W_mc);
    MC_chainKaonP.SetBranchAddress("kaon_px_mc", &kaonp_px_mc);
    MC_chainKaonP.SetBranchAddress("kaon_py_mc", &kaonp_py_mc);
    MC_chainKaonP.SetBranchAddress("kaon_pz_mc", &kaonp_pz_mc);
    MC_chainKaonP.SetBranchAddress("kaon_mom_mc", &kaonp_Ph_mc);
    MC_chainKaonP.SetBranchAddress("helicity", &helicity_mc);
    MC_chainKaonP.SetBranchAddress("Ptarget", &Ptarget_mc);
    MC_chainKaonP.SetBranchAddress("Pbeam", &Pbeam_mc);
    MC_chainKaonP.SetBranchAddress("gamma_mc", &kaonp_gamma_mc);
    MC_chainKaonP.SetBranchAddress("epsilon_mc", &kaonp_epsilon_mc);
    MC_chainKaonP.SetBranchAddress("Q2_mc", &kaonp_Q2_mc);
    MC_chainKaonP.SetBranchAddress("xF_mc", &kaonp_xF_mc);
    MC_chainKaonP.SetBranchAddress("xB_mc", &kaonp_xB_mc);
    MC_chainKaonP.SetBranchAddress("y_mc", &kaonp_y_mc);
    MC_chainKaonP.SetBranchAddress("z_mc", &kaonp_z_mc);
    MC_chainKaonP.SetBranchAddress("kaon_Pt_mc", &kaonp_PhT_mc);
    MC_chainKaonP.SetBranchAddress("Pt_over_zQ_mc", &kaonp_Pt_zQ_mc);
    MC_chainKaonP.SetBranchAddress("kaon_phi_lab_mc", &kaonp_Phi_mc);
    MC_chainKaonP.SetBranchAddress("kaon_theta_mc", &kaonp_Theta_mc);
    MC_chainKaonP.SetBranchAddress("eta_mc", &kaonp_eta_mc);
    MC_chainKaonP.SetBranchAddress("kaon_phi_h_mc", &kaonp_Phi_h_mc);
    MC_chainKaonP.SetBranchAddress("Mx_mc", &kaonp_Mx_mc);
    // tree
    // Kaon +
    treeKaonP.Branch("el_px", &electron_px, "el_px/D");
    treeKaonP.Branch("el_py", &electron_py, "el_py/D");
    treeKaonP.Branch("el_pz", &electron_pz, "el_pz/D");
    treeKaonP.Branch("el_mom", &electron_mom, "el_mom/D");
    treeKaonP.Branch("el_theta", &electron_Theta, "el_theta/D");
    treeKaonP.Branch("el_phi", &electron_Phi, "el_phi/D");
    treeKaonP.Branch("el_W", &electron_W, "el_W/D");
    treeKaonP.Branch("kaon_px", &kaonp_px, "kaon_px/D");
    treeKaonP.Branch("kaon_py", &kaonp_py, "kaon_py/D");
    treeKaonP.Branch("kaon_pz", &kaonp_pz, "kaon_pz/D");
    treeKaonP.Branch("kaon_mom", &kaonp_Ph, "kaon_mom/D");
    treeKaonP.Branch("Polarization", &kaonp_Pol, "Polarization/D");
    treeKaonP.Branch("Helicity", &helicity, "Helicity/D");
    treeKaonP.Branch("spin", &proton_spin, "spin/D");
    treeKaonP.Branch("beta", &kaonp_beta, "beta/D"); 
    treeKaonP.Branch("gamma", &kaonp_gamma, "gamma/D");
    treeKaonP.Branch("epsilon", &kaonp_epsilon, "epsilon/D");
    treeKaonP.Branch("sintheta", &kaonp_sintheta, "sintheta/D");
    treeKaonP.Branch("C_LSA_TSA", &C_LSA_TSA, "C_LSA_TSA/D");
    treeKaonP.Branch("W", &kaonp_W, "W/D");
    treeKaonP.Branch("Q2", &kaonp_Q2, "Q2/D");
    treeKaonP.Branch("xF", &kaonp_xF, "xF/D");
    treeKaonP.Branch("xB", &kaonp_xB, "xB/D");
    treeKaonP.Branch("y", &kaonp_y, "y/D");
    treeKaonP.Branch("z", &kaonp_z, "z/D");
    treeKaonP.Branch("Pt", &kaonp_PhT, "kaon_Pt/D");
    treeKaonP.Branch("Pt_over_zQ", &kaonp_Pt_zQ, "Pt_over_zQ/D");
    treeKaonP.Branch("kaon_phi_lab", &kaonp_Phi, "kaon_phi_lab/D");
    treeKaonP.Branch("kaon_theta", &kaonp_Theta, "kaon_theta/D");
    treeKaonP.Branch("eta", &kaonp_eta, "eta/D");
    treeKaonP.Branch("kaon_phi_h", &kaonp_Phi_h, "kaon_phi_h/D");
    treeKaonP.Branch("Mx", &kaonp_Mx, "Mx/D");
    treeKaonP.Branch("clas_chi2", &kaonp_chi2pid, "clas_chi2/D");
    treeKaonP.Branch("mc_pid", &mc_PID, "mc_pid/D");
    treeKaonP.Branch("clas_pid", &kaonp_track_pid, "clas_pid/D");
    treeKaonP.Branch("Rich_Id", &kaonp_rich_id, "Rich_Id/D");
    treeKaonP.Branch("Rich_PID", &kaonp_rich_pid, "Rich_PID/D");
    treeKaonP.Branch("Rich_RQ", &kaonp_rich_RQ, "Rich_RQ/D");
    treeKaonP.Branch("Rich_mass", &kaonp_bestMass, "Rich_mass/D");
    treeKaonP.Branch("Rich_RL", &kaonp_rich_RL, "Rich_RL/D");
    treeKaonP.Branch("Rich_nTot", &kaonp_rich_nTot, "Rich_nTot/D");
    treeKaonP.Branch("Rich_ch", &kaonp_rich_ch, "Rich_ch/D");
    treeKaonP.Branch("Rich_chi2", &kaonp_rich_chi2, "Rich_chi2/D");
    treeKaonP.Branch("Rich_Mchi2", &kaonp_rich_Mchi2, "Rich_Mchi2/D");
    treeKaonP.Branch("Rich_PMT_edge", &kaonp_rich_tr1_edge, "Rich_PMT_edge/D");
    treeKaonP.Branch("Rich_PMT_x", &kaonp_rich_tr1_x, "Rich_PMT_x/D");
    treeKaonP.Branch("Rich_PMT_y", &kaonp_rich_tr1_y, "Rich_PMT_y/D");
    treeKaonP.Branch("Rich_aerogel1_edge", &kaonp_rich_tr2_edge, "Rich_aerogel1_edge/D");
    treeKaonP.Branch("Rich_aerogel1_x", &kaonp_rich_tr2_x, "Rich_aerogel1_x/D");
    treeKaonP.Branch("Rich_aerogel1_y", &kaonp_rich_tr2_y, "Rich_aerogel1_y/D");
    treeKaonP.Branch("Rich_aerogel2_edge", &kaonp_rich_tr3_edge, "Rich_aerogel2_edge/D");
    treeKaonP.Branch("Rich_aerogel2_x", &kaonp_rich_tr3_x, "Rich_aerogel2_x/D");
    treeKaonP.Branch("Rich_aerogel2_y", &kaonp_rich_tr3_y, "Rich_aerogel2_y/D");
    treeKaonP.Branch("Rich_aerogel3_edge", &kaonp_rich_tr4_edge, "Rich_aerogel3_edge/D");
    treeKaonP.Branch("Rich_aerogel3_x", &kaonp_rich_tr4_x, "Rich_aerogel3_x/D");
    treeKaonP.Branch("Rich_aerogel3_y", &kaonp_rich_tr4_y, "Rich_aerogel3_y/D");

    //
    // KAON+ PLOT
    //dirKaonp->cd();
    double bin = 300;
    double chi_min = 0.1;
    double chi_max = 2000;
    auto make_bins = [](int bins, double min, double max) {
        return CreateLogBinning(bins, min, max);
      };
    const auto log_chi2 = make_bins(bin, chi_min, chi_max);
    double nbin = 150;
    // Mom
    TH1D kp_evnt_chi2 ("_evnt_chi2", "#chi^{2} EventBuilder PID | only EventBuilder | 1.2 < Mom < 8 GeV ; #chi^{2}; count", 300, -8, 8);
    TH1D kp_mc_PID ("_mc_PID", "#chi^{2} MC PID | only MC | 1.2 < Mom < 8 GeV ; #chi^{2}; count", 200, -400, 2300);
    TH1D kp_m ("_best_mass", "m extracted from #beta | 1.2 < Mom < 8 GeV ; m [GeV]; count", 300, 0, 1);
    TH1D kp_deltaB ("_delta_beta", "#beta_{meas} - #beta_{th} | 1.2 < Mom < 8 GeV ; #Delta_{#beta}; count", 300, -0.05, 0.05);
    TH1D kp_Mx ("_missing_mass", "missing mass | 1.2 < Mom < 8 GeV ; M_{x} [GeV]; count", 300, 0.4, 2);
    TH1D target_pol ("_target_polarization", "NH_{3} polarization | 1.2 < Mom < 8 GeV ; polarization (z axis); count", 100, -1, 1);
    TH1D target_spin ("_target_spin", "NH_{3} spin | 1.2 < Mom < 8 GeV ; spin (z axis); count", 10, -2, 2);
    TH1D beam_helicity ("_beam_helicity", "helicity | 1.2 < Mom < 8 GeV ; helicity (z axis); count", 10, -2, 2);
    TH1D kp_phi_plus ("_kp_phi_plus", "#Phi_{h} when spin = 1 | 1.2 < Mom < 8 GeV ; #Phi_{h} [Rad]; count", 200, -M_PI, M_PI);
    TH1D kp_phi_minus ("_kp_phi_minus", "#Phi_{h} when spin = -1 | 1.2 < Mom < 8 GeV ; #Phi_{h} [Rad]; count", 200, -M_PI, M_PI);
    //TH1D kp_rich_chi2 ("_rich_chi2", "#chi^{2} RICH PID ; #chi^{2}; count", 300, -3, 3);
    TH2D kp_MomVsPhT ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  K+  | with EventBuilder + RICH ; P_{hT} [GeV]; Mom [GeV]", nbin, 0, 1.2, nbin, 1, 8);
    TH2D kp_MomVsXb ("_MomVsXb", "Correlation Mom vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; Mom [GeV]", nbin, 0, 0.8, nbin, 1, 8);
    TH2D kp_MomVsXf ("_MomVsXf", "Correlation Mom vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; Mom [GeV]", nbin, 0, 0.5, nbin, 1, 8);
    TH2D kp_MomVsZ ("_MomVsZ", "Correlation Mom vs z  |  K+  | with EventBuilder + RICH ; z; Mom [GeV]", nbin, 0.2, 1.0, nbin, 1, 8);
    TH2D kp_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  K+  | with EventBuilder + RICH ; y; Mom [GeV]", nbin, 0.2, 0.8, nbin, 1, 8);
    TH2D kp_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  K+  | with EventBuilder + RICH ; Eta; Mom [GeV]", nbin, 1.5, 3.0, nbin, 1, 8);
    TH2D kp_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  K+  | with EventBuilder + RICH ; Mom [GeV]; Theta [Rad]", nbin, 1, 8, nbin, 0.09, 0.4);
    TH2D kp_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Mom [GeV]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 1, 8);
    TH2D kp_MomVsPhi_lab ("_MomVsPhi_lab", "Correlation Mom vs #Phi_{lab}  |  K+  | with EventBuilder + RICH ; #Phi_{lab} [Rad]; Mom [GeV]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 1, 8);
    TH2D kp_MomVsMx ("_MomVsMx", "correlation Mom vs M_{x} | K+ |; Mom [GeV]; M_{x} [GeV]", nbin, 1, 8, nbin, 0.6, 3.5);
    TH2D kp_MomVsBeta ("_MomVsBeta", "correlation Mom vs #beta | K+ |; Mom [GeV]; #beta", nbin, 1, 8, nbin, 0.88, 1.02);
    TH2D kp_MomVsCh ("_MomVsCh", "correlation Mom vs #theta_{ch} | K+ |; Mom [GeV]; #theta_{ch} [rad]", nbin, 1, 8, nbin, 0.0, 0.35);
    TH2D kp_MomVsMass_Rich ("_MomVsMass_RICH", "correlation Mom vs Mass | K+ |; Mom [GeV]; Mass [GeV]", nbin, 1, 8, nbin, -0.2, 1.2);
    TH2D kp_MomvsSinTheta ("_MomvsSinTheta", "Correlation sin(#theta_{#gamma}) vs Mom  |  K+  | with EventBuilder + RICH ; Mom [GeV]; sin#theta", nbin, 1, 8, nbin, 0.0, 0.6);
    TH2D kp_MomVsEps ("_MomVsEps", "Correlation Mom vs #epsilon  |  K+  | with EventBuilder + RICH ; #epsilon; Mom [GeV]", nbin, 0.3, 1.0, nbin, 1, 8);
    // Q2
    TH2D kp_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; Q^{2} [GeV^{2}]", nbin, 0, 0.8, nbin, 1, 11);
    TH2D kp_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; Q^{2} [GeV^{2}]", nbin, 0, 0.6, nbin, 1, 10);
    TH2D kp_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  K+  | with EventBuilder + RICH ; Mom [GeV]; Q^{2} [GeV^{2}]", nbin, 0.9, 6, nbin, 1, 10);
    TH2D kp_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  K+  | with EventBuilder + RICH ; P_{hT} [GeV]; Q^{2} [GeV^{2}]", nbin, 0, 1.2, nbin, 1, 10);
    TH2D kp_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs z  |  K+  | with EventBuilder + RICH ; z; Q^{2} [GeV^{2}]", nbin, 0.2, 1.0, nbin, 1, 10);
    TH2D kp_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  K+  | with EventBuilder + RICH ; y; Q^{2} [GeV^{2}]", nbin, 0.2, 0.8, nbin, 0.9, 10);
    TH2D kp_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  K+  | with EventBuilder + RICH ; Eta; Q^{2} [GeV^{2}]", nbin, 1.5, 3.0, nbin, 1, 10);
    TH2D kp_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 1, 10);
    TH2D kp_Q2vsSinTheta ("_Q2vsSinTheta", "Correlation sin(#theta_{#gamma}) vs Q^{2}  |  K+  | with EventBuilder + RICH ; Q^{2} [GeV^{2}]; sin#theta", nbin, 1, 10, nbin, 0.0, 0.6);
    TH2D kp_Q2VsEps ("_Q2VsEps", "Correlation Q^{2} vs #epsilon  |  K+  | with EventBuilder + RICH ; #epsilon; Q^{2} [GeV^{2}]", nbin, 0.3, 1.0, nbin, 1, 10);
    TH2D kp_Q2VsMx ("_Q2VsMx", "Correlation Q^{2} vs M_{x}  |  K+  | with EventBuilder + RICH ; M_{x} [GeV]; Q^{2} [GeV^{2}]", nbin, 0.3, 3.5, nbin, 1 , 10);
    // PhT
    TH2D kp_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs z  |  K+  | with EventBuilder + RICH ; z; P_{hT} [GeV]", nbin, 0.2, 1, nbin, 0, 1.4);
    TH2D kp_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; P_{hT} [GeV]", nbin, 0, 0.8, nbin, 0, 1.2);
    TH2D kp_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  K+  | with EventBuilder + RICH ; Eta; P_{hT} [GeV]", nbin, 1.5, 3.0, nbin, 0, 1.2);
    TH2D kp_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; P_{hT} [GeV]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0, 1.2);
    TH2D kp_PhTvsSinTheta ("_PhTvsSinTheta", "Correlation sin(#theta_{#gamma}) vs P_{hT}  |  K+  | with EventBuilder + RICH ; P_{hT} [GeV]; sin#theta", nbin, 0.0, 1.2, nbin, 0.0, 0.6);
    TH2D kp_PhTVsEps ("_PhTVsEps", "Correlation P_{hT} vs #epsilon  |  K+  | with EventBuilder + RICH ; #epsilon; P_{hT} [GeV]", nbin, 0.3, 1.0, nbin, 0, 1.2);
    TH2D kp_PhTVsMx ("_PhTVsMx", "Correlation P_{hT} vs M_{x}  |  K+  | with EventBuilder + RICH ; P_{hT} [GeV]; M_{x} [GeV]", nbin, 0.0, 1.2, nbin, 0.3 , 3.5);
    // Z
    TH2D kp_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; z", nbin, 0, 0.8, nbin, 0.2, 1.0);
    TH2D kp_zVsXf ("_zVsXf", "Correlation Z vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; z", nbin, 0, 0.6, nbin, 0.2, 1.0);
    TH2D kp_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  K+  | with EventBuilder + RICH ; Eta; z", nbin, 1.5, 3.0, nbin, 0.2, 1.0);
    TH2D kp_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; z", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0.2, 1.0);
    TH2D kp_zvsSinTheta ("_zvsSinTheta", "Correlation sin(#theta_{#gamma}) vs z  |  K+  | with EventBuilder + RICH ; z; sin#theta", nbin, 0.2, 1.0, nbin, 0.0, 0.6);
    TH2D kp_zvsCosPhi ("_zvsCosPhi", "Correlation cos(#Phi_{h}) vs z  |  K+  | with EventBuilder + RICH ; z; cos_{#Phi_{h}}", nbin, 0.2, 1.0, nbin, -1, 1);
    TH2D kp_zVsEps ("_zVsEps", "Correlation z vs #epsilon  |  K+  | with EventBuilder + RICH ; #epsilon; z", nbin, 0.3, 1.0, nbin, 0.2 , 1.0);
    TH2D kp_zVsEps2 ("_zVsEps_2", "Correlation z vs #epsilon  |  K+  | with EventBuilder + RICH ; z; #epsilon", nbin, 0.2, 1.0, nbin, 0.3 , 1.0);
    TH2D kp_zVsMx ("_zVsMx", "Correlation z vs M_{x}  |  K+  | with EventBuilder + RICH ; z; M_{x} [GeV]", nbin, 0.2, 1.0, nbin, 0.3 , 3.5);
    //
    TH2D kp_xBvsY ("_xBvsY", "Correlation y vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; y", nbin, 0, 0.8, nbin, 0.2, 0.8);
    TH2D kp_xBvsEta ("_xBvsEta", "Correlation #eta vs x_{B}  |  K+  | with EventBuilder + RICH ; #eta; x_{B}", nbin, 1.5, 3.0, nbin, 0.0, 0.8);
    TH2D kp_xBvsSinTheta ("_xBvsSinTheta", "Correlation sin(#theta_{#gamma}) vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; sin#theta", nbin, 0, 0.8, nbin, 0.0, 0.6);
    TH2D kp_xBvsCosPhi ("_xBvsCosPhi", "Correlation cos(#Phi_{h}) vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; cos_{#Phi_{h}}", nbin, 0.0, 0.8, nbin, -1, 1);
    TH2D kp_xBVsEps ("_xBVsEps", "Correlation x_{B} vs #epsilon  |  K+  | with EventBuilder + RICH ; #epsilon; x_{B}", nbin, 0.3, 1.0, nbin, 0, 0.8);
    TH2D kp_xBVsMx ("_xBVsMx", "Correlation x_{B} vs M_{x}  |  K+  | with EventBuilder + RICH ; x_{B}; M_{x} [GeV]", nbin, 0.0, 0.8, nbin, 0.3 , 3.5);
    // Angles
    TH2D kp_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0.1, 0.4);
    TH2D kp_ThetaVsPhi_Lab ("_ThetaVsPhi_Lab", "Correlation #theta vs #Phi_{Lab} | K+ | with EventBuilder + RICH ; #Phi_{Lab} [Rad]; #theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0.05, 0.4);
    TH2D kp_PxVsPy ("_PxVsPy", "Correlation P_{x} vs P_{y}  |  K+  | with EventBuilder + RICH ; P_{x} [GeV]; P_{y} [GeV]", nbin, -2, 2, nbin, -2, 2);
    TH2D el_PxVsPy ("el_PxVsPy", "Correlation P_{x} vs P_{y}  |  e-  | with EventBuilder + RICH ; P_{x} [GeV]; P_{y} [GeV]", nbin, -2, 2, nbin, -2, 2);
    
    TH2D kp_rich_pmt_xy_4 ("rich_pmt_sect4_xy", "PMT RICH 4th sector, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, -170, -70, nbin, -70, 70);
    TH2D kp_rich_pmt_xy_1 ("rich_pmt_sect1_xy", "PMT RICH 1th sector, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, 70, 230, nbin, -70, 70);
    TH2D kp_rich_aerogel_xy_4 ("rich_aerogel_sect4_xy", "Aerogel RICH 4th sector, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, -230, -70, nbin, -80, 80);
    TH2D kp_rich_aerogel_xy_1 ("rich_aerogel_sect1_xy", "Aerogel RICH 1th sector, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, 70, 230, nbin, -90, 90);
    TH2D kp_rich_aerogel_l1_xy_4 ("rich_aerogel_sect4_l1_xy", "Aerogel RICH 4th sector, layer 1, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, -230, -70, nbin, -80, 80);
    TH2D kp_rich_aerogel_l2_xy_4 ("rich_aerogel_sect4_l2_xy", "Aerogel RICH 4th sector, layer 2, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, -230, -70, nbin, -80, 80);
    TH2D kp_rich_aerogel_l3_xy_4 ("rich_aerogel_sect4_l3_xy", "Aerogel RICH 4th sector, layer 3, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, -230, -70, nbin, -80, 80);
    TH2D kp_rich_aerogel_l1_xy_1 ("rich_aerogel_sect1_l1_xy", "Aerogel RICH 1th sector, layer 1, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, 70, 230, nbin, -90, 90);
    TH2D kp_rich_aerogel_l2_xy_1 ("rich_aerogel_sect1_l2_xy", "Aerogel RICH 1th sector, layer 2, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, 70, 230, nbin, -90, 90);
    TH2D kp_rich_aerogel_l3_xy_1 ("rich_aerogel_sect1_l3_xy", "Aerogel RICH 1th sector, layer 3, xy plane | t-1 & K+; x [cm]; y [cm]", nbin, 70, 230, nbin, -90, 90);

    TH1D kp_Mx_tof ("_missing_mass_ToF", "missing mass | 1.2 < Mom < 8 GeV | ToF ; M_{x} [GeV]; count", 300, 0.4, 2);
    TH1D kp_Mx_rich ("_missing_mass_RICH", "missing mass | 1.2 < Mom < 8 GeV | RICH ; M_{x} [GeV]; count", 300, 0.4, 2);


    // multidim

    double nbin_xQ2 = 14;
    double nbin_zPt = 25;
    // 2D
    vector<TH2D*> hist_xBvsQ2_binned(nbin_xQ2);
    vector<TH2D*> hist_zvsPt_binned(nbin_zPt);
    // xQ2 bins
    vector<vector<double>> vec_kaonp_phih_2d(nbin_xQ2);
    vector<vector<double>> vec_helicity_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_2phih_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_z_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_Pt_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_xB_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_Q2_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_y_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_pol_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_eps_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_spin_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_sintheta_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_phih_plus_2d(nbin_xQ2);
    vector<vector<double>> vec_kaonp_phih_minus_2d(nbin_xQ2);
    // zPt bins
    vector<vector<double>> vec_kaonp_phih_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_zPt(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_zPt(nbin_zPt);
    //
    vector<vector<vector<double>>> vec_kaonp_phih_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_2phih_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_z_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_Pt_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_xB_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_Q2_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_y_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_pol_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_eps_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_spin_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_sintheta_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_phih_plus_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_kaonp_phih_minus_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    //
    vector<vector<vector<double>>> vec_z_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_Pt_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_xB_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_Q2_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    // LSA
    vector<double> AUL_sin_2d(nbin_xQ2);
    vector<double> AUL_sin_err_2d(nbin_xQ2);
    vector<double> AUL_2sin_2d(nbin_xQ2);
    vector<double> AUL_2sin_err_2d(nbin_xQ2);
    // LL
    vector<double> ALL_0_2d(nbin_xQ2);
    vector<double> ALL_0_err_2d(nbin_xQ2);
    vector<double> ALL_cos_2d(nbin_xQ2);
    vector<double> ALL_cos_err_2d(nbin_xQ2);
    // UL vs z
    vector<double> AUL_sin_2d_zPt(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt(nbin_zPt);
    // LL vs z
    vector<double> ALL_0_2d_zPt(nbin_zPt);
    vector<double> ALL_0_err_2d_zPt(nbin_zPt);
    vector<double> ALL_cos_2d_zPt(nbin_zPt);
    vector<double> ALL_cos_err_2d_zPt(nbin_zPt);
    // LU
    vector<double> ALU_sin_2d(nbin_xQ2);
    vector<double> ALU_sin_err_2d(nbin_xQ2);
    vector<double> ALU_sin_2d_zPt(nbin_zPt);
    vector<double> ALU_sin_err_2d_zPt(nbin_zPt);
    //
    vector<vector<double>> AUL_sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    // TSA
    vector<double> AUL_gamma_sin_2d(nbin_xQ2);
    vector<double> AUL_gamma_sin_err_2d(nbin_xQ2);
    vector<double> AUL_gamma_2sin_2d(nbin_xQ2);
    vector<double> AUL_gamma_2sin_err_2d(nbin_xQ2);
    vector<double> AUL_gamma_sin_2d_zPt(nbin_zPt);
    vector<double> AUL_gamma_sin_err_2d_zPt(nbin_zPt);
    vector<double> AUL_gamma_2sin_2d_zPt(nbin_zPt);
    vector<double> AUL_gamma_2sin_err_2d_zPt(nbin_zPt);
    vector<vector<double>> AUL_gamma_sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_gamma_sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_gamma_2sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_gamma_2sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    // correction
    vector<double> Corr_sin_2d(nbin_xQ2);
    vector<double> Corr_2sin_2d(nbin_xQ2);
    vector<double> Corr_sin_2d_zPt(nbin_zPt);
    vector<double> Corr_2sin_2d_zPt(nbin_zPt);
    vector<vector<double>> Corr_sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Corr_2sin_4d(nbin_xQ2, vector<double> (nbin_zPt));

    // systematic
    vector<double> AUL_sin_sys_2d(nbin_xQ2);
    vector<double> AUL_sin_sys_err_2d(nbin_xQ2);
    vector<double> AUL_2sin_sys_2d(nbin_xQ2);
    vector<double> AUL_2sin_sys_err_2d(nbin_xQ2);
    vector<double> AUL_sin_sys_2d_zPt(nbin_zPt);
    vector<double> AUL_sin_sys_err_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_sys_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_sys_err_2d_zPt(nbin_zPt);
    vector<vector<double>> AUL_sin_sys_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_sin_sys_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_sys_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_sys_err_4d(nbin_xQ2, vector<double> (nbin_zPt));

    vector<TH1D*> hist_Phi_h_plus_binxQ2(nbin_xQ2);
    vector<TH1D*> hist_Phi_h_minus_binxQ2(nbin_xQ2);
    vector<TH1D*> hist_Phi_h_plus_binzPt(nbin_zPt);
    vector<TH1D*> hist_Phi_h_minus_binzPt(nbin_zPt);
    for(int i = 0; i < nbin_xQ2; i++){
        hist_xBvsQ2_binned[i] = new TH2D(Form("hist_xBvsQ2_bin%d", i+1), Form("x_{B} vs Q^{2} for bin (%d, x_{B}-Q^{2}) | pion+; x_{B}; Q^{2} [GeV^{2}]", i+1), 120, 0, 0.8, 120, 1, 11);
        hist_Phi_h_plus_binxQ2[i] = new TH1D(Form("hist_Phi_h_plus_binxQ2_%d", i+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) | spin +| pion+; #Phi_{h} [Rad]; count", i+1), 50, -TMath::Pi(), TMath::Pi());
        hist_Phi_h_minus_binxQ2[i] = new TH1D(Form("hist_Phi_h_minus_binxQ2_%d", i+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) | spin - | pion+; #Phi_{h} [Rad]; count", i+1), 50, -TMath::Pi(), TMath::Pi());
    }
    for(int i = 0; i < nbin_zPt; i++){
        hist_zvsPt_binned[i] = new TH2D(Form("hist_zvsPt_bin%d", i+1), Form("z vs P_{hT} for bin (%d, z-P_{hT}) | pion+; z; P_{hT} [GeV^]", i+1), 120, 0.2, 1.0, 120, 0, 1.4);
        hist_Phi_h_plus_binzPt[i] = new TH1D(Form("hist_Phi_h_plus_binzPt_%d", i+1), Form("#Phi_{h} for bin (%d, z-P_{hT}) | spin +| pion+; #Phi_{h} [Rad]; count", i+1), 50, -TMath::Pi(), TMath::Pi());
        hist_Phi_h_minus_binzPt[i] = new TH1D(Form("hist_Phi_h_minus_binzPt_%d", i+1), Form("#Phi_{h} for bin (%d, z-P_{hT}) | spin - | pion+; #Phi_{h} [Rad]; count", i+1), 50, -TMath::Pi(), TMath::Pi());
    }

    //
    double bin_z_plot[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    double bin_Pt_plot[] = {0, 0.25, 0.5, 0.8, 1.4};
    vector<TH2D*> hist_AUL_sin_xQ2_zPt(nbin_xQ2);
    vector<TH2D*> hist_AUL_sin2_xQ2_zPt(nbin_xQ2);
    vector<vector<TH1D*>> hist_Phi_h_plus_reco(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_minus_reco(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_plus_mc(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_minus_mc(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_plus_true(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_minus_true(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_plus_preIDcut(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_minus_preIDcut(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_reco(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    for (int ix = 0; ix < nbin_xQ2; ix++){
        hist_AUL_sin_xQ2_zPt[ix] = new TH2D(Form("hist_AUL_sin_xQ2_%d", ix+1), Form("AUL_sinPhi error bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
        hist_AUL_sin2_xQ2_zPt[ix] = new TH2D(Form("hist_AUL_sin2_xQ2_%d", ix+1), Form("AUL_sin2Phi error bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
        for (int iz = 0; iz < nbin_zPt; iz++){
            hist_Phi_h_plus_reco[ix][iz] = new TH1D(Form("hist_Phi_h_plus_reco_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H+ | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_minus_reco[ix][iz] = new TH1D(Form("hist_Phi_h_minus_reco_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H- | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_plus_mc[ix][iz] = new TH1D(Form("hist_Phi_h_plus_mc_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H+ | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_minus_mc[ix][iz] = new TH1D(Form("hist_Phi_h_minus_mc_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H- | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_plus_true[ix][iz] = new TH1D(Form("hist_Phi_h_plus_true_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H+ | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_minus_true[ix][iz] = new TH1D(Form("hist_Phi_h_minus_true_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H- | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_plus_preIDcut[ix][iz] = new TH1D(Form("hist_Phi_h_plus_preIDcut_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H+ | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_minus_preIDcut[ix][iz] = new TH1D(Form("hist_Phi_h_minus_preIDcut_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H- | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_reco[ix][iz] = new TH1D(Form("hist_Phi_h_reco_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | H+ | K+; #Phi_{h} [Rad]; count", ix+1, iz+1), 18, -TMath::Pi(), TMath::Pi());
        }
    }
    

    vector<vector<vector<double>>> vec_kaonp_xB_4d_mc(nbin_xQ2, vector<vector<double>> (nbin_zPt));

    TH1D kp_mc_mom ("_MC_momentum", " MC momentum | K+ ; mom [GeV]; count", 150, 1, 8);
    TH1D kp_mom ("_momentum", " momentum | K+ ; mom [GeV]; count", 150, 1, 8);
    TH2D kp_mc_mom_vs_theta ("_MC_mom_vs_theta", " MC momentum vs theta | K+ ; mom [GeV]; theta [deg]", 150, 1, 8, 150, 5, 25);
    TH2D kp_mom_vs_theta ("_mom_vs_theta", " MC momentum vs theta | K+ ; mom [GeV]; theta [deg]", 150, 1, 8, 150, 5, 25);

    TH2D kp_binmig_xB ("bin_migration_xB", "bin migration matrix for xB | K+ ; x_{B} (RECO); x_{B} (MC)", nbin, 0, 0.8, nbin, 0, 0.8);
    TH2D kp_binmig_Q2 ("bin_migration_Q2", "bin migration matrix for Q^{2} | K+ ; Q^{2} (RECO); Q^{2} (MC)", nbin, 1, 10, nbin, 1, 10);
    TH2D kp_binmig_z ("bin_migration_z", "bin migration matrix for z | K+ ; z (RECO); z (MC)", nbin, 0.2, 1.0, nbin, 0.2, 1.0);
    TH2D kp_binmig_Pt ("bin_migration_Pt", "bin migration matrix for P_{hT} | K+ ; P_{hT} (RECO); P_{hT} (MC)", nbin, 0, 1.4, nbin, 0, 1.4);
    TH2D kp_bin_xQ2 ("bin_mig_xQ2", "bin migration matrix for xB | K+ ; bin xQ2 (MC); bin xQ2 (RECO)", 14, 1, 15, 14, 1, 15);
    TH2D kp_bin_zPt ("bin_mig_zPt", "bin migration matrix for z | K+ ; bin zPt (MC); bin zPt (RECO)", 24, 1, 25, 24, 1, 25);



    vector<vector<vector<double>>> vec_Ass_true(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_pos_true(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_neg_true(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_pos_reco(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_neg_reco(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_pos_mc(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_neg_mc(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_pos_preIDcut(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_neg_preIDcut(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_pos_raw(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_neg_raw(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_pos_all_ID(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_neg_all_ID(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_pos_all_preIDcut(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    vector<vector<vector<double>>> vec_helicity_neg_all_preIDcut(nbin_xQ2, vector<vector<double>> (nbin_zPt));
    // MC
    Long64_t nEntries_mc = MC_chainKaonP.GetEntries();
    for (Long64_t i = 0; i < nEntries_mc; i++) {
        MC_chainKaonP.GetEntry(i);
        if(std::fabs(kaonp_Phi_mc) < 2.5 && std::fabs(kaonp_Phi_mc) > 0.65 && kaonp_Ph_mc >= 3) continue; 
        double index_xQ2 = getBinIndex_xQ2(kaonp_xB_mc, kaonp_Q2_mc); 
        double index_zPt = getBinIndex_zPt(kaonp_z_mc, kaonp_PhT_mc);
        kp_mc_mom.Fill(kaonp_Ph_mc);
        double kp_mc_theta_deg = kaonp_Theta_mc * 180.0 / TMath::Pi();
        kp_mc_mom_vs_theta.Fill(kaonp_Ph_mc, kp_mc_theta_deg);
        if(index_xQ2 >= 0 && index_zPt >= 0) {
            vec_kaonp_xB_4d_mc[index_xQ2-1][index_zPt-1].push_back(kaonp_xB_mc);
            if(helicity_mc == 1) {
                vec_helicity_pos_mc[index_xQ2-1][index_zPt-1].push_back(helicity_mc);
                hist_Phi_h_plus_mc[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h_mc);

            } else if(helicity_mc == -1){
                vec_helicity_neg_mc[index_xQ2-1][index_zPt-1].push_back(helicity_mc);
                hist_Phi_h_minus_mc[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h_mc);
            }
        }
    }

    // RECO
    // Kaon+
    Long64_t nEntries_kp = chainKaonP.GetEntries();
    for (Long64_t i = 0; i < nEntries_kp; i++) {
        // RICH cut
        //if (electron_ass != 321) continue;
        chainKaonP.GetEntry(i);
        if(std::fabs(kaonp_Phi) < 2.5 && std::fabs(kaonp_Phi) > 0.65 && kaonp_Ph >= 3) continue;
        // omega all
        double mc_index_xQ2_preCut_allID = getBinIndex_xQ2(kaonp_xB_recoMC, kaonp_Q2_recoMC);
        double mc_index_zPt_preCut_allID = getBinIndex_zPt(kaonp_z_recoMC, kaonp_Pt_recoMC);
        if(mc_index_xQ2_preCut_allID >= 0 && mc_index_zPt_preCut_allID >= 0) {
            if(helicity > 0) vec_helicity_pos_all_preIDcut[mc_index_xQ2_preCut_allID-1][mc_index_zPt_preCut_allID-1].push_back(helicity);
            else vec_helicity_neg_all_preIDcut[mc_index_xQ2_preCut_allID-1][mc_index_zPt_preCut_allID-1].push_back(helicity);
        }
        if (passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph, kaonp_Phi)){
            double mc_index_xQ2_all_ID = getBinIndex_xQ2(kaonp_xB_recoMC, kaonp_Q2_recoMC);
            double mc_index_zPt_all_ID = getBinIndex_zPt(kaonp_z_recoMC, kaonp_Pt_recoMC);
            if(mc_index_xQ2_all_ID >= 0 && mc_index_zPt_all_ID >= 0) {
                if(helicity > 0) vec_helicity_pos_all_ID[mc_index_xQ2_all_ID-1][mc_index_zPt_all_ID-1].push_back(helicity);
                else vec_helicity_neg_all_ID[mc_index_xQ2_all_ID-1][mc_index_zPt_all_ID-1].push_back(helicity);
            }
        }

        if(mc_PID != 321) continue;
        // raw for cavariance
        double mc_index_xQ2_raw = getBinIndex_xQ2(kaonp_xB_recoMC, kaonp_Q2_recoMC);
        double mc_index_zPt_raw = getBinIndex_zPt(kaonp_z_recoMC, kaonp_Pt_recoMC);
        if(mc_index_xQ2_raw >= 0 && mc_index_zPt_raw >= 0) {
            if(helicity > 0) vec_helicity_pos_raw[mc_index_xQ2_raw-1][mc_index_zPt_raw-1].push_back(helicity);
            else vec_helicity_neg_raw[mc_index_xQ2_raw-1][mc_index_zPt_raw-1].push_back(helicity);
        }

        double radius = sqrt(kaonp_px*kaonp_px + kaonp_py*kaonp_py);
        //
        double mc_index_xQ2_preCut = getBinIndex_xQ2(kaonp_xB_recoMC, kaonp_Q2_recoMC);
        double mc_index_zPt_preCut = getBinIndex_zPt(kaonp_z_recoMC, kaonp_Pt_recoMC);
        if(mc_index_xQ2_preCut >= 0 && mc_index_zPt_preCut >= 0) {
            if(helicity > 0){
                vec_helicity_pos_preIDcut[mc_index_xQ2_preCut-1][mc_index_zPt_preCut-1].push_back(helicity);
                hist_Phi_h_plus_preIDcut[mc_index_xQ2_preCut-1][mc_index_zPt_preCut-1]->Fill(kaonp_Phi_h_recoMC);
            } else {
                vec_helicity_neg_preIDcut[mc_index_xQ2_preCut-1][mc_index_zPt_preCut-1].push_back(helicity);
                hist_Phi_h_minus_preIDcut[mc_index_xQ2_preCut-1][mc_index_zPt_preCut-1]->Fill(kaonp_Phi_h_recoMC);
            }
        }

        // ID
        if(1==1){
            kp_Mx_tof.Fill(kaonp_Mx);
            // RICH selection passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph)
            if(passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph, kaonp_Phi)){
            //if(1 == 1){
                double ph2 = kaonp_Ph*kaonp_Ph;
                double beta2 = kaonp_beta*kaonp_beta;
                kaonp_m = sqrt((ph2)*(1-beta2)/beta2);
                double index_xQ2 = getBinIndex_xQ2(kaonp_xB, kaonp_Q2); 
                double index_zPt = getBinIndex_zPt(kaonp_z, kaonp_PhT);
                double mc_index_xQ2 = getBinIndex_xQ2(kaonp_xB_recoMC, kaonp_Q2_recoMC);
                double mc_index_zPt = getBinIndex_zPt(kaonp_z_recoMC, kaonp_Pt_recoMC);
                proton_spin = +1;
                if(kaonp_Pol < 0) proton_spin = -1;
                kaonp_sintheta = kaonp_gamma*sqrt((1-kaonp_y-0.25*kaonp_y*kaonp_y*kaonp_gamma*kaonp_gamma)/(1+kaonp_y*kaonp_y));
                C_LSA_TSA = kaonp_sintheta * (sqrt(2*kaonp_epsilon*(1+kaonp_epsilon))/kaonp_epsilon);
                //if (index_xQ2 != 15) continue;
                //
                if(mc_index_xQ2 >= 0 && mc_index_zPt >= 0) {
                    if(helicity > 0) vec_helicity_pos_true[mc_index_xQ2-1][mc_index_zPt-1].push_back(helicity);
                    else vec_helicity_neg_true[mc_index_xQ2-1][mc_index_zPt-1].push_back(helicity);
                }
                if(index_xQ2 >= 0) {
                    hist_xBvsQ2_binned[index_xQ2-1]->Fill(kaonp_xB, kaonp_Q2);
                    vec_kaonp_phih_2d[index_xQ2-1].push_back(kaonp_Phi_h);
                    vec_helicity_2d[index_xQ2-1].push_back(helicity);
                    vec_kaonp_2phih_2d[index_xQ2-1].push_back(2*kaonp_Phi_h);
                    vec_kaonp_z_2d[index_xQ2-1].push_back(kaonp_z);
                    vec_kaonp_Pt_2d[index_xQ2-1].push_back(kaonp_PhT);
                    vec_kaonp_xB_2d[index_xQ2-1].push_back(kaonp_xB);
                    vec_kaonp_Q2_2d[index_xQ2-1].push_back(kaonp_Q2);
                    vec_kaonp_y_2d[index_xQ2-1].push_back(kaonp_y);
                    vec_kaonp_pol_2d[index_xQ2-1].push_back(kaonp_Pol);
                    vec_kaonp_eps_2d[index_xQ2-1].push_back(kaonp_epsilon);
                    vec_kaonp_spin_2d[index_xQ2-1].push_back(proton_spin);
                    vec_kaonp_sintheta_2d[index_xQ2-1].push_back(kaonp_sintheta);
                    if(proton_spin == 1) hist_Phi_h_plus_binxQ2[index_xQ2-1]->Fill(kaonp_Phi_h);
                    else if(proton_spin == -1) hist_Phi_h_minus_binxQ2[index_xQ2-1]->Fill(kaonp_Phi_h);
                }
                if(index_zPt >= 0){ 
                    hist_zvsPt_binned[index_zPt-1]->Fill(kaonp_z, kaonp_PhT);
                    vec_kaonp_phih_2d_zPt[index_zPt-1].push_back(kaonp_Phi_h);
                    vec_helicity_2d_zPt[index_zPt-1].push_back(helicity);
                    vec_kaonp_2phih_2d_zPt[index_zPt-1].push_back(2*kaonp_Phi_h);
                    vec_kaonp_z_2d_zPt[index_zPt-1].push_back(kaonp_z);
                    vec_kaonp_Pt_2d_zPt[index_zPt-1].push_back(kaonp_PhT);
                    vec_kaonp_xB_2d_zPt[index_zPt-1].push_back(kaonp_xB);
                    vec_kaonp_Q2_2d_zPt[index_zPt-1].push_back(kaonp_Q2);
                    vec_kaonp_y_2d_zPt[index_zPt-1].push_back(kaonp_y);
                    vec_kaonp_pol_2d_zPt[index_zPt-1].push_back(kaonp_Pol);
                    vec_kaonp_eps_2d_zPt[index_zPt-1].push_back(kaonp_epsilon);
                    //if (index_zPt == 16) cout << kaonp_epsilon << endl;
                    vec_kaonp_spin_2d_zPt[index_zPt-1].push_back(proton_spin);
                    vec_kaonp_sintheta_2d_zPt[index_zPt-1].push_back(kaonp_sintheta);
                    if(proton_spin == 1) hist_Phi_h_plus_binzPt[index_zPt-1]->Fill(kaonp_Phi_h);
                    else if(proton_spin == -1) hist_Phi_h_minus_binzPt[index_zPt-1]->Fill(kaonp_Phi_h);
                    if(index_xQ2 >= 0 && index_zPt >= 0){
                        hist_Phi_h_reco[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h);
                        vec_kaonp_phih_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_Phi_h);
                        vec_helicity_4d[index_xQ2-1][index_zPt-1].push_back(helicity);
                        vec_kaonp_2phih_4d[index_xQ2-1][index_zPt-1].push_back(2*kaonp_Phi_h);
                        vec_kaonp_z_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_z);
                        vec_kaonp_Pt_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_PhT);
                        vec_kaonp_xB_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_xB);
                        vec_kaonp_Q2_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_Q2);
                        vec_kaonp_y_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_y);
                        vec_kaonp_pol_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_Pol);
                        vec_kaonp_eps_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_epsilon);
                        vec_kaonp_spin_4d[index_xQ2-1][index_zPt-1].push_back(proton_spin);
                        vec_kaonp_sintheta_4d[index_xQ2-1][index_zPt-1].push_back(kaonp_sintheta);
                        if(helicity > 0){
                            vec_helicity_pos_reco[index_xQ2-1][index_zPt-1].push_back(helicity);
                            hist_Phi_h_plus_reco[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h);
                            hist_Phi_h_plus_true[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h_recoMC);
                        } else {
                            vec_helicity_neg_reco[index_xQ2-1][index_zPt-1].push_back(helicity);
                            hist_Phi_h_minus_reco[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h);
                            hist_Phi_h_minus_true[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h_recoMC);
                        }
                    }
                }
                //
                kp_evnt_chi2.Fill(kaonp_chi2pid);
                kp_m.Fill(kaonp_bestMass);
                kp_mc_PID.Fill(mc_PID);
                //
                kp_mom.Fill(kaonp_Ph);
                double kp_theta_deg = kaonp_Theta * 180.0 / TMath::Pi();
                kp_mom_vs_theta.Fill(kaonp_Ph, kp_theta_deg);
                //
                kp_MomVsMass_Rich.Fill(kaonp_Ph, kaonp_bestMass);
                double beta_th = kaonp_Ph/(sqrt(ph2+0.01948816));
                double delta_beta = kaonp_beta - beta_th;
                kp_deltaB.Fill(delta_beta);
                kp_Mx.Fill(kaonp_Mx);
                kp_Mx_rich.Fill(kaonp_Mx);
                target_pol.Fill(kaonp_Pol);
                target_spin.Fill(proton_spin);
                beam_helicity.Fill(helicity);
                if(proton_spin == 1) kp_phi_plus.Fill(kaonp_Phi_h);
                if(proton_spin == -1) kp_phi_minus.Fill(kaonp_Phi_h);
                // Mom
                kp_MomVsPhT.Fill(kaonp_PhT, kaonp_Ph);
                kp_MomVsEta.Fill(kaonp_eta, kaonp_Ph);
                kp_MomVsMx.Fill(kaonp_Ph, kaonp_Mx);
                kp_MomVsPhi_h.Fill(kaonp_Phi_h, kaonp_Ph);
                kp_MomVsPhi_lab.Fill(kaonp_Phi, kaonp_Ph);
                kp_MomVsTheta.Fill(kaonp_Ph, kaonp_Theta);
                kp_MomVsXb.Fill(kaonp_xB, kaonp_Ph);
                kp_MomVsXf.Fill(kaonp_xF, kaonp_Ph);
                kp_MomVsY.Fill(kaonp_y, kaonp_Ph);
                kp_MomVsZ.Fill(kaonp_z, kaonp_Ph);
                kp_MomVsBeta.Fill(kaonp_Ph, kaonp_beta);
                kp_MomVsCh.Fill(kaonp_Ph, kaonp_rich_ch);
                kp_MomvsSinTheta.Fill(kaonp_Ph, kaonp_sintheta);
                kp_MomVsEps.Fill(kaonp_epsilon, kaonp_Ph);
                // Q2
                kp_Q2VsEta.Fill(kaonp_eta, kaonp_Q2);
                kp_Q2VsMom.Fill(kaonp_Ph, kaonp_Q2);
                kp_Q2VsPhi_h.Fill(kaonp_Phi_h, kaonp_Q2);
                kp_Q2VsPhT.Fill(kaonp_PhT, kaonp_Q2);
                kp_Q2VsXb.Fill(kaonp_xB, kaonp_Q2);
                kp_Q2VsXf.Fill(kaonp_xF, kaonp_Q2);
                kp_Q2VsY.Fill(kaonp_y, kaonp_Q2);
                kp_Q2VsZ.Fill(kaonp_z, kaonp_Q2);
                kp_Q2vsSinTheta.Fill(kaonp_Q2, kaonp_sintheta);
                kp_Q2VsEps.Fill(kaonp_epsilon, kaonp_Q2);
                kp_Q2VsMx.Fill(kaonp_Mx, kaonp_Q2);
                // PhT
                kp_PhTvsEta.Fill(kaonp_eta, kaonp_PhT);
                kp_PhTvsXb.Fill(kaonp_xB, kaonp_PhT);
                kp_PhTvsZ.Fill(kaonp_z, kaonp_PhT);
                kp_PhTvsPhi_h.Fill(kaonp_Phi_h, kaonp_PhT);
                kp_PhTvsSinTheta.Fill(kaonp_PhT, kaonp_sintheta);
                kp_PhTVsEps.Fill(kaonp_epsilon, kaonp_PhT);
                kp_PhTVsMx.Fill(kaonp_PhT, kaonp_Mx);
                // z
                kp_zVsEta.Fill(kaonp_eta, kaonp_z);
                kp_zVsPhi_h.Fill(kaonp_Phi_h, kaonp_z);
                kp_zVsXb.Fill(kaonp_xB, kaonp_z);
                kp_zVsXf.Fill(kaonp_xF, kaonp_z);
                kp_zvsSinTheta.Fill(kaonp_z, kaonp_sintheta);
                kp_zvsCosPhi.Fill(kaonp_z, cos(kaonp_Phi_h));
                kp_zVsEps.Fill(kaonp_epsilon, kaonp_z);
                kp_zVsEps2.Fill(kaonp_z, kaonp_epsilon);
                kp_zVsMx.Fill(kaonp_z, kaonp_Mx);
                //
                kp_xBvsY.Fill(kaonp_xB, kaonp_y);
                kp_xBvsEta.Fill(kaonp_eta, kaonp_xB);
                kp_xBvsSinTheta.Fill(kaonp_xB, kaonp_sintheta);
                kp_xBvsCosPhi.Fill(kaonp_xB, cos(kaonp_Phi_h));
                kp_xBVsEps.Fill(kaonp_epsilon, kaonp_xB);
                kp_xBVsMx.Fill(kaonp_xB, kaonp_Mx);
                //
                kp_ThetaVsPhi_h.Fill(kaonp_Phi_h, kaonp_Theta);
                kp_ThetaVsPhi_Lab.Fill(kaonp_Phi, kaonp_Theta);
                kp_PxVsPy.Fill(kaonp_px, kaonp_py);
                el_PxVsPy.Fill(electron_px, electron_py);
                treeKaonP.Fill();

                // bin migration
                kp_binmig_xB.Fill(kaonp_xB, kaonp_xB_recoMC);
                kp_binmig_Q2.Fill(kaonp_Q2, kaonp_Q2_recoMC);
                kp_binmig_z.Fill(kaonp_z, kaonp_z_recoMC);
                kp_binmig_Pt.Fill(kaonp_PhT, kaonp_Pt_recoMC);
                kp_bin_xQ2.Fill(mc_index_xQ2, index_xQ2);
                kp_bin_zPt.Fill(mc_index_zPt, index_zPt);

                // plot del RICH\
                if(kaonp_rich_nTot > 0){
                    kp_rich_pmt_xy_4.Fill(kaonp_rich_tr1_x, kaonp_rich_tr1_y);
                    kp_rich_pmt_xy_1.Fill(kaonp_rich_tr1_x, kaonp_rich_tr1_y);
                    kp_rich_aerogel_xy_4.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_y);
                    kp_rich_aerogel_xy_4.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_y);
                    kp_rich_aerogel_xy_4.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_y);
                    kp_rich_aerogel_xy_1.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_y);
                    kp_rich_aerogel_xy_1.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_y);
                    kp_rich_aerogel_xy_1.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_y);
                    kp_rich_aerogel_l1_xy_4.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_y);
                    kp_rich_aerogel_l2_xy_4.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_y);
                    kp_rich_aerogel_l3_xy_4.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_y);
                    kp_rich_aerogel_l1_xy_1.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_y);
                    kp_rich_aerogel_l2_xy_1.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_y);
                    kp_rich_aerogel_l3_xy_1.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_y);
                }

            }
    }
    
    
    
    treeKaonP.Write();

    //
    vector<double> x_points;
    vector<double> y_points;

    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {

            const auto& z_vec  = vec_kaonp_z_4d[ix][iz];
            const auto& pt_vec = vec_kaonp_Pt_4d[ix][iz];
            const auto& q2_vec = vec_kaonp_Q2_4d[ix][iz];

            if (z_vec.empty()) continue;  // evita bin vuoti

            double z_mean  = MeanVect(z_vec);
            double pt_mean = MeanVect(pt_vec);
            double q2_mean = MeanVect(q2_vec);
            double q_mean  = sqrt(q2_mean);

            if (z_mean <= 0 || pt_mean <= 0) continue;

            double ratio = pt_mean / (z_mean * q_mean);

            // Asse x: ad esempio <xB> o indice del bin
            double xB_mean = MeanVect(vec_kaonp_xB_4d[ix][iz]);

            x_points.push_back(xB_mean);
            y_points.push_back(ratio);
        }
    }

    TGraph* gr = new TGraph(x_points.size(),x_points.data(),y_points.data());
    gr->SetTitle("TMD factorization parameter; <x_{B}>; <P_{hT}>/(<z><Q>)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1);
    gr->SetMarkerColor(kAzure-4);
    //gr->SetLineWidth(2);
    TCanvas* c_ratio = new TCanvas("c_ratio","P_{hT}/(zQ)",800,600);
    gr->Draw("AP");
    c_ratio->Write();



    kp_binmig_xB.Write();
    kp_binmig_Q2.Write();
    kp_binmig_z.Write();
    kp_binmig_Pt.Write();
    kp_bin_xQ2.Write();
    kp_bin_zPt.Write();
    //
    kp_evnt_chi2.Write();
    kp_mc_PID.Write();
    kp_m.Write();
    kp_deltaB.Write();
    kp_Mx.Write();
    kp_Mx_tof.Write();
    kp_Mx_rich.Write();
    TCanvas* c_kp_Mx = new TCanvas("_missing_mass_overlay","Missing Mass comparison",800,600);
    kp_Mx_tof.SetLineColor(kRed-3); kp_Mx_tof.SetTitle("Missing Mass; M_{X} [GeV]; Counts"); kp_Mx_tof.SetStats(0);
    kp_Mx_rich.SetLineColor(kBlue-3);
    kp_Mx_tof.Draw();
    kp_Mx_rich.Draw("SAME");
    auto leg = new TLegend(0.15,0.78,0.28,0.88);
    leg->AddEntry(&kp_Mx_tof,"CLAS","l");
    leg->AddEntry(&kp_Mx_rich,"RICH","l");
    leg->Draw();
    c_kp_Mx->Write();
    target_pol.Write();
    target_spin.Write();
    beam_helicity.Write();
    kp_phi_plus.Write();
    kp_phi_minus.Write();

    //

    std::vector<TH2D*> hists_kp = {
        &kp_MomVsPhT, &kp_MomVsXb, &kp_MomVsXf, &kp_MomVsZ, &kp_MomVsY, &kp_MomVsEta,
        &kp_MomVsTheta, &kp_MomVsPhi_h, &kp_MomVsPhi_lab, &kp_MomVsMx, &kp_MomVsBeta, &kp_MomVsCh, &kp_MomVsMass_Rich, &kp_MomvsSinTheta, &kp_MomVsEps,
        &kp_Q2VsXb, &kp_Q2VsXf, &kp_Q2VsMom, &kp_Q2VsPhT, &kp_Q2VsZ, &kp_Q2VsY, &kp_Q2VsEta, &kp_Q2VsPhi_h, &kp_Q2vsSinTheta, &kp_Q2VsEps, &kp_Q2VsMx,
        &kp_PhTvsZ, &kp_PhTvsXb, &kp_PhTvsEta, &kp_PhTvsPhi_h, &kp_PhTvsSinTheta, &kp_PhTVsEps, &kp_PhTVsMx,
        &kp_zVsXb, &kp_zVsXf, &kp_zVsEta, &kp_zVsPhi_h, &kp_zvsSinTheta, & kp_zvsCosPhi, &kp_zVsEps, &kp_zVsEps2, &kp_zVsMx,
        &kp_xBvsY, &kp_xBvsEta, &kp_xBvsSinTheta, &kp_xBvsCosPhi, &kp_xBVsEps, &kp_xBVsMx,
        &kp_ThetaVsPhi_h, &kp_ThetaVsPhi_Lab, &kp_PxVsPy, &el_PxVsPy,
        &kp_rich_pmt_xy_4, &kp_rich_pmt_xy_1, &kp_rich_aerogel_xy_4, &kp_rich_aerogel_xy_1, 
        &kp_rich_aerogel_l1_xy_4, &kp_rich_aerogel_l2_xy_4,  &kp_rich_aerogel_l3_xy_4,
        &kp_rich_aerogel_l1_xy_1, &kp_rich_aerogel_l2_xy_1, &kp_rich_aerogel_l3_xy_1,
    };

    // list to set the statbox2 on the canvas
    std::set<std::string> id_box2 = { "_MomVsXf", "_MomVsZ", "_MomVsY", "_MomVsEta", "_Q2VsXb", "_Q2VsY", "_zVsEta", "_zVsPhi_h",
    "_MomVsMass_RICH", "_rich_MomVsXf", "_rich_MomVsZ", "_rich_MomVsY", "_rich_MomVsEta", "_rich_MomVsPhi_h", "_rich_Q2VsXb", "_rich_Q2VsY",
    "_rich_Q2VsEta", "_rich_Q2VsPhi_h", "_rich_PhTvsEta", "_rich_zVsEta", "_rich_zVsPhi_h"};
    for (size_t i = 0; i < hists_kp.size(); ++i) {
        const char* histName = hists_kp[i]->GetName();
        //std::string cname_kp = std::string("c_kp") + histName;
        TCanvas *c = new TCanvas(histName, histName, 800, 600);
        c->SetLogz();   
        // if(i == ...) c->SetLogx();
        // if(i == ...) c->SetLogy();
        hists_kp[i]->Draw("COLZ");
        gPad->Update();
        if (id_box2.count(histName)) SetStatsBox2(hists_kp[i]);
        else SetStatsBox(hists_kp[i]);

        c->Write();
    }

    dir_xQ2->cd();
    TCanvas *c_bin_xQ2 = new TCanvas("Q2_vs_xB_Bin", "Q^{2} vs x_{B} bin", 800, 700);
    c_bin_xQ2->SetLogz();
    kp_Q2VsXb.Draw("COLZ");

    vector<TPolyLine*> rectangles;
    vector<TText*> labels;
    int bin_index = 1;
    // --- Loop sui bin ---
    for (int ix = 0; ix < 4; ++ix) {
        for (size_t iq = 0; iq < binning_xB_for_Q2[ix].size(); ++iq) {
            double Q2[5] = {bin_Q2[ix][0], bin_Q2[ix][1], bin_Q2[ix][1], bin_Q2[ix][0], bin_Q2[ix][0]};
            double xB[5] = {binning_xB_for_Q2[ix][iq][0], binning_xB_for_Q2[ix][iq][0], binning_xB_for_Q2[ix][iq][1], binning_xB_for_Q2[ix][iq][1], binning_xB_for_Q2[ix][iq][0]};

            TPolyLine *rect = new TPolyLine(5, xB, Q2);
            rect->SetLineColor(kBlack);
            rect->SetLineWidth(2);
            rect->Draw("same");
            rectangles.push_back(rect);

            // posizione indici
            double Q2_center = 0.5 *(bin_Q2[ix][0] + bin_Q2[ix][1]);
            double xB_center = 0.5 * (binning_xB_for_Q2[ix][iq][0] + binning_xB_for_Q2[ix][iq][1]);

            TText *label = new TText(xB_center, Q2_center, Form("%d", bin_index++));
            label->SetTextAlign(22);
            label->SetTextSize(0.03);
            label->SetTextColor(kRed+1); //forse rosso si vede di più
            label->Draw("same");
            labels.push_back(label);
        }
    }

    c_bin_xQ2->Update();
    c_bin_xQ2->Write();


    for(int i = 0; i < nbin_xQ2; i++) hist_xBvsQ2_binned[i]->Write();


    dir_zPt->cd();
    TCanvas *c_bin_zPt = new TCanvas("z_vs_Pt_Bin", "z vs P_{hT} bin", 800, 700);
    c_bin_zPt->SetLogz();
    kp_PhTvsZ.Draw("COLZ");

    vector<TPolyLine*> gridLines_zp;
    vector<TText*> labels_zp;

    int bin_index_zp = 1;

    // Loop generale
    for (int iz = 0; iz < 4; ++iz) {
        for (size_t ip = 0; ip < binning_z_for_Pt[iz].size(); ++ip) {

            double Pt[5] = {bin_Pt[iz][0], bin_Pt[iz][1], bin_Pt[iz][1], bin_Pt[iz][0], bin_Pt[iz][0]};
            double z[5] = {binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][0]};

            // Rettangolo
            TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
            rect_zp->SetLineWidth(2);
            rect_zp->SetLineColor(kBlack);
            rect_zp->Draw("same");
            gridLines_zp.push_back(rect_zp);

            // Centro geometrico (utile se log-scale in Pt)
            double Pt_center = 0.5 * (bin_Pt[iz][0] + bin_Pt[iz][1]);
            double z_center = 0.5 * (binning_z_for_Pt[iz][ip][0] + binning_z_for_Pt[iz][ip][1]);

            TText *label = new TText(z_center, Pt_center, Form("%d", bin_index_zp++));
            label->SetTextAlign(22);
            label->SetTextSize(0.03);
            label->SetTextColor(kRed+1);
            label->Draw("same");
            labels_zp.push_back(label);
        }
    }

    c_bin_zPt->Update();
    c_bin_zPt->Write();

    for(int i = 0; i < nbin_zPt; i++) hist_zvsPt_binned[i]->Write();



    // ========================================== EFFICIENCY STUDIES ====================================================

    auto quadCorrect = [](double sys, double sys_err) -> double {
        double diff2 = sys * sys - sys_err * sys_err;
        return (diff2 > 0.0) ? sqrt(diff2) : 0.0;
    };
    //
    dir_eff->cd();
    //
    kp_mc_mom.Write();
    kp_mom.Write();
    TH1D* kp_mom_ratio = (TH1D*)kp_mom.Clone("kp_mom_ratio");
    kp_mom_ratio->Divide(&kp_mom, &kp_mc_mom, 1.0, 1.0, "B");
    //kp_mom_ratio->Divide(&kp_mc_mom);
    kp_mom_ratio->SetTitle("Kaon+ efficiency ; Momentum [GeV]; Ratio (Data/MC)");
    kp_mom_ratio->GetYaxis()->SetRangeUser(0.0, 1.0);
    kp_mom_ratio->Write();

    kp_mc_mom_vs_theta.Write();
    kp_mom_vs_theta.Write();
    TH2D* kp_mom_theta_ratio = (TH2D*)kp_mom_vs_theta.Clone("kp_mom_theta_ratio");
    kp_mom_theta_ratio->Divide(&kp_mom_vs_theta, &kp_mc_mom_vs_theta, 1.0, 1.0, "B");
    kp_mom_theta_ratio->SetTitle("Kaon+ efficiency ; Momentum [GeV]; Theta [deg]");
    //kp_mom_theta_ratio->Write("COLZ");
    TCanvas* c_kp_mom_theta_ratio = new TCanvas("c_kp_mom_theta_ratio","Kaon+ efficiency vs Momentum and Theta",800,600);
    kp_mom_theta_ratio->SetTitle("Kaon+ efficiency ; Momentum [GeV]; Theta [deg]");
    kp_mom_theta_ratio->SetStats(0);
    kp_mom_theta_ratio->SetMinimum(0.0);
    kp_mom_theta_ratio->SetMaximum(1.0);
    kp_mom_theta_ratio->Draw("COLZ");
    c_kp_mom_theta_ratio->Write();


    vector<TH2D*> hist_efficiency_xQ2_zPt(nbin_xQ2);
    vector<vector<double>> efficiency_values(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> efficiency_err(nbin_xQ2, vector<double> (nbin_zPt));
    for (int ix = 0; ix < nbin_xQ2; ix++){
      hist_efficiency_xQ2_zPt[ix] = new TH2D(Form("hist_efficiency4D_xQ2_%d", ix+1), Form("Kaon+ efficiency as P_{hT} vs z for bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
    }

    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {
            double den = vec_kaonp_xB_4d_mc[ix][iz].size();
            double num = vec_kaonp_xB_4d[ix][iz].size();
            double eff = (den > 0) ? num / den : 0.0;
            if (eff > 1) eff = 1.0;
            efficiency_values[ix][iz] = eff;
            // Calcolo dell'errore binomiale
            double eff_err = 0.0;
            if (den > 0) {
                eff_err = sqrt(eff * (1 - eff) / den);
            }
            efficiency_err[ix][iz] = eff_err;
            //cout << "Bin (" << ix+1 << ", " << iz+1 << "): Efficiency = " << eff << " (" << num << "/" << den << ")" << endl;
            if(iz < 7) hist_efficiency_xQ2_zPt[ix]->SetBinContent(iz+1, 1, eff);
            else if (iz < 14) hist_efficiency_xQ2_zPt[ix]->SetBinContent(iz-7+1, 2, eff);
            else if (iz < 21) hist_efficiency_xQ2_zPt[ix]->SetBinContent(iz-14+1, 3, eff);
            else if (iz == 21){
                hist_efficiency_xQ2_zPt[ix]->SetBinContent(1, 4, eff);
                hist_efficiency_xQ2_zPt[ix]->SetBinContent(2, 4, eff);
                hist_efficiency_xQ2_zPt[ix]->SetBinContent(3, 4, eff);
            } else if (iz == 24) {
                hist_efficiency_xQ2_zPt[ix]->SetBinContent(6, 4, eff);
                hist_efficiency_xQ2_zPt[ix]->SetBinContent(7, 4, eff);
            }
            else hist_efficiency_xQ2_zPt[ix]->SetBinContent(iz-22+4, 4, eff);
        }

        TCanvas *c_bin_zPt = new TCanvas(Form("efficiency_zPt_xQ2_%d", ix+1), "z vs P_{hT} bin", 800, 700);
        //c_bin_zPt->SetLogz();
        hist_efficiency_xQ2_zPt[ix]->SetMinimum(0.0);
        hist_efficiency_xQ2_zPt[ix]->SetMaximum(1.0);
        hist_efficiency_xQ2_zPt[ix]->SetStats(0);
        hist_efficiency_xQ2_zPt[ix]->Draw("colz");

        vector<TPolyLine*> gridLines_zp;
        vector<TText*> labels_zp;

        int bin_index_zp = 1;

        // Loop generale
        for (int iz = 0; iz < 4; ++iz) {
            for (size_t ip = 0; ip < binning_z_for_Pt[iz].size(); ++ip) {

                double Pt[5] = {bin_Pt[iz][0], bin_Pt[iz][1], bin_Pt[iz][1], bin_Pt[iz][0], bin_Pt[iz][0]};
                double z[5] = {binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][0]};

                // Rettangolo
                TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
                rect_zp->SetLineWidth(2);
                rect_zp->SetLineColor(kBlack);
                rect_zp->Draw("same");
                gridLines_zp.push_back(rect_zp);

                // Centro geometrico (utile se log-scale in Pt)
                double Pt_center = 0.5 * (bin_Pt[iz][0] + bin_Pt[iz][1]);
                double z_center = 0.5 * (binning_z_for_Pt[iz][ip][0] + binning_z_for_Pt[iz][ip][1]);

                int binx = hist_efficiency_xQ2_zPt[ix]->GetXaxis()->FindBin(z_center);
                int biny = hist_efficiency_xQ2_zPt[ix]->GetYaxis()->FindBin(Pt_center);

                double eff_val = hist_efficiency_xQ2_zPt[ix]->GetBinContent(binx, biny);

                TText *label = new TText(z_center, Pt_center, Form("%.3f", eff_val));
                label->SetTextAlign(22);
                label->SetTextSize(0.03);
                //label->SetTextColor(kWhite);
                label->Draw("same");
                labels_zp.push_back(label);
            }
        }

        c_bin_zPt->Update();
        c_bin_zPt->Write();
        //hist_efficiency_xQ2_zPt[ix]->Write();
    }


    // systematic
    double r = 0.727688;
    double r_inv = 1.0 / r;
    vector<TH2D*> hist_eff_sys_4D(nbin_xQ2);
    vector<double> A_eff_sys_2d_xQ2(nbin_xQ2);
    vector<double> A_eff_sys_2d_zPt(nbin_zPt);
    vector<double> A_eff_sys_2d_xQ2_err(nbin_xQ2);
    vector<double> A_eff_sys_2d_zPt_err(nbin_zPt);
    vector<vector<double>> A_eff_sys_4d(nbin, vector<double> (nbin));
    vector<vector<double>> A_eff_sys_err_4d(nbin, vector<double> (nbin));
    TH2D* h_eff_sys = new TH2D("h_eff_sys", "Systematic; bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 1, nbin_xQ2+1, nbin_zPt, 1, nbin_zPt+1);
    TH2D* h_eff_sys_rel = new TH2D("h_eff_sys_rel", "RelativeSystematic; bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 0, nbin_xQ2, nbin_zPt, 0, nbin_zPt);
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        hist_eff_sys_4D[ix] = new TH2D(Form("hist_eff_sys_xQ2_%d", ix+1), Form("Efficiency systematic | bin:(%d, x_{B}-Q^{2}) | #sigma_{sys}^{K+} = |A_{true}^{K+} - A_{reco}^{K+}|; z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
        for (int iz = 0; iz < nbin_zPt; iz++) {
            double N_pos_true = vec_helicity_pos_preIDcut[ix][iz].size();
            double N_neg_true = vec_helicity_neg_preIDcut[ix][iz].size();
            double N_pos_reco = vec_helicity_pos_true[ix][iz].size();
            double N_neg_reco = vec_helicity_neg_true[ix][iz].size();

            double N_pos_raw = vec_helicity_pos_raw[ix][iz].size();       // "Omega"
            double N_neg_raw = vec_helicity_neg_raw[ix][iz].size();

            if(N_pos_true*N_neg_true*N_pos_reco*N_neg_reco == 0){ 
                cout << "Efficiency systematic: Bin (" << ix+1 << ", " << iz+1 << "): N_pos_true=" << N_pos_true << ", N_neg_true=" << N_neg_true 
                 << ", N_pos_reco=" << N_pos_reco << ", N_neg_reco=" << N_neg_reco << endl;
            }
            if ((N_pos_true + r_inv*N_neg_true) == 0 || (N_pos_reco + r_inv*N_neg_reco) == 0) continue;
            /*
            double A_true = (N_pos_true - r_inv*N_neg_true) / (N_pos_true + r_inv*N_neg_true);
            double A_reco = (N_pos_reco - r_inv*N_neg_reco) / (N_pos_reco + r_inv*N_neg_reco);
            double A_eff_sys = std::abs(A_true - A_reco);
            double sigma_A_true = sqrt( (1 - A_true*A_true) / (N_pos_true + r_inv*r_inv*N_neg_true) );
            double sigma_A_reco = sqrt( (1 - A_reco*A_reco) / (N_pos_reco + r_inv*r_inv*N_neg_reco) );
            double A_eff_sys_err = sqrt(sigma_A_true*sigma_A_true + sigma_A_reco*sigma_A_reco);
            double relative_sys = A_eff_sys / std::abs(A_true);
            double A_eff_sys_corr = quadCorrect(A_eff_sys, A_eff_sys_err);
            A_eff_sys_4d[ix][iz] = A_eff_sys_corr;  // here attention to reuse the non 'corr' one
            A_eff_sys_err_4d[ix][iz] = A_eff_sys_err;
            h_eff_sys->SetBinContent(ix+1, iz+1, A_eff_sys);
            h_eff_sys_rel->SetBinContent(ix+1, iz+1, relative_sys);
            */

            SysWithCovResult r_cov = ComputeSysWithCov(N_pos_true, N_neg_true,N_pos_reco, N_neg_reco,N_pos_raw, N_neg_raw,r_inv);
 
            double A_true = r_cov.A_A;
            double A_eff_sys = r_cov.delta;
            double A_eff_sys_err = r_cov.sigma_delta_corrected;  // <-- ora con -2*Cov
            double relative_sys = A_eff_sys / std::abs(A_true);
            double A_eff_sys_corr = r_cov.delta_corrected;       // <-- gia' quadCorrect-ato con sigma corretto
    
            A_eff_sys_4d[ix][iz] = A_eff_sys_corr;
            A_eff_sys_err_4d[ix][iz] = A_eff_sys_err;
            h_eff_sys->SetBinContent(ix + 1, iz + 1, A_eff_sys);
            h_eff_sys_rel->SetBinContent(ix + 1, iz + 1, relative_sys);

            //
            if(iz < 7) hist_eff_sys_4D[ix]->SetBinContent(iz+1, 1, A_eff_sys);
            else if (iz < 14) hist_eff_sys_4D[ix]->SetBinContent(iz-7+1, 2, A_eff_sys);
            else if (iz < 21) hist_eff_sys_4D[ix]->SetBinContent(iz-14+1, 3, A_eff_sys);
            else if (iz == 21){
                hist_eff_sys_4D[ix]->SetBinContent(1, 4, A_eff_sys);
                hist_eff_sys_4D[ix]->SetBinContent(2, 4, A_eff_sys);
                hist_eff_sys_4D[ix]->SetBinContent(3, 4, A_eff_sys);
            } else if (iz == 24) {
                hist_eff_sys_4D[ix]->SetBinContent(6, 4, A_eff_sys);
                hist_eff_sys_4D[ix]->SetBinContent(7, 4, A_eff_sys);
            }
            else hist_eff_sys_4D[ix]->SetBinContent(iz-22+4, 4, A_eff_sys);
        }

        TCanvas *c_bin_zPt = new TCanvas(Form("eff_sys_zPt_xQ2_%d", ix+1), "z vs P_{hT} bin", 800, 700);
        //c_bin_zPt->SetLogz();
        hist_eff_sys_4D[ix]->SetMinimum(0.0);
        hist_eff_sys_4D[ix]->SetMaximum(1.0);
        hist_eff_sys_4D[ix]->SetStats(0);
        hist_eff_sys_4D[ix]->Draw("colz");

        vector<TPolyLine*> gridLines_zp;
        vector<TText*> labels_zp;

        int bin_index_zp = 1;

        // Loop generale
        for (int iz = 0; iz < 4; ++iz) {
            for (size_t ip = 0; ip < binning_z_for_Pt[iz].size(); ++ip) {

                double Pt[5] = {bin_Pt[iz][0], bin_Pt[iz][1], bin_Pt[iz][1], bin_Pt[iz][0], bin_Pt[iz][0]};
                double z[5] = {binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][0]};

                // Rettangolo
                TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
                rect_zp->SetLineWidth(2);
                rect_zp->SetLineColor(kBlack);
                rect_zp->Draw("same");
                gridLines_zp.push_back(rect_zp);

                // Centro geometrico (utile se log-scale in Pt)
                double Pt_center = 0.5 * (bin_Pt[iz][0] + bin_Pt[iz][1]);
                double z_center = 0.5 * (binning_z_for_Pt[iz][ip][0] + binning_z_for_Pt[iz][ip][1]);

                int binx = hist_eff_sys_4D[ix]->GetXaxis()->FindBin(z_center);
                int biny = hist_eff_sys_4D[ix]->GetYaxis()->FindBin(Pt_center);

                double eff_val = hist_eff_sys_4D[ix]->GetBinContent(binx, biny);

                TText *label = new TText(z_center, Pt_center, Form("%.4f", eff_val));
                label->SetTextAlign(22);
                label->SetTextSize(0.03);
                //label->SetTextColor(kWhite);
                label->Draw("same");
                labels_zp.push_back(label);
            }
        }

        c_bin_zPt->Update();
        c_bin_zPt->Write();

    }


    // 2D systematic
    cout << " " << endl;
    cout << "================ EFFICIENCY SYSTEMATIC ====================" << endl;
    cout << " " << endl;
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        double N_pos_true_sum = 0;
        double N_neg_true_sum = 0;
        double N_pos_reco_sum = 0;
        double N_neg_reco_sum = 0;
        double N_pos_raw_sum = 0, N_neg_raw_sum = 0;   
        for (int iz = 0; iz < nbin_zPt; iz++) {
            N_pos_true_sum += vec_helicity_pos_preIDcut[ix][iz].size();
            N_neg_true_sum += vec_helicity_neg_preIDcut[ix][iz].size();
            N_pos_reco_sum += vec_helicity_pos_true[ix][iz].size();
            N_neg_reco_sum += vec_helicity_neg_true[ix][iz].size();
            N_pos_raw_sum  += vec_helicity_pos_raw[ix][iz].size();   
            N_neg_raw_sum  += vec_helicity_neg_raw[ix][iz].size();
        }
        if ((N_pos_true_sum + r_inv*N_neg_true_sum) == 0) continue;
        if ((N_pos_reco_sum + r_inv*N_neg_reco_sum) == 0) continue;

        /*
        double A_true = (N_pos_true_sum - r_inv*N_neg_true_sum) / (N_pos_true_sum + r_inv*N_neg_true_sum);
        double A_reco = (N_pos_reco_sum - r_inv*N_neg_reco_sum) / (N_pos_reco_sum + r_inv*N_neg_reco_sum);
        double sys = std::abs(A_true - A_reco);
        double sigma_A_true = sqrt( (1 - A_true*A_true) / (N_pos_true_sum + r_inv*r_inv*N_neg_true_sum) );
        double sigma_A_reco = sqrt( (1 - A_reco*A_reco) / (N_pos_reco_sum + r_inv*r_inv*N_neg_reco_sum) );
        double sys_err = sqrt(sigma_A_true*sigma_A_true + sigma_A_reco*sigma_A_reco);
        double sys_corr = quadCorrect(sys, sys_err);
        A_eff_sys_2d_xQ2[ix] = sys_corr;
        A_eff_sys_2d_xQ2_err[ix] = sys_err;
        */
       SysWithCovResult r_cov = ComputeSysWithCov(N_pos_true_sum, N_neg_true_sum,N_pos_reco_sum, N_neg_reco_sum,N_pos_raw_sum, N_neg_raw_sum,r_inv);
 
        A_eff_sys_2d_xQ2[ix] = r_cov.delta_corrected;
        A_eff_sys_2d_xQ2_err[ix] = r_cov.sigma_delta_corrected;

        cout << "Efficiency systematic in bin xQ2 " << ix + 1
         << ": raw_delta=" << r_cov.delta
         << "  sigma_naive=" << r_cov.sigma_delta_naive
         << "  sigma_corrected=" << r_cov.sigma_delta_corrected
         << "  delta_corrected=" << r_cov.delta_corrected << endl;

    }
    cout << "===========================================" << endl;
    for (int iz = 0; iz < nbin_zPt; iz++) {
        double N_pos_true_sum = 0;
        double N_neg_true_sum = 0;
        double N_pos_reco_sum = 0;
        double N_neg_reco_sum = 0;
        double N_pos_raw_sum = 0, N_neg_raw_sum = 0;
        for (int ix = 0; ix < nbin_xQ2; ix++) {
            N_pos_true_sum += vec_helicity_pos_preIDcut[ix][iz].size();
            N_neg_true_sum += vec_helicity_neg_preIDcut[ix][iz].size();
            N_pos_reco_sum += vec_helicity_pos_true[ix][iz].size();
            N_neg_reco_sum += vec_helicity_neg_true[ix][iz].size();
            N_pos_raw_sum  += vec_helicity_pos_raw[ix][iz].size();   
            N_neg_raw_sum  += vec_helicity_neg_raw[ix][iz].size();   
        }
        if ((N_pos_true_sum + r_inv*N_neg_true_sum) == 0) continue;
        if ((N_pos_reco_sum + r_inv*N_neg_reco_sum) == 0) continue;
        /*
        double A_true = (N_pos_true_sum - r_inv*N_neg_true_sum) / (N_pos_true_sum + r_inv*N_neg_true_sum);
        double A_reco = (N_pos_reco_sum - r_inv*N_neg_reco_sum) / (N_pos_reco_sum + r_inv*N_neg_reco_sum);
        double sys = std::abs(A_true - A_reco);
        double sigma_A_true = sqrt( (1 - A_true*A_true) / (N_pos_true_sum + r_inv*r_inv*N_neg_true_sum) );
        double sigma_A_reco = sqrt( (1 - A_reco*A_reco) / (N_pos_reco_sum + r_inv*r_inv*N_neg_reco_sum) );
        double sys_err = sqrt(sigma_A_true*sigma_A_true + sigma_A_reco*sigma_A_reco);
        double sys_corr = quadCorrect(sys, sys_err);
        A_eff_sys_2d_zPt[iz] = sys_corr;
        A_eff_sys_2d_zPt_err[iz] = sys_err;
        cout << "Efficiency systematic in bin zPt " << iz+1 << ": " << sys << endl;
        */

        SysWithCovResult r_cov = ComputeSysWithCov(N_pos_true_sum, N_neg_true_sum,N_pos_reco_sum, N_neg_reco_sum,N_pos_raw_sum, N_neg_raw_sum,r_inv);
 
        A_eff_sys_2d_zPt[iz] = r_cov.delta_corrected;
        A_eff_sys_2d_zPt_err[iz] = r_cov.sigma_delta_corrected;
    
        cout << "Efficiency systematic in bin zPt " << iz + 1
            << ": raw_delta=" << r_cov.delta
            << "  sigma_naive=" << r_cov.sigma_delta_naive
            << "  sigma_corrected=" << r_cov.sigma_delta_corrected
            << "  delta_corrected=" << r_cov.delta_corrected << endl;
    }
    cout << "===========================================" << endl;






    // =========================================== BIN MIGRATION SYSTEMATIC ===========================================
    // 
    vector<TH2D*> hist_binMig_4D(nbin_xQ2);
    dir_bin_mig->cd();
    vector<vector<double>> A_binMig_sys_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> A_binMig_sys_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<double> A_bm_sys_2d_xQ2(nbin_xQ2);
    vector<double> A_bm_sys_2d_zPt(nbin_zPt);
    vector<double> A_bm_sys_2d_xQ2_err(nbin_xQ2);
    vector<double> A_bm_sys_2d_zPt_err(nbin_zPt);
    TH2D* h_sys = new TH2D("h_sys", "Systematic; bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 1, nbin_xQ2+1, nbin_zPt, 1, nbin_zPt+1);
    TH2D* h_sys_rel = new TH2D("h_sys_rel", "RelativeSystematic; bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 0, nbin_xQ2, nbin_zPt, 0, nbin_zPt);
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        hist_binMig_4D[ix] = new TH2D(Form("hist_binMig_sys_xQ2_%d", ix+1), Form("Bin migration systematic | bin:(%d, x_{B}-Q^{2}) | #sigma_{sys}^{K+} = |A_{true}^{K+} - A_{reco}^{K+}|; z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
        for (int iz = 0; iz < nbin_zPt; iz++) {
            double N_pos_true = vec_helicity_pos_true[ix][iz].size();
            double N_neg_true = vec_helicity_neg_true[ix][iz].size();
            double N_pos_reco = vec_helicity_pos_reco[ix][iz].size();
            double N_neg_reco = vec_helicity_neg_reco[ix][iz].size();

            if(N_pos_true*N_neg_true*N_pos_reco*N_neg_reco == 0){ 
                cout << "Bin migration systematic: Bin (" << ix+1 << ", " << iz+1 << "): N_pos_true=" << N_pos_true << ", N_neg_true=" << N_neg_true 
                 << ", N_pos_reco=" << N_pos_reco << ", N_neg_reco=" << N_neg_reco << endl;
            }
            double r = 0.727688;
            double r_inv = 1.0/r;
            if ((N_pos_true + r_inv*N_neg_true) == 0 || (N_pos_reco + r_inv*N_neg_reco) == 0) continue;
            double A_true = (N_pos_true - r_inv*N_neg_true) / (N_pos_true + r_inv*N_neg_true);
            double A_reco = (N_pos_reco - r_inv*N_neg_reco) / (N_pos_reco + r_inv*N_neg_reco);
            double A_binMig_sys = std::abs(A_true - A_reco);
            double sigma_A_true = sqrt( (1 - A_true*A_true) / (N_pos_true + r_inv*r_inv*N_neg_true) );
            double sigma_A_reco = sqrt( (1 - A_reco*A_reco) / (N_pos_reco + r_inv*r_inv*N_neg_reco) );
            double A_binMig_sys_err = sqrt(sigma_A_true*sigma_A_true + sigma_A_reco*sigma_A_reco);
            double relative_sys = A_binMig_sys / std::abs(A_true);
            A_binMig_sys_4d[ix][iz] = A_binMig_sys;
            A_binMig_sys_err_4d[ix][iz] = A_binMig_sys_err;
            h_sys->SetBinContent(ix+1, iz+1, A_binMig_sys);
            h_sys_rel->SetBinContent(ix+1, iz+1, relative_sys);

            //
            if(iz < 7) hist_binMig_4D[ix]->SetBinContent(iz+1, 1, A_binMig_sys);
            else if (iz < 14) hist_binMig_4D[ix]->SetBinContent(iz-7+1, 2, A_binMig_sys);
            else if (iz < 21) hist_binMig_4D[ix]->SetBinContent(iz-14+1, 3, A_binMig_sys);
            else if (iz == 21){
                hist_binMig_4D[ix]->SetBinContent(1, 4, A_binMig_sys);
                hist_binMig_4D[ix]->SetBinContent(2, 4, A_binMig_sys);
                hist_binMig_4D[ix]->SetBinContent(3, 4, A_binMig_sys);
            } else if (iz == 24) {
                hist_binMig_4D[ix]->SetBinContent(6, 4, A_binMig_sys);
                hist_binMig_4D[ix]->SetBinContent(7, 4, A_binMig_sys);
            }
            else hist_binMig_4D[ix]->SetBinContent(iz-22+4, 4, A_binMig_sys);
        }


        TCanvas *c_bin_zPt = new TCanvas(Form("binMig_sys_zPt_xQ2_%d", ix+1), "z vs P_{hT} bin", 800, 700);
        //c_bin_zPt->SetLogz();
        hist_binMig_4D[ix]->SetMinimum(0.0);
        hist_binMig_4D[ix]->SetMaximum(1.0);
        hist_binMig_4D[ix]->SetStats(0);
        hist_binMig_4D[ix]->Draw("colz");

        vector<TPolyLine*> gridLines_zp;
        vector<TText*> labels_zp;

        int bin_index_zp = 1;

        // Loop generale
        for (int iz = 0; iz < 4; ++iz) {
            for (size_t ip = 0; ip < binning_z_for_Pt[iz].size(); ++ip) {

                double Pt[5] = {bin_Pt[iz][0], bin_Pt[iz][1], bin_Pt[iz][1], bin_Pt[iz][0], bin_Pt[iz][0]};
                double z[5] = {binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][0], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][1], binning_z_for_Pt[iz][ip][0]};

                // Rettangolo
                TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
                rect_zp->SetLineWidth(2);
                rect_zp->SetLineColor(kBlack);
                rect_zp->Draw("same");
                gridLines_zp.push_back(rect_zp);

                // Centro geometrico (utile se log-scale in Pt)
                double Pt_center = 0.5 * (bin_Pt[iz][0] + bin_Pt[iz][1]);
                double z_center = 0.5 * (binning_z_for_Pt[iz][ip][0] + binning_z_for_Pt[iz][ip][1]);

                int binx = hist_binMig_4D[ix]->GetXaxis()->FindBin(z_center);
                int biny = hist_binMig_4D[ix]->GetYaxis()->FindBin(Pt_center);

                double eff_val = hist_binMig_4D[ix]->GetBinContent(binx, biny);

                TText *label = new TText(z_center, Pt_center, Form("%.4f", eff_val));
                label->SetTextAlign(22);
                label->SetTextSize(0.03);
                //label->SetTextColor(kWhite);
                label->Draw("same");
                labels_zp.push_back(label);
            }
        }

        c_bin_zPt->Update();
        c_bin_zPt->Write();

    }

    TCanvas* c_binMig_sys = new TCanvas("c_bin_mig_sys", "Systematic uncertainty from bin migration", 800, 600);
    h_sys->SetStats(0);
    h_sys->Draw("colz");
    c_binMig_sys->Write();
    TCanvas* c_binMig_sys_rel = new TCanvas("c_bin_mig_relsys", "Systematic uncertainty from bin migration", 800, 600);
    h_sys_rel->SetStats(0);
    h_sys_rel->SetMaximum(1.0);
    h_sys_rel->Draw("colz");
    c_binMig_sys_rel->Write();

    // 2D systematic
    cout << " " << endl;
    cout << "================ BIN MIGRATION SYSTEMATIC ====================" << endl;
    cout << " " << endl;
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        double N_pos_true_sum = 0;
        double N_neg_true_sum = 0;
        double N_pos_reco_sum = 0;
        double N_neg_reco_sum = 0;
        for (int iz = 0; iz < nbin_zPt; iz++) {
            N_pos_true_sum += vec_helicity_pos_true[ix][iz].size();
            N_neg_true_sum += vec_helicity_neg_true[ix][iz].size();
            N_pos_reco_sum += vec_helicity_pos_reco[ix][iz].size();
            N_neg_reco_sum += vec_helicity_neg_reco[ix][iz].size();
        }
        if ((N_pos_true_sum + r_inv*N_neg_true_sum) == 0) continue;
        if ((N_pos_reco_sum + r_inv*N_neg_reco_sum) == 0) continue;
        double A_true = (N_pos_true_sum - r_inv*N_neg_true_sum) / (N_pos_true_sum + r_inv*N_neg_true_sum);
        double A_reco = (N_pos_reco_sum - r_inv*N_neg_reco_sum) / (N_pos_reco_sum + r_inv*N_neg_reco_sum);
        double sys = std::abs(A_true - A_reco);
        double sigma_A_true = sqrt( (1 - A_true*A_true) / (N_pos_true_sum + r_inv*r_inv*N_neg_true_sum) );
        double sigma_A_reco = sqrt( (1 - A_reco*A_reco) / (N_pos_reco_sum + r_inv*r_inv*N_neg_reco_sum) );
        double sys_err = sqrt(sigma_A_true*sigma_A_true + sigma_A_reco*sigma_A_reco);
        A_bm_sys_2d_xQ2[ix] = sys;
        A_bm_sys_2d_xQ2_err[ix] = sys_err;
        cout << "Bin migration systematic in bin xQ2 " << ix+1 << ": " << sys << endl;
    }
    cout << "===========================================" << endl;
    for (int iz = 0; iz < nbin_zPt; iz++) {
        double N_pos_true_sum = 0;
        double N_neg_true_sum = 0;
        double N_pos_reco_sum = 0;
        double N_neg_reco_sum = 0;
        for (int ix = 0; ix < nbin_xQ2; ix++) {
            N_pos_true_sum += vec_helicity_pos_true[ix][iz].size();
            N_neg_true_sum += vec_helicity_neg_true[ix][iz].size();
            N_pos_reco_sum += vec_helicity_pos_reco[ix][iz].size();
            N_neg_reco_sum += vec_helicity_neg_reco[ix][iz].size();
        }
        if ((N_pos_true_sum + r_inv*N_neg_true_sum) == 0) continue;
        if ((N_pos_reco_sum + r_inv*N_neg_reco_sum) == 0) continue;
        double A_true = (N_pos_true_sum - r_inv*N_neg_true_sum) / (N_pos_true_sum + r_inv*N_neg_true_sum);
        double A_reco = (N_pos_reco_sum - r_inv*N_neg_reco_sum) / (N_pos_reco_sum + r_inv*N_neg_reco_sum);
        double sys = std::abs(A_true - A_reco);
        double sigma_A_true = sqrt( (1 - A_true*A_true) / (N_pos_true_sum + r_inv*r_inv*N_neg_true_sum) );
        double sigma_A_reco = sqrt( (1 - A_reco*A_reco) / (N_pos_reco_sum + r_inv*r_inv*N_neg_reco_sum) );
        double sys_err = sqrt(sigma_A_true*sigma_A_true + sigma_A_reco*sigma_A_reco);
        A_bm_sys_2d_zPt[iz] = sys;
        A_bm_sys_2d_zPt_err[iz] = sys_err;
        cout << "Bin migration systematic in bin zPt " << iz+1 << ": " << sys << endl;
    }
    cout << "===========================================" << endl;





    // =========================================== PHI_H SYSTEMATIC ===========================================
    //
    vector<vector<TH1D*>> hist_A_mc(nbin_xQ2, vector<TH1D*>(nbin_zPt));
    //vector<vector<TH1D*>> hist_A_mc_small(nbin_xQ2, vector<TH1D*>(nbin_zPt));
    vector<vector<TH1D*>> hist_A_reco(nbin_xQ2, vector<TH1D*>(nbin_zPt));
    vector<vector<TH1D*>> hist_A_reco_small(nbin_xQ2, vector<TH1D*>(nbin_zPt));
    vector<vector<TH1D*>> hist_A_reco_Auu(nbin_xQ2, vector<TH1D*>(nbin_zPt));
    vector<vector<double>> Sys_Phi_sin(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_Phi_sin2(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_Phi_cos(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_Phi_cos2(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_Phi_sin_err(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_Phi_sin2_err(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_Phi_cos_err(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_Phi_cos2_err(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_uu_cont_sin(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Sys_uu_cont_sin2(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Auu_cos(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Auu_cos_err(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Auu_cos2(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Auu_cos2_err(nbin_xQ2, vector<double> (nbin_zPt));
    TGraphErrors* g_A_phi_h_sys_mc = new TGraphErrors();
    TH2D* h_sys_phi_A0 = new TH2D("h_sys_phi_sinx", "Phi_{h} systematic on sin(#Phi_{h}); bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 1, nbin_xQ2+1, nbin_zPt, 1, nbin_zPt+1);
    TH2D* h_sys_phi_A1 = new TH2D("h_sys_phi_sin2x", "Phi_{h} systematic on sin(2#Phi_{h}); bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 1, nbin_xQ2+1, nbin_zPt, 1, nbin_zPt+1);
    TH2D* h_sys_phi_A2 = new TH2D("h_sys_phi_cosx", "Phi_{h} systematic on cos(#Phi_{h}); bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 1, nbin_xQ2+1, nbin_zPt, 1, nbin_zPt+1);
    TH2D* h_sys_phi_A3 = new TH2D("h_sys_phi_cos2x", "Phi_{h} systematic on cos(2#Phi_{h}); bin: x_{B}-Q^{2}; bin: z-P_{hT}", nbin_xQ2, 1, nbin_xQ2+1, nbin_zPt, 1, nbin_zPt+1);
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {
            TH1D* h_plus_mc  = hist_Phi_h_plus_true[ix][iz]; // non più mc ma true
            TH1D* h_minus_mc = hist_Phi_h_minus_true[ix][iz];
            TH1D* h_plus_reco  = hist_Phi_h_plus_reco[ix][iz];
            TH1D* h_minus_reco = hist_Phi_h_minus_reco[ix][iz];
            int nBins = h_plus_mc->GetNbinsX();
            // crea istogramma A(φ)
            TH1D* h_A_mc = (TH1D*)h_plus_mc->Clone(Form("h_A_mc_%d_%d", ix, iz));
            //TH1D* h_A_mc_small = (TH1D*)h_plus_mc->Clone(Form("h_A_mc_small_%d_%d", ix, iz));
            TH1D* h_A_reco = (TH1D*)h_plus_reco->Clone(Form("h_A_reco_%d_%d", ix, iz));
            TH1D* h_A_reco_small = (TH1D*)h_plus_reco->Clone(Form("h_A_reco_small_%d_%d", ix, iz));
            TH1D* h_A_reco_Auu = (TH1D*)h_plus_reco->Clone(Form("h_A_reco_Auu_%d_%d", ix, iz));
            h_A_mc->Reset();
            //h_A_mc_small->Reset();
            h_A_reco->Reset();
            h_A_reco_small->Reset();
            h_A_reco_Auu->Reset();
            int validBins_mc = 0;
            int validBins_reco = 0;
            for (int ib = 1; ib <= nBins; ib++) {
                // MC
                double Npos_mc = h_plus_mc->GetBinContent(ib);
                double Nneg_mc = h_minus_mc->GetBinContent(ib);
                double denom_mc = (Npos_mc + r_inv * Nneg_mc);
                if (denom_mc == 0) continue;
                validBins_mc++;
                double A_mc = (Npos_mc - r_inv * Nneg_mc) / denom_mc;
                // errore (opzionale ma consigliato)
                double err_mc = sqrt((1 - A_mc*A_mc) / denom_mc);
                h_A_mc->SetBinContent(ib, A_mc);
                h_A_mc->SetBinError(ib, err_mc);
                //h_A_mc_small->SetBinContent(ib, A_mc);
                //h_A_mc_small->SetBinError(ib, err_mc);
                // RECO
                double Npos_reco = h_plus_reco->GetBinContent(ib);
                double Nneg_reco = h_minus_reco->GetBinContent(ib);
                double denom_reco = (Npos_reco + r_inv * Nneg_reco);
                if (denom_reco == 0) continue;
                validBins_reco++;
                double A_reco = (Npos_reco - r_inv * Nneg_reco) / denom_reco;
                double err_reco = sqrt((1 - A_reco*A_reco) / denom_reco);
                h_A_reco->SetBinContent(ib, A_reco);
                h_A_reco->SetBinError(ib, err_reco);
                h_A_reco_small->SetBinContent(ib, A_reco);
                h_A_reco_small->SetBinError(ib, err_reco);
                h_A_reco_Auu->SetBinContent(ib, A_reco);
                h_A_reco_Auu->SetBinError(ib, err_reco);
            }
            hist_A_mc[ix][iz] = h_A_mc;
            //hist_A_mc_small[ix][iz] = h_A_mc_small;
            hist_A_reco[ix][iz] = h_A_reco;
            hist_A_reco_small[ix][iz] = h_A_reco_small;
            hist_A_reco_Auu[ix][iz] = h_A_reco_Auu;
            //fit
            if (validBins_mc > 6 && validBins_reco > 6) {
                TF1* f_mc = new TF1("f_mc", "[4] + [0]*sin(x) + [1]*sin(2*x) + [2]*cos(x) + [3]*cos(2*x)", -TMath::Pi(), TMath::Pi());
                //TF1* f_mc_small = new TF1("f_mc_small", "[2] + [0]*sin(x) + [1]*sin(2*x)", -TMath::Pi(), TMath::Pi());
                f_mc->SetParLimits(0, -1.0, 1.0);
                f_mc->SetParLimits(1, -1.0, 1.0);
                f_mc->SetParLimits(2, -1.0, 1.0);
                f_mc->SetParLimits(3, -1.0, 1.0);
                hist_A_mc[ix][iz]->Fit(f_mc, "Q");
                double A0_mc = f_mc->GetParameter(0);
                double A1_mc = f_mc->GetParameter(1);
                double A2_mc = f_mc->GetParameter(2);
                double A3_mc = f_mc->GetParameter(3);
                //hist_A_mc_small[ix][iz]->Fit(f_mc_small, "Q");
            
                TF1* f_reco = new TF1("f_reco", "[4] + [0]*sin(x) + [1]*sin(2*x) + [2]*cos(x) + [3]*cos(2*x)", -TMath::Pi(), TMath::Pi());
                TF1* f_reco_small = new TF1("f_reco_small", "[2] + [0]*sin(x) + [1]*sin(2*x)", -TMath::Pi(), TMath::Pi());
                TF1* f_reco_Auu = new TF1("f_reco_Auu", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", -TMath::Pi(), TMath::Pi());
                if (hist_Phi_h_reco[ix][iz]->Integral() > 0) hist_Phi_h_reco[ix][iz]->Scale(1.0 / hist_Phi_h_reco[ix][iz]->Integral());
                f_reco->SetParLimits(0, -1.0, 1.0);
                f_reco->SetParLimits(1, -1.0, 1.0);
                f_reco->SetParLimits(2, -1.0, 1.0);
                f_reco->SetParLimits(3, -1.0, 1.0);
                f_reco_small->SetParLimits(0, -1.0, 1.0);
                f_reco_small->SetParLimits(1, -1.0, 1.0);
                f_reco_Auu->SetParLimits(0, 0.0, 1000.0);
                f_reco_Auu->SetParLimits(1, -1.0, 1.0);
                f_reco_Auu->SetParLimits(2, -1.0, 1.0);
                f_reco_Auu->SetParameters(1.0, 0.0, 0.0);
                hist_A_reco[ix][iz]->Fit(f_reco, "Q");
                hist_A_reco_small[ix][iz]->Fit(f_reco_small, "Q");
                hist_Phi_h_reco[ix][iz]->Fit(f_reco_Auu, "Q");
                double A0_reco = f_reco->GetParameter(0);
                double A1_reco = f_reco->GetParameter(1);
                double A2_reco = f_reco->GetParameter(2);
                double A3_reco = f_reco->GetParameter(3);
            
                double A0_sys = std::abs(A0_mc - A0_reco);
                double A1_sys = std::abs(A1_mc - A1_reco);
                double A2_sys = std::abs(A2_mc - A2_reco);
                double A3_sys = std::abs(A3_mc - A3_reco);
                Sys_Phi_sin[ix][iz] = A0_sys;
                Sys_Phi_sin2[ix][iz] = A1_sys;
                Sys_Phi_cos[ix][iz] = A2_sys;
                Sys_Phi_cos2[ix][iz] = A3_sys;
                // cross-talk
                Sys_uu_cont_sin[ix][iz] = std::abs(f_reco->GetParameter(0) - f_reco_small->GetParameter(0));
                Sys_uu_cont_sin2[ix][iz] = std::abs(f_reco->GetParameter(1) - f_reco_small->GetParameter(1));
                //
                Sys_Phi_sin_err[ix][iz] = sqrt( pow(f_mc->GetParError(0), 2) + pow(f_reco->GetParError(0), 2) );
                Sys_Phi_sin2_err[ix][iz] = sqrt( pow(f_mc->GetParError(1), 2) + pow(f_reco->GetParError(1), 2) );
                Sys_Phi_cos_err[ix][iz] = sqrt( pow(f_mc->GetParError(2), 2) + pow(f_reco->GetParError(2), 2) );
                Sys_Phi_cos2_err[ix][iz] = sqrt( pow(f_mc->GetParError(3), 2) + pow(f_reco->GetParError(3), 2) );
                h_sys_phi_A0->SetBinContent(ix+1, iz+1, A0_sys);
                h_sys_phi_A1->SetBinContent(ix+1, iz+1, A1_sys);
                h_sys_phi_A2->SetBinContent(ix+1, iz+1, A2_sys);
                h_sys_phi_A3->SetBinContent(ix+1, iz+1, A3_sys);


                // for Auu
                double p0   = f_reco_Auu->GetParameter(0);
                double p1   = f_reco_Auu->GetParameter(1);
                double p2   = f_reco_Auu->GetParameter(2);
                double e0   = f_reco_Auu->GetParError(0);
                double e1   = f_reco_Auu->GetParError(1);
                double e2   = f_reco_Auu->GetParError(2);
                /*
                double cov02 = f_reco_Auu->GetCovarianceMatrixElement(0,2);
                double cov12 = f_reco_Auu->GetCovarianceMatrixElement(1,2);
                double AUU_cos = p0 / p2;
                double errAUU_cos    = sqrt((e0*e0)/(p2*p2) + (p0*p0*e2*e2)/(p2*p2*p2*p2) - 2.0*p0*cov02/(p2*p2*p2));
                double AUU_cos2 = p1 / p2;
                double errAUU_cos2   = sqrt((e1*e1)/(p2*p2) + (p1*p1*e2*e2)/(p2*p2*p2*p2) - 2.0*p1*cov12/(p2*p2*p2));
                */
                Auu_cos[ix][iz] = p1;
                Auu_cos_err[ix][iz] = e1;
                Auu_cos2[ix][iz] = p2;
                Auu_cos2_err[ix][iz] = e2;
            }
        }
    }
    dir_Phih_sin->cd();
    h_sys_phi_A0->SetStats(0);
    TCanvas* c_sys_phi_A0 = new TCanvas("c_sys_phi_A0", "Systematic uncertainty on A_{0} from phi_h fit", 800, 600);
    h_sys_phi_A0->SetMaximum(1.0);
    h_sys_phi_A0->Draw("colz");
    c_sys_phi_A0->Write();

    hist_A_reco[0][6]->Write();
    hist_A_reco[0][5]->Write();
    hist_Phi_h_reco[0][2]->Write();
    hist_Phi_h_reco[0][5]->Write();
    hist_Phi_h_reco[2][2]->Write();
    hist_Phi_h_reco[2][5]->Write();

    TH1D hist_sys_sin ("hist_sys_sin", "Systematic uncertainty in sin(x); diff; count", 100, 0.0, 0.5);
    TH1D hist_diff_asym_sin ("hist_diff_asym_sin", "Difference in asymmetry |f(sin+cos) - f(sin)| in sin(x); diff; count", 100, 0.0, 0.1);
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {
            hist_sys_sin.Fill(Sys_Phi_sin[ix][iz]);
            hist_diff_asym_sin.Fill(Sys_uu_cont_sin[ix][iz]);
        }
    }

    hist_sys_sin.Write();
    hist_diff_asym_sin.Write();

    dir_Phih_sin2->cd();
    h_sys_phi_A1->SetStats(0);
    TCanvas* c_sys_phi_A1 = new TCanvas("c_sys_phi_A1", "Systematic uncertainty on A_{1} from phi_h fit", 800, 600);
    h_sys_phi_A1->SetMaximum(1.0);
    h_sys_phi_A1->Draw("colz");
    c_sys_phi_A1->Write();

    TH1D hist_sys_sin2 ("hist_sys_sin2", "Systematic uncertainty in sin(2x); diff; count", 100, 0.0, 0.5);
    TH1D hist_diff_asym_sin2 ("hist_diff_asym_sin2", "Difference in asymmetry |f(sin+cos) - f(sin)| in sin(2x); diff; count", 100, 0.0, 0.1);
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {
            hist_sys_sin2.Fill(Sys_Phi_sin2[ix][iz]);
            hist_diff_asym_sin2.Fill(Sys_uu_cont_sin2[ix][iz]);
        }
    }
    hist_sys_sin2.Write();
    hist_diff_asym_sin2.Write();


    dir_Phih_cos->cd();
    h_sys_phi_A2->SetStats(0);
    TCanvas* c_sys_phi_A2 = new TCanvas("c_sys_phi_A2", "Systematic uncertainty on A_{2} from phi_h fit", 800, 600);
    h_sys_phi_A2->SetMaximum(1.0);
    h_sys_phi_A2->Draw("colz");
    c_sys_phi_A2->Write();

    TH1D hist_sys_cos ("hist_sys_cos", "Systematic uncertainty in cos(x); diff; count", 100, 0.0, 0.5);
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {
            hist_sys_cos.Fill(Sys_Phi_cos[ix][iz]);
        }
    }
    hist_sys_cos.Write();

    dir_Phih_cos2->cd();
    h_sys_phi_A3->SetStats(0);
    TCanvas* c_sys_phi_A3 = new TCanvas("c_sys_phi_A3", "Systematic uncertainty on A_{3} from phi_h fit", 800, 600);
    h_sys_phi_A3->SetMaximum(1.0);
    h_sys_phi_A3->Draw("colz");
    c_sys_phi_A3->Write();

    TH1D hist_sys_cos2 ("hist_sys_cos2", "Systematic uncertainty in cos(2x); diff; count", 100, 0.0, 0.5);
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {
            hist_sys_cos2.Fill(Sys_Phi_cos2[ix][iz]);
        }
    }
    hist_sys_cos2.Write();


    // ===========================================================================
    // CASO 2D integra su tutti gli ix per ogni iz
    //
    // ===========================================================================
    vector<double> Sys_Phi_sin_zPt(nbin_zPt);
    vector<double> Sys_Phi_sin2_zPt(nbin_zPt);
    vector<double> Sys_Phi_cos_zPt(nbin_zPt);
    vector<double> Sys_Phi_cos2_zPt(nbin_zPt);
    vector<double> Sys_Phi_sin_err_zPt(nbin_zPt);
    vector<double> Sys_Phi_sin2_err_zPt(nbin_zPt);
    vector<double> Sys_Phi_cos_err_zPt(nbin_zPt);
    vector<double> Sys_Phi_cos2_err_zPt(nbin_zPt);
    vector<double> Auu_cos_zPt_mc(nbin_zPt);
    vector<double> Auu_cos_zPt_err_mc(nbin_zPt);
    vector<double> Auu_cos2_zPt_mc(nbin_zPt);
    vector<double> Auu_cos2_zPt_err_mc(nbin_zPt);
    vector<double> Auu_cos_zPt_reco(nbin_zPt);
    vector<double> Auu_cos_zPt_err_reco(nbin_zPt);
    vector<double> Auu_cos2_zPt_reco(nbin_zPt);
    vector<double> Auu_cos2_zPt_err_reco(nbin_zPt);
    vector<double> Auu_cos_zPt_true(nbin_zPt);
    vector<double> Auu_cos_zPt_err_true(nbin_zPt);
    vector<double> Auu_cos2_zPt_true(nbin_zPt);
    vector<double> Auu_cos2_zPt_err_true(nbin_zPt);
    vector<double> Auu_cos_zPt_preIDcut(nbin_zPt);
    vector<double> Auu_cos_zPt_err_preIDcut(nbin_zPt);
    vector<double> Auu_cos2_zPt_preIDcut(nbin_zPt);
    vector<double> Auu_cos2_zPt_err_preIDcut(nbin_zPt);
    
    for (int iz = 0; iz < nbin_zPt; iz++) {
    
        // --- somma (integra) gli istogrammi su tutti i bin xB-Q2 ---
        TH1D* h_plus_true_sum   = (TH1D*)hist_Phi_h_plus_true[0][iz]->Clone(Form("h_plus_true_sum_%d", iz));
        TH1D* h_minus_true_sum  = (TH1D*)hist_Phi_h_minus_true[0][iz]->Clone(Form("h_minus_true_sum_%d", iz));
        TH1D* h_plus_reco_sum = (TH1D*)hist_Phi_h_plus_reco[0][iz]->Clone(Form("h_plus_reco_sum_%d", iz));
        TH1D* h_minus_reco_sum= (TH1D*)hist_Phi_h_minus_reco[0][iz]->Clone(Form("h_minus_reco_sum_%d", iz));
        // per preIDcut e MC
        TH1D* h_plus_mc_sum   = (TH1D*)hist_Phi_h_plus_mc[0][iz]->Clone(Form("h_plus_mc_sum_%d", iz));
        TH1D* h_minus_mc_sum  = (TH1D*)hist_Phi_h_minus_mc[0][iz]->Clone(Form("h_minus_mc_sum_%d", iz));
        TH1D* h_plus_preIDcut_sum = (TH1D*)hist_Phi_h_plus_preIDcut[0][iz]->Clone(Form("h_plus_preIDcut_sum_%d", iz));
        TH1D* h_minus_preIDcut_sum= (TH1D*)hist_Phi_h_minus_preIDcut[0][iz]->Clone(Form("h_minus_preIDcut_sum_%d", iz));
    
        for (int ix = 1; ix < nbin_xQ2; ix++) {
            h_plus_true_sum->Add(hist_Phi_h_plus_true[ix][iz]);
            h_minus_true_sum->Add(hist_Phi_h_minus_true[ix][iz]);
            h_plus_reco_sum->Add(hist_Phi_h_plus_reco[ix][iz]);
            h_minus_reco_sum->Add(hist_Phi_h_minus_reco[ix][iz]);
        }
    
        // --- stessa identica funzione di prima, ma sugli istogrammi sommati ---
        PhiHSysResult r = ComputePhiHSystematic(h_plus_true_sum, h_minus_true_sum,h_plus_reco_sum, h_minus_reco_sum,r_inv, Form("zPt_%d", iz));
        Auu_cos_zPt_true[iz] = r.auu_cos_mc;
        Auu_cos_zPt_err_true[iz] = r.auu_cos_mc_err;
        Auu_cos2_zPt_true[iz] = r.auu_cos2_mc;
        Auu_cos2_zPt_err_true[iz] = r.auu_cos2_mc_err;
        Auu_cos_zPt_reco[iz] = r.auu_cos_reco;
        Auu_cos_zPt_err_reco[iz] = r.auu_cos_reco_err;
        Auu_cos2_zPt_reco[iz] = r.auu_cos2_reco;
        Auu_cos2_zPt_err_reco[iz] = r.auu_cos2_reco_err;

        if (r.ok) {
            Sys_Phi_sin_zPt[iz]  = r.sys_sin;
            Sys_Phi_sin2_zPt[iz] = r.sys_sin2;
            Sys_Phi_cos_zPt[iz]  = r.sys_cos;
            Sys_Phi_cos2_zPt[iz] = r.sys_cos2;
            Sys_Phi_sin_err_zPt[iz]  = r.sys_sin_err;
            Sys_Phi_sin2_err_zPt[iz] = r.sys_sin2_err;
            Sys_Phi_cos_err_zPt[iz]  = r.sys_cos_err;
            Sys_Phi_cos2_err_zPt[iz] = r.sys_cos2_err;
        }

        // per mc e preIDcut
        PhiHSysResult r2 = ComputePhiHSystematic(h_plus_mc_sum, h_minus_mc_sum,h_plus_preIDcut_sum, h_minus_preIDcut_sum,r_inv, Form("zPt2_%d", iz));
        Auu_cos_zPt_mc[iz] = r2.auu_cos_mc;
        Auu_cos_zPt_err_mc[iz] = r2.auu_cos_mc_err;
        Auu_cos2_zPt_mc[iz] = r2.auu_cos2_mc;
        Auu_cos2_zPt_err_mc[iz] = r2.auu_cos2_mc_err;
        Auu_cos_zPt_preIDcut[iz] = r2.auu_cos_reco;
        Auu_cos_zPt_err_preIDcut[iz] = r2.auu_cos_reco_err;
        Auu_cos2_zPt_preIDcut[iz] = r2.auu_cos2_reco;
        Auu_cos2_zPt_err_preIDcut[iz] = r2.auu_cos2_reco_err;
    }


    /*
    // plot in bin of z-Pt to highlight the effect of phi_h systematic
    for (int iz = 0; iz < nbin_zPt; iz++){
        // sin
        TGraphErrors* Asys_phih_sin = new TGraphErrors();
        TGraphErrors* Asys_phih_sin_xB_1 = new TGraphErrors();
        TGraphErrors* Asys_phih_sin_xB_2 = new TGraphErrors();
        TGraphErrors* Asys_phih_sin_xB_3 = new TGraphErrors();
        TGraphErrors* Asys_phih_sin_xB_4 = new TGraphErrors();
        //
        TGraphErrors* Asys_phih_sin2 = new TGraphErrors();
        TGraphErrors* Asys_phih_cos = new TGraphErrors();
        TGraphErrors* Asys_phih_cos2 = new TGraphErrors();
        int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0;
        for (int ix = 0; ix < nbin_xQ2; ix++){
            if (vec_kaonp_xB_4d[ix][iz].empty()) continue;
            double sum_xB = 0.0;
            double min_x = 1;
            double max_x = 0;
            for (double val : vec_kaonp_xB_4d[ix][iz]) sum_xB += val;
            double mean_xB = sum_xB / vec_kaonp_xB_4d[ix][iz].size();
            double sys_sin = Sys_Phi_sin[ix][iz];
            double sys_sin_err = Sys_Phi_sin_err[ix][iz];
            Asys_phih_sin->SetPoint(p_idx, mean_xB, sys_sin);
            p_idx++;
            if(ix < 5){
                Asys_phih_sin_xB_1->SetPoint(p_idx_1, mean_xB, sys_sin);
                Asys_phih_sin_xB_1->SetPointError(p_idx_1, 0.0, sys_sin_err); // No x error
                p_idx_1++;
            } else if (ix < 10){
                Asys_phih_sin_xB_2->SetPoint(p_idx_2, mean_xB, sys_sin);
                Asys_phih_sin_xB_2->SetPointError(p_idx_2, 0.0, sys_sin_err);
                p_idx_2++;
            } else if (ix < 12){
                Asys_phih_sin_xB_3->SetPoint(p_idx_3, mean_xB, sys_sin);
                Asys_phih_sin_xB_3->SetPointError(p_idx_3, 0.0, sys_sin_err);
                p_idx_3++;
            } else if (ix < 18){
                Asys_phih_sin_xB_4->SetPoint(p_idx_4, mean_xB, sys_sin);
                Asys_phih_sin_xB_4->SetPointError(p_idx_4, 0.0, sys_sin_err);
                p_idx_4++;
            }
        }

        Asys_phih_sin_xB_1->SetMarkerStyle(20), Asys_phih_sin_xB_2->SetMarkerStyle(20), Asys_phih_sin_xB_3->SetMarkerStyle(20), Asys_phih_sin_xB_4->SetMarkerStyle(20);
        Asys_phih_sin_xB_1->SetLineColor(kAzure-5), Asys_phih_sin_xB_1->SetMarkerColor(kAzure-5);
        Asys_phih_sin_xB_2->SetLineColor(kViolet-5), Asys_phih_sin_xB_2->SetMarkerColor(kViolet-5);
        Asys_phih_sin_xB_3->SetLineColor(kPink-5), Asys_phih_sin_xB_3->SetMarkerColor(kPink-5);
        Asys_phih_sin_xB_4->SetLineColor(kOrange-5), Asys_phih_sin_xB_4->SetMarkerColor(kOrange-5);


        TCanvas* c_Aut_xB_2 = new TCanvas(Form("PhiSys_sin_vs_xB_zPt_bin%d", iz+1), "sin(#Phi_{h}) longitudinal asymmetry vs x_{B} for bin (n, z-P_{hT})", 800, 600);
        //Asys_phih_sin_xB->Draw("A");
        TMultiGraph *mg_Aut4D_xB_2 = new TMultiGraph();
        mg_Aut4D_xB_2->SetTitle(Form("Phi_{h} systematic for sin(x) vs x_{B} for bin (%d, z-P_{hT}) ; x_{B}; #sigma_{Phi_{h}}^{sinx}", iz+1)); // #sqrt{2#epsilon(1+#epsilon)}
        mg_Aut4D_xB_2->Add(Asys_phih_sin_xB_1, "P");
        mg_Aut4D_xB_2->Add(Asys_phih_sin_xB_2, "P");
        mg_Aut4D_xB_2->Add(Asys_phih_sin_xB_3, "P");
        mg_Aut4D_xB_2->Add(Asys_phih_sin_xB_4, "P");
        mg_Aut4D_xB_2->Draw("A");
        mg_Aut4D_xB_2->GetYaxis()->SetRangeUser(-1, 1);
        TLine* guideLine2 = new TLine(mg_Aut4D_xB_2->GetXaxis()->GetXmin(), 0, mg_Aut4D_xB_2->GetXaxis()->GetXmax(), 0);
        guideLine2->SetLineStyle(2);  
        guideLine2->SetLineColor(kGray+1);
        guideLine2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(Asys_phih_sin_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
        legend->AddEntry(Asys_phih_sin_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
        legend->AddEntry(Asys_phih_sin_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
        legend->AddEntry(Asys_phih_sin_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Aut_xB_2->Update();
        c_Aut_xB_2->Write();
    }
    */



    // ============================ ACCEPTANCE SYSTEMATIC ============================
    // SISTEMATICA ACCETTANZA + TRACKING: mc (A) vs preIDcut (B)
    const double N_Omega_pos_total = 57.88e6;  // eventi generati, elicita' +1
    const double N_Omega_neg_total = 42.12e6;  // eventi generati, elicita' -1

    // ---------------------------------------------------------------------
    // CASO 4D
    // ---------------------------------------------------------------------
    vector<vector<double>> A_acc_sys_4d(nbin, vector<double>(nbin));
    vector<vector<double>> A_acc_sys_err_4d(nbin, vector<double>(nbin));

    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {

            double N_pos_mc = vec_helicity_pos_mc[ix][iz].size();        // "A" = mc
            double N_neg_mc = vec_helicity_neg_mc[ix][iz].size();
            double N_pos_preIDcut = vec_helicity_pos_preIDcut[ix][iz].size();  // "B" = preIDcut
            double N_neg_preIDcut = vec_helicity_neg_preIDcut[ix][iz].size();

            if ((N_pos_mc + r_inv * N_neg_mc) == 0) continue;
            if ((N_pos_preIDcut + r_inv * N_neg_preIDcut) == 0) continue;

            // NOTA: qui Omega = TOTALE generato (fisso, stesso per tutti i bin),
            // NON i vettori "raw" per bin (che sarebbero concettualmente sbagliati
            // qui, essendo raw un SOTTOinsieme di mc, non un sovrainsieme).
            SysWithCovResult r_cov = ComputeSysWithCov(N_pos_mc, N_neg_mc,
                                                        N_pos_preIDcut, N_neg_preIDcut,
                                                        N_Omega_pos_total, N_Omega_neg_total,
                                                        r_inv);

            A_acc_sys_4d[ix][iz] = r_cov.delta_corrected;
            A_acc_sys_err_4d[ix][iz] = r_cov.sigma_delta_corrected;
        }
    }

    // ---------------------------------------------------------------------
    // CASO 2D -- integrato su iz per ogni ix
    // ---------------------------------------------------------------------
    vector<double> A_acc_sys_2d_xQ2(nbin_xQ2);
    vector<double> A_acc_sys_2d_xQ2_err(nbin_xQ2);

    cout << " " << endl;
    cout << "================ ACCEPTANCE+TRACKING SYSTEMATIC (xQ2) ====================" << endl;

    for (int ix = 0; ix < nbin_xQ2; ix++) {
        double N_pos_mc_sum = 0, N_neg_mc_sum = 0;
        double N_pos_preIDcut_sum = 0, N_neg_preIDcut_sum = 0;

        for (int iz = 0; iz < nbin_zPt; iz++) {
            N_pos_mc_sum += vec_helicity_pos_mc[ix][iz].size();
            N_neg_mc_sum += vec_helicity_neg_mc[ix][iz].size();
            N_pos_preIDcut_sum += vec_helicity_pos_preIDcut[ix][iz].size();
            N_neg_preIDcut_sum += vec_helicity_neg_preIDcut[ix][iz].size();
        }
        if ((N_pos_mc_sum + r_inv * N_neg_mc_sum) == 0) continue;
        if ((N_pos_preIDcut_sum + r_inv * N_neg_preIDcut_sum) == 0) continue;

        // NOTA: anche qui integrato, Omega resta il totale generato GLOBALE fisso
        // (non va sommato/moltiplicato per il numero di bin: e' un singolo numero,
        // lo stesso usato sopra).
        SysWithCovResult r_cov = ComputeSysWithCov(N_pos_mc_sum, N_neg_mc_sum,
                                                    N_pos_preIDcut_sum, N_neg_preIDcut_sum,
                                                    N_Omega_pos_total, N_Omega_neg_total,
                                                    r_inv);

        A_acc_sys_2d_xQ2[ix] = r_cov.delta_corrected;
        A_acc_sys_2d_xQ2_err[ix] = r_cov.sigma_delta_corrected;

        cout << "Acceptance+tracking systematic in bin xQ2 " << ix + 1
            << ": raw_delta=" << r_cov.delta
            << "  sigma_naive=" << r_cov.sigma_delta_naive
            << "  sigma_corrected=" << r_cov.sigma_delta_corrected
            << "  delta_corrected=" << r_cov.delta_corrected << endl;
    }
    cout << "===========================================" << endl;

    // ---------------------------------------------------------------------
    // CASO 2D -- integrato su ix per ogni iz
    // ---------------------------------------------------------------------
    vector<double> A_acc_sys_2d_zPt(nbin_zPt);
    vector<double> A_acc_sys_2d_zPt_err(nbin_zPt);

    for (int iz = 0; iz < nbin_zPt; iz++) {
        double N_pos_mc_sum = 0, N_neg_mc_sum = 0;
        double N_pos_preIDcut_sum = 0, N_neg_preIDcut_sum = 0;

        for (int ix = 0; ix < nbin_xQ2; ix++) {
            N_pos_mc_sum += vec_helicity_pos_mc[ix][iz].size();
            N_neg_mc_sum += vec_helicity_neg_mc[ix][iz].size();
            N_pos_preIDcut_sum += vec_helicity_pos_preIDcut[ix][iz].size();
            N_neg_preIDcut_sum += vec_helicity_neg_preIDcut[ix][iz].size();
        }
        if ((N_pos_mc_sum + r_inv * N_neg_mc_sum) == 0) continue;
        if ((N_pos_preIDcut_sum + r_inv * N_neg_preIDcut_sum) == 0) continue;

        SysWithCovResult r_cov = ComputeSysWithCov(N_pos_mc_sum, N_neg_mc_sum,
                                                    N_pos_preIDcut_sum, N_neg_preIDcut_sum,
                                                    N_Omega_pos_total, N_Omega_neg_total,
                                                    r_inv);

        A_acc_sys_2d_zPt[iz] = r_cov.delta_corrected;
        A_acc_sys_2d_zPt_err[iz] = r_cov.sigma_delta_corrected;

        cout << "Acceptance+tracking systematic in bin zPt " << iz + 1
            << ": raw_delta=" << r_cov.delta
            << "  sigma_naive=" << r_cov.sigma_delta_naive
            << "  sigma_corrected=" << r_cov.sigma_delta_corrected
            << "  delta_corrected=" << r_cov.delta_corrected << endl;
    }
    cout << "===========================================" << endl;





    // ===========================================================================
    // Calcolo sistematica di purezza/contaminazione:
    
    vector<vector<double>> A_purity_sys_4d(nbin, vector<double>(nbin));
    vector<vector<double>> A_purity_sys_err_4d(nbin, vector<double>(nbin));
    
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {
    
            double N_pos_allID = vec_helicity_pos_all_ID[ix][iz].size();      // "A"
            double N_neg_allID = vec_helicity_neg_all_ID[ix][iz].size();
            double N_pos_true  = vec_helicity_pos_true[ix][iz].size();        // "B" (sottoinsieme di A)
            double N_neg_true  = vec_helicity_neg_true[ix][iz].size();

            double N_pos_omega = vec_helicity_pos_all_preIDcut[ix][iz].size(); 
            double N_neg_omega = vec_helicity_neg_all_preIDcut[ix][iz].size();
    
            if ((N_pos_true + r_inv * N_neg_true) == 0) continue;
            if ((N_pos_allID + r_inv * N_neg_allID) == 0) continue;
    
            SysWithCovResult r_cov = ComputeSysWithCov(N_pos_allID, N_neg_allID,N_pos_true, N_neg_true,N_pos_omega, N_neg_omega,r_inv);
    
            A_purity_sys_4d[ix][iz] = r_cov.delta_corrected;
            A_purity_sys_err_4d[ix][iz] = r_cov.sigma_delta_corrected;
        }
    }
    
    // ===========================================================================
    // Caso 2D (integrato), stessa logica, esattamente come per gli altri sys
    // ===========================================================================
    vector<double> A_purity_sys_2d_zPt(nbin_zPt);
    vector<double> A_purity_sys_2d_zPt_err(nbin_zPt);
    
    cout << " " << endl;
    cout << "================ PID PURITY/CONTAMINATION SYSTEMATIC (zPt) ====================" << endl;
    
    for (int iz = 0; iz < nbin_zPt; iz++) {
        double N_pos_true_sum = 0, N_neg_true_sum = 0;
        double N_pos_allID_sum = 0, N_neg_allID_sum = 0;
        double N_pos_omega_sum = 0, N_neg_omega_sum = 0;
    
        for (int ix = 0; ix < nbin_xQ2; ix++) {
            N_pos_true_sum += vec_helicity_pos_true[ix][iz].size();
            N_neg_true_sum += vec_helicity_neg_true[ix][iz].size();
            N_pos_allID_sum += vec_helicity_pos_all_ID[ix][iz].size();
            N_neg_allID_sum += vec_helicity_neg_all_ID[ix][iz].size();
            N_pos_omega_sum += vec_helicity_pos_all_preIDcut[ix][iz].size();
            N_neg_omega_sum += vec_helicity_neg_all_preIDcut[ix][iz].size();
        }
        if ((N_pos_true_sum + r_inv * N_neg_true_sum) == 0) continue;
        if ((N_pos_allID_sum + r_inv * N_neg_allID_sum) == 0) continue;
    
        SysWithCovResult r_cov = ComputeSysWithCov(N_pos_allID_sum, N_neg_allID_sum,N_pos_true_sum, N_neg_true_sum,N_pos_omega_sum, N_neg_omega_sum,r_inv);
    
        A_purity_sys_2d_zPt[iz] = r_cov.delta_corrected;
        A_purity_sys_2d_zPt_err[iz] = r_cov.sigma_delta_corrected;
    
        cout << "Purity systematic in bin zPt " << iz + 1
            << ": raw_delta=" << r_cov.delta
            << "  sigma_corrected=" << r_cov.sigma_delta_corrected
            << "  delta_corrected=" << r_cov.delta_corrected << endl;
    }
    cout << "===========================================" << endl;

    vector<double> A_purity_sys_2d_xQ2(nbin_xQ2);
    vector<double> A_purity_sys_2d_xQ2_err(nbin_xQ2);
    
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        double N_pos_true_sum = 0, N_neg_true_sum = 0;
        double N_pos_allID_sum = 0, N_neg_allID_sum = 0;
        double N_pos_omega_sum = 0, N_neg_omega_sum = 0;
    
        for (int iz = 0; iz < nbin_zPt; iz++) {
            N_pos_true_sum += vec_helicity_pos_true[ix][iz].size();
            N_neg_true_sum += vec_helicity_neg_true[ix][iz].size();
            N_pos_allID_sum += vec_helicity_pos_all_ID[ix][iz].size();
            N_neg_allID_sum += vec_helicity_neg_all_ID[ix][iz].size();
            N_pos_omega_sum += vec_helicity_pos_all_preIDcut[ix][iz].size();
            N_neg_omega_sum += vec_helicity_neg_all_preIDcut[ix][iz].size();
        }
        if ((N_pos_true_sum + r_inv * N_neg_true_sum) == 0) continue;
        if ((N_pos_allID_sum + r_inv * N_neg_allID_sum) == 0) continue;
    
        SysWithCovResult r_cov = ComputeSysWithCov(N_pos_allID_sum, N_neg_allID_sum,N_pos_true_sum, N_neg_true_sum,N_pos_omega_sum, N_neg_omega_sum,r_inv);
    
        A_purity_sys_2d_xQ2[ix] = r_cov.delta_corrected;
        A_purity_sys_2d_xQ2_err[ix] = r_cov.sigma_delta_corrected;
    }

    // ===========================================================================


    // =========================================== CSV ===========================================
    for (int iz = 0; iz < nbin_zPt; iz++) {
        double sum_z = 0.0, sum_xB = 0.0, sum_Q2 = 0.0, sum_Pt = 0.0, sum_eps = 0.0;
        for(double val_z : vec_kaonp_z_2d_zPt[iz]) sum_z += val_z;
        double mean_z = sum_z / vec_kaonp_z_2d_zPt[iz].size();
        for(double val_xB : vec_kaonp_xB_2d_zPt[iz]) sum_xB += val_xB;
        double mean_xB = sum_xB / vec_kaonp_xB_2d_zPt[iz].size();
        for(double val_Q2 : vec_kaonp_Q2_2d_zPt[iz]) sum_Q2 += val_Q2;
        double mean_Q2 = sum_Q2 / vec_kaonp_Q2_2d_zPt[iz].size();
        for(double val_Pt : vec_kaonp_Pt_2d_zPt[iz]) sum_Pt += val_Pt;
        double mean_Pt = sum_Pt / vec_kaonp_Pt_2d_zPt[iz].size();
        for(double valeps : vec_kaonp_eps_2d_zPt[iz]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_zPt[iz].size();
        csvFile_zPt << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                << 0 << "," << 0 << "," << A_eff_sys_2d_zPt[iz] << "," << A_eff_sys_2d_zPt_err[iz] << ","
                << A_purity_sys_2d_zPt[iz] << "," << A_purity_sys_2d_zPt_err[iz] << ","
                << A_acc_sys_2d_zPt[iz] << "," << A_acc_sys_2d_zPt_err[iz] << ","
                << A_bm_sys_2d_zPt[iz] << "," << A_bm_sys_2d_zPt_err[iz] << ","
                << Sys_Phi_sin_zPt[iz] << "," << Sys_Phi_sin_err_zPt[iz] << ","
                << Sys_Phi_sin2_zPt[iz] << "," << Sys_Phi_sin2_err_zPt[iz] << ","
                << Sys_Phi_cos_zPt[iz] << "," << Sys_Phi_cos_err_zPt[iz] << ","
                << Sys_Phi_cos2_zPt[iz] << "," << Sys_Phi_cos2_err_zPt[iz] << ","
                << Auu_cos_zPt_mc[iz] << "," << Auu_cos_zPt_err_mc[iz] << ","
                << Auu_cos2_zPt_mc[iz] << "," << Auu_cos2_zPt_err_mc[iz] << ","
                << Auu_cos_zPt_preIDcut[iz] << "," << Auu_cos_zPt_err_preIDcut[iz] << ","
                << Auu_cos2_zPt_preIDcut[iz] << "," << Auu_cos2_zPt_err_preIDcut[iz] << ","
                << Auu_cos_zPt_true[iz] << "," << Auu_cos_zPt_err_true[iz] << ","
                << Auu_cos2_zPt_true[iz] << "," << Auu_cos2_zPt_err_true[iz] << ","
                << Auu_cos_zPt_reco[iz] << "," << Auu_cos_zPt_err_reco[iz] << ","
                << Auu_cos2_zPt_reco[iz] << "," << Auu_cos2_zPt_err_reco[iz] << endl;
    }
    for (int ix = 0; ix < nbin_xQ2; ix++) {
        double sum_z = 0.0, sum_xB = 0.0, sum_Q2 = 0.0, sum_Pt = 0.0, sum_eps = 0.0;
        for(double val_z : vec_kaonp_z_2d[ix]) sum_z += val_z;
        double mean_z = sum_z / vec_kaonp_z_2d[ix].size();
        for(double val_xB : vec_kaonp_xB_2d[ix]) sum_xB += val_xB;
        double mean_xB = sum_xB / vec_kaonp_xB_2d[ix].size();
        for(double val_Q2 : vec_kaonp_Q2_2d[ix]) sum_Q2 += val_Q2;
        double mean_Q2 = sum_Q2 / vec_kaonp_Q2_2d[ix].size();
        for(double val_Pt : vec_kaonp_Pt_2d[ix]) sum_Pt += val_Pt;
        double mean_Pt = sum_Pt / vec_kaonp_Pt_2d[ix].size();
        for(double valeps : vec_kaonp_eps_2d[ix]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d[ix].size();
        csvFile_xQ2 << vec_kaonp_y_2d[ix].size() << "," << ix+1 << "," << "NaN" << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                << 0 << "," << 0 << "," << A_eff_sys_2d_xQ2[ix] << "," << A_eff_sys_2d_xQ2_err[ix] << ","
                << A_purity_sys_2d_xQ2[ix] << "," << A_purity_sys_2d_xQ2_err[ix] << ","
                << A_acc_sys_2d_xQ2[ix] << "," << A_acc_sys_2d_xQ2_err[ix] << ","
                << A_bm_sys_2d_xQ2[ix] << "," << A_bm_sys_2d_xQ2_err[ix] << endl;
        for (int iz = 0; iz < nbin_zPt; iz++) {
            // mean calculations
            double sum_z = 0.0, sum_xB = 0.0, sum_Q2 = 0.0, sum_Pt = 0.0, sum_eps = 0.0;
            for (double val_z : vec_kaonp_z_4d[ix][iz]) sum_z += val_z;
            double mean_z = sum_z / vec_kaonp_z_4d[ix][iz].size();
            for (double val_xB : vec_kaonp_xB_4d[ix][iz]) sum_xB += val_xB;
            double mean_xB = sum_xB / vec_kaonp_xB_4d[ix][iz].size();
            for (double val_Q2 : vec_kaonp_Q2_4d[ix][iz]) sum_Q2 += val_Q2;
            double mean_Q2 = sum_Q2 / vec_kaonp_Q2_4d[ix][iz].size();
            for (double val_Pt : vec_kaonp_Pt_4d[ix][iz]) sum_Pt += val_Pt;
            double mean_Pt = sum_Pt / vec_kaonp_Pt_4d[ix][iz].size();
            for (double valeps : vec_kaonp_eps_4d[ix][iz]) sum_eps += valeps;
            double mean_eps = sum_eps / vec_kaonp_eps_4d[ix][iz].size();
            csvFile << vec_kaonp_y_4d[ix][iz].size() << "," << ix+1 << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                    << efficiency_values[ix][iz] << "," << efficiency_err[ix][iz] << "," << A_eff_sys_4d[ix][iz] << "," << A_eff_sys_err_4d[ix][iz] << ","
                    << A_purity_sys_4d[ix][iz] << "," << A_purity_sys_err_4d[ix][iz] << ","
                    << A_acc_sys_4d[ix][iz] << "," << A_acc_sys_err_4d[ix][iz] << ","
                    << A_binMig_sys_4d[ix][iz] << "," << A_binMig_sys_err_4d[ix][iz] << ","
                    << Sys_Phi_sin[ix][iz] << "," << Sys_Phi_sin_err[ix][iz] << ","
                    << Sys_Phi_sin2[ix][iz] << "," << Sys_Phi_sin2_err[ix][iz] << ","
                    << Sys_Phi_cos[ix][iz] << "," << Sys_Phi_cos_err[ix][iz] << ","
                    << Sys_Phi_cos2[ix][iz] << "," << Sys_Phi_cos2_err[ix][iz] << ","
                    << Auu_cos[ix][iz] << "," << Auu_cos_err[ix][iz] << "," << Auu_cos2[ix][iz] << "," << Auu_cos2_err[ix][iz] << ","
                    << Sys_uu_cont_sin[ix][iz] << "," << Sys_uu_cont_sin2[ix][iz] << endl;
        }
    }

    //outFile.Write();
    //treeKaonP.Write("", TObject::kOverwrite);
    csvFile.close();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << " and csv: " << csv_filename << endl;
}
