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
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rgc/mc_output/fall22_neg/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/MC_fall22_neg
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/work/clas12/lpolizzi/rgc/sum22_mc/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/MC_fall22_neg
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
    int validBins_mc, validBins_reco;
    bool ok; // false se non ci sono abbastanza bin validi per il fit
};
 
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
                          const vector<double>& y_vec,           // inelasticity
                          const vector<double>& depol,       // eps
                          const vector<double>& Ptarget,     // target polarization
                          const vector<double>& helicity,    // beam helicity (+1 / -1)
                          bool use_spinstate_for_helicity = true,
                          const char* campaign = "fall22")
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

    double pBeam = 1.0;
    double dilution = 1.0;  

    // -------- Event loop --------
    for (size_t i = 0; i < phi_h.size(); ++i) {

        const double eps = depol[i];
        const int lam = helicity[i];
        double y = y_vec[i];

        double A = (y*y)/(2*(1-eps));
        double B = ((y*y)/(2*(1-eps))) * eps;
        double C = ((y*y)/(2*(1-eps))) * sqrt(1.0 - eps * eps);
        double V = ((y*y)/(2*(1-eps))) * sqrt(2.0 * eps * (1.0 + eps));
        double W = ((y*y)/(2*(1-eps))) * sqrt(2.0 * eps * (1.0 - eps));

       // --- UU modulation ---
        const double UU_mod = (V/A) * A_UU_cos1 * cos(phi_h[i]) + (B/A) * A_UU_cos2 * cos(2.0 * phi_h[i]);
        // --- LU modulation ---
        const double LU_mod = (W/A) * A_LU_sin * sin(phi_h[i]);
        // --- UL modulation ---
        const double UL_mod = (V/A) * A_sin * sin(phi_h[i]) + (B/A) * A_sin2 * sin(2.0 * phi_h[i]);
        // --- LL modulation ---
        const double LL_mod = (C/A) * A_LL_0 + (W/A) * A_LL_cos1 * cos(phi_h[i]);

        // --- Effective polarization ---
        const double pol = dilution * Ptarget[i];

        // Full spin-dependent term
        const double Seff = pol * (UL_mod + lam * pBeam *LL_mod);

        // --- Likelihood argument ---
        const double arg = 1.0 + UU_mod + lam*pBeam*(LU_mod) + pol*(UL_mod) + lam*pBeam*pol*(LL_mod);

        // Physical protection
        if (arg <= 0.0) return 1e12;
        

        logLike += log(arg);
    }

    return -logLike;
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


inline double safeSigmaDiff(double err_big_sample, double err_small_sample) {
    double diff2 = err_small_sample*err_small_sample + err_big_sample*err_big_sample;
    return (diff2 > 0.0) ? sqrt(diff2) : 0.0;
}

double BarlowSys(double delta, double err){
    double diff = delta * delta - err*err;
    double delta_corrected = (diff > 0.0) ? sqrt(diff) : 0.0;
    return delta_corrected;
}

void rgc_new_mc_analysis(const char* period) {
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
    TString outputFile = Form("RGC_MC_PLOT/studies_mc_%s_kaonp.root", period);
    double torus = -1;

    //const char* csv_filename = "table_RGC_MC_summer22.csv";
    //TString csv_filename = Form("RGC_MC_PLOT/table_RGC_MC_%s.csv", period);
    //TString csv_filename_xQ2 = Form("RGC_MC_PLOT/table_RGC_MC_%s_xQ2.csv", period);
    TString csv_filename_zPt = Form("RGC_MC_PLOT/table_RGC_MC_%s_zPt.csv", period);
    TString csv_filename_zPt_test = Form("RGC_MC_PLOT/table_RGC_MC_%s_zPt_test.csv", period);
    //std::ofstream csvFile(csv_filename.Data());
    //std::ofstream csvFile_xQ2(csv_filename_xQ2.Data());
    std::ofstream csvFile_zPt(csv_filename_zPt.Data());
    std::ofstream csvFile_zPt_test(csv_filename_zPt_test.Data());
    //csvFile << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err,Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err, Phih_sys_sinx, Phih_sys_sinx_err, Phih_sys_sin2x, Phih_sys_sin2x_err, Phih_sys_cosx, Phih_sys_cosx_err, Phih_sys_cos2x, Phih_sys_cos2x_err, Auu_cos, Auu_cos_err, Auu_cos2, Auu_cos2_err, Auu_cont_sinx, Auu_cont_sin2x\n";
    //csvFile_xQ2 << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err, Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err, Phih_sys_sinx, Phih_sys_sinx_err, Phih_sys_sin2x, Phih_sys_sin2x_err, Phih_sys_cosx, Phih_sys_cosx_err, Phih_sys_cos2x, Phih_sys_cos2x_err, Auu_cos, Auu_cos_err, Auu_cos2, Auu_cos2_err, Auu_cont_sinx, Auu_cont_sin2x\n";
    csvFile_zPt << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err, Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err\n";
    csvFile_zPt_test << "n_events,binxBQ2,bin_zPt,mean_xB,mean_Q2,mean_z,mean_PhT,epsilon,"
            << "modulation,"
            << "Acc_sys,Acc_sys_err,PID_sys,PID_sys_err,Purity_sys,Purity_sys_err,Bin_mig_sys,Bin_mig_sys_err\n";

    TFile outFile(outputFile.Data(), "RECREATE");  // File di output ROOT
    TTree treeKaonP("Kaon+", "");
    TTree MC_treeKaonP("MC_Kaon+", "");

    TDirectory* dir_xQ2 = outFile.mkdir("Binning xB-Q2");
    TDirectory* dir_zPt = outFile.mkdir("Binning z-Pt");
    TDirectory* dir_aul_sin = outFile.mkdir("AUL sin(Phi)");
    TDirectory* dir_aul_sin2 = outFile.mkdir("AUL sin(2Phi)");
    TDirectory* dir_all_0 = outFile.mkdir("ALL");
    TDirectory* dir_all_cos = outFile.mkdir("ALL cos(Phi)");
    TDirectory* dir_alu_sin = outFile.mkdir("ALU sin(Phi)");
    TDirectory* dir_eff = outFile.mkdir("Efficiency");
    TDirectory* dir_bin_mig = outFile.mkdir("Bin migration");
    TDirectory* dir_Phih_sin = outFile.mkdir("Phi_{h} systematics of sin(x)");
    TDirectory* dir_Phih_sin2 = outFile.mkdir("Phi_{h} systematics of sin(2x)");
    TDirectory* dir_Phih_cos = outFile.mkdir("Phi_{h} systematics of cos(x)");
    TDirectory* dir_Phih_cos2 = outFile.mkdir("Phi_{h} systematics of cos(2x)");

    double kaonp_xB_recoMC, kaonp_Q2_recoMC, kaonp_z_recoMC, kaonp_Pt_recoMC;
    double kaonp_Phi_h_recoMC, kaonp_Theta_recoMC, kaonp_Phi_recoMC;
    double kaonp_y_recoMC, kaonp_epsilon_recoMC;

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
    chainKaonP.SetBranchAddress("epsilon_mc", &kaonp_epsilon_recoMC);
    chainKaonP.SetBranchAddress("W", &kaonp_W);
    chainKaonP.SetBranchAddress("Q2", &kaonp_Q2);
    chainKaonP.SetBranchAddress("xF", &kaonp_xF);
    chainKaonP.SetBranchAddress("xB", &kaonp_xB);
    chainKaonP.SetBranchAddress("y", &kaonp_y);
    chainKaonP.SetBranchAddress("y_mc", &kaonp_y_recoMC);
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
    TH1D kp_fmax ("_fmax", "f_{max} injection function ; f_{max}; count", 300, 0, 2);
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
    // RECO
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
    // MC
    vector<vector<double>> vec_kaonp_phih_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_helicity_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_2phih_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_z_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_Pt_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_xB_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_Q2_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_y_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_pol_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_eps_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_spin_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_sintheta_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_phih_plus_2d_mc(nbin_xQ2);
    vector<vector<double>> vec_kaonp_phih_minus_2d_mc(nbin_xQ2);
    // allID
    vector<vector<double>> vec_kaonp_phih_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_all_ID(nbin_zPt);
    // preID
    vector<vector<double>> vec_kaonp_phih_2d_preID(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_preID(nbin_zPt);
    // TRUE
    vector<vector<double>> vec_kaonp_phih_2d_true(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_true(nbin_zPt);
    // zPt bins
    // RECO
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
    // MC
    vector<vector<double>> vec_kaonp_phih_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_zP_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_zPt_mc(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_zPt_mc(nbin_zPt);
    // allID
    vector<vector<double>> vec_kaonp_phih_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_zP_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_zPt_all_ID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_zPt_all_ID(nbin_zPt);
    // preID
    vector<vector<double>> vec_kaonp_phih_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_zP_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_zPt_preID(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_zPt_preID(nbin_zPt);
    // TRUE
    vector<vector<double>> vec_kaonp_phih_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_helicity_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_2phih_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_z_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_Pt_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_xB_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_Q2_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_y_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_pol_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_eps_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_spin_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_sintheta_2d_zP_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_plus_2d_zPt_true(nbin_zPt);
    vector<vector<double>> vec_kaonp_phih_minus_2d_zPt_true(nbin_zPt);
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
        // mc
    vector<double> AUL_sin_2d_mc(nbin_xQ2);
    vector<double> AUL_sin_err_2d_mc(nbin_xQ2);
    vector<double> AUL_2sin_2d_mc(nbin_xQ2);
    vector<double> AUL_2sin_err_2d_mc(nbin_xQ2);
        // true
    vector<double> AUL_sin_2d_true(nbin_xQ2);
    vector<double> AUL_sin_err_2d_true(nbin_xQ2);
    vector<double> AUL_2sin_2d_true(nbin_xQ2);
    vector<double> AUL_2sin_err_2d_true(nbin_xQ2);
        // all ID
    vector<double> AUL_sin_2d_all_ID(nbin_xQ2);
    vector<double> AUL_sin_err_2d_all_ID(nbin_xQ2);
    vector<double> AUL_2sin_2d_all_ID(nbin_xQ2);
    vector<double> AUL_2sin_err_2d_all_ID(nbin_xQ2);
        // pre ID
    vector<double> AUL_sin_2d_preID(nbin_xQ2);
    vector<double> AUL_sin_err_2d_preID(nbin_xQ2);
    vector<double> AUL_2sin_2d_preID(nbin_xQ2);
    vector<double> AUL_2sin_err_2d_preID(nbin_xQ2);
        // reco
    vector<double> AUL_sin_2d_reco(nbin_xQ2);
    vector<double> AUL_sin_err_2d_reco(nbin_xQ2);
    vector<double> AUL_2sin_2d_reco(nbin_xQ2);
    vector<double> AUL_2sin_err_2d_reco(nbin_xQ2);
    // LL
    vector<double> ALL_0_2d(nbin_xQ2);
    vector<double> ALL_0_err_2d(nbin_xQ2);
    vector<double> ALL_cos_2d(nbin_xQ2);
    vector<double> ALL_cos_err_2d(nbin_xQ2);
        // mc
    vector<double> ALL_0_2d_mc(nbin_xQ2);
    vector<double> ALL_0_err_2d_mc(nbin_xQ2);
    vector<double> ALL_cos_2d_mc(nbin_xQ2);
    vector<double> ALL_cos_err_2d_mc(nbin_xQ2);
        // true
    vector<double> ALL_0_2d_true(nbin_xQ2);
    vector<double> ALL_0_err_2d_true(nbin_xQ2);
    vector<double> ALL_cos_2d_true(nbin_xQ2);
    vector<double> ALL_cos_err_2d_true(nbin_xQ2);
        // all ID
    vector<double> ALL_0_2d_all_ID(nbin_xQ2);
    vector<double> ALL_0_err_2d_all_ID(nbin_xQ2);
    vector<double> ALL_cos_2d_all_ID(nbin_xQ2);
    vector<double> ALL_cos_err_2d_all_ID(nbin_xQ2);
        // preID
    vector<double> ALL_0_2d_preID(nbin_xQ2);
    vector<double> ALL_0_err_2d_preID(nbin_xQ2);
    vector<double> ALL_cos_2d_preID(nbin_xQ2);
    vector<double> ALL_cos_err_2d_preID(nbin_xQ2);
        // reco
    vector<double> ALL_0_2d_reco(nbin_xQ2);
    vector<double> ALL_0_err_2d_reco(nbin_xQ2);
    vector<double> ALL_cos_2d_reco(nbin_xQ2);
    vector<double> ALL_cos_err_2d_reco(nbin_xQ2);
    // UL vs z
    vector<double> AUL_sin_2d_zPt(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt(nbin_zPt);
        // mc 
    vector<double> AUL_sin_2d_zPt_mc(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt_mc(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt_mc(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt_mc(nbin_zPt);
        //true
    vector<double> AUL_sin_2d_zPt_true(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt_true(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt_true(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt_true(nbin_zPt);
        // all ID
    vector<double> AUL_sin_2d_zPt_all_ID(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt_all_ID(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt_all_ID(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt_all_ID(nbin_zPt);
        // preID
    vector<double> AUL_sin_2d_zPt_preID(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt_preID(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt_preID(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt_preID(nbin_zPt);
        // reco
    vector<double> AUL_sin_2d_zPt_reco(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt_reco(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt_reco(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt_reco(nbin_zPt);
    // LL vs z
    vector<double> ALL_0_2d_zPt(nbin_zPt);
    vector<double> ALL_0_err_2d_zPt(nbin_zPt);
    vector<double> ALL_cos_2d_zPt(nbin_zPt);
    vector<double> ALL_cos_err_2d_zPt(nbin_zPt);
      // mc
    vector<double> ALL_0_2d_zPt_mc(nbin_zPt);
    vector<double> ALL_0_err_2d_zPt_mc(nbin_zPt);
    vector<double> ALL_cos_2d_zPt_mc(nbin_zPt);
    vector<double> ALL_cos_err_2d_zPt_mc(nbin_zPt);
      // true
    vector<double> ALL_0_2d_zPt_true(nbin_zPt);
    vector<double> ALL_0_err_2d_zPt_true(nbin_zPt);
    vector<double> ALL_cos_2d_zPt_true(nbin_zPt);
    vector<double> ALL_cos_err_2d_zPt_true(nbin_zPt);
        // all ID
    vector<double> ALL_0_2d_zPt_all_ID(nbin_zPt);
    vector<double> ALL_0_err_2d_zPt_all_ID(nbin_zPt);
    vector<double> ALL_cos_2d_zPt_all_ID(nbin_zPt);
    vector<double> ALL_cos_err_2d_zPt_all_ID(nbin_zPt);
        // pre ID
    vector<double> ALL_0_2d_zPt_preID(nbin_zPt);
    vector<double> ALL_0_err_2d_zPt_preID(nbin_zPt);
    vector<double> ALL_cos_2d_zPt_preID(nbin_zPt);
    vector<double> ALL_cos_err_2d_zPt_preID(nbin_zPt);
        // reco
    vector<double> ALL_0_2d_zPt_reco(nbin_zPt);
    vector<double> ALL_0_err_2d_zPt_reco(nbin_zPt);
    vector<double> ALL_cos_2d_zPt_reco(nbin_zPt);
    vector<double> ALL_cos_err_2d_zPt_reco(nbin_zPt);
    // LU
    vector<double> ALU_sin_2d(nbin_xQ2);
    vector<double> ALU_sin_err_2d(nbin_xQ2);
    vector<double> ALU_sin_2d_zPt(nbin_zPt);
    vector<double> ALU_sin_err_2d_zPt(nbin_zPt);
      // mc
    vector<double> ALU_sin_2d_mc(nbin_xQ2);
    vector<double> ALU_sin_err_2d_mc(nbin_xQ2);
    vector<double> ALU_sin_2d_zPt_mc(nbin_zPt);
    vector<double> ALU_sin_err_2d_zPt_mc(nbin_zPt);
      // true
    vector<double> ALU_sin_2d_true(nbin_xQ2);
    vector<double> ALU_sin_err_2d_true(nbin_xQ2);
    vector<double> ALU_sin_2d_zPt_true(nbin_zPt);
    vector<double> ALU_sin_err_2d_zPt_true(nbin_zPt);
        // all ID
    vector<double> ALU_sin_2d_all_ID(nbin_xQ2);
    vector<double> ALU_sin_err_2d_all_ID(nbin_xQ2);
    vector<double> ALU_sin_2d_zPt_all_ID(nbin_zPt);
    vector<double> ALU_sin_err_2d_zPt_all_ID(nbin_zPt);
        // pre ID
    vector<double> ALU_sin_2d_preID(nbin_xQ2);
    vector<double> ALU_sin_err_2d_preID(nbin_xQ2);
    vector<double> ALU_sin_2d_zPt_preID(nbin_zPt);
    vector<double> ALU_sin_err_2d_zPt_preID(nbin_zPt);
        // reco
    vector<double> ALU_sin_2d_reco(nbin_xQ2);
    vector<double> ALU_sin_err_2d_reco(nbin_xQ2);
    vector<double> ALU_sin_2d_zPt_reco(nbin_zPt);
    vector<double> ALU_sin_err_2d_zPt_reco(nbin_zPt);
    //
    vector<vector<double>> AUL_sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    
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


    double A_UL_sin_inject  = 0.1;
    double A_UL_2sin_inject = 0.1;
    double A_LL_0_inject    = 0.0;
    double A_LL_cos_inject  = 0.10;
    double A_LU_sin_inject = 0.05;
    double dilution = 0.25;
    double beam_pol = 0.90;

    double f_max_inject = 1.0 + dilution*0.85*(std::abs(A_UL_sin_inject)+std::abs(A_UL_2sin_inject))
                            + dilution*0.85*beam_pol*(std::abs(A_LL_0_inject)+std::abs(A_LL_cos_inject))
                            + beam_pol*(std::abs(A_LU_sin_inject));

    TRandom3 rng(12345);
    // MC
    Long64_t nEntries_mc = MC_chainKaonP.GetEntries();
    vector<bool>   keep_injected(MC_chainKaonP.GetEntries());
    vector<double> fake_Ptarget_vec(MC_chainKaonP.GetEntries());   // salvato per riuso nel fit
    for (Long64_t i = 0; i < nEntries_mc; i++) {
        MC_chainKaonP.GetEntry(i);
        if(std::fabs(kaonp_Phi_mc) < 2.5 && std::fabs(kaonp_Phi_mc) > 0.65 && kaonp_Ph_mc >= 3) continue; 
        double index_xQ2 = getBinIndex_xQ2(kaonp_xB_mc, kaonp_Q2_mc); 
        double index_zPt = getBinIndex_zPt(kaonp_z_mc, kaonp_PhT_mc);
        double kp_mc_theta_deg = kaonp_Theta_mc * 180.0 / TMath::Pi();

        // spin injection 
        int fake_spinstate = (rng.Rndm() < 0.5) ? +1 : -1;   // deciso UNA VOLTA per questo evento
        double fake_Pt = 1.0 * fake_spinstate;
        fake_Ptarget_vec[i] = fake_Pt;

        double eps = kaonp_epsilon_mc, y = kaonp_y_mc;
        double A = (y*y)/(2*(1-eps));
        double B = A*eps, C = A*sqrt(1-eps*eps);
        double V = A*sqrt(2*eps*(1+eps)), W = A*sqrt(2*eps*(1-eps));

        double UL_mod = (V/A)*A_UL_sin_inject*sin(kaonp_Phi_h_mc) + (B/A)*A_UL_2sin_inject*sin(2*kaonp_Phi_h_mc);
        double LL_mod = (C/A)*A_LL_0_inject + (W/A)*A_LL_cos_inject*cos(kaonp_Phi_h_mc);
        double LU_mod = (W/A) * A_LU_sin_inject * sin(kaonp_Phi_h_mc);

        double f = 1.0 + helicity_mc*LU_mod + fake_Pt*UL_mod + helicity_mc*fake_Pt*LL_mod;  // SOLO UL/LL iniettate

        double u = f_max_inject * rng.Rndm();   // un numero casuale DIVERSO per OGNI evento
        keep_injected[i] = (u < f);
        // 
        if(index_xQ2 >= 0 && index_zPt >= 0) {
            vec_kaonp_xB_4d_mc[index_xQ2-1][index_zPt-1].push_back(kaonp_xB_mc);
            if(keep_injected[i]){
                // xBQ2
                vec_kaonp_phih_2d_mc[index_xQ2-1].push_back(kaonp_Phi_h_mc);
                vec_helicity_2d_mc[index_xQ2-1].push_back(helicity_mc);
                vec_kaonp_2phih_2d_mc[index_xQ2-1].push_back(2*kaonp_Phi_h_mc);
                vec_kaonp_z_2d_mc[index_xQ2-1].push_back(kaonp_z_mc);
                vec_kaonp_Pt_2d_mc[index_xQ2-1].push_back(kaonp_PhT_mc);
                vec_kaonp_xB_2d_mc[index_xQ2-1].push_back(kaonp_xB_mc);
                vec_kaonp_Q2_2d_mc[index_xQ2-1].push_back(kaonp_Q2_mc);
                vec_kaonp_y_2d_mc[index_xQ2-1].push_back(kaonp_y_mc);
                vec_kaonp_pol_2d_mc[index_xQ2-1].push_back(fake_Pt);
                vec_kaonp_eps_2d_mc[index_xQ2-1].push_back(kaonp_epsilon_mc);
                vec_kaonp_spin_2d_mc[index_xQ2-1].push_back(fake_spinstate);
                // zPt
                vec_kaonp_phih_2d_zPt_mc[index_zPt-1].push_back(kaonp_Phi_h_mc);
                vec_helicity_2d_zPt_mc[index_zPt-1].push_back(helicity_mc);
                vec_kaonp_2phih_2d_zPt_mc[index_zPt-1].push_back(2*kaonp_Phi_h_mc);
                vec_kaonp_z_2d_zPt_mc[index_zPt-1].push_back(kaonp_z_mc);
                vec_kaonp_Pt_2d_zPt_mc[index_zPt-1].push_back(kaonp_PhT_mc);
                vec_kaonp_xB_2d_zPt_mc[index_zPt-1].push_back(kaonp_xB_mc);
                vec_kaonp_Q2_2d_zPt_mc[index_zPt-1].push_back(kaonp_Q2_mc);
                vec_kaonp_y_2d_zPt_mc[index_zPt-1].push_back(kaonp_y_mc);
                vec_kaonp_pol_2d_zPt_mc[index_zPt-1].push_back(fake_Pt);
                vec_kaonp_eps_2d_zPt_mc[index_zPt-1].push_back(kaonp_epsilon_mc);
                vec_kaonp_spin_2d_zPt_mc[index_zPt-1].push_back(fake_spinstate);
            }
            if(helicity_mc == 1) {
                vec_helicity_pos_mc[index_xQ2-1][index_zPt-1].push_back(helicity_mc);

            } else if(helicity_mc == -1){
                vec_helicity_neg_mc[index_xQ2-1][index_zPt-1].push_back(helicity_mc);
            }
        }
    }

    // RECO
    // Kaon+
    Long64_t nEntries_kp = chainKaonP.GetEntries();
    vector<bool>   keep_injected_reco(chainKaonP.GetEntries());
    vector<double> fake_Ptarget_vec_reco(chainKaonP.GetEntries());
    for (Long64_t i = 0; i < nEntries_kp; i++) {
        // RICH cut
        //if (electron_ass != 321) continue;
        chainKaonP.GetEntry(i);
        if(std::fabs(kaonp_Phi) < 2.5 && std::fabs(kaonp_Phi) > 0.65 && kaonp_Ph >= 3) continue;

        // ======================================================================================================================
        // spin injection 
        int fake_spinstate = (rng.Rndm() < 0.5) ? +1 : -1;   // deciso UNA VOLTA per questo evento
        double fake_Pt = 1.0 * fake_spinstate;
        fake_Ptarget_vec_reco[i] = fake_Pt;

        double eps = kaonp_epsilon_recoMC, y = kaonp_y_recoMC;
        double A = (y*y)/(2*(1-eps));
        double B = A*eps, C = A*sqrt(1-eps*eps);
        double V = A*sqrt(2*eps*(1+eps)), W = A*sqrt(2*eps*(1-eps));

        double UL_mod = (V/A)*A_UL_sin_inject*sin(kaonp_Phi_h_recoMC) + (B/A)*A_UL_2sin_inject*sin(2*kaonp_Phi_h_recoMC);
        double LL_mod = (C/A)*A_LL_0_inject + (W/A)*A_LL_cos_inject*cos(kaonp_Phi_h_recoMC);
        double LU_mod = (W/A) * A_LU_sin_inject * sin(kaonp_Phi_h_recoMC);

        double f = 1.0 + helicity*LU_mod + fake_Pt*UL_mod + helicity*fake_Pt*LL_mod;  // SOLO UL/LL iniettate

        double u = f_max_inject * rng.Rndm();   // un numero casuale DIVERSO per OGNI evento
        keep_injected_reco[i] = (u < f);
        // ======================================================================================================================

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

                if(keep_injected_reco[i]){
                    // xBQ2
                    vec_kaonp_phih_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(kaonp_Phi_h_recoMC);
                    vec_helicity_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(helicity);
                    vec_kaonp_2phih_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(2*kaonp_Phi_h_recoMC);
                    vec_kaonp_z_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(kaonp_z_recoMC);
                    vec_kaonp_Pt_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(kaonp_Pt_recoMC);
                    vec_kaonp_xB_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(kaonp_xB_recoMC);
                    vec_kaonp_Q2_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(kaonp_Q2_recoMC);
                    vec_kaonp_y_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(kaonp_y_recoMC);
                    vec_kaonp_pol_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(fake_Pt);
                    vec_kaonp_eps_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(kaonp_epsilon_recoMC);
                    vec_kaonp_spin_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(fake_spinstate);
                    // zPt
                    vec_kaonp_phih_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(kaonp_Phi_h_recoMC);
                    vec_helicity_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(helicity);
                    vec_kaonp_2phih_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(2*kaonp_Phi_h_recoMC);
                    vec_kaonp_z_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(kaonp_z_recoMC);
                    vec_kaonp_Pt_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(kaonp_Pt_recoMC);
                    vec_kaonp_xB_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(kaonp_xB_recoMC);
                    vec_kaonp_Q2_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(kaonp_Q2_recoMC);
                    vec_kaonp_y_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(kaonp_y_recoMC);
                    vec_kaonp_pol_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(fake_Pt);
                    vec_kaonp_eps_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(kaonp_epsilon_recoMC);
                    vec_kaonp_spin_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(fake_spinstate);
                }
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
            if(helicity > 0) vec_helicity_pos_preIDcut[mc_index_xQ2_preCut-1][mc_index_zPt_preCut-1].push_back(helicity);
            else vec_helicity_neg_preIDcut[mc_index_xQ2_preCut-1][mc_index_zPt_preCut-1].push_back(helicity);
            if(keep_injected_reco[i]){
                // xBQ2
                vec_kaonp_phih_2d_preID[mc_index_xQ2_preCut-1].push_back(kaonp_Phi_h_recoMC);
                vec_helicity_2d_preID[mc_index_xQ2_preCut-1].push_back(helicity);
                vec_kaonp_2phih_2d_preID[mc_index_xQ2_preCut-1].push_back(2*kaonp_Phi_h_recoMC);
                vec_kaonp_z_2d_preID[mc_index_xQ2_preCut-1].push_back(kaonp_z_recoMC);
                vec_kaonp_Pt_2d_preID[mc_index_xQ2_preCut-1].push_back(kaonp_Pt_recoMC);
                vec_kaonp_xB_2d_preID[mc_index_xQ2_preCut-1].push_back(kaonp_xB_recoMC);
                vec_kaonp_Q2_2d_preID[mc_index_xQ2_preCut-1].push_back(kaonp_Q2_recoMC);
                vec_kaonp_y_2d_preID[mc_index_xQ2_preCut-1].push_back(kaonp_y_recoMC);
                vec_kaonp_pol_2d_preID[mc_index_xQ2_preCut-1].push_back(fake_Pt);
                vec_kaonp_eps_2d_preID[mc_index_xQ2_preCut-1].push_back(kaonp_epsilon_recoMC);
                vec_kaonp_spin_2d_preID[mc_index_xQ2_preCut-1].push_back(fake_spinstate);
                // zPt
                vec_kaonp_phih_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(kaonp_Phi_h_recoMC);
                vec_helicity_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(helicity);
                vec_kaonp_2phih_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(2*kaonp_Phi_h_recoMC);
                vec_kaonp_z_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(kaonp_z_recoMC);
                vec_kaonp_Pt_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(kaonp_Pt_recoMC);
                vec_kaonp_xB_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(kaonp_xB_recoMC);
                vec_kaonp_Q2_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(kaonp_Q2_recoMC);
                vec_kaonp_y_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(kaonp_y_recoMC);
                vec_kaonp_pol_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(fake_Pt);
                vec_kaonp_eps_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(kaonp_epsilon_recoMC);
                vec_kaonp_spin_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(fake_spinstate);
            }
        }

        // ID
        if(1==1){
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
                    if(keep_injected_reco[i]){
                    // xBQ2
                    vec_kaonp_phih_2d_true[mc_index_xQ2-1].push_back(kaonp_Phi_h_recoMC);
                    vec_helicity_2d_true[mc_index_xQ2-1].push_back(helicity);
                    vec_kaonp_2phih_2d_true[mc_index_xQ2-1].push_back(2*kaonp_Phi_h_recoMC);
                    vec_kaonp_z_2d_true[mc_index_xQ2-1].push_back(kaonp_z_recoMC);
                    vec_kaonp_Pt_2d_true[mc_index_xQ2-1].push_back(kaonp_Pt_recoMC);
                    vec_kaonp_xB_2d_true[mc_index_xQ2-1].push_back(kaonp_xB_recoMC);
                    vec_kaonp_Q2_2d_true[mc_index_xQ2-1].push_back(kaonp_Q2_recoMC);
                    vec_kaonp_y_2d_true[mc_index_xQ2-1].push_back(kaonp_y_recoMC);
                    vec_kaonp_pol_2d_true[mc_index_xQ2-1].push_back(fake_Pt);
                    vec_kaonp_eps_2d_true[mc_index_xQ2-1].push_back(kaonp_epsilon_recoMC);
                    vec_kaonp_spin_2d_true[mc_index_xQ2-1].push_back(fake_spinstate);
                    // zPt
                    vec_kaonp_phih_2d_zPt_true[mc_index_zPt-1].push_back(kaonp_Phi_h_recoMC);
                    vec_helicity_2d_zPt_true[mc_index_zPt-1].push_back(helicity);
                    vec_kaonp_2phih_2d_zPt_true[mc_index_zPt-1].push_back(2*kaonp_Phi_h_recoMC);
                    vec_kaonp_z_2d_zPt_true[mc_index_zPt-1].push_back(kaonp_z_recoMC);
                    vec_kaonp_Pt_2d_zPt_true[mc_index_zPt-1].push_back(kaonp_Pt_recoMC);
                    vec_kaonp_xB_2d_zPt_true[mc_index_zPt-1].push_back(kaonp_xB_recoMC);
                    vec_kaonp_Q2_2d_zPt_true[mc_index_zPt-1].push_back(kaonp_Q2_recoMC);
                    vec_kaonp_y_2d_zPt_true[mc_index_zPt-1].push_back(kaonp_y_recoMC);
                    vec_kaonp_pol_2d_zPt_true[mc_index_zPt-1].push_back(fake_Pt);
                    vec_kaonp_eps_2d_zPt_true[mc_index_zPt-1].push_back(kaonp_epsilon_recoMC);
                    vec_kaonp_spin_2d_zPt_true[mc_index_zPt-1].push_back(fake_spinstate);
                    }
                }
                if(index_xQ2 >= 0) {
                    if(keep_injected_reco[i]){
                        vec_kaonp_phih_2d[index_xQ2-1].push_back(kaonp_Phi_h);
                        vec_helicity_2d[index_xQ2-1].push_back(helicity);
                        vec_kaonp_2phih_2d[index_xQ2-1].push_back(2*kaonp_Phi_h);
                        vec_kaonp_z_2d[index_xQ2-1].push_back(kaonp_z);
                        vec_kaonp_Pt_2d[index_xQ2-1].push_back(kaonp_Pt);
                        vec_kaonp_xB_2d[index_xQ2-1].push_back(kaonp_xB);
                        vec_kaonp_Q2_2d[index_xQ2-1].push_back(kaonp_Q2);
                        vec_kaonp_y_2d[index_xQ2-1].push_back(kaonp_y);
                        vec_kaonp_pol_2d[index_xQ2-1].push_back(fake_Pt);
                        vec_kaonp_eps_2d[index_xQ2-1].push_back(kaonp_epsilon);
                        vec_kaonp_spin_2d[index_xQ2-1].push_back(fake_spinstate);
                        vec_kaonp_sintheta_2d[index_xQ2-1].push_back(kaonp_sintheta);
                    }
                }
                if(index_zPt >= 0){ 
                    if(keep_injected_reco[i]){
                        vec_kaonp_phih_2d_zPt[index_zPt-1].push_back(kaonp_Phi_h);
                        vec_helicity_2d_zPt[index_zPt-1].push_back(helicity);
                        vec_kaonp_2phih_2d_zPt[index_zPt-1].push_back(2*kaonp_Phi_h);
                        vec_kaonp_z_2d_zPt[index_zPt-1].push_back(kaonp_z);
                        vec_kaonp_Pt_2d_zPt[index_zPt-1].push_back(kaonp_Pt);
                        vec_kaonp_xB_2d_zPt[index_zPt-1].push_back(kaonp_xB);
                        vec_kaonp_Q2_2d_zPt[index_zPt-1].push_back(kaonp_Q2);
                        vec_kaonp_y_2d_zPt[index_zPt-1].push_back(kaonp_y);
                        vec_kaonp_pol_2d_zPt[index_zPt-1].push_back(fake_Pt);
                        vec_kaonp_eps_2d_zPt[index_zPt-1].push_back(kaonp_epsilon);
                        //if (index_zPt == 16) cout << kaonp_epsilon << endl;
                        vec_kaonp_spin_2d_zPt[index_zPt-1].push_back(fake_spinstate);
                        vec_kaonp_sintheta_2d_zPt[index_zPt-1].push_back(kaonp_sintheta);
                    }
                    if(index_xQ2 >= 0 && index_zPt >= 0){
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
                        } else {
                            vec_helicity_neg_reco[index_xQ2-1][index_zPt-1].push_back(helicity);
                        }
                    }
                }
                double kp_theta_deg = kaonp_Theta * 180.0 / TMath::Pi();
                //
                double beta_th = kaonp_Ph/(sqrt(ph2+0.01948816));
                double delta_beta = kaonp_beta - beta_th;

                
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
    
    std::vector<double> Auu_cos_mc = {
        0.616385, 0.799272, 0.688099, 0.647754, 0.694137,
        0.672366, 0.717932, -0.415006, 0.27936, 0.603572,
        0.807656, 0.852119, 0.917368, 1, -1,
        -0.532044, -0.0623458, 0.328036, 0.579722, 0.795063,
        0.853295, -0.975872, -0.294294, 0.0660673, 0.347074
    };

    std::vector<double> Auu_cos_err_mc = {
        0.00814018, 0.00880552, 0.0117463, 0.0172805, 0.0239639,
        0.0388183, 0.0742917, 0.00837456, 0.00745188, 0.00824765,
        0.0104343, 0.0422769, 0.0225419, 0.0387944, 0.00111135,
        0.0110645, 0.0118706, 0.0147915, 0.018933, 0.0294045,
        0.155625, 0.0215164, 0.0378033, 0.0507369, 0.0715543
    };

    std::vector<double> Auu_cos2_mc = {
        -0.405209, -0.234859, -0.29282, -0.263625, -0.201495,
        -0.182861, -0.105578, -0.0688026, -0.361847, -0.420302,
        -0.516952, -0.412946, -0.389459, -0.480769, 0.223943,
        0.00590098, -0.193538, -0.323644, -0.435759, -0.474997,
        -0.261821, 0.142675, -0.122084, -0.0263655, -0.177311
    };

    std::vector<double> Auu_cos2_err_mc = {
        0.0088949, 0.0102285, 0.0133053, 0.0192356, 0.0271033,
        0.0433261, 0.0839511, 0.00865943, 0.00789726, 0.00872331,
        0.0165172, 0.041117, 0.0399073, 0.106548, 0.0133754,
        0.0114417, 0.0122231, 0.0154746, 0.0200113, 0.0442504,
        0.154216, 0.0234382, 0.0398114, 0.0511441, 0.0731267
    };

    std::vector<double> Auu_cos_preID = {
        1, 1, 0.999999, -0.872938, 0.855044,
        0.999994, 0, 0.166777, 0.518392, 0.848694,
        0.721374, 1, -0.424213, -0.0176376, -0.726464,
        -0.0151637, 0.346015, 0.792786, -0.194847, 0.999995,
        0.160007, -0.429706, 0.163811, 0.373369, -0.197517
    };

    std::vector<double> Auu_cos_err_preID = {
        0.0237486, 0.0378679, 0.17867, 1.02525, 1.57807,
        1.97756, 1.25362, 0.0440704, 0.0495596, 0.115206,
        0.423672, 0.327941, 1.23846, 1.37397, 0.0775626,
        0.0706797, 0.0919822, 0.153389, 0.286238, 1.61251,
        0.705899, 0.116161, 0.249656, 0.418152, 1.37638
    };

    std::vector<double> Auu_cos2_preID = {
        -0.36093, 0.114696, -0.379932, 0.532023, -1,
        -0.549307, 0, -0.294929, -0.31704, -0.259763,
        -0.120545, -0.365729, 0.522757, 0.00684932, -0.0660886,
        0.085462, -0.374534, -0.00797941, 0.0527759, -0.769554,
        0.324263, -0.251983, -0.0309429, -0.0848608, 0.0358296
    };

    std::vector<double> Auu_cos2_err_preID = {
        0.105018, 0.118073, 0.262468, 0.242087, 1.97763,
        1.21967, 1.45857, 0.0462542, 0.0536927, 0.123893,
        0.295568, 0.290059, 0.323423, 1.48276, 0.0839862,
        0.0656973, 0.0977291, 0.162618, 0.312981, 0.495981,
        0.609089, 0.127859, 0.225743, 0.403173, 1.45042
    };

    std::vector<double> Auu_cos_true = {
        0.333875, 0.148239, -0.0225515, -0.144546, -0.188514,
        -0.381285, -0.402697, -0.140435, 0.371473, 0.220551,
        0.024784, -0.279038, -0.680618, -0.506308, -1,
        -0.288023, -0.0392597, 0.0941313, -0.223707, -0.716473,
        -0.193431, -1, -0.501112, -0.103997, -0.198347
    };

    std::vector<double> Auu_cos_err_true = {
        0.0122332, 0.0129728, 0.0164847, 0.0226721, 0.0313557,
        0.0399235, 0.035392, 0.0116768, 0.00932344, 0.0112533,
        0.0152335, 0.0205348, 0.0229342, 0.0294867, 0.00411767,
        0.0152045, 0.0172895, 0.0207423, 0.0255907, 0.0277308,
        0.0727163, 0.0126026, 0.0631858, 0.0631888, 0.134961
    };

    std::vector<double> Auu_cos2_true = {
        -0.541371, -0.410531, -0.193605, -0.0405717, 0.014279,
        0.144683, 0.105171, -0.341814, -0.510568, -0.534217,
        -0.468811, -0.29916, -0.00899866, -0.161152, 0.143266,
        -0.186371, -0.445203, -0.617095, -0.428775, -0.123521,
        -0.381913, 0.117672, 0.0404123, -0.509206, -0.396153
    };

    std::vector<double> Auu_cos2_err_true = {
        0.0133934, 0.0139092, 0.0171701, 0.0235466, 0.0323847,
        0.0394398, 0.0356696, 0.0121204, 0.00982446, 0.0119817,
        0.0162429, 0.0223289, 0.0247766, 0.0326195, 0.0152641,
        0.0155019, 0.0179767, 0.0215692, 0.0267114, 0.0307573,
        0.07776, 0.024321, 0.061185, 0.0641711, 0.144453
    };

    std::vector<double> Auu_cos_reco = {
        0.445845, 0.27578, 0.136039, 0.0318084, -0.0101453,
        -0.182544, -0.311684, -0.110251, 0.443308, 0.293468,
        0.102209, -0.181095, -0.619729, -0.476427, -1,
        -0.257957, 0.00644669, 0.156743, -0.170998, -0.696744,
        -0.206659, -1, -0.499426, -0.0901324, -0.2262
    };

    std::vector<double> Auu_cos_err_reco = {
        0.0118678, 0.0126565, 0.0162117, 0.0225625, 0.0313852,
        0.0405663, 0.0360168, 0.0117766, 0.00932016, 0.0111982,
        0.0151809, 0.0208212, 0.0235976, 0.0295602, 0.00291723,
        0.0153678, 0.0173856, 0.0209748, 0.0267522, 0.0280809,
        0.0708397, 0.00714447, 0.0682617, 0.0686683, 0.134746
    };

    std::vector<double> Auu_cos2_reco = {
        -0.51314, -0.414621, -0.246236, -0.0920178, -0.0193428,
        0.0812403, 0.101458, -0.333557, -0.466727, -0.513419,
        -0.462745, -0.315971, -0.0162759, -0.167516, 0.135431,
        -0.17559, -0.427221, -0.59631, -0.380095, -0.134244,
        -0.394723, 0.101731, 0.190188, -0.387252, -0.394463
    };

    std::vector<double> Auu_cos2_err_reco = {
        0.0130561, 0.0136554, 0.0170336, 0.0236306, 0.032197,
        0.0402596, 0.0360069, 0.01217, 0.00982769, 0.0119169,
        0.0161806, 0.0225788, 0.0253259, 0.0327441, 0.0148564,
        0.0155902, 0.0180031, 0.0217514, 0.0273939, 0.0311547,
        0.0764327, 0.02386, 0.0620767, 0.0684223, 0.147174
    };

    double Auu_cos_fixed = 0.0;
    double Auu_cos2_fixed = 0.0;

    // ASYMMETRIES xQ2
    for (int x = 0; x < nbin_xQ2; x++){ 
        ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_preID[x], vec_kaonp_spin_2d_preID[x], vec_kaonp_y_2d_preID[x],vec_kaonp_eps_2d_preID[x], vec_kaonp_pol_2d_preID[x], vec_helicity_2d_preID[x], true, period);}, 7);
        minimizer.SetFunction(MLE);
        minimizer.SetMaxFunctionCalls(50000);
        minimizer.SetMaxIterations(10000);
        minimizer.SetTolerance(0.001);
        minimizer.SetStrategy(1);
        minimizer.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        minimizer.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer.Minimize();
        // LSA
        AUL_sin_2d_preID[x] = minimizer.X()[0];
        AUL_sin_err_2d_preID[x] = minimizer.Errors()[0];
        AUL_2sin_2d_preID[x] = minimizer.X()[1];
        AUL_2sin_err_2d_preID[x] = minimizer.Errors()[1];
        ALL_0_2d_preID[x] = minimizer.X()[2];
        ALL_0_err_2d_preID[x] = minimizer.Errors()[2];
        ALL_cos_2d_preID[x] = minimizer.X()[3];
        ALL_cos_err_2d_preID[x] = minimizer.Errors()[3];
        ALU_sin_2d_preID[x] = minimizer.X()[6];
        ALU_sin_err_2d_preID[x] = minimizer.Errors()[6];
        //cout << "bin: " << x <<  "  -->    Aul_sin = " << minimizer.X()[0]  << "    Aul_sin2 = " << minimizer.X()[1]  <<"    All = " << minimizer.X()[2] << "    All_cos = " << minimizer.X()[3]  << "    Auu_cos = " << minimizer.X()[4] << "    Auu_cos2 = " << minimizer.X()[5] <<  "    Alu_sin = " << minimizer.X()[6] << endl;
        
        // true
        ROOT::Minuit2::Minuit2Minimizer minimizer2(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE2([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_true[x], vec_kaonp_spin_2d_true[x], vec_kaonp_y_2d_true[x],vec_kaonp_eps_2d_true[x], vec_kaonp_pol_2d_true[x], vec_helicity_2d_true[x], true, period);}, 7);
        minimizer2.SetFunction(MLE2);
        minimizer2.SetMaxFunctionCalls(50000);
        minimizer2.SetMaxIterations(10000);
        minimizer2.SetTolerance(0.001);
        minimizer2.SetStrategy(1);
        minimizer2.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer2.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer2.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer2.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer2.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer2.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        minimizer2.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer2.Minimize();
        // LSA
        AUL_sin_2d_true[x] = minimizer2.X()[0];
        AUL_sin_err_2d_true[x] = minimizer2.Errors()[0];
        AUL_2sin_2d_true[x] = minimizer2.X()[1];
        AUL_2sin_err_2d_true[x] = minimizer2.Errors()[1];
        ALL_0_2d_true[x] = minimizer2.X()[2];
        ALL_0_err_2d_true[x] = minimizer2.Errors()[2];
        ALL_cos_2d_true[x] = minimizer2.X()[3];
        ALL_cos_err_2d_true[x] = minimizer2.Errors()[3];
        ALU_sin_2d_true[x] = minimizer2.X()[6];
        ALU_sin_err_2d_true[x] = minimizer2.Errors()[6];

        // allID
        ROOT::Minuit2::Minuit2Minimizer minimizer3(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE3([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_all_ID[x], vec_kaonp_spin_2d_all_ID[x], vec_kaonp_y_2d_all_ID[x],vec_kaonp_eps_2d_all_ID[x], vec_kaonp_pol_2d_all_ID[x], vec_helicity_2d_all_ID[x], true, period);}, 7);
        minimizer3.SetFunction(MLE3);
        minimizer3.SetMaxFunctionCalls(50000);
        minimizer3.SetMaxIterations(10000);
        minimizer3.SetTolerance(0.001);
        minimizer3.SetStrategy(1);
        minimizer3.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer3.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer3.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer3.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer3.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer3.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        minimizer3.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer3.Minimize();
        // LSA
        AUL_sin_2d_all_ID[x] = minimizer3.X()[0];
        AUL_sin_err_2d_all_ID[x] = minimizer3.Errors()[0];
        AUL_2sin_2d_all_ID[x] = minimizer3.X()[1];
        AUL_2sin_err_2d_all_ID[x] = minimizer3.Errors()[1];
        ALL_0_2d_all_ID[x] = minimizer3.X()[2];
        ALL_0_err_2d_all_ID[x] = minimizer3.Errors()[2];
        ALL_cos_2d_all_ID[x] = minimizer3.X()[3];
        ALL_cos_err_2d_all_ID[x] = minimizer3.Errors()[3];
        ALU_sin_2d_all_ID[x] = minimizer3.X()[6];
        ALU_sin_err_2d_all_ID[x] = minimizer3.Errors()[6];

        // preID
        ROOT::Minuit2::Minuit2Minimizer minimizer4(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE4([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_preID[x], vec_kaonp_spin_2d_preID[x], vec_kaonp_y_2d_preID[x],vec_kaonp_eps_2d_preID[x], vec_kaonp_pol_2d_preID[x], vec_helicity_2d_preID[x], true, period);}, 7);
        minimizer4.SetFunction(MLE4);
        minimizer4.SetMaxFunctionCalls(50000);
        minimizer4.SetMaxIterations(10000);
        minimizer4.SetTolerance(0.001);
        minimizer4.SetStrategy(1);
        minimizer4.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer4.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer4.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer4.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer4.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer4.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        minimizer4.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer4.Minimize();
        // LSA
        AUL_sin_2d_preID[x] = minimizer4.X()[0];
        AUL_sin_err_2d_preID[x] = minimizer4.Errors()[0];
        AUL_2sin_2d_preID[x] = minimizer4.X()[1];
        AUL_2sin_err_2d_preID[x] = minimizer4.Errors()[1];
        ALL_0_2d_preID[x] = minimizer4.X()[2];
        ALL_0_err_2d_preID[x] = minimizer4.Errors()[2];
        ALL_cos_2d_preID[x] = minimizer4.X()[3];
        ALL_cos_err_2d_preID[x] = minimizer4.Errors()[3];
        ALU_sin_2d_preID[x] = minimizer4.X()[6];
        ALU_sin_err_2d_preID[x] = minimizer4.Errors()[6];

        // reco
        ROOT::Minuit2::Minuit2Minimizer minimizer5(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE5([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d[x], vec_kaonp_spin_2d[x], vec_kaonp_y_2d[x],vec_kaonp_eps_2d[x], vec_kaonp_pol_2d[x], vec_helicity_2d[x], true, period);}, 7);
        minimizer5.SetFunction(MLE5);
        minimizer5.SetMaxFunctionCalls(50000);
        minimizer5.SetMaxIterations(10000);
        minimizer5.SetTolerance(0.001);
        minimizer5.SetStrategy(1);
        minimizer5.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer5.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer5.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer5.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer5.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer5.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        minimizer5.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer5.Minimize();
        // LSA
        AUL_sin_2d_reco[x] = minimizer5.X()[0];
        AUL_sin_err_2d_reco[x] = minimizer5.Errors()[0];
        AUL_2sin_2d_reco[x] = minimizer5.X()[1];
        AUL_2sin_err_2d_reco[x] = minimizer5.Errors()[1];
        ALL_0_2d_reco[x] = minimizer5.X()[2];
        ALL_0_err_2d_reco[x] = minimizer5.Errors()[2];
        ALL_cos_2d_reco[x] = minimizer5.X()[3];
        ALL_cos_err_2d_reco[x] = minimizer5.Errors()[3];
        ALU_sin_2d_reco[x] = minimizer5.X()[6];
        ALU_sin_err_2d_reco[x] = minimizer5.Errors()[6];

        // mc
        ROOT::Minuit2::Minuit2Minimizer minimizer6(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE6([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_mc[x], vec_kaonp_spin_2d_mc[x], vec_kaonp_y_2d_mc[x],vec_kaonp_eps_2d_mc[x], vec_kaonp_pol_2d_mc[x], vec_helicity_2d_mc[x], true, period);}, 7);
        minimizer6.SetFunction(MLE6);
        minimizer6.SetMaxFunctionCalls(50000);
        minimizer6.SetMaxIterations(10000);
        minimizer6.SetTolerance(0.001);
        minimizer6.SetStrategy(1);
        minimizer6.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer6.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer6.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer6.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        //minimizer6.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        //minimizer6.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        minimizer6.SetLimitedVariable(4, "Auu", 0.00, 0.000, -0.8, 0.8);
        minimizer6.SetLimitedVariable(5, "Auu2", 0.00, 0.000, -0.8, 0.8);
        minimizer6.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer6.Minimize();
        // LSA
        AUL_sin_2d_mc[x] = minimizer6.X()[0];
        AUL_sin_err_2d_mc[x] = minimizer6.Errors()[0];
        AUL_2sin_2d_mc[x] = minimizer6.X()[1];
        AUL_2sin_err_2d_mc[x] = minimizer6.Errors()[1];
        ALL_0_2d_mc[x] = minimizer6.X()[2];
        ALL_0_err_2d_mc[x] = minimizer6.Errors()[2];
        ALL_cos_2d_mc[x] = minimizer6.X()[3];
        ALL_cos_err_2d_mc[x] = minimizer6.Errors()[3];
        ALU_sin_2d_mc[x] = minimizer6.X()[6];
        ALU_sin_err_2d_mc[x] = minimizer6.Errors()[6];
        cout << "bin: " << x <<  "  -->  Aul_sin = " << minimizer6.X()[0] << "    Aul_sin2 = " << minimizer6.X()[1] << "    All = " << minimizer6.X()[2] << "    All_cos = " << minimizer6.X()[3]  << "    Auu_cos = " << minimizer6.X()[4] << "    Auu_cos2 = " << minimizer6.X()[5] <<  "    Alu_sin = " << minimizer6.X()[6] << endl;
    }
    cout << " -------------------------------------------------------------------------------------------------------------------------------------------- " << endl;
    


    // Asymmetries zPt
    for (int z = 0; z < nbin_zPt; z++){ 
        ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_preID[z], vec_kaonp_spin_2d_zPt_preID[z], vec_kaonp_y_2d_zPt_preID[z],vec_kaonp_eps_2d_zPt_preID[z], vec_kaonp_pol_2d_zPt_preID[z], vec_helicity_2d_zPt_preID[z], true, period);}, 7);
        minimizer.SetFunction(MLE);
        minimizer.SetMaxFunctionCalls(50000);
        minimizer.SetMaxIterations(10000);
        minimizer.SetTolerance(0.001);
        minimizer.SetStrategy(1);
        minimizer.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        //minimizer.SetLimitedVariable(4, "Auu", 0.00, 0.000, -0.8, 0.8);
        //minimizer.SetLimitedVariable(5, "Auu2", 0.00, 0.000, -0.8, 0.8);
        minimizer.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer.Minimize();
        // LSA
        AUL_sin_2d_zPt_preID[z] = minimizer.X()[0];
        AUL_sin_err_2d_zPt_preID[z] = minimizer.Errors()[0];
        AUL_2sin_2d_zPt_preID[z] = minimizer.X()[1];
        AUL_2sin_err_2d_zPt_preID[z] = minimizer.Errors()[1];
        ALL_0_2d_zPt_preID[z] = minimizer.X()[2];
        ALL_0_err_2d_zPt_preID[z] = minimizer.Errors()[2];
        ALL_cos_2d_zPt_preID[z] = minimizer.X()[3];
        ALL_cos_err_2d_zPt_preID[z] = minimizer.Errors()[3];
        ALU_sin_2d_zPt_preID[z] = minimizer.X()[6];
        ALU_sin_err_2d_zPt_preID[z] = minimizer.Errors()[6];
        //cout << "bin: " << z <<  "  -->  Aul_sin = " << minimizer.X()[0] << "    Aul_sin2 = " << minimizer.X()[1] << "    All = " << minimizer.X()[2] << "    All_cos = " << minimizer.X()[3]  << "    Auu_cos = " << minimizer.X()[4] << "    Auu_cos2 = " << minimizer.X()[5] <<  "    Alu_sin = " << minimizer.X()[6] << endl;

        // true
        ROOT::Minuit2::Minuit2Minimizer minimizer2(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE2([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_true[z], vec_kaonp_spin_2d_zPt_true[z], vec_kaonp_y_2d_zPt_true[z],vec_kaonp_eps_2d_zPt_true[z], vec_kaonp_pol_2d_zPt_true[z], vec_helicity_2d_zPt_true[z], true, period);}, 7);
        minimizer2.SetFunction(MLE2);
        minimizer2.SetMaxFunctionCalls(50000);
        minimizer2.SetMaxIterations(10000);
        minimizer2.SetTolerance(0.001);
        minimizer2.SetStrategy(1);
        minimizer2.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer2.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer2.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer2.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer2.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer2.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        //minimizer2.SetLimitedVariable(4, "Auu", 0.00, 0.000, -0.8, 0.8);
        //minimizer2.SetLimitedVariable(5, "Auu2", 0.00, 0.000, -0.8, 0.8);
        minimizer2.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer2.Minimize();
        // LSA
        AUL_sin_2d_zPt_true[z] = minimizer2.X()[0];
        AUL_sin_err_2d_zPt_true[z] = minimizer2.Errors()[0];
        AUL_2sin_2d_zPt_true[z] = minimizer2.X()[1];
        AUL_2sin_err_2d_zPt_true[z] = minimizer2.Errors()[1];
        ALL_0_2d_zPt_true[z] = minimizer2.X()[2];
        ALL_0_err_2d_zPt_true[z] = minimizer2.Errors()[2];
        ALL_cos_2d_zPt_true[z] = minimizer2.X()[3];
        ALL_cos_err_2d_zPt_true[z] = minimizer2.Errors()[3];
        ALU_sin_2d_zPt_true[z] = minimizer2.X()[6];
        ALU_sin_err_2d_zPt_true[z] = minimizer2.Errors()[6];

        // allID
        ROOT::Minuit2::Minuit2Minimizer minimizer3(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE3([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_all_ID[z], vec_kaonp_spin_2d_zPt_all_ID[z], vec_kaonp_y_2d_zPt_all_ID[z],vec_kaonp_eps_2d_zPt_all_ID[z], vec_kaonp_pol_2d_zPt_all_ID[z], vec_helicity_2d_zPt_all_ID[z], true, period);}, 7);
        minimizer3.SetFunction(MLE3);
        minimizer3.SetMaxFunctionCalls(50000);
        minimizer3.SetMaxIterations(10000);
        minimizer3.SetTolerance(0.001);
        minimizer3.SetStrategy(1);
        minimizer3.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer3.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer3.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer3.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer3.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer3.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        //minimizer3.SetLimitedVariable(4, "Auu", 0.00, 0.000, -0.8, 0.8);
        //minimizer3.SetLimitedVariable(5, "Auu2", 0.00, 0.000, -0.8, 0.8);
        minimizer3.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer3.Minimize();
        // LSA
        AUL_sin_2d_zPt_all_ID[z] = minimizer3.X()[0];
        AUL_sin_err_2d_zPt_all_ID[z] = minimizer3.Errors()[0];
        AUL_2sin_2d_zPt_all_ID[z] = minimizer3.X()[1];
        AUL_2sin_err_2d_zPt_all_ID[z] = minimizer3.Errors()[1];
        ALL_0_2d_zPt_all_ID[z] = minimizer3.X()[2];
        ALL_0_err_2d_zPt_all_ID[z] = minimizer3.Errors()[2];
        ALL_cos_2d_zPt_all_ID[z] = minimizer3.X()[3];
        ALL_cos_err_2d_zPt_all_ID[z] = minimizer3.Errors()[3];
        ALU_sin_2d_zPt_all_ID[z] = minimizer3.X()[6];
        ALU_sin_err_2d_zPt_all_ID[z] = minimizer3.Errors()[6];

        // pre ID
        ROOT::Minuit2::Minuit2Minimizer minimizer4(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE4([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_preID[z], vec_kaonp_spin_2d_zPt_preID[z], vec_kaonp_y_2d_zPt_preID[z],vec_kaonp_eps_2d_zPt_preID[z], vec_kaonp_pol_2d_zPt_preID[z], vec_helicity_2d_zPt_preID[z], true, period);}, 7);
        minimizer4.SetFunction(MLE4);
        minimizer4.SetMaxFunctionCalls(50000);
        minimizer4.SetMaxIterations(10000);
        minimizer4.SetTolerance(0.001);
        minimizer4.SetStrategy(1);
        minimizer4.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer4.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer4.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer4.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer4.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer4.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        //minimizer4.SetLimitedVariable(4, "Auu", 0.00, 0.000, -0.8, 0.8);
        //minimizer4.SetLimitedVariable(5, "Auu2", 0.00, 0.000, -0.8, 0.8);
        minimizer4.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer4.Minimize();
        // LSA
        AUL_sin_2d_zPt_preID[z] = minimizer4.X()[0];
        AUL_sin_err_2d_zPt_preID[z] = minimizer4.Errors()[0];
        AUL_2sin_2d_zPt_preID[z] = minimizer4.X()[1];
        AUL_2sin_err_2d_zPt_preID[z] = minimizer4.Errors()[1];
        ALL_0_2d_zPt_preID[z] = minimizer4.X()[2];
        ALL_0_err_2d_zPt_preID[z] = minimizer4.Errors()[2];
        ALL_cos_2d_zPt_preID[z] = minimizer4.X()[3];
        ALL_cos_err_2d_zPt_preID[z] = minimizer4.Errors()[3];
        ALU_sin_2d_zPt_preID[z] = minimizer4.X()[6];
        ALU_sin_err_2d_zPt_preID[z] = minimizer4.Errors()[6];

        // reco
        ROOT::Minuit2::Minuit2Minimizer minimizer5(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE5([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt[z], vec_kaonp_spin_2d_zPt[z], vec_kaonp_y_2d_zPt[z],vec_kaonp_eps_2d_zPt[z], vec_kaonp_pol_2d_zPt[z], vec_helicity_2d_zPt[z], true, period);}, 7);
        minimizer5.SetFunction(MLE5);
        minimizer5.SetMaxFunctionCalls(50000);
        minimizer5.SetMaxIterations(10000);
        minimizer5.SetTolerance(0.001);
        minimizer5.SetStrategy(1);
        minimizer5.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer5.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer5.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer5.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        minimizer5.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        minimizer5.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        //minimizer5.SetLimitedVariable(4, "Auu", 0.00, 0.000, -0.8, 0.8);
        //minimizer5.SetLimitedVariable(5, "Auu2", 0.00, 0.000, -0.8, 0.8);
        minimizer5.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer5.Minimize();
        // LSA
        AUL_sin_2d_zPt_reco[z] = minimizer5.X()[0];
        AUL_sin_err_2d_zPt_reco[z] = minimizer5.Errors()[0];
        AUL_2sin_2d_zPt_reco[z] = minimizer5.X()[1];
        AUL_2sin_err_2d_zPt_reco[z] = minimizer5.Errors()[1];
        ALL_0_2d_zPt_reco[z] = minimizer5.X()[2];
        ALL_0_err_2d_zPt_reco[z] = minimizer5.Errors()[2];
        ALL_cos_2d_zPt_reco[z] = minimizer5.X()[3];
        ALL_cos_err_2d_zPt_reco[z] = minimizer5.Errors()[3];
        ALU_sin_2d_zPt_reco[z] = minimizer5.X()[6];
        ALU_sin_err_2d_zPt_reco[z] = minimizer5.Errors()[6];

        ROOT::Minuit2::Minuit2Minimizer minimizer6(ROOT::Minuit2::kMigrad);
        ROOT::Math::Functor MLE6([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_mc[z], vec_kaonp_spin_2d_zPt_mc[z], vec_kaonp_y_2d_zPt_mc[z],vec_kaonp_eps_2d_zPt_mc[z], vec_kaonp_pol_2d_zPt_mc[z], vec_helicity_2d_zPt_mc[z], true, period);}, 7);
        minimizer6.SetFunction(MLE6);
        minimizer6.SetMaxFunctionCalls(50000);
        minimizer6.SetMaxIterations(10000);
        minimizer6.SetTolerance(0.001);
        minimizer6.SetStrategy(1);
        minimizer6.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
        minimizer6.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
        minimizer6.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
        minimizer6.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
        //minimizer6.SetFixedVariable(4, "Auu", Auu_cos_fixed);
        //minimizer6.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
        minimizer6.SetLimitedVariable(4, "Auu", 0.00, 0.000, -0.8, 0.8);
        minimizer6.SetLimitedVariable(5, "Auu2", 0.00, 0.000, -0.8, 0.8);
        minimizer6.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
        minimizer6.Minimize();
        // LSA
        AUL_sin_2d_zPt_mc[z] = minimizer6.X()[0];
        AUL_sin_err_2d_zPt_mc[z] = minimizer6.Errors()[0];
        AUL_2sin_2d_zPt_mc[z] = minimizer6.X()[1];
        AUL_2sin_err_2d_zPt_mc[z] = minimizer6.Errors()[1];
        ALL_0_2d_zPt_mc[z] = minimizer6.X()[2];
        ALL_0_err_2d_zPt_mc[z] = minimizer6.Errors()[2];
        ALL_cos_2d_zPt_mc[z] = minimizer6.X()[3];
        ALL_cos_err_2d_zPt_mc[z] = minimizer6.Errors()[3];
        ALU_sin_2d_zPt_mc[z] = minimizer6.X()[6];
        ALU_sin_err_2d_zPt_mc[z] = minimizer6.Errors()[6];
        cout << "bin: " << z <<  "  -->  Aul_sin = " << minimizer6.X()[0] << "    Aul_sin2 = " << minimizer6.X()[1] << "    All = " << minimizer6.X()[2] << "    All_cos = " << minimizer6.X()[3]  << "    Auu_cos = " << minimizer6.X()[4] << "    Auu_cos2 = " << minimizer6.X()[5] <<  "    Alu_sin = " << minimizer6.X()[6] << endl;
        
    }
    cout << " -------------------------------------------------------------------------------------------------------------------------------------------- " << endl;
    

    //

    vector<double> x_points_AUL;
    vector<double> y_points_AUL;

    for (int ix = 0; ix < nbin_xQ2; ix++) {
        for (int iz = 0; iz < nbin_zPt; iz++) {

            const auto& z_vec  = vec_z_4d[ix][iz];
            const auto& pt_vec = vec_Pt_4d[ix][iz];
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

            x_points_AUL.push_back(xB_mean);
            y_points_AUL.push_back(ratio);
        }
    }

    TGraph* gr_AUL = new TGraph(x_points_AUL.size(),x_points_AUL.data(),y_points_AUL.data());
    gr_AUL->SetTitle("TMD factorization parameter | K+; <x_{B}>; <P_{hT}>/(<z><Q>)");
    gr_AUL->SetMarkerStyle(20);
    gr_AUL->SetMarkerSize(1);
    gr_AUL->SetMarkerColor(kOrange-3);
    //gr->SetLineWidth(2);
    TCanvas* c_ratio_AUL = new TCanvas("c_ratio_AUL","P_{hT}/(zQ)",800,600);
    gr->Draw("AP");
    gr_AUL->Draw("P");
    TLegend* leg0 = new TLegend(0.70, 0.75, 0.88, 0.88); // Adjust position (x1,y1,x2,y2)
    leg0->AddEntry(gr, "ALL bin coverage", "p");
    leg0->AddEntry(gr_AUL, "A_{UL}^{sin(2#Phi_{h})} err < 5% bin", "p");
    leg0->SetFillStyle(0);  // Transparent background
    leg0->Draw();
    c_ratio_AUL->Write();





    //
    dir_aul_sin->cd();
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_AUL_vs_xB = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_xB_4 = new TGraphErrors();
    int p_idx1 = 0, p_idx1_1 = 0, p_idx1_2 = 0, p_idx1_3 = 0, p_idx1_4 = 0;
    for(int x = 0; x < nbin_xQ2; x++){
        if (vec_kaonp_xB_2d_mc[x].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_mc[x].empty()) continue;
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double val : vec_kaonp_xB_2d_mc[x]){
            sum_xB += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_xB = sum_xB / vec_kaonp_xB_2d_mc[x].size();
        double aul_sin = AUL_sin_2d_mc[x]; 
        double aul_sin_err = AUL_sin_err_2d_mc[x];
        graph2D_AUL_vs_xB->SetPoint(p_idx1, mean_xB, aul_sin);
        graph2D_AUL_vs_xB->SetPointError(p_idx1, 0.0, aul_sin_err); // No x error
        p_idx1++;
        if(x < 5){
            graph2D_AUL_vs_xB_1->SetPoint(p_idx1_1, mean_xB, aul_sin);
            graph2D_AUL_vs_xB_1->SetPointError(p_idx1_1, 0.0, aul_sin_err); 
            p_idx1_1++;
        } else if (x < 10){
            graph2D_AUL_vs_xB_2->SetPoint(p_idx1_2, mean_xB, aul_sin);
            graph2D_AUL_vs_xB_2->SetPointError(p_idx1_2, 0.0, aul_sin_err);
            p_idx1_2++;
        } else if (x < 12){
            graph2D_AUL_vs_xB_3->SetPoint(p_idx1_3, mean_xB, aul_sin);
            graph2D_AUL_vs_xB_3->SetPointError(p_idx1_3, 0.0, aul_sin_err);
            p_idx1_3++;
        } else if (x < 18){
            graph2D_AUL_vs_xB_4->SetPoint(p_idx1_4, mean_xB, aul_sin);
            graph2D_AUL_vs_xB_4->SetPointError(p_idx1_4, 0.0, aul_sin_err);
            p_idx1_4++;
        }
    }
    graph2D_AUL_vs_xB_1->SetMarkerStyle(20), graph2D_AUL_vs_xB_2->SetMarkerStyle(20), graph2D_AUL_vs_xB_3->SetMarkerStyle(20), graph2D_AUL_vs_xB_4->SetMarkerStyle(20);
    graph2D_AUL_vs_xB_1->SetLineColor(kAzure-5), graph2D_AUL_vs_xB_1->SetMarkerColor(kAzure-5);
    graph2D_AUL_vs_xB_2->SetLineColor(kViolet-5), graph2D_AUL_vs_xB_2->SetMarkerColor(kViolet-5);
    graph2D_AUL_vs_xB_3->SetLineColor(kPink-5), graph2D_AUL_vs_xB_3->SetMarkerColor(kPink-5);
    graph2D_AUL_vs_xB_4->SetLineColor(kOrange-5), graph2D_AUL_vs_xB_4->SetMarkerColor(kOrange-5);


    TCanvas* c_Aut2D_xB_2 = new TCanvas("Aul_sin_vs_xB_2d", "sin(#Phi_{h}) longitudinal asymmetry vs x_{B}", 800, 600);
    //graph2D_AUL_vs_xB->Draw("A");
    TMultiGraph *mg_Aut2D_xB_2 = new TMultiGraph();
    mg_Aut2D_xB_2->SetTitle("A_{UL}^{sin#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin) | lepton frame; x_{B}; F_{UL}^{sin(#Phi_{h})}/F_{UU}");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_1, "P");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_2, "P");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_3, "P");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_4, "P");
    mg_Aut2D_xB_2->Draw("A");
    mg_Aut2D_xB_2->GetYaxis()->SetRangeUser(-0.3, 0.3);
    TLine* guideLine11 = new TLine(mg_Aut2D_xB_2->GetXaxis()->GetXmin(), 0, mg_Aut2D_xB_2->GetXaxis()->GetXmax(), 0);
    guideLine11->SetLineStyle(2);  
    guideLine11->SetLineColor(kGray+1);
    guideLine11->Draw();
    TLine* guideLine11_val = new TLine(mg_Aut2D_xB_2->GetXaxis()->GetXmin(), A_UL_sin_inject, mg_Aut2D_xB_2->GetXaxis()->GetXmax(), A_UL_sin_inject);
    guideLine11_val->SetLineStyle(2);  
    guideLine11_val->SetLineColor(kRed+1);
    guideLine11_val->Draw();

    TLegend* legend = new TLegend(0.13, 0.7, 0.33, 0.88); 
    legend->AddEntry(graph2D_AUL_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend->AddEntry(graph2D_AUL_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend->AddEntry(graph2D_AUL_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend->AddEntry(graph2D_AUL_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
    legend->SetFillStyle(0);  // Transparent background
    legend->Draw();
    c_Aut2D_xB_2->Update();
    c_Aut2D_xB_2->Write();

    
    // zPt
    TGraphErrors* graph2D_AUL_vs_z = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph2D_AUL_vs_z_4 = new TGraphErrors();
    int pz_idx1 = 0, pz_idx1_1 = 0, pz_idx1_2 = 0, pz_idx1_3 = 0, pz_idx1_4 = 0;
    for(int z = 0; z < nbin_zPt; z++){
        if (vec_kaonp_xB_2d_zPt_mc[z].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_zPt_mc[z].empty()) continue;
        // Compute mean xB for this bin
        double sum_z = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_zPt_mc[z]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_zPt_mc[z].size();
        for (double val : vec_kaonp_z_2d_zPt_mc[z]){
            sum_z += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_z = sum_z / vec_kaonp_z_2d_zPt_mc[z].size();
        // Get AUL and error from vectors
        double aul_sin = AUL_sin_2d_zPt_mc[z]; 
        double aul_sin_err = AUL_sin_err_2d_zPt_mc[z];
        graph2D_AUL_vs_z->SetPoint(pz_idx1, mean_z, aul_sin);
        graph2D_AUL_vs_z->SetPointError(pz_idx1, 0.0, aul_sin_err); 
        pz_idx1++;
        if(z < 7){
            graph2D_AUL_vs_z_1->SetPoint(pz_idx1_1, mean_z, aul_sin);
            graph2D_AUL_vs_z_1->SetPointError(pz_idx1_1, 0.0, aul_sin_err); 
            pz_idx1_1++;
        } else if (z < 14){
            graph2D_AUL_vs_z_2->SetPoint(pz_idx1_2, mean_z, aul_sin);
            graph2D_AUL_vs_z_2->SetPointError(pz_idx1_2, 0.0, aul_sin_err);
            pz_idx1_2++;
        } else if (z < 21){
            graph2D_AUL_vs_z_3->SetPoint(pz_idx1_3, mean_z, aul_sin);
            graph2D_AUL_vs_z_3->SetPointError(pz_idx1_3, 0.0, aul_sin_err);
            pz_idx1_3++;
        } else if (z < 26){
            graph2D_AUL_vs_z_4->SetPoint(pz_idx1_4, mean_z, aul_sin);
            graph2D_AUL_vs_z_4->SetPointError(pz_idx1_4, 0.0, aul_sin_err);
            pz_idx1_4++;
        }
    }
    graph2D_AUL_vs_z_1->SetMarkerStyle(20), graph2D_AUL_vs_z_2->SetMarkerStyle(20), graph2D_AUL_vs_z_3->SetMarkerStyle(20), graph2D_AUL_vs_z_4->SetMarkerStyle(20);
    graph2D_AUL_vs_z_1->SetLineColor(kAzure-5), graph2D_AUL_vs_z_1->SetMarkerColor(kAzure-5);
    graph2D_AUL_vs_z_2->SetLineColor(kViolet-5), graph2D_AUL_vs_z_2->SetMarkerColor(kViolet-5);
    graph2D_AUL_vs_z_3->SetLineColor(kPink-5), graph2D_AUL_vs_z_3->SetMarkerColor(kPink-5);
    graph2D_AUL_vs_z_4->SetLineColor(kOrange-5), graph2D_AUL_vs_z_4->SetMarkerColor(kOrange-5);



    vector<TGraphErrors*> graphs_AUL_vs_z = {graph2D_AUL_vs_z_1,graph2D_AUL_vs_z_2,graph2D_AUL_vs_z_3,graph2D_AUL_vs_z_4};
    vector<string> titles_AUL_vs_z = {"0.0 < P_{hT} < 0.25 GeV","0.25 < P_{hT} < 0.5 GeV","0.5 < P_{hT} < 0.8 GeV","0.8 < P_{hT} < 1.4 GeV"};

    for (size_t i = 0; i < graphs_AUL_vs_z.size(); ++i) {
        TString cname = Form("c_AUL_sin_vs_z_PtBin_%zu", i);
        TString ctitle = Form("A_{UL}^{sin#Phi_{h}} vs z (%s)", titles_AUL_vs_z[i].c_str());

        TCanvas* c = new TCanvas(cname, ctitle, 800, 600);
        graphs_AUL_vs_z[i]->SetTitle(Form("A_{UL}^{sin#Phi_{h}} vs z | %s | lepton frame; z; F_{UL}^{sin(#Phi_{h})}/F_{UU}",titles_AUL_vs_z[i].c_str()));
        graphs_AUL_vs_z[i]->Draw("AP");
        graphs_AUL_vs_z[i]->GetYaxis()->SetRangeUser(-0.3, 0.3);


        TLine* zeroLine = new TLine(graphs_AUL_vs_z[i]->GetXaxis()->GetXmin(), 0,graphs_AUL_vs_z[i]->GetXaxis()->GetXmax(), 0);
        zeroLine->SetLineStyle(2);
        zeroLine->SetLineColor(kGray+1);
        zeroLine->Draw();
        TLine* zeroLine_v = new TLine(graphs_AUL_vs_z[i]->GetXaxis()->GetXmin(), A_UL_sin_inject,graphs_AUL_vs_z[i]->GetXaxis()->GetXmax(), A_UL_sin_inject);
        zeroLine_v->SetLineStyle(2);
        zeroLine_v->SetLineColor(kRed+1);
        zeroLine_v->Draw();

        c->Write();
    }

    


    

    // AUL sin2 ------------------------------------------------------------------------------------


    dir_aul_sin2->cd();
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_AUL2_vs_xB = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_4 = new TGraphErrors();
    int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0;
    for(int x = 0; x < nbin_xQ2; x++){
        if (vec_kaonp_xB_2d_mc[x].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_mc[x].empty()) continue;
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double val : vec_kaonp_xB_2d_mc[x]){
            sum_xB += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_xB = sum_xB / vec_kaonp_xB_2d_mc[x].size();
        double aul_sin = AUL_2sin_2d_mc[x];
        double aul_sin_err = AUL_2sin_err_2d_mc[x];
        if(aul_sin > 1 || aul_sin < -1) continue;
        graph2D_AUL2_vs_xB->SetPoint(p_idx, mean_xB, aul_sin);
        graph2D_AUL2_vs_xB->SetPointError(p_idx, 0.0, aul_sin_err); // No x error
        p_idx++;
        if(x < 5){
            graph2D_AUL2_vs_xB_1->SetPoint(p_idx_1, mean_xB, aul_sin);
            graph2D_AUL2_vs_xB_1->SetPointError(p_idx_1, 0.0, aul_sin_err); // No x error
            p_idx_1++;
        } else if (x < 10){
            graph2D_AUL2_vs_xB_2->SetPoint(p_idx_2, mean_xB, aul_sin);
            graph2D_AUL2_vs_xB_2->SetPointError(p_idx_2, 0.0, aul_sin_err);
            p_idx_2++;
        } else if (x < 12){
            graph2D_AUL2_vs_xB_3->SetPoint(p_idx_3, mean_xB, aul_sin);
            graph2D_AUL2_vs_xB_3->SetPointError(p_idx_3, 0.0, aul_sin_err);
            p_idx_3++;
        } else if (x < 18){
            graph2D_AUL2_vs_xB_4->SetPoint(p_idx_4, mean_xB, aul_sin);
            graph2D_AUL2_vs_xB_4->SetPointError(p_idx_4, 0.0, aul_sin_err);
            p_idx_4++;
        }
    }

    graph2D_AUL2_vs_xB->SetTitle("A_{UL}^{sin2#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin); x_{B}; F_{UL}^{sin(2#Phi_{h})}/F_{UU}"); // #sqrt{2#epsilon(1+#epsilon)}
    graph2D_AUL2_vs_xB_1->SetMarkerStyle(20), graph2D_AUL2_vs_xB_2->SetMarkerStyle(20), graph2D_AUL2_vs_xB_3->SetMarkerStyle(20), graph2D_AUL2_vs_xB_4->SetMarkerStyle(20);
    graph2D_AUL2_vs_xB_1->SetLineColor(kAzure-5), graph2D_AUL2_vs_xB_1->SetMarkerColor(kAzure-5);
    graph2D_AUL2_vs_xB_2->SetLineColor(kViolet-5), graph2D_AUL2_vs_xB_2->SetMarkerColor(kViolet-5);
    graph2D_AUL2_vs_xB_3->SetLineColor(kPink-5), graph2D_AUL2_vs_xB_3->SetMarkerColor(kPink-5);
    graph2D_AUL2_vs_xB_4->SetLineColor(kOrange-5), graph2D_AUL2_vs_xB_4->SetMarkerColor(kOrange-5);

    
    TCanvas* c_Aut2D2_xB_2 = new TCanvas("Aul_sin2_vs_xB_2d", "sin(2#Phi_{h}) longitudinal asymmetry vs x_{B}", 800, 600);
    TMultiGraph *mg_Aut2D2_xB_2 = new TMultiGraph();
    mg_Aut2D2_xB_2->SetTitle("A_{UL}^{sin2#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin) | lepton frame; x_{B}; F_{UL}^{sin(2#Phi_{h})}/F_{UU}"); // #sqrt{2#epsilon(1+#epsilon)}
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_1, "P");
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_2, "P");
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_3, "P");
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_4, "P");
    //mg_Aut2D2_xB_2->Add(graph2D_corr_2sin, "P");
    mg_Aut2D2_xB_2->Draw("A");
    mg_Aut2D2_xB_2->GetYaxis()->SetRangeUser(-0.3, 0.3);
    TLine* guideLine111 = new TLine(mg_Aut2D2_xB_2->GetXaxis()->GetXmin(), 0, mg_Aut2D2_xB_2->GetXaxis()->GetXmax(), 0);
    guideLine111->SetLineStyle(2);  
    guideLine111->SetLineColor(kGray+1);
    guideLine111->Draw();
    TLine* guideLine111_val = new TLine(mg_Aut2D2_xB_2->GetXaxis()->GetXmin(), A_UL_2sin_inject, mg_Aut2D2_xB_2->GetXaxis()->GetXmax(), A_UL_2sin_inject);
    guideLine111_val->SetLineStyle(2);  
    guideLine111_val->SetLineColor(kRed+1);
    guideLine111_val->Draw();
    TLegend* legend2 = new TLegend(0.13, 0.7, 0.33, 0.88); // Adjust position (x1,y1,x2,y2)
    legend2->AddEntry(graph2D_AUL2_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend2->AddEntry(graph2D_AUL2_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend2->AddEntry(graph2D_AUL2_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend2->AddEntry(graph2D_AUL2_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
    legend2->SetFillStyle(0);  // Transparent background
    legend2->Draw();
    c_Aut2D2_xB_2->Update();
    c_Aut2D2_xB_2->Write();


    // zPt
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_AUL2_vs_z = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_z_4 = new TGraphErrors();
    int pz_idx = 0, pz_idx_1 = 0, pz_idx_2 = 0, pz_idx_3 = 0, pz_idx_4 = 0;
    for(int z = 0; z < nbin_zPt; z++){
        if (vec_kaonp_xB_2d_zPt_mc[z].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_zPt_mc[z].empty()) continue;
        // Compute mean xB for this bin
        double sum_z = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_zPt_mc[z]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_zPt_mc[z].size();
        for (double val : vec_kaonp_z_2d_zPt_mc[z]){
            sum_z += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_z = sum_z / vec_kaonp_z_2d_zPt_mc[z].size();
        // Get AUL2 and error from vectors
        double aul_sin = AUL_2sin_2d_zPt_mc[z];
        double aul_sin_err = AUL_2sin_err_2d_zPt_mc[z];
        graph2D_AUL2_vs_z->SetPoint(pz_idx, mean_z, aul_sin);
        graph2D_AUL2_vs_z->SetPointError(pz_idx, 0.0, aul_sin_err); 
        pz_idx++;
        if(z < 7){
            graph2D_AUL2_vs_z_1->SetPoint(pz_idx_1, mean_z, aul_sin);
            graph2D_AUL2_vs_z_1->SetPointError(pz_idx_1, 0.0, aul_sin_err); 
            pz_idx_1++;
        } else if (z < 14){
            graph2D_AUL2_vs_z_2->SetPoint(pz_idx_2, mean_z, aul_sin);
            graph2D_AUL2_vs_z_2->SetPointError(pz_idx_2, 0.0, aul_sin_err);
            pz_idx_2++;
        } else if (z < 21){
            graph2D_AUL2_vs_z_3->SetPoint(pz_idx_3, mean_z, aul_sin);
            graph2D_AUL2_vs_z_3->SetPointError(pz_idx_3, 0.0, aul_sin_err);
            pz_idx_3++;
        } else if (z < 26){
            graph2D_AUL2_vs_z_4->SetPoint(pz_idx_4, mean_z, aul_sin);
            graph2D_AUL2_vs_z_4->SetPointError(pz_idx_4, 0.0, aul_sin_err);
            pz_idx_4++;
        }
    }
    graph2D_AUL2_vs_z_1->SetMarkerStyle(20), graph2D_AUL2_vs_z_2->SetMarkerStyle(20), graph2D_AUL2_vs_z_3->SetMarkerStyle(20), graph2D_AUL2_vs_z_4->SetMarkerStyle(20);
    graph2D_AUL2_vs_z_1->SetLineColor(kAzure-5), graph2D_AUL2_vs_z_1->SetMarkerColor(kAzure-5);
    graph2D_AUL2_vs_z_2->SetLineColor(kViolet-5), graph2D_AUL2_vs_z_2->SetMarkerColor(kViolet-5);
    graph2D_AUL2_vs_z_3->SetLineColor(kPink-5), graph2D_AUL2_vs_z_3->SetMarkerColor(kPink-5);
    graph2D_AUL2_vs_z_4->SetLineColor(kOrange-5), graph2D_AUL2_vs_z_4->SetMarkerColor(kOrange-5);


    vector<TGraphErrors*> graphs_AUL2_vs_z = {graph2D_AUL2_vs_z_1,graph2D_AUL2_vs_z_2,graph2D_AUL2_vs_z_3,graph2D_AUL2_vs_z_4};
    vector<string> titles_AUL2_vs_z = {"0.0 < P_{hT} < 0.25 GeV","0.25 < P_{hT} < 0.5 GeV","0.5 < P_{hT} < 0.8 GeV","0.8 < P_{hT} < 1.4 GeV"};

    for (size_t i = 0; i < graphs_AUL2_vs_z.size(); ++i) {
        TString cname = Form("c_AUL_sin2_vs_z_PtBin_%zu", i);
        TString ctitle = Form("A_{UL}^{sin2#Phi_{h}} vs z (%s)", titles_AUL2_vs_z[i].c_str());

        TCanvas* c = new TCanvas(cname, ctitle, 800, 600);
        graphs_AUL2_vs_z[i]->SetTitle(Form("A_{UL}^{sin2#Phi_{h}} vs z | %s | lepton frame; z; F_{UL}^{sin(2#Phi_{h})}/F_{UU}",titles_AUL2_vs_z[i].c_str()));
        graphs_AUL2_vs_z[i]->Draw("AP");
        graphs_AUL2_vs_z[i]->GetYaxis()->SetRangeUser(-0.3, 0.3);

        TLine* zeroLine = new TLine(graphs_AUL2_vs_z[i]->GetXaxis()->GetXmin(), 0,graphs_AUL2_vs_z[i]->GetXaxis()->GetXmax(), 0);
        zeroLine->SetLineStyle(2);
        zeroLine->SetLineColor(kGray+1);
        zeroLine->Draw();
        TLine* zeroLine_v = new TLine(graphs_AUL2_vs_z[i]->GetXaxis()->GetXmin(), A_UL_2sin_inject,graphs_AUL2_vs_z[i]->GetXaxis()->GetXmax(), A_UL_2sin_inject);
        zeroLine_v->SetLineStyle(2);
        zeroLine_v->SetLineColor(kRed+1);
        zeroLine_v->Draw();

        c->Write();
    }





    dir_all_0->cd();
    // NO GAMMA correction
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_ALL_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_vs_xB_4 = new TGraphErrors();
    //TGraphErrors* graph2D_corr_sin = new TGraphErrors();
    int pll_idx1 = 0, pll_idx1_1 = 0, pll_idx1_2 = 0, pll_idx1_3 = 0, pll_idx1_4 = 0;
    for(int x = 0; x < nbin_xQ2; x++){
        if (vec_kaonp_xB_2d_mc[x].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_mc[x].empty()) continue;
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_mc[x]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_mc[x].size();
        for (double val : vec_kaonp_xB_2d_mc[x]){
            sum_xB += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_xB = sum_xB / vec_kaonp_xB_2d_mc[x].size();
        // Get AUL and error from vectors
        double all_0 = ALL_0_2d_mc[x]; // correction for TSA
        double all_0_err = ALL_0_err_2d_mc[x];
        //double all_0 = 0;
        //if(all_0 > 1 || all_0 < -1) continue;
        //graph2D_corr_sin->SetPoint(pll_idx1, mean_xB, Corr_sin_2d[x]);
        pll_idx1++;
        if(x < 5){
            graph2D_ALL_vs_xB_1->SetPoint(pll_idx1_1, mean_xB, all_0);
            graph2D_ALL_vs_xB_1->SetPointError(pll_idx1_1, 0.0, all_0_err); // No x error
            pll_idx1_1++;
        } else if (x < 10){
            graph2D_ALL_vs_xB_2->SetPoint(pll_idx1_2, mean_xB, all_0);
            graph2D_ALL_vs_xB_2->SetPointError(pll_idx1_2, 0.0, all_0_err);
            pll_idx1_2++;
        } else if (x < 12){
            graph2D_ALL_vs_xB_3->SetPoint(pll_idx1_3, mean_xB, all_0);
            graph2D_ALL_vs_xB_3->SetPointError(pll_idx1_3, 0.0, all_0_err);
            pll_idx1_3++;
        } else if (x < 18){
            graph2D_ALL_vs_xB_4->SetPoint(pll_idx1_4, mean_xB, all_0);
            graph2D_ALL_vs_xB_4->SetPointError(pll_idx1_4, 0.0, all_0_err);
            pll_idx1_4++;
        }
    }
    graph2D_ALL_vs_xB_1->SetMarkerStyle(20), graph2D_ALL_vs_xB_2->SetMarkerStyle(20), graph2D_ALL_vs_xB_3->SetMarkerStyle(20), graph2D_ALL_vs_xB_4->SetMarkerStyle(20);
    graph2D_ALL_vs_xB_1->SetLineColor(kAzure-5), graph2D_ALL_vs_xB_1->SetMarkerColor(kAzure-5);
    graph2D_ALL_vs_xB_2->SetLineColor(kViolet-5), graph2D_ALL_vs_xB_2->SetMarkerColor(kViolet-5);
    graph2D_ALL_vs_xB_3->SetLineColor(kPink-5), graph2D_ALL_vs_xB_3->SetMarkerColor(kPink-5);
    graph2D_ALL_vs_xB_4->SetLineColor(kOrange-5), graph2D_ALL_vs_xB_4->SetMarkerColor(kOrange-5);


    TCanvas* c_All2D_xB_0 = new TCanvas("ALL_0_vs_xB_2d", "sin(#Phi_{h}) longitudinal asymmetry vs x_{B}", 800, 600);
    //graph2D_ALL_vs_xB->Draw("A");
    TMultiGraph *mg_All2D_xB_0 = new TMultiGraph();
    mg_All2D_xB_0->SetTitle("A_{LL} vs x_{B} (x_{B}-Q^{2} bin) | lepton frame; x_{B}; F_{LL}/F_{UU}");
    mg_All2D_xB_0->Add(graph2D_ALL_vs_xB_1, "P");
    mg_All2D_xB_0->Add(graph2D_ALL_vs_xB_2, "P");
    mg_All2D_xB_0->Add(graph2D_ALL_vs_xB_3, "P");
    mg_All2D_xB_0->Add(graph2D_ALL_vs_xB_4, "P");
    //mg_All2D_xB_0->Add(graph2D_corr_sin, "P");
    mg_All2D_xB_0->Draw("A");
    mg_All2D_xB_0->GetYaxis()->SetRangeUser(-0.1, 0.8);
    TLine* guideLine_ll0 = new TLine(mg_All2D_xB_0->GetXaxis()->GetXmin(), 0, mg_All2D_xB_0->GetXaxis()->GetXmax(), 0);
    guideLine_ll0->SetLineStyle(2);  
    guideLine_ll0->SetLineColor(kGray+1);
    guideLine_ll0->Draw();
    TLine* guideLine_ll0_val = new TLine(mg_All2D_xB_0->GetXaxis()->GetXmin(), A_LL_0_inject, mg_All2D_xB_0->GetXaxis()->GetXmax(), A_LL_0_inject);
    guideLine_ll0_val->SetLineStyle(2);  
    guideLine_ll0_val->SetLineColor(kRed+1);
    guideLine_ll0_val->Draw();

    TLegend* legend_ll0 = new TLegend(0.13, 0.7, 0.33, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_ll0->AddEntry(graph2D_ALL_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend_ll0->AddEntry(graph2D_ALL_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend_ll0->AddEntry(graph2D_ALL_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend_ll0->AddEntry(graph2D_ALL_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
    legend_ll0->SetFillStyle(0);  // Transparent background
    legend_ll0->Draw();
    c_All2D_xB_0->Update();
    c_All2D_xB_0->Write();






    // zPt
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_ALL0_vs_z = new TGraphErrors();
    TGraphErrors* graph2D_ALL0_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph2D_ALL0_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph2D_ALL0_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph2D_ALL0_vs_z_4 = new TGraphErrors();
    //TGraphErrors* graph2D_corr_sin = new TGraphErrors();
    int pzll0_idx = 0, pzll0_idx_1 = 0, pzll0_idx_2 = 0, pzll0_idx_3 = 0, pzll0_idx_4 = 0;
    for(int z = 0; z < nbin_zPt; z++){
        if (vec_kaonp_xB_2d_zPt_mc[z].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_zPt_mc[z].empty()) continue;
        // Compute mean xB for this bin
        double sum_z = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_zPt_mc[z]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_zPt_mc[z].size();
        for (double val : vec_kaonp_z_2d_zPt_mc[z]){
            sum_z += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_z = sum_z / vec_kaonp_z_2d_zPt_mc[z].size();
        double all_0 = ALL_0_2d_zPt_mc[z];
        double all_0_err = ALL_0_err_2d_zPt_mc[z];
        //double all_0 = 0;
        //if(all_0 > 1 || all_0 < -1) continue;
        graph2D_ALL0_vs_z->SetPoint(pzll0_idx, mean_z, all_0);
        graph2D_ALL0_vs_z->SetPointError(pzll0_idx, 0.0, all_0_err); // No x error
        pzll0_idx++;
        if(z < 7){
            graph2D_ALL0_vs_z_1->SetPoint(pzll0_idx_1, mean_z, all_0);
            graph2D_ALL0_vs_z_1->SetPointError(pzll0_idx_1, 0.0, all_0_err); // No x error
            pzll0_idx_1++;
        } else if (z < 14){
            graph2D_ALL0_vs_z_2->SetPoint(pzll0_idx_2, mean_z, all_0);
            graph2D_ALL0_vs_z_2->SetPointError(pzll0_idx_2, 0.0, all_0_err);
            pzll0_idx_2++;
        } else if (z < 21){
            graph2D_ALL0_vs_z_3->SetPoint(pzll0_idx_3, mean_z, all_0);
            graph2D_ALL0_vs_z_3->SetPointError(pzll0_idx_3, 0.0, all_0_err);
            pzll0_idx_3++;
        } else if (z < 26){
            graph2D_ALL0_vs_z_4->SetPoint(pzll0_idx_4, mean_z, all_0);
            graph2D_ALL0_vs_z_4->SetPointError(pzll0_idx_4, 0.0, all_0_err);
            pzll0_idx_4++;
        }
    }
    graph2D_ALL0_vs_z_1->SetMarkerStyle(20), graph2D_ALL0_vs_z_2->SetMarkerStyle(20), graph2D_ALL0_vs_z_3->SetMarkerStyle(20), graph2D_ALL0_vs_z_4->SetMarkerStyle(20);
    graph2D_ALL0_vs_z_1->SetLineColor(kAzure-5), graph2D_ALL0_vs_z_1->SetMarkerColor(kAzure-5);
    graph2D_ALL0_vs_z_2->SetLineColor(kViolet-5), graph2D_ALL0_vs_z_2->SetMarkerColor(kViolet-5);
    graph2D_ALL0_vs_z_3->SetLineColor(kPink-5), graph2D_ALL0_vs_z_3->SetMarkerColor(kPink-5);
    graph2D_ALL0_vs_z_4->SetLineColor(kOrange-5), graph2D_ALL0_vs_z_4->SetMarkerColor(kOrange-5);



    vector<TGraphErrors*> graphs_ALL0_vs_z = {graph2D_ALL0_vs_z_1,graph2D_ALL0_vs_z_2,graph2D_ALL0_vs_z_3,graph2D_ALL0_vs_z_4};
    vector<string> titles_ALL0_vs_z = {"0.0 < P_{hT} < 0.25 GeV","0.25 < P_{hT} < 0.5 GeV","0.5 < P_{hT} < 0.8 GeV","0.8 < P_{hT} < 1.4 GeV"};

    for (size_t i = 0; i < graphs_ALL0_vs_z.size(); ++i) {
        TString cname = Form("c_ALL_0_vs_z_PtBin_%zu", i);
        TString ctitle = Form("A_{LL} vs z (%s)", titles_ALL0_vs_z[i].c_str());

        TCanvas* c = new TCanvas(cname, ctitle, 800, 600);
        graphs_ALL0_vs_z[i]->SetTitle(Form("A_{LL} vs z | %s | lepton frame; z; F_{LL}/F_{UU}",titles_ALL0_vs_z[i].c_str()));
        graphs_ALL0_vs_z[i]->Draw("AP");
        graphs_ALL0_vs_z[i]->GetYaxis()->SetRangeUser(-0.1, 0.8);

        TLine* zeroLine = new TLine(graphs_ALL0_vs_z[i]->GetXaxis()->GetXmin(), 0,graphs_ALL0_vs_z[i]->GetXaxis()->GetXmax(), 0);
        zeroLine->SetLineStyle(2);
        zeroLine->SetLineColor(kGray+1);
        zeroLine->Draw();
        TLine* zeroLine_v = new TLine(graphs_ALL0_vs_z[i]->GetXaxis()->GetXmin(), A_LL_0_inject,graphs_ALL0_vs_z[i]->GetXaxis()->GetXmax(), A_LL_0_inject);
        zeroLine_v->SetLineStyle(2);
        zeroLine_v->SetLineColor(kRed+1);
        zeroLine_v->Draw();

        c->Write();
    }




    dir_all_cos->cd();
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_ALL_cos_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_cos_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_cos_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_cos_vs_xB_4 = new TGraphErrors();
    //TGraphErrors* graph2D_corr_sin = new TGraphErrors();
    int pll2_idx1 = 0, pll2_idx1_1 = 0, pll2_idx1_2 = 0, pll2_idx1_3 = 0, pll2_idx1_4 = 0;
    for(int x = 0; x < nbin_xQ2; x++){
        if (vec_kaonp_xB_2d_mc[x].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_mc[x].empty()) continue;
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_mc[x]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_mc[x].size();
        for (double val : vec_kaonp_xB_2d_mc[x]){
            sum_xB += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_xB = sum_xB / vec_kaonp_xB_2d_mc[x].size();
        // Get AUL and error from vectors
        double all_cos = ALL_cos_2d_mc[x]; // correction for TSA
        double all_cos_err = ALL_cos_err_2d_mc[x];
        pll2_idx1++;
        if(x < 5){
            graph2D_ALL_cos_vs_xB_1->SetPoint(pll2_idx1_1, mean_xB, all_cos);
            graph2D_ALL_cos_vs_xB_1->SetPointError(pll2_idx1_1, 0.0, all_cos_err); // No x error
            pll2_idx1_1++;
        } else if (x < 10){
            graph2D_ALL_cos_vs_xB_2->SetPoint(pll2_idx1_2, mean_xB, all_cos);
            graph2D_ALL_cos_vs_xB_2->SetPointError(pll2_idx1_2, 0.0, all_cos_err);
            pll2_idx1_2++;
        } else if (x < 12){
            graph2D_ALL_cos_vs_xB_3->SetPoint(pll2_idx1_3, mean_xB, all_cos);
            graph2D_ALL_cos_vs_xB_3->SetPointError(pll2_idx1_3, 0.0, all_cos_err);
            pll2_idx1_3++;
        } else if (x < 18){
            graph2D_ALL_cos_vs_xB_4->SetPoint(pll2_idx1_4, mean_xB, all_cos);
            graph2D_ALL_cos_vs_xB_4->SetPointError(pll2_idx1_4, 0.0, all_cos_err);
            pll2_idx1_4++;
        }
    }
    graph2D_ALL_cos_vs_xB_1->SetMarkerStyle(20), graph2D_ALL_cos_vs_xB_2->SetMarkerStyle(20), graph2D_ALL_cos_vs_xB_3->SetMarkerStyle(20), graph2D_ALL_cos_vs_xB_4->SetMarkerStyle(20);
    graph2D_ALL_cos_vs_xB_1->SetLineColor(kAzure-5), graph2D_ALL_cos_vs_xB_1->SetMarkerColor(kAzure-5);
    graph2D_ALL_cos_vs_xB_2->SetLineColor(kViolet-5), graph2D_ALL_cos_vs_xB_2->SetMarkerColor(kViolet-5);
    graph2D_ALL_cos_vs_xB_3->SetLineColor(kPink-5), graph2D_ALL_cos_vs_xB_3->SetMarkerColor(kPink-5);
    graph2D_ALL_cos_vs_xB_4->SetLineColor(kOrange-5), graph2D_ALL_cos_vs_xB_4->SetMarkerColor(kOrange-5);


    TCanvas* c_All2D_xB_cos = new TCanvas("ALL_cos_vs_xB_2d", "sin(#Phi_{h}) longitudinal asymmetry vs x_{B}", 800, 600);
    //graph2D_ALL_cos_vs_xB->Draw("A");
    TMultiGraph *mg_All2D_xB_cos = new TMultiGraph();
    mg_All2D_xB_cos->SetTitle("A_{LL}^{cos#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin) | lepton frame; x_{B}; F_{LL}^{cos#Phi_{h}}/F_{UU}");
    mg_All2D_xB_cos->Add(graph2D_ALL_cos_vs_xB_1, "P");
    mg_All2D_xB_cos->Add(graph2D_ALL_cos_vs_xB_2, "P");
    mg_All2D_xB_cos->Add(graph2D_ALL_cos_vs_xB_3, "P");
    mg_All2D_xB_cos->Add(graph2D_ALL_cos_vs_xB_4, "P");
    //mg_All2D_xB_cos->Add(graph2D_corr_sin, "P");
    mg_All2D_xB_cos->Draw("A");
    mg_All2D_xB_cos->GetYaxis()->SetRangeUser(-0.4, 0.4);
    TLine* guideLine_llcos = new TLine(mg_All2D_xB_cos->GetXaxis()->GetXmin(), 0, mg_All2D_xB_cos->GetXaxis()->GetXmax(), 0);
    guideLine_llcos->SetLineStyle(2);  
    guideLine_llcos->SetLineColor(kGray+1);
    guideLine_llcos->Draw();
    TLine* guideLine_llcos_val = new TLine(mg_All2D_xB_cos->GetXaxis()->GetXmin(), A_LL_cos_inject, mg_All2D_xB_cos->GetXaxis()->GetXmax(), A_LL_cos_inject);
    guideLine_llcos_val->SetLineStyle(2);  
    guideLine_llcos_val->SetLineColor(kRed+1);
    guideLine_llcos_val->Draw();

    TLegend* legend_llcos = new TLegend(0.13, 0.7, 0.33, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_llcos->AddEntry(graph2D_ALL_cos_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend_llcos->AddEntry(graph2D_ALL_cos_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend_llcos->AddEntry(graph2D_ALL_cos_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend_llcos->AddEntry(graph2D_ALL_cos_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
    legend_llcos->SetFillStyle(0);  // Transparent background
    legend_llcos->Draw();
    c_All2D_xB_cos->Update();
    c_All2D_xB_cos->Write();




    // zPt
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_ALL_cos_vs_z = new TGraphErrors();
    TGraphErrors* graph2D_ALL_cos_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_cos_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_cos_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph2D_ALL_cos_vs_z_4 = new TGraphErrors();
    //TGraphErrors* graph2D_corr_sin = new TGraphErrors();
    int pzllcos_idx = 0, pzllcos_idx_1 = 0, pzllcos_idx_2 = 0, pzllcos_idx_3 = 0, pzllcos_idx_4 = 0;
    for(int z = 0; z < nbin_zPt; z++){
        if (vec_kaonp_xB_2d_zPt_mc[z].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_zPt_mc[z].empty()) continue;
        // Compute mean xB for this bin
        double sum_z = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_zPt_mc[z]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_zPt_mc[z].size();
        for (double val : vec_kaonp_z_2d_zPt_mc[z]){
            sum_z += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_z = sum_z / vec_kaonp_z_2d_zPt_mc[z].size();
        // Get ALL_cos and error from vectors
        double all_cos = ALL_cos_2d_zPt_mc[z];
        double all_cos_err = ALL_cos_err_2d_zPt_mc[z];
        //double all_cos = 0;
        //if(all_cos > 1 || all_cos < -1) continue;
        graph2D_ALL_cos_vs_z->SetPoint(pzllcos_idx, mean_z, all_cos);
        graph2D_ALL_cos_vs_z->SetPointError(pzllcos_idx, 0.0, all_cos_err); // No x error
        //graph2D_corr_sin->SetPoint(pzllcos_idx, mean_z, Corr_sin_2d[z]);
        pzllcos_idx++;
        if(z < 7){
            graph2D_ALL_cos_vs_z_1->SetPoint(pzllcos_idx_1, mean_z, all_cos);
            graph2D_ALL_cos_vs_z_1->SetPointError(pzllcos_idx_1, 0.0, all_cos_err); // No x error
            pzllcos_idx_1++;
        } else if (z < 14){
            graph2D_ALL_cos_vs_z_2->SetPoint(pzllcos_idx_2, mean_z, all_cos);
            graph2D_ALL_cos_vs_z_2->SetPointError(pzllcos_idx_2, 0.0, all_cos_err);
            pzllcos_idx_2++;
        } else if (z < 21){
            graph2D_ALL_cos_vs_z_3->SetPoint(pzllcos_idx_3, mean_z, all_cos);
            graph2D_ALL_cos_vs_z_3->SetPointError(pzllcos_idx_3, 0.0, all_cos_err);
            pzllcos_idx_3++;
        } else if (z < 26){
            graph2D_ALL_cos_vs_z_4->SetPoint(pzllcos_idx_4, mean_z, all_cos);
            graph2D_ALL_cos_vs_z_4->SetPointError(pzllcos_idx_4, 0.0, all_cos_err);
            pzllcos_idx_4++;
        }
    }
    graph2D_ALL_cos_vs_z_1->SetMarkerStyle(20), graph2D_ALL_cos_vs_z_2->SetMarkerStyle(20), graph2D_ALL_cos_vs_z_3->SetMarkerStyle(20), graph2D_ALL_cos_vs_z_4->SetMarkerStyle(20);
    graph2D_ALL_cos_vs_z_1->SetLineColor(kAzure-5), graph2D_ALL_cos_vs_z_1->SetMarkerColor(kAzure-5);
    graph2D_ALL_cos_vs_z_2->SetLineColor(kViolet-5), graph2D_ALL_cos_vs_z_2->SetMarkerColor(kViolet-5);
    graph2D_ALL_cos_vs_z_3->SetLineColor(kPink-5), graph2D_ALL_cos_vs_z_3->SetMarkerColor(kPink-5);
    graph2D_ALL_cos_vs_z_4->SetLineColor(kOrange-5), graph2D_ALL_cos_vs_z_4->SetMarkerColor(kOrange-5);



    vector<TGraphErrors*> graphs_ALL_cos_vs_z = {graph2D_ALL_cos_vs_z_1,graph2D_ALL_cos_vs_z_2,graph2D_ALL_cos_vs_z_3,graph2D_ALL_cos_vs_z_4};
    vector<string> titles_ALL_cos_vs_z = {"0.0 < P_{hT} < 0.25 GeV","0.25 < P_{hT} < 0.5 GeV","0.5 < P_{hT} < 0.8 GeV","0.8 < P_{hT} < 1.4 GeV"};

    for (size_t i = 0; i < graphs_ALL_cos_vs_z.size(); ++i) {
        TString cname = Form("c_ALL_cos_vs_z_PtBin_%zu", i);
        TString ctitle = Form("A_{LL}^{cos#Phi_{h}} vs z (%s)", titles_ALL_cos_vs_z[i].c_str());

        TCanvas* c = new TCanvas(cname, ctitle, 800, 600);
        graphs_ALL_cos_vs_z[i]->SetTitle(Form("A_{LL}^{cos#Phi_{h}} vs z | %s | lepton frame; z; F_{LL}^{cos#Phi_{h}}/F_{UU}",titles_ALL_cos_vs_z[i].c_str()));
        graphs_ALL_cos_vs_z[i]->Draw("AP");
        graphs_ALL_cos_vs_z[i]->GetYaxis()->SetRangeUser(-0.4, 0.4);

        TLine* zeroLine = new TLine(graphs_ALL_cos_vs_z[i]->GetXaxis()->GetXmin(), 0,graphs_ALL_cos_vs_z[i]->GetXaxis()->GetXmax(), 0);
        zeroLine->SetLineStyle(2);
        zeroLine->SetLineColor(kGray+1);
        zeroLine->Draw();
        TLine* zeroLine_v = new TLine(graphs_ALL_cos_vs_z[i]->GetXaxis()->GetXmin(), A_LL_cos_inject,graphs_ALL_cos_vs_z[i]->GetXaxis()->GetXmax(), A_LL_cos_inject);
        zeroLine_v->SetLineStyle(2);
        zeroLine_v->SetLineColor(kRed+1);
        zeroLine_v->Draw();

        c->Write();
    }




    // LU 


    dir_alu_sin->cd();
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_ALU_sin_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph2D_ALU_sin_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph2D_ALU_sin_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph2D_ALU_sin_vs_xB_4 = new TGraphErrors();
    //TGraphErrors* graph2D_corr_sin = new TGraphErrors();
    int plu2_idx1 = 0, plu2_idx1_1 = 0, plu2_idx1_2 = 0, plu2_idx1_3 = 0, plu2_idx1_4 = 0;
    for(int x = 0; x < nbin_xQ2; x++){
        if (vec_kaonp_xB_2d_mc[x].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_mc[x].empty()) continue;
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_mc[x]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_mc[x].size();
        for (double val : vec_kaonp_xB_2d_mc[x]){
            sum_xB += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_xB = sum_xB / vec_kaonp_xB_2d_mc[x].size();
        // Get AUL and error from vectors
        double alu_sin = ALU_sin_2d_mc[x]; // correction for TSA
        double alu_sin_err = ALU_sin_err_2d_mc[x];
        plu2_idx1++;
        if(x < 5){
            graph2D_ALU_sin_vs_xB_1->SetPoint(plu2_idx1_1, mean_xB, alu_sin);
            graph2D_ALU_sin_vs_xB_1->SetPointError(plu2_idx1_1, 0.0, alu_sin_err); // No x error
            plu2_idx1_1++;
        } else if (x < 10){
            graph2D_ALU_sin_vs_xB_2->SetPoint(plu2_idx1_2, mean_xB, alu_sin);
            graph2D_ALU_sin_vs_xB_2->SetPointError(plu2_idx1_2, 0.0, alu_sin_err);
            plu2_idx1_2++;
        } else if (x < 12){
            graph2D_ALU_sin_vs_xB_3->SetPoint(plu2_idx1_3, mean_xB, alu_sin);
            graph2D_ALU_sin_vs_xB_3->SetPointError(plu2_idx1_3, 0.0, alu_sin_err);
            plu2_idx1_3++;
        } else if (x < 18){
            graph2D_ALU_sin_vs_xB_4->SetPoint(plu2_idx1_4, mean_xB, alu_sin);
            graph2D_ALU_sin_vs_xB_4->SetPointError(plu2_idx1_4, 0.0, alu_sin_err);
            plu2_idx1_4++;
        }
    }
    graph2D_ALU_sin_vs_xB_1->SetMarkerStyle(20), graph2D_ALU_sin_vs_xB_2->SetMarkerStyle(20), graph2D_ALU_sin_vs_xB_3->SetMarkerStyle(20), graph2D_ALU_sin_vs_xB_4->SetMarkerStyle(20);
    graph2D_ALU_sin_vs_xB_1->SetLineColor(kAzure-5), graph2D_ALU_sin_vs_xB_1->SetMarkerColor(kAzure-5);
    graph2D_ALU_sin_vs_xB_2->SetLineColor(kViolet-5), graph2D_ALU_sin_vs_xB_2->SetMarkerColor(kViolet-5);
    graph2D_ALU_sin_vs_xB_3->SetLineColor(kPink-5), graph2D_ALU_sin_vs_xB_3->SetMarkerColor(kPink-5);
    graph2D_ALU_sin_vs_xB_4->SetLineColor(kOrange-5), graph2D_ALU_sin_vs_xB_4->SetMarkerColor(kOrange-5);


    TCanvas* c_Alu2D_xB_sin = new TCanvas("ALU_sin_vs_xB_2d", "sin(#Phi_{h}) longitudinal asymmetry vs x_{B}", 800, 600);
    //graph2D_ALU_sin_vs_xB->Draw("A");
    TMultiGraph *mg_Alu2D_xB_sin = new TMultiGraph();
    mg_Alu2D_xB_sin->SetTitle("A_{LU}^{sin#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin) | lepton frame; x_{B}; F_{LU}^{sin#Phi_{h}}/F_{UU}");
    mg_Alu2D_xB_sin->Add(graph2D_ALU_sin_vs_xB_1, "P");
    mg_Alu2D_xB_sin->Add(graph2D_ALU_sin_vs_xB_2, "P");
    mg_Alu2D_xB_sin->Add(graph2D_ALU_sin_vs_xB_3, "P");
    mg_Alu2D_xB_sin->Add(graph2D_ALU_sin_vs_xB_4, "P");
    //mg_Alu2D_xB_sin->Add(graph2D_corr_sin, "P");
    mg_Alu2D_xB_sin->Draw("A");
    mg_Alu2D_xB_sin->GetYaxis()->SetRangeUser(-0.2, 0.2);
    TLine* guideLine_lusin = new TLine(mg_Alu2D_xB_sin->GetXaxis()->GetXmin(), 0, mg_Alu2D_xB_sin->GetXaxis()->GetXmax(), 0);
    guideLine_lusin->SetLineStyle(2);  
    guideLine_lusin->SetLineColor(kGray+1);
    guideLine_lusin->Draw();

    TLegend* legend_lusin = new TLegend(0.13, 0.7, 0.33, 0.88); // Adjust position (x1,y1,x2,y2)
    legend_lusin->AddEntry(graph2D_ALU_sin_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend_lusin->AddEntry(graph2D_ALU_sin_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend_lusin->AddEntry(graph2D_ALU_sin_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend_lusin->AddEntry(graph2D_ALU_sin_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
    legend_lusin->SetFillStyle(0);  // Transparent background
    legend_lusin->Draw();
    c_Alu2D_xB_sin->Update();
    c_Alu2D_xB_sin->Write();




    // zPt
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_ALU_sin_vs_z = new TGraphErrors();
    TGraphErrors* graph2D_ALU_sin_vs_z_1 = new TGraphErrors();
    TGraphErrors* graph2D_ALU_sin_vs_z_2 = new TGraphErrors();
    TGraphErrors* graph2D_ALU_sin_vs_z_3 = new TGraphErrors();
    TGraphErrors* graph2D_ALU_sin_vs_z_4 = new TGraphErrors();
    //TGraphErrors* graph2D_corr_sin = new TGraphErrors();
    int pzlusin_idx = 0, pzlusin_idx_1 = 0, pzlusin_idx_2 = 0, pzlusin_idx_3 = 0, pzlusin_idx_4 = 0;
    for(int z = 0; z < nbin_zPt; z++){
        if (vec_kaonp_xB_2d_zPt_mc[z].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d_zPt_mc[z].empty()) continue;
        // Compute mean xB for this bin
        double sum_z = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        for (double valeps : vec_kaonp_eps_2d_zPt_mc[z]) sum_eps += valeps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d_zPt_mc[z].size();
        for (double val : vec_kaonp_z_2d_zPt_mc[z]){
            sum_z += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_z = sum_z / vec_kaonp_z_2d_zPt_mc[z].size();
        // Get ALU_sin and error from vectors
        double alu_sin = ALU_sin_2d_zPt_mc[z];
        double alu_sin_err = ALU_sin_err_2d_zPt_mc[z];
        //double alu_sin = 0;
        //if(alu_sin > 1 || alu_sin < -1) continue;
        graph2D_ALU_sin_vs_z->SetPoint(pzlusin_idx, mean_z, alu_sin);
        graph2D_ALU_sin_vs_z->SetPointError(pzlusin_idx, 0.0, alu_sin_err); // No x error
        pzlusin_idx++;
        if(z < 7){
            graph2D_ALU_sin_vs_z_1->SetPoint(pzlusin_idx_1, mean_z, alu_sin);
            graph2D_ALU_sin_vs_z_1->SetPointError(pzlusin_idx_1, 0.0, alu_sin_err); // No x error
            pzlusin_idx_1++;
        } else if (z < 14){
            graph2D_ALU_sin_vs_z_2->SetPoint(pzlusin_idx_2, mean_z, alu_sin);
            graph2D_ALU_sin_vs_z_2->SetPointError(pzlusin_idx_2, 0.0, alu_sin_err);
            pzlusin_idx_2++;
        } else if (z < 21){
            graph2D_ALU_sin_vs_z_3->SetPoint(pzlusin_idx_3, mean_z, alu_sin);
            graph2D_ALU_sin_vs_z_3->SetPointError(pzlusin_idx_3, 0.0, alu_sin_err);
            pzlusin_idx_3++;
        } else if (z < 26){
            graph2D_ALU_sin_vs_z_4->SetPoint(pzlusin_idx_4, mean_z, alu_sin);
            graph2D_ALU_sin_vs_z_4->SetPointError(pzlusin_idx_4, 0.0, alu_sin_err);
            pzlusin_idx_4++;
        }
    }
    graph2D_ALU_sin_vs_z_1->SetMarkerStyle(20), graph2D_ALU_sin_vs_z_2->SetMarkerStyle(20), graph2D_ALU_sin_vs_z_3->SetMarkerStyle(20), graph2D_ALU_sin_vs_z_4->SetMarkerStyle(20);
    graph2D_ALU_sin_vs_z_1->SetLineColor(kAzure-5), graph2D_ALU_sin_vs_z_1->SetMarkerColor(kAzure-5);
    graph2D_ALU_sin_vs_z_2->SetLineColor(kViolet-5), graph2D_ALU_sin_vs_z_2->SetMarkerColor(kViolet-5);
    graph2D_ALU_sin_vs_z_3->SetLineColor(kPink-5), graph2D_ALU_sin_vs_z_3->SetMarkerColor(kPink-5);
    graph2D_ALU_sin_vs_z_4->SetLineColor(kOrange-5), graph2D_ALU_sin_vs_z_4->SetMarkerColor(kOrange-5);



    vector<TGraphErrors*> graphs_ALU_sin_vs_z = {graph2D_ALU_sin_vs_z_1,graph2D_ALU_sin_vs_z_2,graph2D_ALU_sin_vs_z_3,graph2D_ALU_sin_vs_z_4};
    vector<string> titles_ALU_sin_vs_z = {"0.0 < P_{hT} < 0.25 GeV","0.25 < P_{hT} < 0.5 GeV","0.5 < P_{hT} < 0.8 GeV","0.8 < P_{hT} < 1.4 GeV"};

    for (size_t i = 0; i < graphs_ALU_sin_vs_z.size(); ++i) {
        TString cname = Form("c_ALU_sin_vs_z_PtBin_%zu", i);
        TString ctitle = Form("A_{LU}^{sin#Phi_{h}} vs z (%s)", titles_ALU_sin_vs_z[i].c_str());

        TCanvas* c = new TCanvas(cname, ctitle, 800, 600);
        graphs_ALU_sin_vs_z[i]->SetTitle(Form("A_{LU}^{sin#Phi_{h}} vs z | %s | lepton frame; z; F_{LU}^{sin#Phi_{h}}/F_{UU}",titles_ALU_sin_vs_z[i].c_str()));
        graphs_ALU_sin_vs_z[i]->Draw("AP");
        graphs_ALU_sin_vs_z[i]->GetYaxis()->SetRangeUser(-0.2, 0.2);

        TLine* zeroLine = new TLine(graphs_ALU_sin_vs_z[i]->GetXaxis()->GetXmin(), 0,graphs_ALU_sin_vs_z[i]->GetXaxis()->GetXmax(), 0);
        zeroLine->SetLineStyle(2);
        zeroLine->SetLineColor(kGray+1);
        zeroLine->Draw();

        c->Write();
    }


    // ------------------------------------------------------------------------------------------------------------------------------


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
                << 0 << "," << 0 << "," 
                << AUL_sin_2d_zPt_preID[iz] - AUL_sin_2d_zPt_true[iz] << "," << sqrt(AUL_sin_err_2d_zPt_preID[iz]*AUL_sin_err_2d_zPt_preID[iz] - AUL_sin_err_2d_zPt_true[iz]*AUL_sin_err_2d_zPt_true[iz]) << ","
                << 0 << "," << 0 << ","
                <<  AUL_sin_2d_zPt_mc[iz] - AUL_sin_2d_zPt_preID[iz]<< "," << sqrt(AUL_sin_err_2d_zPt_mc[iz]*AUL_sin_err_2d_zPt_mc[iz] - AUL_sin_err_2d_zPt_preID[iz]*AUL_sin_err_2d_zPt_preID[iz]) << ","
                << 0 << "," << 0 << endl;
        // riga per AUL_sinPhi
        csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                    << "AUL_sinPhi" << ","
                    <<  std::abs(AUL_sin_2d_zPt_mc[iz] - AUL_sin_2d_zPt_preID[iz]) << "," << safeSigmaDiff(AUL_sin_err_2d_zPt_mc[iz], AUL_sin_err_2d_zPt_preID[iz]) << ","
                    << std::abs(AUL_sin_2d_zPt_preID[iz] - AUL_sin_2d_zPt_true[iz]) << "," << safeSigmaDiff(AUL_sin_err_2d_zPt_preID[iz], AUL_sin_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(AUL_sin_2d_zPt_all_ID[iz] - AUL_sin_2d_zPt_true[iz])<< "," << safeSigmaDiff(AUL_sin_err_2d_zPt_all_ID[iz], AUL_sin_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(AUL_sin_2d_zPt_true[iz] - AUL_sin_2d_zPt_reco[iz])<< "," << safeSigmaDiff(AUL_sin_err_2d_zPt_true[iz], AUL_sin_err_2d_zPt_reco[iz]) << "\n";
        csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                    << "AUL_sin2Phi" << ","
                    <<  std::abs(AUL_2sin_2d_zPt_mc[iz] - AUL_2sin_2d_zPt_preID[iz])<< "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_mc[iz], AUL_2sin_err_2d_zPt_preID[iz]) << ","
                    << std::abs(AUL_2sin_2d_zPt_preID[iz] - AUL_2sin_2d_zPt_true[iz]) << "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_preID[iz], AUL_2sin_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(AUL_2sin_2d_zPt_all_ID[iz] - AUL_2sin_2d_zPt_true[iz])<< "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_all_ID[iz], AUL_2sin_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(AUL_2sin_2d_zPt_true[iz] - AUL_2sin_2d_zPt_reco[iz])<< "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_true[iz], AUL_2sin_err_2d_zPt_reco[iz]) << "\n";
        csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                    << "ALL_const" << ","
                    <<  std::abs(ALL_0_2d_zPt_mc[iz] - ALL_0_2d_zPt_preID[iz])<< "," << safeSigmaDiff(ALL_0_err_2d_zPt_mc[iz], ALL_0_err_2d_zPt_preID[iz]) << ","
                    << std::abs(ALL_0_2d_zPt_preID[iz] - ALL_0_2d_zPt_true[iz]) << "," << safeSigmaDiff(ALL_0_err_2d_zPt_preID[iz], ALL_0_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(ALL_0_2d_zPt_all_ID[iz] - ALL_0_2d_zPt_true[iz])<< "," << safeSigmaDiff(ALL_0_err_2d_zPt_all_ID[iz], ALL_0_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(ALL_0_2d_zPt_true[iz] - ALL_0_2d_zPt_reco[iz])<< "," << safeSigmaDiff(ALL_0_err_2d_zPt_true[iz], ALL_0_err_2d_zPt_reco[iz]) << "\n";
        csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                    << "ALL_cosPhi" << ","
                    <<  std::abs(ALL_cos_2d_zPt_mc[iz] - ALL_cos_2d_zPt_preID[iz])<< "," << safeSigmaDiff(ALL_cos_err_2d_zPt_mc[iz], ALL_cos_err_2d_zPt_preID[iz]) << ","
                    << std::abs(ALL_cos_2d_zPt_preID[iz] - ALL_cos_2d_zPt_true[iz]) << "," << safeSigmaDiff(ALL_cos_err_2d_zPt_preID[iz], ALL_cos_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(ALL_cos_2d_zPt_all_ID[iz] - ALL_cos_2d_zPt_true[iz])<< "," << safeSigmaDiff(ALL_cos_err_2d_zPt_all_ID[iz], ALL_cos_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(ALL_cos_2d_zPt_true[iz] - ALL_cos_2d_zPt_reco[iz])<< "," << safeSigmaDiff(ALL_cos_err_2d_zPt_true[iz], ALL_cos_err_2d_zPt_reco[iz]) << "\n";
        csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                    << "ALU_sinPhi" << ","
                    <<  std::abs(ALU_sin_2d_zPt_mc[iz] - ALU_sin_2d_zPt_preID[iz])<< "," << safeSigmaDiff(ALU_sin_err_2d_zPt_mc[iz], ALU_sin_err_2d_zPt_preID[iz]) << ","
                    << std::abs(ALU_sin_2d_zPt_preID[iz] - ALU_sin_2d_zPt_true[iz]) << "," << safeSigmaDiff(ALU_sin_err_2d_zPt_preID[iz], ALU_sin_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(ALU_sin_2d_zPt_all_ID[iz] - ALU_sin_2d_zPt_true[iz])<< "," << safeSigmaDiff(ALU_sin_err_2d_zPt_all_ID[iz], ALU_sin_err_2d_zPt_true[iz]) << ","
                    <<  std::abs(ALU_sin_2d_zPt_true[iz] - ALU_sin_2d_zPt_reco[iz])<< "," << safeSigmaDiff(ALU_sin_err_2d_zPt_true[iz], ALU_sin_err_2d_zPt_reco[iz]) << "\n";
    }

    //outFile.Write();
    //treeKaonP.Write("", TObject::kOverwrite);
    csvFile_zPt.close();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << " and csv: " << csv_filename_zPt << endl;
}
