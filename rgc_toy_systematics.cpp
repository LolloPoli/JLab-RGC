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
                          const vector<double>& y_vec,       // inelasticity
                          const vector<double>& depol,       // eps
                          const vector<double>& Ptarget,     // target polarization
                          const vector<double>& helicity,    // beam helicity (+1 / -1)
                          const vector<double>& bootw, 
                          double r_local,
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

    double pBeam = 0.84;  // beam polarization
    double dilution = 0.25;  // not anymore 1.0

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

        // Full target spin-dependent term
        const double Seff = pol * (UL_mod + lam * pBeam *LL_mod);

        // --- Likelihood argument ---
        const double arg = 1.0 + UU_mod + lam*pBeam*(LU_mod) + pol*(UL_mod) + lam*pBeam*pol*(LL_mod);

        // Physical protection
        if (arg <= 0.0) return 1e12;
        /*
        double w_hel = 1.0;
        // luminosity beam state contribution
        double r = r_local;
        if (helicity[i] == +1) w_hel = 1.0;
        else if (helicity[i] == -1) w_hel = 1/r;
        */
        //logLike += bootw[i] * w_hel * log(arg); // with beam helicity weight
        logLike += bootw[i] * log(arg);
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

void rgc_toy_systematics(const char* period) {
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
    //TString outputFile = Form("fmax_%s.root", period);
    //TFile outFile(outputFile.Data(), "RECREATE");
    double torus = -1;

    TH1D kp_fmax ("_fmax", "f_{max} injection function ; f_{max}; count", 100, 0, 2);
    TH1D kp_fkeep ("_fkeep", "f_{keep} injection function ; f_{keep}; count", 100, 0, 2);

    //const char* csv_filename = "table_RGC_MC_summer22.csv";
    //TString csv_filename = Form("RGC_MC_PLOT/table_RGC_MC_%s.csv", period);
    //TString csv_filename_xQ2 = Form("RGC_MC_PLOT/table_RGC_MC_%s_xQ2.csv", period);
    //TString csv_filename_zPt = Form("RGC_MC_PLOT/table_RGC_MC_%s_zPt.csv", period);


    // root 'rgc_toy_systematics.cpp("fall22")' -l -b -q
    int toy_step = 500;
    for (int s = 0; s < toy_step; s++){
        TString csv_filename_zPt_test = Form("toy_model2/table_RGC_MC_%s_zPt_test_%d.csv", period,s);
        //std::ofstream csvFile(csv_filename.Data());
        //std::ofstream csvFile_xQ2(csv_filename_xQ2.Data());
        //std::ofstream csvFile_zPt(csv_filename_zPt.Data());
        std::ofstream csvFile_zPt_test(csv_filename_zPt_test.Data());
        //csvFile << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err,Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err, Phih_sys_sinx, Phih_sys_sinx_err, Phih_sys_sin2x, Phih_sys_sin2x_err, Phih_sys_cosx, Phih_sys_cosx_err, Phih_sys_cos2x, Phih_sys_cos2x_err, Auu_cos, Auu_cos_err, Auu_cos2, Auu_cos2_err, Auu_cont_sinx, Auu_cont_sin2x\n";
        //csvFile_xQ2 << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err, Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err, Phih_sys_sinx, Phih_sys_sinx_err, Phih_sys_sin2x, Phih_sys_sin2x_err, Phih_sys_cosx, Phih_sys_cosx_err, Phih_sys_cos2x, Phih_sys_cos2x_err, Auu_cos, Auu_cos_err, Auu_cos2, Auu_cos2_err, Auu_cont_sinx, Auu_cont_sin2x\n";
        //csvFile_zPt << " n_events, binxBQ2, bin_zPt, mean_xB, mean_Q2, mean_z, mean_PhT, epsilon, efficiency, eff_binom_err, PID_sys, PID_sys_err, Purity_sys, Purity_sys_err, Acc_sys, Acc_sys_err, Bin_mig_sys, Bin_mig_sys_err\n";
        csvFile_zPt_test << "n_events,binxBQ2,bin_zPt,mean_xB,mean_Q2,mean_z,mean_PhT,epsilon,"
                << "modulation,"
                << "Acc_sys,Acc_sys_err,PID_sys,PID_sys_err,Purity_sys,Purity_sys_err,Bin_mig_sys,Bin_mig_sys_err\n";


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
        

        // multidim

        double nbin_xQ2 = 14;
        double nbin_zPt = 25;
        // 2D
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
        vector<vector<double>> vec_kaonp_bootw_2d_preID(nbin_xQ2);
        vector<vector<double>> vec_kaonp_bootw_2d_zPt_preID(nbin_zPt);
        vector<vector<double>> vec_kaonp_bootw_2d_all_ID(nbin_xQ2);
        vector<vector<double>> vec_kaonp_bootw_2d_zPt_all_ID(nbin_zPt);
        vector<vector<double>> vec_kaonp_bootw_2d_true(nbin_xQ2);
        vector<vector<double>> vec_kaonp_bootw_2d_zPt_true(nbin_zPt);
        vector<vector<double>> vec_kaonp_bootw_2d_mc(nbin_xQ2);
        vector<vector<double>> vec_kaonp_bootw_2d_zPt_mc(nbin_zPt);
        vector<vector<double>> vec_kaonp_bootw_2d(nbin_xQ2);
        vector<vector<double>> vec_kaonp_bootw_2d_zPt(nbin_zPt);
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

        //
        double bin_z_plot[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
        double bin_Pt_plot[] = {0, 0.25, 0.5, 0.8, 1.4};
        

        vector<vector<vector<double>>> vec_kaonp_xB_4d_mc(nbin_xQ2, vector<vector<double>> (nbin_zPt));



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

        // fall22 ALU = 0
        // summer22 ALU = 0.05 and w
        double A_UL_sin_inject  = 0.0;
        double A_UL_2sin_inject = 0.0;
        double A_LL_0_inject    = 0.0;
        double A_LL_cos_inject  = 0.0;
        double A_LU_sin_inject = 0.0;
        double dilution = 0.25;
        double beam_pol = 0.84;
        double PbPt = 0.71;

        double f_max_inject = 1.0 + dilution*0.85*(std::abs(A_UL_sin_inject)+std::abs(A_UL_2sin_inject))
                                + dilution*0.85*beam_pol*(std::abs(A_LL_0_inject)+std::abs(A_LL_cos_inject))
                                + beam_pol*(std::abs(A_LU_sin_inject));

        TRandom3 rng(12345 + s);
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
            int fake_spinstate = (rng.Rndm() < 0.5) ? +1 : -1;   // one time per event
            //int fake_spinstate = (rng.Rndm() < 0.5788) ? +1 : -1;
            helicity_mc = (rng.Rndm() < 0.5) ? +1 : -1;
            double fake_Pt = 0.85 * fake_spinstate;
            fake_Ptarget_vec[i] = fake_Pt;
            double bootw = rng.PoissonD(1.0);

            double eps = kaonp_epsilon_mc, y = kaonp_y_mc;
            double A = (y*y)/(2*(1-eps));
            double B = A*eps, C = A*sqrt(1-eps*eps);
            double V = A*sqrt(2*eps*(1+eps)), W = A*sqrt(2*eps*(1-eps));

            double UL_mod = (V/A)*A_UL_sin_inject*sin(kaonp_Phi_h_mc) + (B/A)*A_UL_2sin_inject*sin(2*kaonp_Phi_h_mc);
            double LL_mod = (C/A)*A_LL_0_inject + (W/A)*A_LL_cos_inject*cos(kaonp_Phi_h_mc);
            double LU_mod = (W/A) * A_LU_sin_inject * sin(kaonp_Phi_h_mc);

            //double f = 1.0 + helicity_mc*LU_mod + dilution*fake_Pt*beam_pol*UL_mod + dilution*beam_pol*helicity_mc*fake_Pt*LL_mod;  
            double f = 1.0 + helicity_mc*beam_pol*LU_mod + dilution*fake_Pt*UL_mod + dilution*beam_pol*helicity_mc*fake_Pt*LL_mod;  // beampol
            double u = f_max_inject * rng.Rndm();   
            keep_injected[i] = (u < f);
            // like writing rng.Rndm() < f/f_max_inject
            //kp_fmax.Fill(u);
            //kp_fkeep.Fill(f);
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
                    vec_kaonp_bootw_2d_mc[index_xQ2-1].push_back(bootw);
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
                    vec_kaonp_bootw_2d_zPt_mc[index_zPt-1].push_back(bootw);
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
            int fake_spinstate = (rng.Rndm() < 0.5) ? +1 : -1;   // one time per event
            //int fake_spinstate = (rng.Rndm() < 0.5788) ? +1 : -1;
            helicity = (rng.Rndm() < 0.5) ? +1 : -1;
            double fake_Pt = 0.85 * fake_spinstate;
            fake_Ptarget_vec_reco[i] = fake_Pt;
            double bootw = rng.PoissonD(1.0);

            double eps = kaonp_epsilon_recoMC, y = kaonp_y_recoMC;
            double A = (y*y)/(2*(1-eps));
            double B = A*eps, C = A*sqrt(1-eps*eps);
            double V = A*sqrt(2*eps*(1+eps)), W = A*sqrt(2*eps*(1-eps));

            double UL_mod = (V/A)*A_UL_sin_inject*sin(kaonp_Phi_h_recoMC) + (B/A)*A_UL_2sin_inject*sin(2*kaonp_Phi_h_recoMC);
            double LL_mod = (C/A)*A_LL_0_inject + (W/A)*A_LL_cos_inject*cos(kaonp_Phi_h_recoMC);
            double LU_mod = (W/A) * A_LU_sin_inject * sin(kaonp_Phi_h_recoMC);

            //double f = 1.0 + helicity*LU_mod + dilution*fake_Pt*UL_mod + dilution*helicity*fake_Pt*LL_mod;   
            double f = 1.0 + helicity*beam_pol*LU_mod + dilution*fake_Pt*UL_mod + dilution*beam_pol*helicity*fake_Pt*LL_mod;  // beampol 

            double u = f_max_inject * rng.Rndm();   // different for each event
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
                        vec_kaonp_bootw_2d_all_ID[mc_index_xQ2_all_ID-1].push_back(bootw);
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
                        vec_kaonp_bootw_2d_zPt_all_ID[mc_index_zPt_all_ID-1].push_back(bootw);
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
                    vec_kaonp_bootw_2d_preID[mc_index_xQ2_preCut-1].push_back(bootw);
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
                    vec_kaonp_bootw_2d_zPt_preID[mc_index_zPt_preCut-1].push_back(bootw);
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
                            vec_kaonp_bootw_2d_true[mc_index_xQ2-1].push_back(bootw);
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
                            vec_kaonp_bootw_2d_zPt_true[mc_index_zPt-1].push_back(bootw);
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
                            vec_kaonp_bootw_2d[index_xQ2-1].push_back(bootw);
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
                            vec_kaonp_bootw_2d_zPt[index_zPt-1].push_back(bootw);
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
        /*
        for (int x = 0; x < nbin_xQ2; x++){ 
            ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
            ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_mc[x], vec_kaonp_spin_2d_mc[x], vec_kaonp_y_2d_mc[x],vec_kaonp_eps_2d_mc[x], vec_kaonp_pol_2d_mc[x], vec_helicity_2d_mc[x], vec_kaonp_bootw_2d_mc[x], true, period);}, 7);
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
            //minimizer.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
            minimizer.SetLimitedVariable(6, "Alu", 0.00, 0.001, -0.8, 0.8);
            minimizer.Minimize();
            // LSA
            AUL_sin_2d_mc[x] = minimizer.X()[0];
            AUL_sin_err_2d_mc[x] = minimizer.Errors()[0];
            AUL_2sin_2d_mc[x] = minimizer.X()[1];
            AUL_2sin_err_2d_mc[x] = minimizer.Errors()[1];
            ALL_0_2d_mc[x] = minimizer.X()[2];
            ALL_0_err_2d_mc[x] = minimizer.Errors()[2];
            ALL_cos_2d_mc[x] = minimizer.X()[3];
            ALL_cos_err_2d_mc[x] = minimizer.Errors()[3];
            ALU_sin_2d_mc[x] = minimizer.X()[6];
            ALU_sin_err_2d_mc[x] = minimizer.Errors()[6];
            cout << "bin: " << x <<  "  -->    Aul_sin = " << minimizer.X()[0]  << "    Aul_sin2 = " << minimizer.X()[1]  <<"    All = " << minimizer.X()[2] << "    All_cos = " << minimizer.X()[3]  << "    Auu_cos = " << minimizer.X()[4] << "    Auu_cos2 = " << minimizer.X()[5] <<  "    Alu_sin = " << minimizer.X()[6] << endl;
            
            // true
            ROOT::Minuit2::Minuit2Minimizer minimizer2(ROOT::Minuit2::kMigrad);
            ROOT::Math::Functor MLE2([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_true[x], vec_kaonp_spin_2d_true[x], vec_kaonp_y_2d_true[x],vec_kaonp_eps_2d_true[x], vec_kaonp_pol_2d_true[x], vec_helicity_2d_true[x], vec_kaonp_bootw_2d_true[x], true, period);}, 7);
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
            //minimizer2.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer2.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            ROOT::Math::Functor MLE3([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_all_ID[x], vec_kaonp_spin_2d_all_ID[x], vec_kaonp_y_2d_all_ID[x],vec_kaonp_eps_2d_all_ID[x], vec_kaonp_pol_2d_all_ID[x], vec_helicity_2d_all_ID[x], vec_kaonp_bootw_2d_all_ID[x], true, period);}, 7);
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
            //minimizer3.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer3.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            ROOT::Math::Functor MLE4([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_preID[x], vec_kaonp_spin_2d_preID[x], vec_kaonp_y_2d_preID[x],vec_kaonp_eps_2d_preID[x], vec_kaonp_pol_2d_preID[x], vec_helicity_2d_preID[x], vec_kaonp_bootw_2d_preID[x], true, period);}, 7);
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
            //minimizer4.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer4.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            ROOT::Math::Functor MLE5([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d[x], vec_kaonp_spin_2d[x], vec_kaonp_y_2d[x],vec_kaonp_eps_2d[x], vec_kaonp_pol_2d[x], vec_helicity_2d[x], vec_kaonp_bootw_2d[x], true, period);}, 7);
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
            //minimizer5.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer5.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            ROOT::Math::Functor MLE6([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_mc[x], vec_kaonp_spin_2d_mc[x], vec_kaonp_y_2d_mc[x],vec_kaonp_eps_2d_mc[x], vec_kaonp_pol_2d_mc[x], vec_helicity_2d_mc[x], vec_kaonp_bootw_2d_mc[x], true, period);}, 7);
            minimizer6.SetFunction(MLE6);
            minimizer6.SetMaxFunctionCalls(50000);
            minimizer6.SetMaxIterations(10000);
            minimizer6.SetTolerance(0.001);
            minimizer6.SetStrategy(1);
            minimizer6.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
            minimizer6.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
            minimizer6.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
            minimizer6.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
            minimizer6.SetFixedVariable(4, "Auu", Auu_cos_fixed);
            minimizer6.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
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
            
        }
        
        cout << " -------------------------------------------------------------------------------------------------------------------------------------------- " << endl;
        */
        

        // Asymmetries zPt
        for (int z = 0; z < nbin_zPt; z++){ 

            double r_local_mc, r_local_preID, r_local_true, r_local_all_ID, r_local_reco;
            /*
            double N_pos_mc = std::count(vec_helicity_2d_zPt_mc[z].begin(), vec_helicity_2d_zPt_mc[z].end(), +1);
            double N_neg_mc = std::count(vec_helicity_2d_zPt_mc[z].begin(), vec_helicity_2d_zPt_mc[z].end(), -1);
            double r_local_mc = N_neg_mc / N_pos_mc;

            double N_pos_preID = std::count(vec_helicity_2d_zPt_preID[z].begin(), vec_helicity_2d_zPt_preID[z].end(), +1);
            double N_neg_preID = std::count(vec_helicity_2d_zPt_preID[z].begin(), vec_helicity_2d_zPt_preID[z].end(), -1);
            double r_local_preID = N_neg_preID / N_pos_preID;

            double N_pos_true = std::count(vec_helicity_2d_zPt_true[z].begin(), vec_helicity_2d_zPt_true[z].end(), +1);
            double N_neg_true = std::count(vec_helicity_2d_zPt_true[z].begin(), vec_helicity_2d_zPt_true[z].end(), -1);
            double r_local_true = N_neg_true / N_pos_true;

            double N_pos_all_ID = std::count(vec_helicity_2d_zPt_all_ID[z].begin(), vec_helicity_2d_zPt_all_ID[z].end(), +1);
            double N_neg_all_ID = std::count(vec_helicity_2d_zPt_all_ID[z].begin(), vec_helicity_2d_zPt_all_ID[z].end(), -1);
            double r_local_all_ID = N_neg_all_ID / N_pos_all_ID;

            double N_pos_reco = std::count(vec_helicity_2d_zPt[z].begin(), vec_helicity_2d_zPt[z].end(), +1);
            double N_neg_reco = std::count(vec_helicity_2d_zPt[z].begin(), vec_helicity_2d_zPt[z].end(), -1);
            double r_local_reco = N_neg_reco / N_pos_reco;
            */


            ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
            ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_preID[z], vec_kaonp_spin_2d_zPt_preID[z], vec_kaonp_y_2d_zPt_preID[z],vec_kaonp_eps_2d_zPt_preID[z], vec_kaonp_pol_2d_zPt_preID[z], vec_helicity_2d_zPt_preID[z], vec_kaonp_bootw_2d_zPt_preID[z], r_local_preID, true, period);}, 7);
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
            //minimizer.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            cout << "bin: " << z <<  "  -->  Aul_sin = " << minimizer.X()[0] << "    Aul_sin2 = " << minimizer.X()[1] << "    All = " << minimizer.X()[2] << "    All_cos = " << minimizer.X()[3]  << "    Auu_cos = " << minimizer.X()[4] << "    Auu_cos2 = " << minimizer.X()[5] <<  "    Alu_sin = " << minimizer.X()[6] << endl;

            // true
            ROOT::Minuit2::Minuit2Minimizer minimizer2(ROOT::Minuit2::kMigrad);
            ROOT::Math::Functor MLE2([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_true[z], vec_kaonp_spin_2d_zPt_true[z], vec_kaonp_y_2d_zPt_true[z],vec_kaonp_eps_2d_zPt_true[z], vec_kaonp_pol_2d_zPt_true[z], vec_helicity_2d_zPt_true[z], vec_kaonp_bootw_2d_zPt_true[z], r_local_true, true, period);}, 7);
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
            //minimizer2.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer2.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            ROOT::Math::Functor MLE3([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_all_ID[z], vec_kaonp_spin_2d_zPt_all_ID[z], vec_kaonp_y_2d_zPt_all_ID[z],vec_kaonp_eps_2d_zPt_all_ID[z], vec_kaonp_pol_2d_zPt_all_ID[z], vec_helicity_2d_zPt_all_ID[z], vec_kaonp_bootw_2d_zPt_all_ID[z], r_local_all_ID, true, period);}, 7);
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
            //minimizer3.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer3.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            /*
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
            */
            // reco
            ROOT::Minuit2::Minuit2Minimizer minimizer5(ROOT::Minuit2::kMigrad);
            ROOT::Math::Functor MLE5([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt[z], vec_kaonp_spin_2d_zPt[z], vec_kaonp_y_2d_zPt[z],vec_kaonp_eps_2d_zPt[z], vec_kaonp_pol_2d_zPt[z], vec_helicity_2d_zPt[z], vec_kaonp_bootw_2d_zPt[z], r_local_reco, true, period);}, 7);
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
            //minimizer5.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer5.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
            ROOT::Math::Functor MLE6([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d_zPt_mc[z], vec_kaonp_spin_2d_zPt_mc[z], vec_kaonp_y_2d_zPt_mc[z],vec_kaonp_eps_2d_zPt_mc[z], vec_kaonp_pol_2d_zPt_mc[z], vec_helicity_2d_zPt_mc[z], vec_kaonp_bootw_2d_zPt_mc[z], r_local_mc, true, period);}, 7);
            minimizer6.SetFunction(MLE6);
            minimizer6.SetMaxFunctionCalls(50000);
            minimizer6.SetMaxIterations(10000);
            minimizer6.SetTolerance(0.001);
            minimizer6.SetStrategy(1);
            minimizer6.SetLimitedVariable(0, "Aul_sin", 0.00, 0.001, -0.8, 0.8); 
            minimizer6.SetLimitedVariable(1, "Aul_2sin", 0.00, 0.001, -0.8, 0.8);
            minimizer6.SetLimitedVariable(2, "All", 0.00, 0.001, -0.8, 0.8); 
            minimizer6.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.8, 0.8);
            minimizer6.SetFixedVariable(4, "Auu", Auu_cos_fixed);
            minimizer6.SetFixedVariable(5, "Auu2", Auu_cos2_fixed);
            //minimizer6.SetLimitedVariable(4, "Auu", 0.00, 0.001, -0.8, 0.8);
            //minimizer6.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -0.8, 0.8);
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
        
        for (int x = 0; x < nbin_xQ2; x++){
            // HELICITY
            double N_pos_mc = std::count(vec_helicity_2d_mc[x].begin(), vec_helicity_2d_mc[x].end(), +1);
            double N_neg_mc = std::count(vec_helicity_2d_mc[x].begin(), vec_helicity_2d_mc[x].end(), -1);
            double r_local_mc = N_neg_mc / N_pos_mc;

            double N_pos_preID = std::count(vec_helicity_2d_preID[x].begin(), vec_helicity_2d_preID[x].end(), +1);
            double N_neg_preID = std::count(vec_helicity_2d_preID[x].begin(), vec_helicity_2d_preID[x].end(), -1);
            double r_local_preID = N_neg_preID / N_pos_preID;

            double N_pos_true = std::count(vec_helicity_2d_true[x].begin(), vec_helicity_2d_true[x].end(), +1);
            double N_neg_true = std::count(vec_helicity_2d_true[x].begin(), vec_helicity_2d_true[x].end(), -1);
            double r_local_true = N_neg_true / N_pos_true;

            double N_pos_all_ID = std::count(vec_helicity_2d_all_ID[x].begin(), vec_helicity_2d_all_ID[x].end(), +1);
            double N_neg_all_ID = std::count(vec_helicity_2d_all_ID[x].begin(), vec_helicity_2d_all_ID[x].end(), -1);
            double r_local_all_ID = N_neg_all_ID / N_pos_all_ID;

            double N_pos_reco = std::count(vec_helicity_2d[x].begin(), vec_helicity_2d[x].end(), +1);
            double N_neg_reco = std::count(vec_helicity_2d[x].begin(), vec_helicity_2d[x].end(), -1);
            double r_local_reco = N_neg_reco / N_pos_reco;

            //cout << "bin " << x << ": r_local_mc=" << r_local_mc << "  r_local_preID=" << r_local_preID << "r_local_all_ID=" << r_local_all_ID << "  r_local_true=" << r_local_true << "  r_local_reco=" << r_local_reco << "  r_teorico=0.7277" << endl;

            // SPIN
            double S_pos_mc = std::count(vec_kaonp_spin_2d_mc[x].begin(), vec_kaonp_spin_2d_mc[x].end(), +1);
            double S_neg_mc = std::count(vec_kaonp_spin_2d_mc[x].begin(), vec_kaonp_spin_2d_mc[x].end(), -1);
            double Sr_local_mc = S_neg_mc / S_pos_mc;

            double S_pos_preID = std::count(vec_kaonp_spin_2d_preID[x].begin(), vec_kaonp_spin_2d_preID[x].end(), +1);
            double S_neg_preID = std::count(vec_kaonp_spin_2d_preID[x].begin(), vec_kaonp_spin_2d_preID[x].end(), -1);
            double Sr_local_preID = S_neg_preID / S_pos_preID;

            double S_pos_true = std::count(vec_kaonp_spin_2d_true[x].begin(), vec_kaonp_spin_2d_true[x].end(), +1);
            double S_neg_true = std::count(vec_kaonp_spin_2d_true[x].begin(), vec_kaonp_spin_2d_true[x].end(), -1);
            double Sr_local_true = S_neg_true / S_pos_true;

            double S_pos_all_ID = std::count(vec_kaonp_spin_2d_all_ID[x].begin(), vec_kaonp_spin_2d_all_ID[x].end(), +1);
            double S_neg_all_ID = std::count(vec_kaonp_spin_2d_all_ID[x].begin(), vec_kaonp_spin_2d_all_ID[x].end(), -1);
            double Sr_local_all_ID = S_neg_all_ID / S_pos_all_ID;

            double S_pos_reco = std::count(vec_kaonp_spin_2d[x].begin(), vec_kaonp_spin_2d[x].end(), +1);
            double S_neg_reco = std::count(vec_kaonp_spin_2d[x].begin(), vec_kaonp_spin_2d[x].end(), -1);
            double Sr_local_reco = S_neg_reco / S_pos_reco;

            //cout << "bin " << x << ": Sr_local_mc=" << Sr_local_mc << "  Sr_local_preID=" << Sr_local_preID << "Sr_local_all_ID=" << Sr_local_all_ID << "  Sr_local_true=" << Sr_local_true << "  Sr_local_reco=" << Sr_local_reco << "  Sr_teorico=0.5" << endl;

            //cout << "bin " << x << ": H x S mc=" << Sr_local_mc*r_local_mc<< "  H x S preID=" << Sr_local_preID*r_local_preID << "H x S all_ID=" << Sr_local_all_ID*r_local_all_ID << "  H x S true=" << Sr_local_true*r_local_true << "  H x S reco=" << Sr_local_reco*r_local_reco << "  Sr_teorico=0.727" << endl;
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

            csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                        << "AUL_sinPhi" << ","
                        <<  (AUL_sin_2d_zPt_mc[iz] - AUL_sin_2d_zPt_preID[iz]) << "," << safeSigmaDiff(AUL_sin_err_2d_zPt_mc[iz], AUL_sin_err_2d_zPt_preID[iz]) << ","
                        << (AUL_sin_2d_zPt_preID[iz] - AUL_sin_2d_zPt_true[iz]) << "," << safeSigmaDiff(AUL_sin_err_2d_zPt_preID[iz], AUL_sin_err_2d_zPt_true[iz]) << ","
                        <<  (AUL_sin_2d_zPt_all_ID[iz] - AUL_sin_2d_zPt_true[iz])<< "," << safeSigmaDiff(AUL_sin_err_2d_zPt_all_ID[iz], AUL_sin_err_2d_zPt_true[iz]) << ","
                        <<  (AUL_sin_2d_zPt_true[iz] - AUL_sin_2d_zPt_reco[iz])<< "," << safeSigmaDiff(AUL_sin_err_2d_zPt_true[iz], AUL_sin_err_2d_zPt_reco[iz]) << "\n";
            csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                        << "AUL_sin2Phi" << ","
                        <<  (AUL_2sin_2d_zPt_mc[iz] - AUL_2sin_2d_zPt_preID[iz])<< "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_mc[iz], AUL_2sin_err_2d_zPt_preID[iz]) << ","
                        << (AUL_2sin_2d_zPt_preID[iz] - AUL_2sin_2d_zPt_true[iz]) << "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_preID[iz], AUL_2sin_err_2d_zPt_true[iz]) << ","
                        <<  (AUL_2sin_2d_zPt_all_ID[iz] - AUL_2sin_2d_zPt_true[iz])<< "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_all_ID[iz], AUL_2sin_err_2d_zPt_true[iz]) << ","
                        <<  (AUL_2sin_2d_zPt_true[iz] - AUL_2sin_2d_zPt_reco[iz])<< "," << safeSigmaDiff(AUL_2sin_err_2d_zPt_true[iz], AUL_2sin_err_2d_zPt_reco[iz]) << "\n";
            csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                        << "ALL_const" << ","
                        <<  (ALL_0_2d_zPt_mc[iz] - ALL_0_2d_zPt_preID[iz])<< "," << safeSigmaDiff(ALL_0_err_2d_zPt_mc[iz], ALL_0_err_2d_zPt_preID[iz]) << ","
                        << (ALL_0_2d_zPt_preID[iz] - ALL_0_2d_zPt_true[iz]) << "," << safeSigmaDiff(ALL_0_err_2d_zPt_preID[iz], ALL_0_err_2d_zPt_true[iz]) << ","
                        <<  (ALL_0_2d_zPt_all_ID[iz] - ALL_0_2d_zPt_true[iz])<< "," << safeSigmaDiff(ALL_0_err_2d_zPt_all_ID[iz], ALL_0_err_2d_zPt_true[iz]) << ","
                        <<  (ALL_0_2d_zPt_true[iz] - ALL_0_2d_zPt_reco[iz])<< "," << safeSigmaDiff(ALL_0_err_2d_zPt_true[iz], ALL_0_err_2d_zPt_reco[iz]) << "\n";
            csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                        << "ALL_cosPhi" << ","
                        <<  (ALL_cos_2d_zPt_mc[iz] - ALL_cos_2d_zPt_preID[iz])<< "," << safeSigmaDiff(ALL_cos_err_2d_zPt_mc[iz], ALL_cos_err_2d_zPt_preID[iz]) << ","
                        << (ALL_cos_2d_zPt_preID[iz] - ALL_cos_2d_zPt_true[iz]) << "," << safeSigmaDiff(ALL_cos_err_2d_zPt_preID[iz], ALL_cos_err_2d_zPt_true[iz]) << ","
                        <<  (ALL_cos_2d_zPt_all_ID[iz] - ALL_cos_2d_zPt_true[iz])<< "," << safeSigmaDiff(ALL_cos_err_2d_zPt_all_ID[iz], ALL_cos_err_2d_zPt_true[iz]) << ","
                        <<  (ALL_cos_2d_zPt_true[iz] - ALL_cos_2d_zPt_reco[iz])<< "," << safeSigmaDiff(ALL_cos_err_2d_zPt_true[iz], ALL_cos_err_2d_zPt_reco[iz]) << "\n";
            csvFile_zPt_test << vec_kaonp_y_2d_zPt[iz].size() << "," << "NaN" << "," << iz+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << mean_eps << ","
                        << "ALU_sinPhi" << ","
                        <<  (ALU_sin_2d_zPt_mc[iz] - ALU_sin_2d_zPt_preID[iz])<< "," << safeSigmaDiff(ALU_sin_err_2d_zPt_mc[iz], ALU_sin_err_2d_zPt_preID[iz]) << ","
                        << (ALU_sin_2d_zPt_preID[iz] - ALU_sin_2d_zPt_true[iz]) << "," << safeSigmaDiff(ALU_sin_err_2d_zPt_preID[iz], ALU_sin_err_2d_zPt_true[iz]) << ","
                        <<  (ALU_sin_2d_zPt_all_ID[iz] - ALU_sin_2d_zPt_true[iz])<< "," << safeSigmaDiff(ALU_sin_err_2d_zPt_all_ID[iz], ALU_sin_err_2d_zPt_true[iz]) << ","
                        <<  (ALU_sin_2d_zPt_true[iz] - ALU_sin_2d_zPt_reco[iz])<< "," << safeSigmaDiff(ALU_sin_err_2d_zPt_true[iz], ALU_sin_err_2d_zPt_reco[iz]) << "\n";
        }

        csvFile_zPt_test.close();

        cout << "csv: " << csv_filename_zPt_test << endl;
    }

    //outFile.Write();
    //outFile.Close();

}