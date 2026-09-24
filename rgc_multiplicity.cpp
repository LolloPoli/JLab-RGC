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
//#include "LHAPDF/LHAPDF.h"


using namespace std;
namespace fs = std::filesystem;
gROOT->SetBatch(kTRUE);


// to download the data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rgc/summer22/NH3/multi/kaon_plus/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/sum22_multi_NH3_data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/w/hallb-scshelf2102/clas12/lpolizzi/rgc/output_prova3.root /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/fall22_NH3_data

double MeanVect(const vector<double>& v) {
    if (v.empty()) return 0.0;
    double sum = 0.0;
    for (double x : v) sum += x;
    return sum / v.size();
}

double AUL_loglike(const double* Aul,
                   const vector<double>& phi_h,
                   const vector<double>& spinstate, // +1 or -1 for target orientation
                   const vector<double>& depol,  // D(kin) or 1 if not used
                   const vector<double>& Ptarget, // nominal P_target
                   const vector<double>& sintheta,
                   bool use_spinstate_for_helicity = true)
{
    double logLike = 0.0;
    const double A_sin  = Aul[0];
    const double A_sin2 = Aul[1]; 
    //const double A_sin3 = Aul[2];
    double dilution = 1;
    // rapporto di luminosità osservato (N+ / N−)
    double Nplus = 0.0;
    double Nminus = 0.0;
    for (size_t i = 0; i < Ptarget.size(); ++i) {
        if (Ptarget[i] > 0) Nplus++;
        else               Nminus++;
    }
    double r = Nplus / Nminus;

    for (size_t i = 0; i < phi_h.size(); ++i) {
        double eps = depol[i];
        double D_sin  = sqrt(2*eps*(1+eps));
        double D_2sin = eps;

        // modello di asimmetria
        double modulation = D_sin * A_sin * sin(phi_h[i]) + D_2sin * A_sin2 * sin(2.0*phi_h[i]);

        double pol_factor = Ptarget[i] * dilution; 

        // correzione per luminosità asimmetrica
        double Seff = pol_factor * modulation;
        double numer = r * (1.0 + Seff);
        double denom = r * (1.0 + Seff) + (1.0 - Seff);
        double arg = numer / denom; // probabilità per spin = +1 

        // gestione spin state
        double prob;
        if (use_spinstate_for_helicity) {
            if (spinstate[i] == +1) {
                prob = arg;
            } else {
                prob = 1.0 - arg;
            }
        } else {
            prob = arg;
        }

        if (prob <= 1e-8) return 1e-8;
        logLike += log(prob);
    }

    return -logLike;
}

// OLD FUNCTION WITH LL
/*
double AUL_loglike_withLL(const double* Aul,
                   const vector<double>& phi_h,
                   const vector<double>& spinstate,    // +1 or -1 for target orientation
                   const vector<double>& depol,    // eps per-event
                   const vector<double>& Ptarget,  // P_target per-event (aligned to gamma* frame)
                   const vector<double>& helicity,    // beam helicity per-event (+1/-1)
                   const vector<double>& sintheta,
                   bool use_spinstate_for_helicity = true)
{
    double logLike = 0.0;

    // [0]=A_UL_sin, [1]=A_UL_sin2, [2]=A_UL_sin3,
    // [3]=A_LL_0 (phi-independent), [4]=A_LL_cosphi, [5]=A_LL_cos2phi
    const double A_sin      = Aul[0];
    const double A_sin2     = Aul[1];
    const double A_sin3     = Aul[2];
    const double A_LL_0     = Aul[3];
    const double A_LL_cos1  = Aul[4];
    const double A_LL_cos2  = Aul[5];
    //const double A_UU_cos1  = Aul[6];
    //const double A_UU_cos2  = Aul[7];
    //const double Azim_asy   = Aul[6];
    const double dilution = 1.0;

    // rapporto di luminosità 
    double Nplus = 0.0;
    double Nminus = 0.0;
    for (size_t i = 0; i < spinstate.size(); ++i) {
        if (spinstate[i] > 0) Nplus++;
        else Nminus++;
    }
    double r = Nplus / Nminus;

    for (size_t i = 0; i < phi_h.size(); ++i) {
        double eps = depol[i];
        // depolarization factors
        double D_sin   = sqrt(2.0 * eps * (1.0 + eps));      // per sin phi (UL)
        double D_2sin  = eps;                                // per sin 2phi (UL)

        // UL modulation - Add A_sin2eff if you want to switch to TSA
        double UL_mod = D_sin * A_sin * sin(phi_h[i])
                      + D_2sin * A_sin2 * sin(2.0 * phi_h[i])
                      - A_sin3 * sin(3.0 * phi_h[i]);

        // LL depolarization factors
        double D_LL_0   = sqrt(1.0 - eps*eps);               // phi-independent
        double D_LL_cos = sqrt(2.0 * eps * (1.0 - eps));     // cos phi
        double D_LL_cos2 = 1.0;                              // 

        int lam = helicity[i]; 
        double LL_mod = D_LL_0 * A_LL_0
                      + D_LL_cos * A_LL_cos1 * cos(phi_h[i])
                      - D_LL_cos2 * A_LL_cos2 * cos(2.0 * phi_h[i]);


        double D_UU_1 = sqrt(2*eps*(1+eps));
        double D_UU_2 = eps;
        //double UU_mod = 1 + D_UU_1 * A_UU_cos1 * cos(phi_h[i]) + D_UU_2 * A_UU_cos2 + cos(2*phi_h[i]);

        // totale (LL moltiplicato per lambda)
        double modulation = UL_mod + lam * LL_mod;

        // effective polarization (target * dilution)
        double pol_factor = fabs(Ptarget[i]); 
        double SUL =  dilution * pol_factor * UL_mod;
        double SLL = dilution * pol_factor * lam * LL_mod;
        //double Seff =  (dilution * pol_factor * modulation);
        double Seff = SUL + SLL;

        // CORREZIONE luminosità
        double numer = r * (1.0 + Seff);
        double denom = numer + (1.0 - Seff);
        double arg = numer / denom; // prob 
        double prob;
        if (use_spinstate_for_helicity) {
            // attenzione: qui usi convenzione spinstate -1 -> arg
            if (spinstate[i] > 0) prob = arg;
            else    prob = 1.0 - arg;
        } else {
            prob = arg;
        }

        if (prob <= 1e-6) return 1e-6;
        logLike += log(prob);
    }

    return -logLike; // minimize -logL
}
*/

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

    const double dilution = 0.4;

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
    if(nTot == 0){
        if(Ph <= 3.5){
            if(track_pid == 321) return true;
            else return false;
        }
    }
    else if(nTot > 0){
        if(std::abs(phi) >= 2.6){
            if(Ph >= 1.5){
                // RICH: kaon and CLAS: hadron
                if (rich_pid == 321 && track_pid != 321) {
                    if (RL > 6 || RQ < 0.1 || nTot <= 2) return false;

                    double m = 0.0, q = 0.0;
                    if (rich_chi2 < 50 && Ph <= 3) {
                        m = mA; q = qA;
                    } else if (rich_chi2 < 50 && Ph > 3) {
                        m = mB; q = qB;
                    } else if (rich_chi2 >= 50 && Ph <= 3) {
                        m = mC; q = qC;
                    } else { // rich_chi2 >= 50 && Ph > 3
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
            }
        }
        else if(Ph < 1.5){
            if(track_pid == 321) return true;
            else return false;
        }
        // RICH: hadron and CLAS: kaon
        if (rich_pid != 321 && Ph <= 3) { // since we cannot trust CLAS above 3 GeV
            if(std::abs(phi) >= 2.6){
                if (RL > 6.5 && nTot <= 5) return true; // the oppisite, if the RICH is bad, we want these data, we trust CLAS
                if (RQ < 0.1) return true;

                double m = 0.0, q = 0.0;
                if (rich_chi2 < 50 && Ph <= 3) {
                    m = mA; q = qA;
                } else if (rich_chi2 < 50 && Ph > 3) {
                    m = mB; q = qB;
                } else if (rich_chi2 >= 50 && Ph <= 3) {
                    m = mC; q = qC;
                } else { // rich_chi2 >= 50 && Ph > 3
                    m = mD; q = qD;
                }
                double y_limit = m * RL + q;
                double x_limit = (RQ - q)/m;
                if (RL > x_limit) return true; 

                return false;
            }
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

int FindBin1D(double value, double bins[][2], int nBins){
    
    for(int i = 0; i < nBins; i++){
        if(value >= bins[i][0] && value < bins[i][1])
            return i;
    }
    
    return -1; // fuori range
}

int FindBinDependent(double value, const vector<array<double,2>>& bins){
    for(int i = 0; i < bins.size(); i++){
        if(value >= bins[i][0] && value < bins[i][1])
            return i;
    }

    return -1;
}

double GetDelta(double bins[][2], int i){
    return bins[i][1] - bins[i][0];
}

double GetDeltaDep(const vector<array<double,2>>& bins, int i){
    return bins[i][1] - bins[i][0];
}


double chargesq(int flav){
    if(abs(flav)==2) return 4.0/9.0; // u
    if(abs(flav)==1) return 1.0/9.0; // d
    if(abs(flav)==3) return 1.0/9.0; // s
    return 0.0;
}
/*
double multiplicity_theory(LHAPDF::PDF* pdf,LHAPDF::PDF* ff,double x,double z,double Q2){
    double num = 0.0;
    double den = 0.0;

    vector<int> flavors = {2,3,-2,-1,-3};

    for(int f : flavors){

        double fq = pdf->xfxQ2(f, x, Q2)/x;
        double Dq = ff->xfxQ2(f, z, Q2)/z;

        num += chargesq(f)*fq*Dq;
        den += chargesq(f)*fq;
    }

    if(den > 0) return num/den;
    else return 0;
}

*/



void rgc_multiplicity(const char* period, const char* target) {
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
    double C_LSA_TSA;
    double bin_Q2[][2] = {{1,3}, {3, 5}, {5, 7}, {7,11}};
    double bin_Q2_long[][2] = {{1, 3}, {1,3}, {1,3}, {1,3}, {1,3}, {3, 5}, {3, 5}, {3, 5}, {3, 5}, {3, 5}, {5, 7}, {5, 7}, {7, 11}, {7, 11}};
    double bin_xB[][2] = {{0.0, 0.12}, {0.12, 0.18}, {0.18, 0.24}, {0.24, 0.3}, {0.3, 0.8}, {0.0, 0.24}, {0.24, 0.3}, {0.3, 0.36}, {0.36, 0.44}, {0.44, 0.8}, {0.0, 0.44}, {0.44, 0.8}, {0.0, 0.55}, {0.55, 0.8}};
    vector<vector<array<double,2>>> binning_xB_for_Q2 = {
        // in sequenza i bin di Q2 per i rispettivi bin di xB
        {{0.0, 0.12}, {0.12, 0.18}, {0.18, 0.24}, {0.24, 0.3}, {0.3, 0.8}},
        {{0.0, 0.24}, {0.24, 0.3}, {0.3, 0.36}, {0.36, 0.44}, {0.44, 0.8}},
        {{0.0, 0.44}, {0.44, 0.8}},
        {{0.0, 0.55}, {0.55, 0.8}}
    };
    double bin_Pt[][2] = {{0, 0.25}, {0.25, 0.5}, {0.5, 0.8}, {0.8, 1.4}};
    double bin_Pt_long[][2] = {{0, 0.25}, {0, 0.25}, {0, 0.25}, {0, 0.25}, {0, 0.25}, {0, 0.25}, {0, 0.25}, {0.25, 0.5}, {0.25, 0.5}, {0.25, 0.5}, {0.25, 0.5}, {0.25, 0.5}, {0.25, 0.5}, {0.25, 0.5}, {0.5, 0.8}, {0.5, 0.8}, {0.5, 0.8}, {0.5, 0.8}, {0.5, 0.8}, {0.5, 0.8}, {0.5, 0.8}, {0.8, 1.4}, {0.8, 1.4}, {0.8, 1.4}, {0.8, 1.4}};
    double bin_z[][2] = {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}, {0.2, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 1.0}};
    vector<vector<array<double,2>>> binning_z_for_Pt = {
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8}, {0.8, 1.0}},
        {{0.2, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 1.0}}
    };

    string inputDir = string(period) + "_" + string(target) + "_data";
    //vector<string> rootFiles = {"fall22_NH3_data", "sum22_NH3_data"};

    TChain chainKaonP("Kaon+");
    TChain chainElectron("Electron");
    int fileCount = 0;
    
    for (const auto &entry : fs::directory_iterator(inputDir)) {
        if (entry.path().extension() == ".root") {
            string filePath = entry.path().string();
            chainElectron.Add(Form("%s/Electron", filePath.c_str()));
            chainKaonP.Add(Form("%s/Kaon+", filePath.c_str()));


            fileCount++;
        }
    }
    
    /*
    for (const auto& dir : rootFiles) {

        if (!fs::exists(dir)) {
            cerr << "Directory non trovata: " << dir << endl;
            continue;
        }

        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.path().extension() == ".root") {

                string filePath = entry.path().string();
                chainKaonP.Add(Form("%s/Kaon+", filePath.c_str()));
                fileCount++;
            }
        }
    }
    */
    if (fileCount == 0) {
        cerr << "Nessun file .root trovato nelle directory di input" << endl;
        return;
    }


    // creo un output root 
    TString outputFile = Form("plot_%s_1D_rgc_%s_kaonp.root", target, period);
    double torus = -1;
    TString csv_filename = Form("output_RGC_%s_asymmetries_%s.csv", target, period);
    std::ofstream csvFile(csv_filename.Data());
    csvFile << "bin, bin_xB, bin_zPt, mean_xB, mean_Q2 , mean_z, mean_PhT, Multi, Multi_err\n";

    TFile outFile(outputFile, "RECREATE");  // File di output ROOT
    TTree treeKaonP("Kaon+", "");

    TDirectory* dir_xQ2 = outFile.mkdir("Binning xB-Q2");
    TDirectory* dir_zPt = outFile.mkdir("Binning z-Pt");
    TDirectory* dir_aul_sin = outFile.mkdir("AUL sin(Phi)");
    TDirectory* dir_aul_2sin = outFile.mkdir("AUL sin(2Phi)");
    TDirectory* dir_error = outFile.mkdir("Error");
    TDirectory* dir_multi = outFile.mkdir("Multiplicity");

    // To save all the variables
    // Electron
    double electron_xB, electron_y, electron_theta, electron_phi, electron_chi2;
    chainElectron.SetBranchAddress("el_px", &electron_px);
    chainElectron.SetBranchAddress("el_py", &electron_py);
    chainElectron.SetBranchAddress("el_pz", &electron_pz);
    chainElectron.SetBranchAddress("el_mom", &electron_mom);
    chainElectron.SetBranchAddress("el_xB", &electron_xB);
    chainElectron.SetBranchAddress("el_Q2", &electron_Q2);
    chainElectron.SetBranchAddress("el_W", &electron_W);
    chainElectron.SetBranchAddress("el_theta", &electron_theta);
    chainElectron.SetBranchAddress("el_phi", &electron_phi);
    chainElectron.SetBranchAddress("el_y", &electron_y);
    chainElectron.SetBranchAddress("el_chi2", &electron_chi2);
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
    chainKaonP.SetBranchAddress("kaon_Pt", &kaonp_PhT);
    chainKaonP.SetBranchAddress("Pt_over_zQ", &kaonp_Pt_zQ);
    chainKaonP.SetBranchAddress("kaon_phi_lab", &kaonp_Phi);
    chainKaonP.SetBranchAddress("kaon_theta", &kaonp_Theta);
    chainKaonP.SetBranchAddress("eta", &kaonp_eta);
    chainKaonP.SetBranchAddress("kaon_phi_h", &kaonp_Phi_h);
    chainKaonP.SetBranchAddress("Mx", &kaonp_Mx);
    chainKaonP.SetBranchAddress("clas_chi2", &kaonp_chi2pid);
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
    treeKaonP.Branch("kaon_Pt", &kaonp_PhT, "kaon_Pt/D");
    treeKaonP.Branch("Pt_over_zQ", &kaonp_Pt_zQ, "Pt_over_zQ/D");
    treeKaonP.Branch("kaon_phi_lab", &kaonp_Phi, "kaon_phi_lab/D");
    treeKaonP.Branch("kaon_theta", &kaonp_Theta, "kaon_theta/D");
    treeKaonP.Branch("eta", &kaonp_eta, "eta/D");
    treeKaonP.Branch("kaon_phi_h", &kaonp_Phi_h, "kaon_phi_h/D");
    treeKaonP.Branch("Mx", &kaonp_Mx, "Mx/D");
    treeKaonP.Branch("clas_chi2", &kaonp_chi2pid, "clas_chi2/D");
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
    double nbin = 300;
    // Mom
    TH1D kp_evnt_chi2 ("_evnt_chi2", "#chi^{2} EventBuilder PID | only EventBuilder | 1.2 < Mom < 8 GeV ; #chi^{2}; count", 300, -8, 8);
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
    TH2D kp_MomVsMx ("_MomVsMx", "correlation Mom vs M | K+ |; Mom [GeV]; M_{x} [GeV]", nbin, 1, 8, nbin, 0.6, 3.5);
    TH2D kp_MomVsBeta ("_MomVsBeta", "correlation Mom vs #beta | K+ |; Mom [GeV]; #beta", nbin, 1, 8, nbin, 0.88, 1.02);
    TH2D kp_MomVsCh ("_MomVsCh", "correlation Mom vs #theta_{ch} | K+ |; Mom [GeV]; #theta_{ch} [rad]", nbin, 1, 8, nbin, 0.0, 0.35);
    TH2D kp_MomVsMass_Rich ("_MomVsMass_RICH", "correlation Mom vs Mass | K+ |; Mom [GeV]; Mass [GeV]", nbin, 1, 8, nbin, -0.2, 1.2);
    // Q2
    TH2D kp_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; Q^{2} [GeV^{2}]", nbin, 0, 0.8, nbin, 1, 11);
    TH2D kp_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; Q^{2} [GeV^{2}]", nbin, 0, 0.6, nbin, 1, 10);
    TH2D kp_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  K+  | with EventBuilder + RICH ; Mom [GeV]; Q^{2} [GeV^{2}]", nbin, 0.9, 6, nbin, 1, 10);
    TH2D kp_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  K+  | with EventBuilder + RICH ; P_{hT} [GeV]; Q^{2} [GeV^{2}]", nbin, 0, 1.2, nbin, 1, 10);
    TH2D kp_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs z  |  K+  | with EventBuilder + RICH ; z; Q^{2} [GeV^{2}]", nbin, 0.2, 1.0, nbin, 1, 10);
    TH2D kp_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  K+  | with EventBuilder + RICH ; y; Q^{2} [GeV^{2}]", nbin, 0.2, 0.8, nbin, 0.9, 10);
    TH2D kp_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  K+  | with EventBuilder + RICH ; Eta; Q^{2} [GeV^{2}]", nbin, 1.5, 3.0, nbin, 1, 10);
    TH2D kp_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 1, 10);
    // PhT
    TH2D kp_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs z  |  K+  | with EventBuilder + RICH ; z; P_{hT} [GeV]", nbin, 0.2, 1, nbin, 0, 1.4);
    TH2D kp_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; P_{hT} [GeV]", nbin, 0, 0.8, nbin, 0, 1.2);
    TH2D kp_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  K+  | with EventBuilder + RICH ; Eta; P_{hT} [GeV]", nbin, 1.5, 3.0, nbin, 0, 1.2);
    TH2D kp_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; P_{hT} [GeV]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0, 1.2);
    // Z
    TH2D kp_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; z", nbin, 0, 0.8, nbin, 0.2, 1.0);
    TH2D kp_zVsXf ("_zVsXf", "Correlation Z vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; z", nbin, 0, 0.6, nbin, 0.2, 1.0);
    TH2D kp_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  K+  | with EventBuilder + RICH ; Eta; z", nbin, 1.5, 3.0, nbin, 0.2, 1.0);
    TH2D kp_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; z", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0.2, 1.0);
    //
    TH2D kp_xBvsY ("_xBvsY", "Correlation y vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; y", nbin, 0, 0.8, nbin, 0.2, 0.8);
    TH2D kp_xBvsEta ("_xBvsEta", "Correlation #eta vs x_{B}  |  K+  | with EventBuilder + RICH ; #eta; x_{B}", nbin, 1.5, 3.0, nbin, 0.0, 0.8);
    TH2D kp_xBvsSinTheta ("_xBvsSinTheta", "Correlation sin(#theta_{#gamma}) vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; sin#theta", nbin, 0, 0.8, nbin, 0.0, 0.6);
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
    vector<vector<double>> vec_electron_y_2d(nbin_xQ2);
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
    vector<vector<vector<double>>> vec_electron_y_4d(nbin_xQ2, vector<vector<double>> (nbin_zPt));
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
    vector<double> AUL_sin_2d_zPt(nbin_zPt);
    vector<double> AUL_sin_err_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_2d_zPt(nbin_zPt);
    vector<double> AUL_2sin_err_2d_zPt(nbin_zPt);
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
    vector<vector<TH1D*>> hist_Phi_h_plus_4D(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    vector<vector<TH1D*>> hist_Phi_h_minus_4D(nbin_xQ2, vector<TH1D*> (nbin_zPt));
    for (int ix = 0; ix < nbin_xQ2; ix++){
        hist_AUL_sin_xQ2_zPt[ix] = new TH2D(Form("hist_AUL_sin_xQ2_%d", ix+1), Form("AUL_sinPhi error bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
        hist_AUL_sin2_xQ2_zPt[ix] = new TH2D(Form("hist_AUL_sin2_xQ2_%d", ix+1), Form("AUL_sin2Phi error bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
        for (int iz = 0; iz < nbin_zPt; iz++){
            hist_Phi_h_plus_4D[ix][iz] = new TH1D(Form("hist_Phi_h_plus_4D_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | spin +| pion+; #Phi_{h} [Rad]; count", ix+1, iz+1), 40, -TMath::Pi(), TMath::Pi());
            hist_Phi_h_minus_4D[ix][iz] = new TH1D(Form("hist_Phi_h_minus_4D_%d_%d", ix+1, iz+1), Form("#Phi_{h} for bin (%d, x_{B}-Q^{2}) & (%d, z-P_{hT}) | spin - | pion+; #Phi_{h} [Rad]; count", ix+1, iz+1), 40, -TMath::Pi(), TMath::Pi());
        }
    }
    


    //

    // ORA RIEMPI I GRAFICI
    // Kaon+
    Long64_t nEntries_el = chainElectron.GetEntries();
    for (Long64_t i = 0; i < nEntries_el; i++) {
        chainElectron.GetEntry(i);
        //if(electron_track_pid != -11) continue;
        if(std::fabs(electron_chi2) > 3) continue;
        if(electron_y <= 0.2) continue;
        double index_xQ2 = getBinIndex_xQ2(electron_xB, electron_Q2); 
        if(index_xQ2 >= 0) {
            vec_electron_y_2d[index_xQ2-1].push_back(electron_y);
            for (int i = 0; i < nbin_zPt; i++){
                vec_electron_y_4d[index_xQ2-1][i].push_back(electron_y);
            }
        }
    }
    Long64_t nEntries_kp = chainKaonP.GetEntries();
    for (Long64_t i = 0; i < nEntries_kp; i++) {
        // RICH cut
        //if (electron_ass != 321) continue;
        chainKaonP.GetEntry(i);
        if (kaonp_track_pid == -11) continue;
        if (kaonp_y <= 0.2) continue;
        //if (std::abs(kaonp_Pol) < 0.1) continue;
        //
        //if(kaonp_Mx < 2 || kaonp_Mx > 2.5) continue;
        if(std::abs(kaonp_chi2pid) < 3){
            kp_Mx_tof.Fill(kaonp_Mx);
            // passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph)
            if(passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph, kaonp_Phi)){
                double ph2 = kaonp_Ph*kaonp_Ph;
                double beta2 = kaonp_beta*kaonp_beta;
                kaonp_m = sqrt((ph2)*(1-beta2)/beta2);
                double index_xQ2 = getBinIndex_xQ2(kaonp_xB, kaonp_Q2); 
                double index_zPt = getBinIndex_zPt(kaonp_z, kaonp_PhT);
                proton_spin = +1;
                if(kaonp_Pol < 0) proton_spin = -1;
                kaonp_sintheta = kaonp_gamma*sqrt((1-kaonp_y-0.25*kaonp_y*kaonp_y*kaonp_gamma*kaonp_gamma)/(1+kaonp_y*kaonp_y));
                C_LSA_TSA = kaonp_sintheta * (sqrt(2*kaonp_epsilon*(1+kaonp_epsilon))/kaonp_epsilon);
                //if (index_xQ2 != 15) continue;
                //
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
                    vec_kaonp_spin_2d_zPt[index_zPt-1].push_back(proton_spin);
                    vec_kaonp_sintheta_2d_zPt[index_zPt-1].push_back(kaonp_sintheta);
                    if(proton_spin == 1) hist_Phi_h_plus_binzPt[index_zPt-1]->Fill(kaonp_Phi_h);
                    else if(proton_spin == -1) hist_Phi_h_minus_binzPt[index_zPt-1]->Fill(kaonp_Phi_h);
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
                        if(proton_spin == 1) hist_Phi_h_plus_4D[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h);
                        else if(proton_spin == -1) hist_Phi_h_minus_4D[index_xQ2-1][index_zPt-1]->Fill(kaonp_Phi_h);   
                    }
                }
                //
                kp_evnt_chi2.Fill(kaonp_chi2pid);
                kp_m.Fill(kaonp_bestMass);
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
                kp_MomVsTheta.Fill(kaonp_Ph, kaonp_Theta);
                kp_MomVsXb.Fill(kaonp_xB, kaonp_Ph);
                kp_MomVsXf.Fill(kaonp_xF, kaonp_Ph);
                kp_MomVsY.Fill(kaonp_y, kaonp_Ph);
                kp_MomVsZ.Fill(kaonp_z, kaonp_Ph);
                kp_MomVsBeta.Fill(kaonp_Ph, kaonp_beta);
                kp_MomVsCh.Fill(kaonp_Ph, kaonp_rich_ch);
                // Q2
                kp_Q2VsEta.Fill(kaonp_eta, kaonp_Q2);
                kp_Q2VsMom.Fill(kaonp_Ph, kaonp_Q2);
                kp_Q2VsPhi_h.Fill(kaonp_Phi_h, kaonp_Q2);
                kp_Q2VsPhT.Fill(kaonp_PhT, kaonp_Q2);
                kp_Q2VsXb.Fill(kaonp_xB, kaonp_Q2);
                kp_Q2VsXf.Fill(kaonp_xF, kaonp_Q2);
                kp_Q2VsY.Fill(kaonp_y, kaonp_Q2);
                kp_Q2VsZ.Fill(kaonp_z, kaonp_Q2);
                // PhT
                kp_PhTvsEta.Fill(kaonp_eta, kaonp_PhT);
                kp_PhTvsXb.Fill(kaonp_xB, kaonp_PhT);
                kp_PhTvsZ.Fill(kaonp_z, kaonp_PhT);
                kp_PhTvsPhi_h.Fill(kaonp_Phi_h, kaonp_PhT);
                // z
                kp_zVsEta.Fill(kaonp_eta, kaonp_z);
                kp_zVsPhi_h.Fill(kaonp_Phi_h, kaonp_z);
                kp_zVsXb.Fill(kaonp_xB, kaonp_z);
                kp_zVsXf.Fill(kaonp_xF, kaonp_z);
                //
                kp_xBvsY.Fill(kaonp_xB, kaonp_y);
                kp_xBvsEta.Fill(kaonp_eta, kaonp_xB);
                kp_xBvsSinTheta.Fill(kaonp_xB, kaonp_sintheta);
                kp_ThetaVsPhi_h.Fill(kaonp_Phi_h, kaonp_Theta);
                kp_ThetaVsPhi_Lab.Fill(kaonp_Phi, kaonp_Theta);
                kp_PxVsPy.Fill(kaonp_px, kaonp_py);
                el_PxVsPy.Fill(electron_px, electron_py);
                treeKaonP.Fill();

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
    
    

    
    //treeKaonP.Write();


    TH1D* hist_ratio_phi = (TH1D*)kp_phi_plus.Clone("kp_phi_ratio");
    hist_ratio_phi->Scale(1.0/kp_phi_plus.Integral());
    TH1D* hist_ratio_phi_minus = (TH1D*)kp_phi_minus.Clone("kp_phi_minus_ratio");
    hist_ratio_phi_minus->Scale(1.0/kp_phi_minus.Integral());
    hist_ratio_phi->Divide(hist_ratio_phi_minus);
    TCanvas* c_phi_ratio = new TCanvas("c_phi_ratio","#phi_{h} ratio",800,600);
    hist_ratio_phi->SetTitle("K^{+} #Phi_{h} + / #Phi_{h} - ; #Phi_{h} [Rad]; Ratio");
    hist_ratio_phi->SetMarkerStyle(20);
    hist_ratio_phi->SetMarkerColor(kAzure-4);
    hist_ratio_phi->Draw("E1");
    TF1* fit_ratio_phi = new TF1("fit_ratio_phi","[0] + [1]*sin(x) + [2]*sin(2*x) + [3]*sin(3*x) + [4]*cos(x) + [5]*cos(2*x)", -TMath::Pi(), TMath::Pi());
    fit_ratio_phi->SetParameters(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    hist_ratio_phi->Fit("fit_ratio_phi","QER");
    fit_ratio_phi->SetLineColor(kRed);
    c_phi_ratio->Write();


    vector<PhiRatioFitResult> results_xQ2(nbin_xQ2);
    vector<PhiRatioFitResult> results_zPt(nbin_zPt);
    vector<vector<PhiRatioFitResult>> results_4d(nbin_xQ2, vector<PhiRatioFitResult> (nbin_zPt));
    for (int i = 0; i < nbin_xQ2; i++) {
        if(hist_Phi_h_plus_binxQ2[i]->GetEntries() < 20){
            cout << "2D: kipping xQ2 bin " << i+1 << " due to insufficient s=+1 entries: " << hist_Phi_h_plus_binxQ2[i]->GetEntries() << endl;
            continue;
        } if (hist_Phi_h_minus_binxQ2[i]->GetEntries() < 20){
            cout << "2D: kipping xQ2 bin " << i+1 << " due to insufficient s=-1 entries: " << hist_Phi_h_minus_binxQ2[i]->GetEntries() << endl;
            continue;
        }
        results_xQ2[i] = FitPhiRatio(hist_Phi_h_plus_binxQ2[i],hist_Phi_h_minus_binxQ2[i],Form("phi_ratio_xQ2_bin%d", i+1));  
        AUL_sin_sys_2d[i] = results_xQ2[i].p1;
        AUL_sin_sys_err_2d[i] = results_xQ2[i].e1;
        AUL_2sin_sys_2d[i] = results_xQ2[i].p2;
        AUL_2sin_sys_err_2d[i] = results_xQ2[i].e2;
        for(int z = 0; z < nbin_zPt; z++){
            if (hist_Phi_h_plus_4D[i][z]->GetEntries() < 20){
                cout << "4D: kipping xQ2 bin " << i+1 << ", zPt bin " << z+1 << " due to insufficient s=+1 entries: " << hist_Phi_h_plus_4D[i][z]->GetEntries() << endl;
                continue;
            } if (hist_Phi_h_minus_4D[i][z]->GetEntries() < 20){
                cout << "4D: kipping xQ2 bin " << i+1 << ", zPt bin " << z+1 << " due to insufficient s=-1 entries: " << hist_Phi_h_minus_4D[i][z]->GetEntries() << endl;
                continue;
            }
            results_4d[i][z] = FitPhiRatio(hist_Phi_h_plus_4D[i][z], hist_Phi_h_minus_4D[i][z], Form("phi_ratio_4d_xQ2bin%d_zPtbin%d", i+1, z+1));
            if(results_4d[i][z].p1 >= 1 || results_4d[i][z].p2 >=1) continue;
            AUL_sin_sys_4d[i][z] = results_4d[i][z].p1;
            AUL_sin_sys_err_4d[i][z] = results_4d[i][z].e1;
            AUL_2sin_sys_4d[i][z] = results_4d[i][z].p2;
            AUL_2sin_sys_err_4d[i][z] = results_4d[i][z].e2;
        }
    }
    for (int i = 0; i < nbin_zPt; i++) {
        if(hist_Phi_h_plus_binzPt[i]->GetEntries() < 20){
            cout << "2D: skipping zPt bin " << i+1 << " due to insufficient s=+1 entries" << endl;
            continue;
        } if (hist_Phi_h_minus_binzPt[i]->GetEntries() < 20){
            cout << "2D: skipping zPt bin " << i+1 << " due to insufficient s=-1 entries" << endl;
            continue;
        }
        results_zPt[i] = FitPhiRatio(hist_Phi_h_plus_binzPt[i],hist_Phi_h_minus_binzPt[i],Form("phi_ratio_zPt_bin%d", i+1)); 
        AUL_sin_sys_2d_zPt[i] = results_zPt[i].p1;
        AUL_sin_sys_err_2d_zPt[i] = results_zPt[i].e1;
        AUL_2sin_sys_2d_zPt[i] = results_zPt[i].p2; 
        AUL_2sin_sys_err_2d_zPt[i] = results_zPt[i].e2;
    }

    TGraphErrors* g_p1_xQ2 = new TGraphErrors(nbin_xQ2);
    TGraphErrors* g_p2_xQ2 = new TGraphErrors(nbin_xQ2);
    TGraphErrors* g_p1_zPt = new TGraphErrors(nbin_zPt);
    TGraphErrors* g_p2_zPt = new TGraphErrors(nbin_zPt);
    TGraphErrors* g_p1_4d = new TGraphErrors(nbin_xQ2*nbin_zPt);
    TGraphErrors* g_p2_4d = new TGraphErrors(nbin_xQ2*nbin_zPt);

    for (int i = 0; i < nbin_xQ2; i++) {
        g_p1_xQ2->SetPoint(i, i+1, std::fabs(AUL_sin_sys_2d[i]/2));
        g_p1_xQ2->SetPointError(i, 0.0, AUL_sin_sys_err_2d[i]/2);
        g_p2_xQ2->SetPoint(i, i+1, std::fabs(AUL_2sin_sys_2d[i]/2));
        g_p2_xQ2->SetPointError(i, 0.0, AUL_2sin_sys_err_2d[i]/2);
        for (int z = 0; z < nbin_zPt; z++) {
            double sum_xB = 0.0;
            if (vec_kaonp_xB_4d[i][z].empty()) continue;
            for (double val_xB : vec_kaonp_xB_4d[i][z]) sum_xB += val_xB;
            double mean_xB = sum_xB / vec_kaonp_xB_4d[i][z].size();
            int index = i * nbin_zPt + z;
            g_p1_4d->SetPoint(index, mean_xB, std::fabs(AUL_sin_sys_4d[i][z]/2));
            g_p1_4d->SetPointError(index, 0.0, 0);
            g_p2_4d->SetPoint(index, mean_xB, std::fabs(AUL_2sin_sys_4d[i][z]/2));
            g_p2_4d->SetPointError(index, 0.0, 0);
        }
    }
    for (int i = 0; i < nbin_zPt; i++) {
        g_p1_zPt->SetPoint(i, i+1, std::fabs(AUL_sin_sys_2d_zPt[i]/2));
        g_p1_zPt->SetPointError(i, 0.0, AUL_sin_sys_err_2d_zPt[i]/2);
        g_p2_zPt->SetPoint(i, i+1, std::fabs(AUL_2sin_sys_2d_zPt[i]/2));
        g_p2_zPt->SetPointError(i, 0.0, AUL_2sin_sys_err_2d_zPt[i]/2);
    }

    TCanvas* c_p1_xQ2 = new TCanvas("c_p1_sinPhi_xQ2","sinphi systematic vs xQ2 bin",800,600);
    g_p1_xQ2->SetTitle("Residual sin(#Phi) in R(#Phi_{h}); x_{B}-Q^{2} bin; p_{1}");
    g_p1_xQ2->SetMarkerStyle(20);
    g_p1_xQ2->SetMarkerColor(kRed+1), g_p1_xQ2->SetLineColor(kRed+1);
    g_p1_xQ2->Draw("AP");
    g_p1_xQ2->GetYaxis()->SetRangeUser(-0.01, 0.05);
    TLine* l0 = new TLine(0.0, 0.0, nbin_xQ2+1, 0.0);
    l0->SetLineStyle(2);
    l0->SetLineColor(kGray+1);
    l0->Draw();
    c_p1_xQ2->Write();

    TCanvas* c_p2_xQ2 = new TCanvas("c_p2_sin2Phi_xQ2","sin2phi systematic vs xQ2 bin",800,600);
    g_p2_xQ2->SetTitle("Residual sin(2#Phi) in R(#Phi_{h}); x_{B}-Q^{2} bin; p_{2}");
    g_p2_xQ2->SetMarkerStyle(20);
    g_p2_xQ2->SetMarkerColor(kBlue+1), g_p2_xQ2->SetLineColor(kBlue+1);
    g_p2_xQ2->Draw("AP");
    g_p2_xQ2->GetYaxis()->SetRangeUser(-0.01, 0.05);
    l0->Draw();
    c_p2_xQ2->Write();

    TCanvas* c_p1_zPt = new TCanvas("c_p1_sinPhi_zPt","sinphi systematic vs zPt bin",800,600);
    g_p1_zPt->SetTitle("Residual sin(#Phi) in R(#Phi_{h}); z-P_{hT} bin; p_{1}");
    g_p1_zPt->SetMarkerStyle(20);
    g_p1_zPt->SetMarkerColor(kRed+1), g_p1_zPt->SetLineColor(kRed+1);
    g_p1_zPt->Draw("AP");
    g_p1_zPt->GetYaxis()->SetRangeUser(-0.01, 0.06);
    TLine* l1 = new TLine(0.0, 0.0, nbin_zPt+1, 0.0);
    l1->SetLineStyle(2);
    l1->SetLineColor(kGray+1);
    l1->Draw();
    c_p1_zPt->Write();  
    TCanvas* c_p2_zPt = new TCanvas("c_p2_sin2Phi_zPt","sin2phi systematic vs zPt bin",800,600);
    g_p2_zPt->SetTitle("Residual sin(2#Phi) in R(#Phi_{h}); z-P_{hT} bin; p_{2}");
    g_p2_zPt->SetMarkerStyle(20);
    g_p2_zPt->SetMarkerColor(kBlue+1), g_p2_zPt->SetLineColor(kBlue+1);
    g_p2_zPt->Draw("AP");
    g_p2_zPt->GetYaxis()->SetRangeUser(-0.01, 0.06);
    l1->Draw();
    c_p2_zPt->Write();

    TCanvas* c_p1_4d = new TCanvas("c_p1_sinPhi_4d","sinphi systematic vs 4d bin",800,600);
    g_p1_4d->SetTitle("Residual sin(#Phi) in R(#Phi_{h}); <x_{B}>; p_{1}");
    g_p1_4d->SetMarkerStyle(20);
    g_p1_4d->SetMarkerColor(kRed+1), g_p1_4d->SetLineColor(kRed+1);
    g_p1_4d->Draw("AP");
    g_p1_4d->GetYaxis()->SetRangeUser(-0.02, 0.2);
    c_p1_4d->Write();   
    TCanvas* c_p2_4d = new TCanvas("c_p2_sin2Phi_4d","sin2phi systematic vs 4d bin",800,600);
    g_p2_4d->SetTitle("Residual sin(2#Phi) in R(#Phi_{h}); <x_{B}>; p_{2}");
    g_p2_4d->SetMarkerStyle(20);
    g_p2_4d->SetMarkerColor(kBlue+1), g_p2_4d->SetLineColor(kBlue+1);
    g_p2_4d->Draw("AP");
    g_p2_4d->GetYaxis()->SetRangeUser(-0.02, 0.2);
    c_p2_4d->Write();


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


    //auto pdf = LHAPDF::mkPDF("NNPDF40_nlo_as_01180",0);
    //auto ff  = LHAPDF::mkPDF("JAM20-SIDIS_FF_kaon_nlo",0);

    // MULTIPLICITY ANALYSIS
    TGraphErrors* g_multi_xQ2 = new TGraphErrors(nbin_xQ2);
    TGraphErrors* g_multi_xQ2_TH = new TGraphErrors(nbin_xQ2);
    TGraphErrors* g_multi_xQ2zPt = new TGraphErrors(nbin_zPt*nbin_xQ2);
    TGraphErrors* g_multi_xQ2zPt_TH = new TGraphErrors(nbin_zPt*nbin_xQ2);
    for (int i = 0; i < nbin_xQ2; i++) {
        // mean z
        double sum_z = 0.0, sum_xB = 0.0, sum_Q2 = 0.0;
        for(double val_z : vec_kaonp_z_2d[i]) sum_z += val_z;
        double mean_z = sum_z / vec_kaonp_z_2d[i].size();
        // mean x
        for(double val_xB : vec_kaonp_xB_2d[i]) sum_xB += val_xB;
        double mean_xB = sum_xB / vec_kaonp_xB_2d[i].size();
        // mean Q2
        for(double val_Q2 : vec_kaonp_Q2_2d[i]) sum_Q2 += val_Q2;
        double mean_Q2 = sum_Q2 / vec_kaonp_Q2_2d[i].size();
        double delta_Q2 = bin_Q2_long[i][1] - bin_Q2_long[i][0];
        double delta_xB = bin_xB[i][1] - bin_xB[i][0];
        double num_events = vec_kaonp_y_2d[i].size();
        double denom_events = vec_electron_y_2d[i].size()*delta_xB*delta_Q2;
        double multiplicity = num_events / denom_events;
        double err = 0.0;
        double x_centre = 0.5*(bin_xB[i][0] + bin_xB[i][1]);
        double Q2_centre = 0.5*(bin_Q2_long[i][0] + bin_Q2_long[i][1]);
        //double M_th_2d = multiplicity_theory(pdf, ff, mean_xB, mean_z, mean_Q2);
        //cout << "xQ2 bin " << i+1 << " Nh: " << num_events << "   NDIS: "  << vec_electron_y_2d[i].size() << endl;
        for (int j = 0; j < nbin_zPt; j++) {
            double z_centre = 0.5*(bin_z[j][0] + bin_z[j][1]);
            double delta_z = bin_z[j][1] - bin_z[j][0];
            double delta_Pt = bin_Pt_long[j][1] - bin_Pt_long[j][0];
            double multiplicity_4d = vec_kaonp_y_4d[i][j].size() / (denom_events*delta_z*delta_Pt);
            g_multi_xQ2zPt->SetPoint(i*nbin_zPt + j, i+1 + j*0.04, multiplicity_4d);
            //if (multiplicity_4d < 0.01) cout << "entries : " << vec_kaonp_y_4d[i][j].size() << "  high multiplicity - xB-Q2: " << i << "  z-Pt: " << j << "  multiplicity: " << multiplicity_4d << endl; 
            double err_4d = 0.0;
            //double M_th = multiplicity_theory(pdf, ff, x_centre, z_centre, Q2_centre);
            g_multi_xQ2zPt->SetPointError(i*nbin_zPt + j, 0.0, 0.0);
        }
        if (num_events > 0) err = 1/sqrt(num_events);
        g_multi_xQ2->SetPoint(i, i+1, multiplicity);
        g_multi_xQ2->SetPointError(i, 0.0, err);
        //g_multi_xQ2_TH->SetPoint(i, i+1, M_th_2d);
    }

    TMultiGraph *mg_multi2D = new TMultiGraph();
    TCanvas* c_multip = new TCanvas("c_multiplicity","sin2phi systematic vs xQ2 bin",800,600);
    g_multi_xQ2->SetTitle("Multiplicity in x_{B}-Q^{2} bins; x_{B}-Q^{2} bin; Multiplicity: N^{K^{+}}/(N_{DIS}*#Delta x_{B}*#Delta Q^{2})");
    g_multi_xQ2->SetMarkerStyle(20), g_multi_xQ2_TH->SetMarkerStyle(21);
    g_multi_xQ2->SetMarkerColor(kBlue+1), g_multi_xQ2->SetLineColor(kBlue+1);
    g_multi_xQ2_TH->SetMarkerColor(kRed+1), g_multi_xQ2_TH->SetLineColor(kRed+1);
    g_multi_xQ2->Draw("AP");
    g_multi_xQ2->GetYaxis()->SetRangeUser(-0.02, 0.2);
    //mg_multi2D->Add(g_multi_xQ2, "P");
    //mg_multi2D->Add(g_multi_xQ2_TH, "P");
    //mg_multi2D->Draw("AP");
    TLegend* leg_multi = new TLegend(0.75,0.75,0.88,0.85);
    leg_multi->AddEntry(g_multi_xQ2, "RGC", "p");
    leg_multi->AddEntry(g_multi_xQ2_TH, "NNPDF40+JAM20", "p");
    //leg_multi->Draw();
    c_multip->Update();
    c_multip->Write();

    TCanvas* c_multi_4d = new TCanvas("c_multiplicity_4d","Multiplicity in 4D bins",800,600);
    g_multi_xQ2zPt->SetTitle("Multiplicity in 4D bins; x_{B}-Q^{2} bin; Multiplicity: N^{K^{+}}/(N_{DIS}*#Delta4D)");
    g_multi_xQ2zPt->SetMarkerStyle(20);
    g_multi_xQ2zPt->SetMarkerColor(kRed+1), g_multi_xQ2zPt->SetLineColor(kRed+1);
    g_multi_xQ2zPt->Draw("AP");
    c_multi_4d->Write();



    //

    kp_evnt_chi2.Write();
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

    std::vector<TH2D*> hists_kp = {
        &kp_MomVsPhT, &kp_MomVsXb, &kp_MomVsXf, &kp_MomVsZ, &kp_MomVsY, &kp_MomVsEta,
        &kp_MomVsTheta, &kp_MomVsPhi_h, &kp_MomVsMx, &kp_MomVsBeta, &kp_MomVsCh, &kp_MomVsMass_Rich,
        &kp_Q2VsXb, &kp_Q2VsXf, &kp_Q2VsMom, &kp_Q2VsPhT, &kp_Q2VsZ, &kp_Q2VsY, &kp_Q2VsEta, &kp_Q2VsPhi_h,
        &kp_PhTvsZ, &kp_PhTvsXb, &kp_PhTvsEta, &kp_PhTvsPhi_h,
        &kp_zVsXb, &kp_zVsXf, &kp_zVsEta, &kp_zVsPhi_h,
        &kp_xBvsY, &kp_xBvsEta, &kp_xBvsSinTheta, &kp_ThetaVsPhi_h, &kp_ThetaVsPhi_Lab, &kp_PxVsPy, &el_PxVsPy,
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

    dir_multi->cd();
    int bin_csv = 1;
    for (int x = 0; x < nbin_xQ2; x++) {
        TGraphErrors* g_multi_xQ2_binselected_1 = new TGraphErrors(7);
        TGraphErrors* g_multi_xQ2_binselected_2 = new TGraphErrors(7);
        TGraphErrors* g_multi_xQ2_binselected_3 = new TGraphErrors(7);
        TGraphErrors* g_multi_xQ2_binselected_4 = new TGraphErrors(4);
        double delta_Q2 = bin_Q2_long[x][1] - bin_Q2_long[x][0];
        double delta_xB = bin_xB[x][1] - bin_xB[x][0];
        double denom_events = vec_electron_y_2d[x].size()*delta_xB*delta_Q2;
        for (int z = 0; z < nbin_zPt; z++) {
            double sum_z = 0.0, sum_xB = 0.0, sum_Q2 = 0.0, sum_Pt = 0.0;
            if (vec_kaonp_z_4d[x][z].empty()){ 
                csvFile << bin_csv << "," << x+1 << "," << z+1 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << endl;
                bin_csv++;
                continue;
            }
            // -------
            for (double val_z : vec_kaonp_z_4d[x][z]) sum_z += val_z;
            double mean_z = sum_z / vec_kaonp_z_4d[x][z].size();
            for (double val_xB : vec_kaonp_xB_4d[x][z]) sum_xB += val_xB;
            double mean_xB = sum_xB / vec_kaonp_xB_4d[x][z].size();
            for (double val_Q2 : vec_kaonp_Q2_4d[x][z]) sum_Q2 += val_Q2;
            double mean_Q2 = sum_Q2 / vec_kaonp_Q2_4d[x][z].size();
            for (double val_Pt : vec_kaonp_Pt_4d[x][z]) sum_Pt += val_Pt;
            double mean_Pt = sum_Pt / vec_kaonp_Pt_4d[x][z].size();
            // -------
            double delta_z = bin_z[z][1] - bin_z[z][0];
            double delta_Pt = bin_Pt_long[z][1] - bin_Pt_long[z][0];
            double multiplicity_4d = vec_kaonp_y_4d[x][z].size() / (denom_events*delta_z*delta_Pt);
            double err_4d = 0.0;
            if (vec_kaonp_y_4d[x][z].size() > 0) err_4d = 1/sqrt(vec_kaonp_y_4d[x][z].size());
            if (err_4d > 0.3) err_4d = 0.3; 
            csvFile << bin_csv << "," << x+1 << "," << z+1 << "," << mean_xB << "," << mean_Q2 << "," << mean_z << "," << mean_Pt << "," << multiplicity_4d << "," << err_4d << endl;
            bin_csv++;
            if(z < 7){
                g_multi_xQ2_binselected_1->SetPoint(z, mean_z, multiplicity_4d);
                g_multi_xQ2_binselected_1->SetPointError(z, 0.0, err_4d);   
            } else if (z < 14){
                g_multi_xQ2_binselected_2->SetPoint(z-7, mean_z, multiplicity_4d);
                g_multi_xQ2_binselected_2->SetPointError(z-7, 0.0, err_4d); 
            } else if (z < 21){
                g_multi_xQ2_binselected_3->SetPoint(z-14, mean_z, multiplicity_4d);
                g_multi_xQ2_binselected_3->SetPointError(z-14, 0.0, err_4d); 
            } else {
                g_multi_xQ2_binselected_4->SetPoint(z-21, mean_z, multiplicity_4d);
                g_multi_xQ2_binselected_4->SetPointError(z-21, 0.0, err_4d); 
            }
        }
        vector<string> titles_multi_xQ2bin = {"x_{B} < 0.12 & 1 < Q^{2} < 3 GeV^{2}", "0.12 < x_{B} < 0.18 & 1 < Q^{2} < 3 GeV^{2}", "0.18 < x_{B} < 0.24 & 1 < Q^{2} < 3 GeV^{2}", "0.24 < x_{B} < 0.30 & 1 < Q^{2} < 3 GeV^{2}", "x_{B} > 0.30 & 1 < Q^{2} < 3 GeV^{2}",
            "x_{B} < 0.24 & 3 < Q^{2} < 5 GeV^{2}", "0.24 < x_{B} < 0.30 & 3 < Q^{2} < 5 GeV^{2}", "0.30 < x_{B} < 0.36 & 3 < Q^{2} < 5 GeV^{2}", "0.36 < x_{B} < 0.44 & 3 < Q^{2} < 5 GeV^{2}", "x_{B} > 0.44 & 3 < Q^{2} < 5 GeV^{2}",    
            "x_{B} < 0.44 & 5 < Q^{2} < 7 GeV^{2}", "x_{B} > 0.44 & 5 < Q^{2} < 7 GeV^{2}", 
            "x_{B} < 0.55 & 7 < Q^{2} < 10 GeV^{2}", "x_{B} > 0.55 & 7 < Q^{2} < 10 GeV^{2}",
        };
        TCanvas* c_multi_xQ2_bin = new TCanvas(Form("c_multi4D_xQ2_bin%d", x+1), Form("Multiplicity in xQ2 bin %d", x+1), 800, 600);
        //g_multi_xQ2_binselected_1->SetTitle(Form("Multiplicity in x_{B}-Q^{2} bin %d; z; Multiplicity", x+1));
        g_multi_xQ2_binselected_1->SetMarkerStyle(20),  g_multi_xQ2_binselected_2->SetMarkerStyle(20),  g_multi_xQ2_binselected_3->SetMarkerStyle(20),  g_multi_xQ2_binselected_4->SetMarkerStyle(20);
        g_multi_xQ2_binselected_1->SetMarkerColor(kAzure-5), g_multi_xQ2_binselected_1->SetLineColor(kAzure-5);
        g_multi_xQ2_binselected_2->SetMarkerColor(kViolet-5), g_multi_xQ2_binselected_2->SetLineColor(kViolet-5);
        g_multi_xQ2_binselected_3->SetMarkerColor(kPink-5), g_multi_xQ2_binselected_3->SetLineColor(kPink-5);
        g_multi_xQ2_binselected_4->SetMarkerColor(kOrange-5), g_multi_xQ2_binselected_4->SetLineColor(kOrange-5);
        //
        TMultiGraph *mg_multiplicity = new TMultiGraph();
        mg_multiplicity->SetTitle(Form("Multiplicity in x_{B}-Q^{2} bin %d | %s; z; M^{K^{+}}(x_{B},Q^{2},z,P_{hT})", x+1, titles_multi_xQ2bin[x].c_str())); // #sqrt{2#epsilon(1+#epsilon)}
        mg_multiplicity->Add(g_multi_xQ2_binselected_1, "PL");
        mg_multiplicity->Add(g_multi_xQ2_binselected_2, "PL");
        mg_multiplicity->Add(g_multi_xQ2_binselected_3, "PL");
        mg_multiplicity->Add(g_multi_xQ2_binselected_4, "PL");
        mg_multiplicity->Draw("A");

        TLegend* legend = new TLegend(0.67, 0.7, 0.87, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(g_multi_xQ2_binselected_1, "0.00 < P_{hT} < 0.25 GeV", "ep");
        legend->AddEntry(g_multi_xQ2_binselected_2, "0.25 < P_{hT} < 0.50 GeV", "ep");
        legend->AddEntry(g_multi_xQ2_binselected_3, "0.50 < P_{hT} < 0.80 GeV", "ep");
        legend->AddEntry(g_multi_xQ2_binselected_4, "0.80 < P_{hT} < 1.40 GeV", "ep");
        //legend2->AddEntry(graph2D_corr_2sin, "TSA contribution", "p");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_multi_xQ2_bin->Update();
        c_multi_xQ2_bin->Write();
    }



    //

    outFile.cd();

    //outFile.Write();
    csvFile.close();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << " and csv: " << csv_filename << endl;
}
