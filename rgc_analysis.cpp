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
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rgc/fall22/NH3/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/fall22_NH3_data
// 
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
    const double dilution = 17.0/3.0;

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

    const double dilution = 0.24;

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



void rgc_analysis() {
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

    string inputDir = "fall22_NH3_data_new";
    vector<string> rootFiles = {"fall22_NH3_data_new", "sum22_NH3_data"};

    TTree treeKaonP("Kaon+", "");
    TChain chainKaonP("Kaon+");
    int fileCount = 0;
    
    for (const auto &entry : fs::directory_iterator(inputDir)) {
        if (entry.path().extension() == ".root") {
            string filePath = entry.path().string();
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
        cerr << "Nessun file .root trovato in " << inputDir << endl;
        return;
    }
    // creo un output root 
    const char* outputFile = "plot_rgc_fall22_kaonp.root"; 
    double torus = -1;

    TFile outFile(outputFile, "RECREATE");  // File di output ROOT

    TDirectory* dir_xQ2 = outFile.mkdir("Binning xB-Q2");
    TDirectory* dir_zPt = outFile.mkdir("Binning z-Pt");
    TDirectory* dir_aul_sin = outFile.mkdir("AUL sin(Phi)");
    TDirectory* dir_aul_2sin = outFile.mkdir("AUL sin(2Phi)");
    TDirectory* dir_error = outFile.mkdir("Error");

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
    TH1D kp_phi_plus ("_kp_phi_plus", "#Phi_{h} when spin = 1 | 1.2 < Mom < 8 GeV ; #Phi_{h} [Rad]; count", nbin, -M_PI, M_PI);
    TH1D kp_phi_minus ("_kp_phi_minus", "#Phi_{h} when spin = -1 | 1.2 < Mom < 8 GeV ; #Phi_{h} [Rad]; count", nbin, -M_PI, M_PI);
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
    // Angles
    TH2D kp_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0.1, 0.4);
    TH2D kp_ThetaVsPhi_Lab ("_ThetaVsPhi_Lab", "Correlation #theta vs #Phi_{Lab} | K+ | with EventBuilder + RICH ; #Phi_{Lab} [Rad]; #theta [Rad]", nbin, -TMath::Pi(), TMath::Pi(), nbin, 0.05, 0.4);
    
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
    vector<TH2D*> hist_xBvsQ2_binned(nbin_xQ2);
    vector<TH2D*> hist_zvsPt_binned(nbin_zPt);
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
    vector<vector<double>> AUL_sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_2sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    // TSA
    vector<double> AUL_gamma_sin_2d(nbin_xQ2);
    vector<double> AUL_gamma_sin_err_2d(nbin_xQ2);
    vector<double> AUL_gamma_2sin_2d(nbin_xQ2);
    vector<double> AUL_gamma_2sin_err_2d(nbin_xQ2);
    vector<vector<double>> AUL_gamma_sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_gamma_sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_gamma_2sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> AUL_gamma_2sin_err_4d(nbin_xQ2, vector<double> (nbin_zPt));
    // correction
    vector<double> Corr_sin_2d(nbin_xQ2);
    vector<double> Corr_2sin_2d(nbin_xQ2);
    vector<vector<double>> Corr_sin_4d(nbin_xQ2, vector<double> (nbin_zPt));
    vector<vector<double>> Corr_2sin_4d(nbin_xQ2, vector<double> (nbin_zPt));


    for(int i = 0; i < nbin_xQ2; i++){
        hist_xBvsQ2_binned[i] = new TH2D(Form("hist_xBvsQ2_bin%d", i+1),
                                Form("x_{B} vs Q^{2} for bin (%d, x_{B}-Q^{2}) | pion+; x_{B}; Q^{2} [GeV^{2}]", i+1), 120, 0, 0.8, 120, 1, 11);
    }
    for(int i = 0; i < nbin_zPt; i++){
        hist_zvsPt_binned[i] = new TH2D(Form("hist_zvsPt_bin%d", i+1),
                                Form("z vs P_{hT} for bin (%d, z-P_{hT}) | pion+; z; P_{hT} [GeV^]", i+1), 120, 0.2, 1.0, 120, 0, 1.4);
    }

    //
    double bin_z_plot[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    double bin_Pt_plot[] = {0, 0.25, 0.5, 0.8, 1.4};
    vector<TH2D*> hist_AUL_sin_xQ2_zPt(nbin_xQ2);
    vector<TH2D*> hist_AUL_sin2_xQ2_zPt(nbin_xQ2);
    for (int ix = 0; ix < nbin_xQ2; ix++){
        hist_AUL_sin_xQ2_zPt[ix] = new TH2D(Form("hist_AUL_sin_xQ2_%d", ix+1),
            Form("AUL_sinPhi error bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
        hist_AUL_sin2_xQ2_zPt[ix] = new TH2D(Form("hist_AUL_sin2_xQ2_%d", ix+1),
            Form("AUL_sin2Phi error bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 7, bin_z_plot, 4, bin_Pt_plot);
    }
    


    //

    // ORA RIEMPI I GRAFICI
    // Kaon+
    Long64_t nEntries_kp = chainKaonP.GetEntries();
    for (Long64_t i = 0; i < nEntries_kp; i++) {
        // RICH cut
        //if (electron_ass != 321) continue;
        chainKaonP.GetEntry(i);
        if (kaonp_track_pid == -11) continue;
        if (kaonp_y <= 0.2) continue;
        if (std::abs(kaonp_Pol) < 0.1) continue;
        //
        //if(kaonp_Mx < 1.5) continue;
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
                }
                if(index_zPt >= 0) hist_zvsPt_binned[index_zPt-1]->Fill(kaonp_z, kaonp_PhT);
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
                kp_ThetaVsPhi_h.Fill(kaonp_Phi_h, kaonp_Theta);
                kp_ThetaVsPhi_Lab.Fill(kaonp_Phi, kaonp_Theta);
                treeKaonP.Fill();

                // plot del RICH
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
        &kp_xBvsY, &kp_ThetaVsPhi_h, &kp_ThetaVsPhi_Lab,
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




    // ASYMMETRIES
    for (int x = 0; x < nbin_xQ2; x++){ 
        ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
        //ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike(p, vec_kaonp_phih_2d[x], vec_kaonp_spin_2d[x], vec_kaonp_eps_2d[x], vec_kaonp_pol_2d[x], vec_kaonp_sintheta_2d[x], false);}, 2);
        ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_2d[x], vec_kaonp_spin_2d[x], vec_kaonp_eps_2d[x], vec_kaonp_pol_2d[x], vec_helicity_2d[x], vec_kaonp_sintheta_2d[x], true);}, 7);
        minimizer.SetFunction(MLE);
        minimizer.SetMaxFunctionCalls(10000);
        minimizer.SetMaxIterations(10000);
        minimizer.SetTolerance(0.005);
        minimizer.SetVariable(0, "Aul_sin", 0.00, 0.001); 
        minimizer.SetVariable(1, "Aul_2sin", 0.00, 0.001);
        minimizer.SetLimitedVariable(2, "All", 0.00, 0.001, -0.5, 0.5); 
        minimizer.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.5, 0.5);
        minimizer.SetLimitedVariable(4, "Auu", 0.00, 0.001, -1, 1);
        minimizer.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -1, 1);
        minimizer.SetLimitedVariable(6, "Alu", 0.00, 0.001, -1, 1);
        //minimizer.SetLimitedVariable(7, "Aul_3sin", 0.00, 0.000, -1, 1);
        //minimizer.SetLimitedVariable(8, "All_cos2", 0.00, 0.000, -0.03, 0.03);
        minimizer.Minimize();
        // LSA
        AUL_sin_2d[x] = minimizer.X()[0];
        AUL_sin_err_2d[x] = minimizer.Errors()[0];
        AUL_2sin_2d[x] = minimizer.X()[1];
        AUL_2sin_err_2d[x] = minimizer.Errors()[1];
        cout <<  "    All = " << minimizer.X()[2] << "    All_cos = " << minimizer.X()[3]  << "    Auu_cos = " << minimizer.X()[4] << "    Auu_cos2 = " << minimizer.X()[5] <<  "    Alu_sin = " << minimizer.X()[6] << endl;
        // TSA
        // mean value
        double sum_sint = 0, sum_eps = 0;
        for (double val_sint : vec_kaonp_sintheta_2d[x]) sum_sint += val_sint;
        double mean_sint = sum_sint / vec_kaonp_sintheta_2d[x].size();
        for (double val_eps : vec_kaonp_eps_2d[x]) sum_eps += val_eps;
        double mean_eps = sum_eps / vec_kaonp_eps_2d[x].size();
        // correction
        double C_UL2 = mean_sint * (sqrt(2*mean_eps*(1+mean_eps))/mean_eps);
        double C_UL1 = mean_sint * 1/(sqrt(2*mean_eps*(1+mean_eps)));
        double Aut_2d = (hist_xBvsQ2_binned[x]->GetMean(1)*hist_xBvsQ2_binned[x]->GetMean(1)*hist_xBvsQ2_binned[x]->GetMean(1)*hist_xBvsQ2_binned[x]->GetMean(1))/4;
        double correction_1_2d = Aut_2d*C_UL1;
        double correction_2_2d = Aut_2d*C_UL2;
        Corr_sin_2d[x] =  correction_1_2d;
        Corr_2sin_2d[x] =  correction_2_2d;
        //double C_LL2 = mean_sint * (sqrt(2*mean_eps*(1-mean_eps))/sqrt(1-mean_eps*mean_eps));
        //double C_LL1 = mean_sint * sqrt(1-mean_eps*mean_eps)/(sqrt(2*mean_eps*(1-mean_eps)));
        // effective
        AUL_gamma_sin_2d[x] = minimizer.X()[0] + 0;
        AUL_gamma_sin_err_2d[x] = minimizer.Errors()[0];
        AUL_gamma_2sin_2d[x] = minimizer.X()[1] + 0;
        AUL_gamma_2sin_err_2d[x] = minimizer.Errors()[1];
    }
    cout << " -------------------------------------------------------------------------------------------------------------------------------------------- " << endl;
    for (int x = 0; x < nbin_xQ2; x++){ 
        for(int z = 0; z < nbin_zPt; z++){
            ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
            //ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike(p, vec_kaonp_phih_4d[x][z], vec_kaonp_spin_4d[x][z], vec_kaonp_eps_4d[x][z], vec_kaonp_pol_4d[x][z], vec_kaonp_sintheta_4d[x][z], false);}, 2);
            ROOT::Math::Functor MLE([&](const double* p) {return AUL_loglike_withLL(p, vec_kaonp_phih_4d[x][z], vec_kaonp_spin_4d[x][z], vec_kaonp_eps_4d[x][z], vec_kaonp_pol_4d[x][z], vec_helicity_4d[x][z], vec_kaonp_sintheta_4d[x][z], true);}, 7);
            minimizer.SetFunction(MLE);
            minimizer.SetMaxFunctionCalls(10000);
            minimizer.SetMaxIterations(10000);
            minimizer.SetTolerance(0.005);
            minimizer.SetVariable(0, "Aul_sin", 0.00, 0.001); 
            minimizer.SetVariable(1, "Aul_2sin", 0.00, 0.001);
            minimizer.SetLimitedVariable(2, "All", 0.00, 0.001, -0.5, 0.5); 
            minimizer.SetLimitedVariable(3, "All_cos", 0.00, 0.001, -0.5, 0.5);
            minimizer.SetLimitedVariable(4, "Auu", 0.00, 0.001, -1, 1);
            minimizer.SetLimitedVariable(5, "Auu2", 0.00, 0.001, -1, 1);
            minimizer.SetLimitedVariable(6, "Alu", 0.00, 0.001, -1, 1);
            //minimizer.SetLimitedVariable(7, "Aul_3sin", 0.00, 0.000, -1, 1);
            //minimizer.SetLimitedVariable(8, "All_cos2", 0.00, 0.000, -0.03, 0.03);
            minimizer.Minimize();
            // LSA
            AUL_sin_4d[x][z] = minimizer.X()[0];
            AUL_sin_err_4d[x][z] = minimizer.Errors()[0];
            AUL_2sin_4d[x][z] = minimizer.X()[1];
            AUL_2sin_err_4d[x][z] = minimizer.Errors()[1];
            int status = minimizer.Status();
            if (status != 0) {
                cout << "BAD FIT at x=" << x << " z=" << z
                    << "  status=" << status
                    << "  N=" << vec_kaonp_phih_4d[x][z].size()
                    << endl;
            }
            cout <<  "    All = " << minimizer.X()[2] << "    All_cos = " << minimizer.X()[3]  << "    Auu_cos = " << minimizer.X()[4] << "    Auu_cos2 = " << minimizer.X()[5] <<  "    Alu_sin = " << minimizer.X()[6] << endl;
            // TSA
            // mean value
            double sum_sint = 0, sum_eps = 0, sum_xB = 0;
            for (double val_sint : vec_kaonp_sintheta_4d[x][z]) sum_sint += val_sint;
            double mean_sint = sum_sint / vec_kaonp_sintheta_4d[x][z].size();
            for (double val_eps : vec_kaonp_eps_4d[x][z]) sum_eps += val_eps;
            double mean_eps = sum_eps / vec_kaonp_eps_4d[x][z].size();
            for (double val_xB : vec_kaonp_xB_4d[x][z]) sum_xB += val_xB;
            double mean_xB = sum_xB / vec_kaonp_xB_4d[x][z].size();
            // correction
            double C_UL1 = mean_sint * (sqrt(2*mean_eps*(1+mean_eps))/mean_eps);
            double C_UL2 = mean_sint * mean_eps/(sqrt(2*mean_eps*(1+mean_eps)));
            double Aut_4d = (mean_xB*mean_xB*mean_xB*mean_xB)/4;
            double correction_1_4d = Aut_4d*C_UL1;
            double correction_2_4d = Aut_4d*C_UL2;
            Corr_sin_4d[x][z] =  correction_1_4d;
            Corr_2sin_4d[x][z] =  correction_2_4d;
            //double C_LL2 = mean_sint * (sqrt(2*mean_eps*(1-mean_eps))/sqrt(1-mean_eps*mean_eps));
            //double C_LL1 = mean_sint * sqrt(1-mean_eps*mean_eps)/(sqrt(2*mean_eps*(1-mean_eps)));
            // effective
            AUL_gamma_sin_4d[x][z] = minimizer.X()[0] + correction_1_4d;
            AUL_gamma_sin_err_4d[x][z] = minimizer.Errors()[0];
            AUL_gamma_2sin_4d[x][z] = minimizer.X()[1] + correction_2_4d;
            AUL_gamma_2sin_err_4d[x][z] = minimizer.Errors()[1];

            if(AUL_gamma_sin_err_4d[x][z] <= 0.05){
                if(z < 7) {
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(z+1, 1, AUL_gamma_sin_err_4d[x][z]);
                } else if (z < 14) {
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(z+1-7, 2, AUL_gamma_sin_err_4d[x][z]);
                } else if (z < 21) {
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(z+1-14, 3, AUL_gamma_sin_err_4d[x][z]);
                } else if (z == 21) {
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(1, 4, AUL_gamma_sin_err_4d[x][z]);
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(2, 4, AUL_gamma_sin_err_4d[x][z]);
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(3, 4, AUL_gamma_sin_err_4d[x][z]);
                } else if (z < 24){
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(z-18, 4, AUL_gamma_sin_err_4d[x][z]);
                } else if (z == 24){
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(6, 4, AUL_gamma_sin_err_4d[x][z]);
                    hist_AUL_sin_xQ2_zPt[x]->SetBinContent(7, 4, AUL_gamma_sin_err_4d[x][z]);
                }
            }
            if(AUL_gamma_2sin_err_4d[x][z] <= 0.05){
                vec_z_4d[x][z] = vec_kaonp_z_4d[x][z];
                vec_xB_4d[x][z] = vec_kaonp_xB_4d[x][z];
                vec_Q2_4d[x][z] = vec_kaonp_Q2_4d[x][z];
                vec_Pt_4d[x][z] = vec_kaonp_Pt_4d[x][z];
                if(z < 7) {
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(z+1, 1, AUL_gamma_2sin_err_4d[x][z]);
                } else if (z < 14) {
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(z+1-7, 2, AUL_gamma_2sin_err_4d[x][z]);
                } else if (z < 21) {
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(z+1-14, 3, AUL_gamma_2sin_err_4d[x][z]);
                } else if (z == 21) {
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(1, 4, AUL_gamma_2sin_err_4d[x][z]);
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(2, 4, AUL_gamma_2sin_err_4d[x][z]);
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(3, 4, AUL_gamma_2sin_err_4d[x][z]);
                } else if (z < 24){
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(z-18, 4, AUL_gamma_2sin_err_4d[x][z]);
                } else if (z == 24){
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(6, 4, AUL_gamma_2sin_err_4d[x][z]);
                    hist_AUL_sin2_xQ2_zPt[x]->SetBinContent(7, 4, AUL_gamma_2sin_err_4d[x][z]);
                }
            }

            hist_AUL_sin2_xQ2_zPt[x]->SetMinimum(0);
            hist_AUL_sin2_xQ2_zPt[x]->SetMaximum(0.05);
            hist_AUL_sin_xQ2_zPt[x]->SetMinimum(0);
            hist_AUL_sin_xQ2_zPt[x]->SetMaximum(0.05);


            //
        }
        cout << " -------------------------------------------------------------------------------------------------------------------------------------------- " << endl;
    }


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
    TGraphErrors* graph2D_corr_sin = new TGraphErrors();
    int p_idx1 = 0, p_idx1_1 = 0, p_idx1_2 = 0, p_idx1_3 = 0, p_idx1_4 = 0;
    for(int x = 0; x < nbin_xQ2; x++){
        if (vec_kaonp_xB_2d[x].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d[x].empty()) continue;
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        //for (double valeps : vec_kaonp_eps_2d[x]) sum_eps += valeps;
        //double mean_eps = sum_eps / vec_kaonp_eps_2d[x].size();
        for (double val : vec_kaonp_xB_2d[x]){
            sum_xB += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_xB = sum_xB / vec_kaonp_xB_2d[x].size();
        // Get AUL and error from vectors
        //double D = sqrt(2*mean_eps*(1+mean_eps));
        double aul_sin = AUL_gamma_sin_2d[x];
        double aul_sin_err = AUL_gamma_sin_err_2d[x];
        //double aul_sin = 0;
        //if(aul_sin > 1 || aul_sin < -1) continue;
        graph2D_AUL_vs_xB->SetPoint(p_idx1, mean_xB, aul_sin);
        graph2D_AUL_vs_xB->SetPointError(p_idx1, 0.0, aul_sin_err); // No x error
        graph2D_corr_sin->SetPoint(p_idx1, mean_xB, Corr_sin_2d[x]);
        p_idx1++;
        if(x < 5){
            graph2D_AUL_vs_xB_1->SetPoint(p_idx1_1, mean_xB, aul_sin);
            graph2D_AUL_vs_xB_1->SetPointError(p_idx1_1, 0.0, aul_sin_err); // No x error
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
    graph2D_corr_sin->SetMarkerStyle(26);
    graph2D_corr_sin->SetMarkerColor(kRed);


    TCanvas* c_Aut2D_xB_2 = new TCanvas("Aul_sin_vs_xB_2d", "sin(#Phi_{h}) longitudinal asymmetry vs x_{B}", 800, 600);
    //graph2D_AUL_vs_xB->Draw("A");
    TMultiGraph *mg_Aut2D_xB_2 = new TMultiGraph();
    mg_Aut2D_xB_2->SetTitle("A_{UL}^{sin#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin); x_{B}; A_{UL}^{sin(#Phi_{h})}");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_1, "P");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_2, "P");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_3, "P");
    mg_Aut2D_xB_2->Add(graph2D_AUL_vs_xB_4, "P");
    mg_Aut2D_xB_2->Add(graph2D_corr_sin, "P");
    mg_Aut2D_xB_2->Draw("A");
    //mg_Aut2D_xB_2->GetYaxis()->SetRangeUser(-0.05, 0.05);
    TLine* guideLine11 = new TLine(mg_Aut2D_xB_2->GetXaxis()->GetXmin(), 0, mg_Aut2D_xB_2->GetXaxis()->GetXmax(), 0);
    guideLine11->SetLineStyle(2);  
    guideLine11->SetLineColor(kGray+1);
    guideLine11->Draw();

    TLegend* legend = new TLegend(0.13, 0.7, 0.33, 0.88); // Adjust position (x1,y1,x2,y2)
    legend->AddEntry(graph2D_AUL_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend->AddEntry(graph2D_AUL_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend->AddEntry(graph2D_AUL_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend->AddEntry(graph2D_AUL_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
    legend->AddEntry(graph2D_corr_sin, "TSA contribution", "p");
    legend->SetFillStyle(0);  // Transparent background
    legend->Draw();
    c_Aut2D_xB_2->Update();
    c_Aut2D_xB_2->Write();


    // 4D

    for (int z = 0; z < nbin_zPt; z++) {
        //if(z == 0 || z == 7) continue;
        TGraphErrors* graph_AUL_vs_xB = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_1 = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_2 = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_3 = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_4 = new TGraphErrors();
        TGraphErrors* graph4D_corr_sin = new TGraphErrors();
        int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0;
        for(int x = 0; x < nbin_xQ2; x++){
            if (vec_kaonp_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_kaonp_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_xB = 0.0;
            double min_x = 1;
            double max_x = 0;
            double sum_eps = 0;
            //for (double valeps : vec_kaonp_eps_4d[x][z]) sum_eps += valeps;
            //double mean_eps = sum_eps / vec_kaonp_eps_4d[x][z].size();
            for (double val : vec_kaonp_xB_4d[x][z]){
                sum_xB += val;
                if(val > max_x) max_x = val;
                if(val < min_x) min_x = val;
            }
            double mean_xB = sum_xB / vec_kaonp_xB_4d[x][z].size();
            // Get AUL and error from vectors
            //double D = sqrt(2*mean_eps*(1+mean_eps));
            double aul_sin = AUL_gamma_sin_4d[x][z];
            //double aul_sin = 0;
            double aul_sin_err = AUL_gamma_sin_err_4d[x][z];
            if(aul_sin > 1 || aul_sin < -1) continue;
            if(aul_sin_err > 0.25) aul_sin_err = 0.25;
            graph_AUL_vs_xB->SetPoint(p_idx, mean_xB, aul_sin);
            graph_AUL_vs_xB->SetPointError(p_idx, 0.0, aul_sin_err); // No x error
            graph4D_corr_sin->SetPoint(p_idx, mean_xB, Corr_sin_4d[x][z]);
            p_idx++;
            if(x < 5){
                graph_AUL_vs_xB_1->SetPoint(p_idx_1, mean_xB, aul_sin);
                graph_AUL_vs_xB_1->SetPointError(p_idx_1, 0.0, aul_sin_err); // No x error
                p_idx_1++;
            } else if (x < 10){
                graph_AUL_vs_xB_2->SetPoint(p_idx_2, mean_xB, aul_sin);
                graph_AUL_vs_xB_2->SetPointError(p_idx_2, 0.0, aul_sin_err);
                p_idx_2++;
            } else if (x < 12){
                graph_AUL_vs_xB_3->SetPoint(p_idx_3, mean_xB, aul_sin);
                graph_AUL_vs_xB_3->SetPointError(p_idx_3, 0.0, aul_sin_err);
                p_idx_3++;
            } else if (x < 18){
                graph_AUL_vs_xB_4->SetPoint(p_idx_4, mean_xB, aul_sin);
                graph_AUL_vs_xB_4->SetPointError(p_idx_4, 0.0, aul_sin_err);
                p_idx_4++;
            }
        }

        graph_AUL_vs_xB_1->SetMarkerStyle(20), graph_AUL_vs_xB_2->SetMarkerStyle(20), graph_AUL_vs_xB_3->SetMarkerStyle(20), graph_AUL_vs_xB_4->SetMarkerStyle(20);
        graph_AUL_vs_xB_1->SetLineColor(kAzure-5), graph_AUL_vs_xB_1->SetMarkerColor(kAzure-5);
        graph_AUL_vs_xB_2->SetLineColor(kViolet-5), graph_AUL_vs_xB_2->SetMarkerColor(kViolet-5);
        graph_AUL_vs_xB_3->SetLineColor(kPink-5), graph_AUL_vs_xB_3->SetMarkerColor(kPink-5);
        graph_AUL_vs_xB_4->SetLineColor(kOrange-5), graph_AUL_vs_xB_4->SetMarkerColor(kOrange-5);


        TCanvas* c_Aut_xB_2 = new TCanvas(Form("Aul_sin_vs_xB_4d_zPt_bin%d", z+1), "sin(#Phi_{h}) longitudinal asymmetry vs x_{B} for bin (n, z-P_{hT})", 800, 600);
        //graph_AUL_vs_xB->Draw("A");
        TMultiGraph *mg_Aut4D_xB_2 = new TMultiGraph();
        mg_Aut4D_xB_2->SetTitle(Form("A_{UL}^{sin#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin) for bin (%d, z-P_{hT}) ; x_{B}; A_{UL}^{sin(#Phi_{h})}", z+1)); // #sqrt{2#epsilon(1+#epsilon)}
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_1, "P");
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_2, "P");
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_3, "P");
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_4, "P");
        mg_Aut4D_xB_2->Draw("A");
        //mg_Aut4D_xB_2->GetYaxis()->SetRangeUser(-0.1, 0.1);
        TLine* guideLine2 = new TLine(mg_Aut4D_xB_2->GetXaxis()->GetXmin(), 0, mg_Aut4D_xB_2->GetXaxis()->GetXmax(), 0);
        guideLine2->SetLineStyle(2);  
        guideLine2->SetLineColor(kGray+1);
        guideLine2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(graph_AUL_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
        legend->AddEntry(graph_AUL_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
        legend->AddEntry(graph_AUL_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
        legend->AddEntry(graph_AUL_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Aut_xB_2->Update();
        c_Aut_xB_2->Write();
    }


    dir_aul_2sin->cd();
    // Q2 color scale vs xB 2D
    TGraphErrors* graph2D_AUL2_vs_xB = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_1 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_2 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_3 = new TGraphErrors();
    TGraphErrors* graph2D_AUL2_vs_xB_4 = new TGraphErrors();
    TGraphErrors* graph2D_corr_2sin = new TGraphErrors();
    int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0;
    for(int x = 0; x < nbin_xQ2; x++){
        if (vec_kaonp_xB_2d[x].empty()) continue; // Skip empty bins
        if (vec_kaonp_z_2d[x].empty()) continue;
        // Compute mean xB for this bin
        double sum_xB = 0.0;
        double min_x = 1;
        double max_x = 0;
        double sum_eps = 0;
        //for (double valeps : vec_kaonp_eps_2d[x]) sum_eps += valeps;
        //double mean_eps = sum_eps / vec_kaonp_eps_2d[x].size();
        for (double val : vec_kaonp_xB_2d[x]){
            sum_xB += val;
            if(val > max_x) max_x = val;
            if(val < min_x) min_x = val;
        }
        double mean_xB = sum_xB / vec_kaonp_xB_2d[x].size();
        // Get AUL2 and error from vectors
        //double D = sqrt(2*mean_eps*(1+mean_eps));
        double aul_sin = AUL_gamma_2sin_2d[x];
        //double aul_sin = 0;
        double aul_sin_err = AUL_gamma_2sin_err_2d[x];
        if(aul_sin > 1 || aul_sin < -1) continue;
        graph2D_AUL2_vs_xB->SetPoint(p_idx, mean_xB, aul_sin);
        graph2D_AUL2_vs_xB->SetPointError(p_idx, 0.0, aul_sin_err); // No x error
        graph2D_corr_2sin->SetPoint(p_idx, mean_xB, Corr_2sin_2d[x]);
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

    graph2D_AUL2_vs_xB->SetTitle("A_{UL}^{sin2#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin); x_{B}; A_{UL}^{sin(2#Phi_{h})}"); // #sqrt{2#epsilon(1+#epsilon)}
    graph2D_AUL2_vs_xB_1->SetMarkerStyle(20), graph2D_AUL2_vs_xB_2->SetMarkerStyle(20), graph2D_AUL2_vs_xB_3->SetMarkerStyle(20), graph2D_AUL2_vs_xB_4->SetMarkerStyle(20);
    graph2D_AUL2_vs_xB_1->SetLineColor(kAzure-5), graph2D_AUL2_vs_xB_1->SetMarkerColor(kAzure-5);
    graph2D_AUL2_vs_xB_2->SetLineColor(kViolet-5), graph2D_AUL2_vs_xB_2->SetMarkerColor(kViolet-5);
    graph2D_AUL2_vs_xB_3->SetLineColor(kPink-5), graph2D_AUL2_vs_xB_3->SetMarkerColor(kPink-5);
    graph2D_AUL2_vs_xB_4->SetLineColor(kOrange-5), graph2D_AUL2_vs_xB_4->SetMarkerColor(kOrange-5);
    graph2D_corr_2sin->SetMarkerStyle(26);
    graph2D_corr_2sin->SetMarkerColor(kRed);
    


    TCanvas* c_Aut2D2_xB_2 = new TCanvas("Aul_sin2_vs_xB_2d", "sin(2#Phi_{h}) longitudinal asymmetry vs x_{B}", 800, 600);
    //graph2D_AUL2_vs_xB->Draw("A");
    TMultiGraph *mg_Aut2D2_xB_2 = new TMultiGraph();
    mg_Aut2D2_xB_2->SetTitle("A_{UL}^{sin2#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin); x_{B}; A_{UL}^{sin(2#Phi_{h})}"); // #sqrt{2#epsilon(1+#epsilon)}
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_1, "P");
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_2, "P");
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_3, "P");
    mg_Aut2D2_xB_2->Add(graph2D_AUL2_vs_xB_4, "P");
    mg_Aut2D2_xB_2->Add(graph2D_corr_2sin, "P");
    mg_Aut2D2_xB_2->Draw("A");
    //mg_Aut2D2_xB_2->GetYaxis()->SetRangeUser(-0.05, 0.05);
    guideLine11->Draw();
    TLegend* legend2 = new TLegend(0.13, 0.7, 0.33, 0.88); // Adjust position (x1,y1,x2,y2)
    legend2->AddEntry(graph2D_AUL2_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
    legend2->AddEntry(graph2D_AUL2_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
    legend2->AddEntry(graph2D_AUL2_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
    legend2->AddEntry(graph2D_AUL2_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
    legend2->AddEntry(graph2D_corr_2sin, "TSA contribution", "p");
    legend2->SetFillStyle(0);  // Transparent background
    legend2->Draw();
    c_Aut2D2_xB_2->Update();
    c_Aut2D2_xB_2->Write();

    // 4D
    for (int z = 0; z < nbin_zPt; z++) {
        //if(z == 0 || z == 7) continue;
        TGraphErrors* graph_AUL_vs_xB = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_1 = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_2 = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_3 = new TGraphErrors();
        TGraphErrors* graph_AUL_vs_xB_4 = new TGraphErrors();
        TGraphErrors* graph4D_corr_2sin = new TGraphErrors();
        int p_idx = 0, p_idx_1 = 0, p_idx_2 = 0, p_idx_3 = 0, p_idx_4 = 0;
        for(int x = 0; x < nbin_xQ2; x++){
            if (vec_kaonp_xB_4d[x][z].empty()) continue; // Skip empty bins
            if (vec_kaonp_z_4d[x][z].empty()) continue;
            // Compute mean xB for this bin
            double sum_eps = 0;
            //for (double valeps : vec_kaonp_eps_4d[x][z]) sum_eps += valeps;
            //double mean_eps = sum_eps / vec_kaonp_eps_4d[x][z].size();
            double sum_xB = 0.0;
            double min_x = 1;
            double max_x = 0;
            for (double val : vec_kaonp_xB_4d[x][z]){
                sum_xB += val;
                if(val > max_x) max_x = val;
                if(val < min_x) min_x = val;
            }
            double mean_xB = sum_xB / vec_kaonp_xB_4d[x][z].size();
            // Get AUL and error from vectors
            double aul_sin = AUL_gamma_2sin_4d[x][z];
            //double aul_sin = 0;
            double aul_sin_err = AUL_gamma_2sin_err_4d[x][z];
            if(aul_sin > 1 || aul_sin < -1) continue;
            if(aul_sin_err > 0.4) aul_sin_err = 0.4;
            graph_AUL_vs_xB->SetPoint(p_idx, mean_xB, aul_sin);
            graph_AUL_vs_xB->SetPointError(p_idx, 0.0, aul_sin_err); // No x error
            graph4D_corr_2sin->SetPoint(p_idx, mean_xB, Corr_2sin_4d[x][z]);
            p_idx++;
            if(x < 5){
                graph_AUL_vs_xB_1->SetPoint(p_idx_1, mean_xB, aul_sin);
                graph_AUL_vs_xB_1->SetPointError(p_idx_1, 0.0, aul_sin_err); // No x error
                p_idx_1++;
            } else if (x < 10){
                graph_AUL_vs_xB_2->SetPoint(p_idx_2, mean_xB, aul_sin);
                graph_AUL_vs_xB_2->SetPointError(p_idx_2, 0.0, aul_sin_err);
                p_idx_2++;
            } else if (x < 12){
                graph_AUL_vs_xB_3->SetPoint(p_idx_3, mean_xB, aul_sin);
                graph_AUL_vs_xB_3->SetPointError(p_idx_3, 0.0, aul_sin_err);
                p_idx_3++;
            } else if (x < 18){
                graph_AUL_vs_xB_4->SetPoint(p_idx_4, mean_xB, aul_sin);
                graph_AUL_vs_xB_4->SetPointError(p_idx_4, 0.0, aul_sin_err);
                p_idx_4++;
            }
        }

        graph_AUL_vs_xB->SetTitle(Form("A_{UL}^{2sin#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin) for bin (%d, z-P_{hT}) ; x_{B}; A_{UL}^{2sin(#Phi_{h})}", z+1));
        graph_AUL_vs_xB_1->SetMarkerStyle(20), graph_AUL_vs_xB_2->SetMarkerStyle(20), graph_AUL_vs_xB_3->SetMarkerStyle(20), graph_AUL_vs_xB_4->SetMarkerStyle(20);
        graph_AUL_vs_xB_1->SetLineColor(kAzure-5);
        graph_AUL_vs_xB_1->SetMarkerColor(kAzure-5);
        graph_AUL_vs_xB_2->SetLineColor(kViolet-5);
        graph_AUL_vs_xB_2->SetMarkerColor(kViolet-5);
        graph_AUL_vs_xB_3->SetLineColor(kPink-5);
        graph_AUL_vs_xB_3->SetMarkerColor(kPink-5);
        graph_AUL_vs_xB_4->SetLineColor(kOrange-5);
        graph_AUL_vs_xB_4->SetMarkerColor(kOrange-5);


        TCanvas* c_Aut_xB_2 = new TCanvas(Form("Aul_2sin_vs_xB_4d_zPt_bin%d", z+1), "sin(#Phi_{h}) longitudinal asymmetry vs x_{B} for bin (n, z-P_{hT})", 800, 600);
        //graph_AUL_vs_xB->Draw("A");
        TMultiGraph *mg_Aut4D_xB_2 = new TMultiGraph();
        mg_Aut4D_xB_2->SetTitle(Form("A_{UL}^{2sin#Phi_{h}} vs x_{B} (x_{B}-Q^{2} bin) for bin (%d, z-P_{hT}) ; x_{B}; A_{UL}^{2sin(#Phi_{h})}", z+1));
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_1, "P");
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_2, "P");
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_3, "P");
        mg_Aut4D_xB_2->Add(graph_AUL_vs_xB_4, "P");
        mg_Aut4D_xB_2->Draw("A");
        //mg_Aut4D_xB_2->GetYaxis()->SetRangeUser(-0.1, 0.1);
        TLine* guideLine2 = new TLine(mg_Aut4D_xB_2->GetXaxis()->GetXmin(), 0, mg_Aut4D_xB_2->GetXaxis()->GetXmax(), 0);
        guideLine2->SetLineStyle(2);  
        guideLine2->SetLineColor(kGray+1);
        guideLine2->Draw();

        TLegend* legend = new TLegend(0.13, 0.7, 0.35, 0.88); // Adjust position (x1,y1,x2,y2)
        legend->AddEntry(graph_AUL_vs_xB_1, "1 < Q^{2} [GeV^{2}] < 3", "ep");
        legend->AddEntry(graph_AUL_vs_xB_2, "3 < Q^{2} [GeV^{2}] < 5", "ep");
        legend->AddEntry(graph_AUL_vs_xB_3, "5 < Q^{2} [GeV^{2}] < 7", "ep");
        legend->AddEntry(graph_AUL_vs_xB_4, "7 < Q^{2} [GeV^{2}] < 11", "ep");
        legend->SetFillStyle(0);  // Transparent background
        legend->Draw();
        c_Aut_xB_2->Update();
        c_Aut_xB_2->Write();
    }












    // errors
    dir_error->cd();
    for(int x = 0; x < nbin_xQ2; x++){
        TCanvas *c_bin_zPt = new TCanvas(Form("AULerr_sin_z_vs_Pt_Bin%d",x), "z vs P_{hT} bin", 800, 700);
        hist_AUL_sin_xQ2_zPt[x]->SetStats(0);
        hist_AUL_sin_xQ2_zPt[x]->Draw("COLZ");

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

            }
        }

        c_bin_zPt->Update();
        c_bin_zPt->Write();
    }
    for(int x = 0; x < nbin_xQ2; x++){
        TCanvas *c_bin_zPt = new TCanvas(Form("AULerr_sin2_z_vs_Pt_Bin%d",x), "z vs P_{hT} bin", 800, 700);
        hist_AUL_sin2_xQ2_zPt[x]->SetStats(0);
        hist_AUL_sin2_xQ2_zPt[x]->Draw("COLZ");

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

            }
        }

        c_bin_zPt->Update();
        c_bin_zPt->Write();
    }





    // multiplot

    // =======================================================================
    // CONFIGURAZIONE
    // =======================================================================

    const int nRows = 4;
    const int nCols = 8;

    int layout[nRows][nCols] = {
        { 0,  0,  0, 0, 0, 0, 13, 14},
        { 0,  0,  0, 0, 0, 11, 12, 0},
        { 0,  0,  6, 7, 8, 9, 10, 0},
        { 1,  2,  3, 4, 5, 0,  0,  0}
    };

    // Nuovi assi lineari
    double x_min = 0.0,  x_max = 0.8;   // xB
    double y_min = 1.0,  y_max = 10.0;  // Q²

    // =======================================================================
    // CANVAS
    // =======================================================================

    TCanvas *c_layout = new TCanvas("c_AUL_sin", "AUL sin multipanel", 2000, 1300);
    c_layout->cd();

    // =======================================================================
    // PAD PER ASSI GLOBALI (lineari)
    // =======================================================================

    TPad *pad_axes = new TPad("pad_axes", "Global Axes", 0, 0, 1, 1);
    pad_axes->SetFrameLineWidth(0);
    pad_axes->Draw();
    pad_axes->cd();

    pad_axes->SetTickx(0);
    pad_axes->SetTicky(0);

    double x_ndc_min = 0.07, x_ndc_max = 0.96;
    double y_ndc_min = 0.07, y_ndc_max = 0.96;

    // -------- Draw axes lines --------
    TLine *xAxisLine = new TLine(x_ndc_min, y_ndc_min, x_ndc_max, y_ndc_min);
    xAxisLine->SetLineWidth(2);
    xAxisLine->SetNDC(true);
    xAxisLine->Draw();

    TLine *yAxisLine = new TLine(x_ndc_min, y_ndc_min, x_ndc_min, y_ndc_max);
    yAxisLine->SetLineWidth(2);
    yAxisLine->SetNDC(true);
    yAxisLine->Draw();


    // =======================================================================
    // TICKS E LABELS — asse x lineare
    // =======================================================================

    int nTicksX = 8;
    for (int i = 0; i <= nTicksX; i++) {

        double xv = x_min + i*(x_max-x_min)/nTicksX;
        double pos = x_ndc_min + (xv - x_min)/(x_max - x_min)*(x_ndc_max - x_ndc_min);

        TLine *tick = new TLine(pos, y_ndc_min, pos, y_ndc_min - 0.015);
        tick->SetNDC(true);
        tick->Draw();

        TLatex *lab = new TLatex(pos, y_ndc_min - 0.035, Form("%.1f", xv));
        lab->SetNDC(true);
        lab->SetTextAlign(22);
        lab->SetTextSize(0.02);
        lab->Draw();
    }

    // Label asse X
    TLatex *xlabel = new TLatex(0.85, y_ndc_min - 0.065, "x_{B}");
    xlabel->SetTextSize(0.03);
    xlabel->SetNDC(true);
    xlabel->Draw();


    // =======================================================================
    // TICKS E LABELS — asse y lineare
    // =======================================================================

    int nTicksY = 9;
    for (int i = 0; i <= nTicksY; i++) {

        double yv = y_min + i*(y_max-y_min)/nTicksY;
        double pos = y_ndc_min + (yv - y_min)/(y_max - y_min)*(y_ndc_max - y_ndc_min);

        TLine *tick = new TLine(x_ndc_min, pos, x_ndc_min - 0.015, pos);
        tick->SetNDC(true);
        tick->Draw();

        TLatex *lab = new TLatex(x_ndc_min - 0.025, pos, Form("%.1f", yv));
        lab->SetNDC(true);
        lab->SetTextAlign(32);
        lab->SetTextSize(0.02);
        lab->Draw();
    }

    // Label asse Y
    TLatex *ylabel = new TLatex(x_ndc_min - 0.05, 0.9, "Q^{2}  (GeV^{2})");
    ylabel->SetTextSize(0.03);
    ylabel->SetTextAngle(90);
    ylabel->SetNDC(true);
    ylabel->Draw();

    pad_axes->Modified();
    pad_axes->Update();


    // =======================================================================
    // PAD MULTIPLI PER I PLOT
    // =======================================================================

    c_layout->cd();

    double xMargin = 0.07, yMargin = 0.07;
    double padW = (1.0 - 2*xMargin) / nCols;
    double padH = (1.0 - 2*yMargin) / nRows;

    for (int iRow = 0; iRow < nRows; ++iRow) {
        for (int iCol = 0; iCol < nCols; ++iCol) {

            int idx = layout[iRow][iCol];
            if (idx == 0) continue;

            double x1 = xMargin + iCol * padW;
            double x2 = x1 + padW;

            double y2 = 1 - yMargin - iRow * padH;
            double y1 = y2 - padH;

            TString padName = Form("pad_r%d_c%d", iRow, iCol);

            TPad *pad = new TPad(padName, padName, x1, y1, x2, y2);
            pad->SetRightMargin(0);
            pad->SetLeftMargin(0);
            pad->SetBottomMargin(0);
            pad->SetTopMargin(0);
            pad->Draw();
            pad->cd();

            hist_AUL_sin_xQ2_zPt[idx - 1]->SetTitle("");
            hist_AUL_sin_xQ2_zPt[idx - 1]->Draw("COL");
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
                    rect_zp->SetLineWidth(1);
                    rect_zp->SetLineColor(kBlack);
                    rect_zp->Draw("same");
                    gridLines_zp.push_back(rect_zp);

                }
            }

            c_layout->cd();
        }
    }


    // =======================================================================
    // PALETTE
    // =======================================================================

    TPad *pad_palette = new TPad("pad_palette", "", 0.95, 0.1, 1, 0.9);
    pad_palette->SetRightMargin(0.5);
    pad_palette->Draw();
    pad_palette->cd();

    TH2F *h_ref = (TH2F*) hist_AUL_sin_xQ2_zPt[2]->Clone("h_ref_for_palette");
    h_ref->SetStats(0);
    h_ref->Draw("COLZ");

    gPad->Update();

    TPaletteAxis *pal = (TPaletteAxis*) h_ref->GetListOfFunctions()->FindObject("palette");
    if (pal) {
        pal->SetLabelSize(0.22);
        pal->SetX1NDC(0.05);
        pal->SetX2NDC(0.6);
        pal->SetY1NDC(0.05);
        pal->SetY2NDC(0.95);
    }

    pad_palette->Modified();
    c_layout->Update();
    c_layout->Write();





    // ------------------------------------------------------------------------------



    TCanvas *c_layout2 = new TCanvas("c_AUL_sin2", "AUL sin2 multipanel", 2000, 1300);
    c_layout2->cd();

    // =======================================================================
    // PAD PER ASSI GLOBALI (lineari)
    // =======================================================================

    TPad *pad_axes2 = new TPad("pad_axes2", "Global Axes", 0, 0, 1, 1);
    pad_axes2->SetFrameLineWidth(0);
    pad_axes2->Draw();
    pad_axes2->cd();

    pad_axes2->SetTickx(0);
    pad_axes2->SetTicky(0);


    // -------- Draw axes lines --------
    TLine *xAxisLine2 = new TLine(x_ndc_min, y_ndc_min, x_ndc_max, y_ndc_min);
    xAxisLine2->SetLineWidth(2);
    xAxisLine2->SetNDC(true);
    xAxisLine2->Draw();

    TLine *yAxisLine2 = new TLine(x_ndc_min, y_ndc_min, x_ndc_min, y_ndc_max);
    yAxisLine2->SetLineWidth(2);
    yAxisLine2->SetNDC(true);
    yAxisLine2->Draw();


    // =======================================================================
    // TICKS E LABELS — asse x lineare
    // =======================================================================

    for (int i = 0; i <= nTicksX; i++) {

        double xv = x_min + i*(x_max-x_min)/nTicksX;
        double pos = x_ndc_min + (xv - x_min)/(x_max - x_min)*(x_ndc_max - x_ndc_min);

        TLine *tick2 = new TLine(pos, y_ndc_min, pos, y_ndc_min - 0.015);
        tick2->SetNDC(true);
        tick2->Draw();

        TLatex *lab2 = new TLatex(pos, y_ndc_min - 0.035, Form("%.1f", xv));
        lab2->SetNDC(true);
        lab2->SetTextAlign(22);
        lab2->SetTextSize(0.02);
        lab2->Draw();
    }

    // Label asse X
    TLatex *xlabel2 = new TLatex(0.85, y_ndc_min - 0.065, "x_{B}");
    xlabel2->SetTextSize(0.03);
    xlabel2->SetNDC(true);
    xlabel2->Draw();


    // =======================================================================
    // TICKS E LABELS — asse y lineare
    // =======================================================================

    for (int i = 0; i <= nTicksY; i++) {

        double yv = y_min + i*(y_max-y_min)/nTicksY;
        double pos = y_ndc_min + (yv - y_min)/(y_max - y_min)*(y_ndc_max - y_ndc_min);

        TLine *tick = new TLine(x_ndc_min, pos, x_ndc_min - 0.015, pos);
        tick->SetNDC(true);
        tick->Draw();

        TLatex *lab = new TLatex(x_ndc_min - 0.025, pos, Form("%.1f", yv));
        lab->SetNDC(true);
        lab->SetTextAlign(32);
        lab->SetTextSize(0.02);
        lab->Draw();
    }

    // Label asse Y
    TLatex *ylabel2 = new TLatex(x_ndc_min - 0.05, 0.9, "Q^{2}  (GeV^{2})");
    ylabel2->SetTextSize(0.03);
    ylabel2->SetTextAngle(90);
    ylabel2->SetNDC(true);
    ylabel2->Draw();

    pad_axes2->Modified();
    pad_axes2->Update();


    // =======================================================================
    // PAD MULTIPLI PER I PLOT
    // =======================================================================

    c_layout2->cd();



    for (int iRow = 0; iRow < nRows; ++iRow) {
        for (int iCol = 0; iCol < nCols; ++iCol) {

            int idx = layout[iRow][iCol];
            if (idx == 0) continue;

            double x1 = xMargin + iCol * padW;
            double x2 = x1 + padW;

            double y2 = 1 - yMargin - iRow * padH;
            double y1 = y2 - padH;

            TString padName = Form("pad_r%d_c%d", iRow, iCol);

            TPad *pad = new TPad(padName, padName, x1, y1, x2, y2);
            pad->SetRightMargin(0);
            pad->SetLeftMargin(0);
            pad->SetBottomMargin(0);
            pad->SetTopMargin(0);
            pad->Draw();
            pad->cd();

            hist_AUL_sin2_xQ2_zPt[idx - 1]->SetTitle("");
            hist_AUL_sin2_xQ2_zPt[idx - 1]->Draw("COL");
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
                    rect_zp->SetLineWidth(1);
                    rect_zp->SetLineColor(kBlack);
                    rect_zp->Draw("same");
                    gridLines_zp.push_back(rect_zp);

                }
            }

            c_layout2->cd();
        }
    }


    // =======================================================================
    // PALETTE
    // =======================================================================

    TPad *pad_palette2 = new TPad("pad_palette2", "", 0.95, 0.1, 1, 0.9);
    pad_palette2->SetRightMargin(0.5);
    pad_palette2->Draw();
    pad_palette2->cd();

    TH2F *h_ref2 = (TH2F*) hist_AUL_sin2_xQ2_zPt[2]->Clone("h_ref_for_palette2");
    h_ref2->SetStats(0);
    h_ref2->Draw("COLZ");

    gPad->Update();

    TPaletteAxis *pal2 = (TPaletteAxis*) h_ref2->GetListOfFunctions()->FindObject("palette");
    if (pal2) {
        pal2->SetLabelSize(0.22);
        pal2->SetX1NDC(0.05);
        pal2->SetX2NDC(0.6);
        pal2->SetY1NDC(0.05);
        pal2->SetY2NDC(0.95);
    }

    pad_palette2->Modified();
    c_layout2->Update();
    c_layout2->Write();

    











    //



    //outFile.Write();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
