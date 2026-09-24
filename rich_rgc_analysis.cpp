#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
namespace fs = std::filesystem;
gROOT->SetBatch(kTRUE);
// to download the data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/RICH_data_generation/RGC/analysis/first_prod/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rgc/RICH_production/first_prod
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rga/torus-1/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rga/output_fall2018_torus-1
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/rga/torus+1/kaonp/ /Users/lorenzopolizzi/Desktop/PhD/JLAB/rga/output_fall2018_t+1_kaonp
// 

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


bool passRichSelection(double RL, double RQ, int nTot, int rich_pid, int track_pid, double rich_chi2, double Ph){

    double mA = 0.04, qA = -0.02;
    double mB = 0.04, qB = -0.02;
    double mC = 0.06, qC = -0.08;
    double mD = 0.06, qD = -0.08;
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
    else if(Ph < 1.5){
        if(track_pid == 321) return true;
        else return false;
    }
    // RICH: hadron and CLAS: kaon
    if (rich_pid != 321 && Ph <= 3) { // since we cannot trust CLAS above 3 GeV
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

    return false; // other cases
}




void rich_rgc_analysis() {
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


    // open file
    //TFile inFile("out_t-1.root", "READ");

    //string inputDir = "output_fall2018_torus+1"; 
    string inputDir = "sum22_NH3_data";
    /*
    if(rich_yes){
        if (torus == -1) {
            inputDir = "output_fall2018_t-1_kaonp";
        } else if (torus == +1){
            inputDir = "output_fall2018_t+1_kaonp";
        } 
    }
    else{
        if (torus == -1) {
            inputDir = "output_fall2018_torus-1"; 
        } else if (torus == +1){
            inputDir = "output_fall2018_torus+1";
        }
    }
    */
    TTree treeKaonP_badRich("not RICH Kaon+", "");
    TTree treeKaonP_goodRich("Kaon+ from RICH", "");
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
    if (fileCount == 0) {
        cerr << "Nessun file .root trovato in " << inputDir << endl;
        return;
    }
    // creo un output root 
    //const char* outputFile = "plot_rga_t+1.root";
    const char* outputFile = "plot_rgc_kaonp_rich.root"; 
    double torus = -1;
    /*
    if(rich_yes){
        if (torus == -1) outputFile = "plot_rich_t-1_kaonp.root";
        else if (torus == +1) outputFile = "plot_rich_t+1_kaonp.root";
    }
    else{
        if (torus == -1){
            outputFile = "plot_rga_t-1_kaonp.root";
        } else if (torus == +1){
            outputFile = "plot_rga_t+1_kaonp.root";
        }
    }
    */  
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT
    // Creiamo un TChain per ogni TTree

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
    //chainKaonP.SetBranchAddress("Beta", &kaonp_beta);
    chainKaonP.SetBranchAddress("Polarization", &kaonp_Pol);
    chainKaonP.SetBranchAddress("gamma", &kaonp_gamma);
    chainKaonP.SetBranchAddress("epsilon", &kaonp_epsilon);
    chainKaonP.SetBranchAddress("W", &kaonp_W);
    //chainKaonP.SetBranchAddress("vz", &kaonp_vz);
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
    //chainKaonP.SetBranchAddress("helicity", &kaonp_helicity);
    chainKaonP.SetBranchAddress("Mx", &kaonp_Mx);
    //chainKaonP.SetBranchAddress("PID_event", &kaonp_PID_event); // 1 se rich, 0 se event builder
    chainKaonP.SetBranchAddress("Rich_Id", &kaonp_rich_id);
    chainKaonP.SetBranchAddress("Rich_PID", &kaonp_rich_pid); // PID_rich == 321 e RQ > 0.1 con 1.2<P<8 GeV sono i tagli da richiedere per il rich
    chainKaonP.SetBranchAddress("Rich_RQ", &kaonp_rich_RQ);   // PID_rich != 321 con 1.2<P<3 GeV + HadronPID_Chi2Pid + KinematicPID_Vertex per event builder
    chainKaonP.SetBranchAddress("Rich_mass", &kaonp_bestMass);
    chainKaonP.SetBranchAddress("Rich_RL", &kaonp_rich_RL);
    chainKaonP.SetBranchAddress("Rich_nTot", &kaonp_rich_nTot);
    chainKaonP.SetBranchAddress("Rich_ch", &kaonp_rich_ch);
    chainKaonP.SetBranchAddress("Rich_chi2", &kaonp_rich_chi2);
    chainKaonP.SetBranchAddress("Rich_Mchi2", &kaonp_rich_Mchi2);
    chainKaonP.SetBranchAddress("clas_chi2", &kaonp_chi2pid);
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
    chainKaonP.SetBranchAddress("clas_pid", &kaonp_track_pid);
    // tree
    // Kaon +
    treeKaonP.Branch("E", &kaonp_E, "E/D");
    treeKaonP.Branch("px", &kaonp_px, "px/D");
    treeKaonP.Branch("py", &kaonp_py, "py/D");
    treeKaonP.Branch("pz", &kaonp_pz, "pz/D");
    treeKaonP.Branch("Mom", &kaonp_Ph, "Mom/D");
    treeKaonP.Branch("Beta", &kaonp_beta, "Beta/D");
    treeKaonP.Branch("W", &kaonp_W, "W/D");
    treeKaonP.Branch("vz", &kaonp_vz, "vz/D");
    treeKaonP.Branch("Best_mass", &kaonp_bestMass, "Best_mass/D");
    treeKaonP.Branch("Q2", &kaonp_Q2, "Q2/D");
    treeKaonP.Branch("xF", &kaonp_xF, "xF/D");
    treeKaonP.Branch("xB", &kaonp_xB, "xB/D");
    treeKaonP.Branch("y", &kaonp_y, "y/D");
    treeKaonP.Branch("z", &kaonp_z, "z/D");
    treeKaonP.Branch("Pt", &kaonp_Pt, "Pt/D");
    treeKaonP.Branch("PhT", &kaonp_PhT, "PhT/D");
    treeKaonP.Branch("Phi_Lab", &kaonp_Phi, "Phi_Lab/D");
    treeKaonP.Branch("Theta_Lab", &kaonp_Theta, "Theta_Lab/D");
    treeKaonP.Branch("Pseudorapidity", &kaonp_eta, "Pseudorapidity/D");
    treeKaonP.Branch("Phi_h", &kaonp_Phi_h, "Phi_h/D");
    treeKaonP.Branch("helicity", &kaonp_helicity, "helicity/D");
    treeKaonP.Branch("Mx", &kaonp_Mx, "Mx/D");
    treeKaonP.Branch("PID_event", &kaonp_PID_event, "PID_event/D");
    treeKaonP.Branch("Rich_Id", &kaonp_rich_id, "Rich_Id/D");
    treeKaonP.Branch("Rich_PID", &kaonp_rich_pid, "Rich_PID/D");
    treeKaonP.Branch("Rich_RL", &kaonp_rich_RL, "Rich_RL/D");
    treeKaonP.Branch("Rich_RQ", &kaonp_rich_RQ, "Rich_RQ/D");
    treeKaonP.Branch("Rich_nTot", &kaonp_rich_nTot, "Rich_nTot/D");
    treeKaonP.Branch("Rich_ch", &kaonp_rich_ch, "Rich_ch/D");
    treeKaonP.Branch("Rich_chi2", &kaonp_rich_chi2, "Rich_chi2/D");
    treeKaonP.Branch("Rich_Mchi2", &kaonp_rich_Mchi2, "Rich_Mchi2/D");
    treeKaonP.Branch("Polarization", &kaonp_Pol, "Polarization/D");
    treeKaonP.Branch("epsilon", &kaonp_epsilon, "epsilon/D");
    treeKaonP.Branch("chi2_pid", &kaonp_chi2pid, "chi2_pid/D");
    treeKaonP.Branch("EventBuilder_PID", &kaonp_track_pid, "EventBuilder_PID/D");
    treeKaonP.Branch("Rich_PMT_x", &kaonp_rich_tr1_x, "Rich_PMT_x/D");
    treeKaonP.Branch("Rich_PMT_y", &kaonp_rich_tr1_y, "Rich_PMT_y/D");
    treeKaonP.Branch("Rich_aerogel1_x", &kaonp_rich_tr2_x, "Rich_aerogel1_x/D");
    treeKaonP.Branch("Rich_aerogel1_y", &kaonp_rich_tr2_y, "Rich_aerogel1_y/D");
    treeKaonP.Branch("Rich_aerogel1_z", &kaonp_rich_tr2_z, "Rich_aerogel1_z/D");
    treeKaonP.Branch("Rich_aerogel1_path", &kaonp_rich_tr2_path, "Rich_aerogel1_path/D");
    treeKaonP.Branch("Rich_aerogel2_x", &kaonp_rich_tr3_x, "Rich_aerogel2_x/D");
    treeKaonP.Branch("Rich_aerogel2_y", &kaonp_rich_tr3_y, "Rich_aerogel2_y/D");
    treeKaonP.Branch("Rich_aerogel2_z", &kaonp_rich_tr3_z, "Rich_aerogel2_z/D");
    treeKaonP.Branch("Rich_aerogel2_path", &kaonp_rich_tr3_path, "Rich_aerogel2_path/D");
    treeKaonP.Branch("Rich_aerogel3_x", &kaonp_rich_tr4_x, "Rich_aerogel3_x/D");
    treeKaonP.Branch("Rich_aerogel3_y", &kaonp_rich_tr4_y, "Rich_aerogel3_y/D");
    treeKaonP.Branch("Rich_aerogel3_z", &kaonp_rich_tr4_z, "Rich_aerogel3_z/D");
    treeKaonP.Branch("Rich_aerogel3_path", &kaonp_rich_tr4_path, "Rich_aerogel3_path/D");
    //
    // evnt 321 & rich != 321
    treeKaonP_badRich.Branch("E", &kaonp_E, "E/D");
    treeKaonP_badRich.Branch("px", &kaonp_px, "px/D");
    treeKaonP_badRich.Branch("py", &kaonp_py, "py/D");
    treeKaonP_badRich.Branch("pz", &kaonp_pz, "pz/D");
    treeKaonP_badRich.Branch("Mom", &kaonp_Ph, "Mom/D");
    treeKaonP_badRich.Branch("Beta", &kaonp_beta, "Beta/D");
    treeKaonP_badRich.Branch("W", &kaonp_W, "W/D");
    treeKaonP_badRich.Branch("vz", &kaonp_vz, "vz/D");
    treeKaonP_badRich.Branch("Q2", &kaonp_Q2, "Q2/D");
    treeKaonP_badRich.Branch("xF", &kaonp_xF, "xF/D");
    treeKaonP_badRich.Branch("xB", &kaonp_xB, "xB/D");
    treeKaonP_badRich.Branch("y", &kaonp_y, "y/D");
    treeKaonP_badRich.Branch("z", &kaonp_z, "z/D");
    treeKaonP_badRich.Branch("Pt", &kaonp_Pt, "Pt/D");
    treeKaonP_badRich.Branch("PhT", &kaonp_PhT, "PhT/D");
    treeKaonP_badRich.Branch("Phi_Lab", &kaonp_Phi, "Phi_Lab/D");
    treeKaonP_badRich.Branch("Theta_Lab", &kaonp_Theta, "Theta_Lab/D");
    treeKaonP_badRich.Branch("Pseudorapidity", &kaonp_eta, "Pseudorapidity/D");
    treeKaonP_badRich.Branch("Phi_h", &kaonp_Phi_h, "Phi_h/D");
    treeKaonP_badRich.Branch("helicity", &kaonp_helicity, "helicity/D");
    treeKaonP_badRich.Branch("Mx", &kaonp_Mx, "Mx/D");
    treeKaonP_badRich.Branch("PID_event", &kaonp_PID_event, "PID_event/D");
    treeKaonP_badRich.Branch("Rich_Id", &kaonp_rich_id, "Rich_Id/D");
    treeKaonP_badRich.Branch("Rich_PID", &kaonp_rich_pid, "Rich_PID/D");
    treeKaonP_badRich.Branch("Rich_RL", &kaonp_rich_RL, "Rich_RL/D");
    treeKaonP_badRich.Branch("Rich_RQ", &kaonp_rich_RQ, "Rich_RQ/D");
    treeKaonP_badRich.Branch("Rich_nTot", &kaonp_rich_nTot, "Rich_nTot/D");
    treeKaonP_badRich.Branch("Rich_ch", &kaonp_rich_ch, "Rich_ch/D");
    treeKaonP_badRich.Branch("Rich_chi2", &kaonp_rich_chi2, "Rich_chi2/D");
    treeKaonP_badRich.Branch("Polarization", &kaonp_Pol, "Polarization/D");
    treeKaonP_badRich.Branch("epsilon", &kaonp_epsilon, "epsilon/D");
    treeKaonP_badRich.Branch("chi2_pid", &kaonp_chi2pid, "chi2_pid/D");
    treeKaonP_badRich.Branch("EventBuilder_PID", &kaonp_track_pid, "EventBuilder_PID/D");
    treeKaonP_badRich.Branch("Rich_PMT_x", &kaonp_rich_tr1_x, "Rich_PMT_x/D");
    treeKaonP_badRich.Branch("Rich_PMT_y", &kaonp_rich_tr1_y, "Rich_PMT_y/D");
    treeKaonP_badRich.Branch("Rich_aerogel1_x", &kaonp_rich_tr2_x, "Rich_aerogel1_x/D");
    treeKaonP_badRich.Branch("Rich_aerogel1_y", &kaonp_rich_tr2_y, "Rich_aerogel1_y/D");
    treeKaonP_badRich.Branch("Rich_aerogel1_z", &kaonp_rich_tr2_z, "Rich_aerogel1_z/D");
    treeKaonP_badRich.Branch("Rich_aerogel1_path", &kaonp_rich_tr2_path, "Rich_aerogel1_path/D");
    treeKaonP_badRich.Branch("Rich_aerogel2_x", &kaonp_rich_tr3_x, "Rich_aerogel2_x/D");
    treeKaonP_badRich.Branch("Rich_aerogel2_y", &kaonp_rich_tr3_y, "Rich_aerogel2_y/D");
    treeKaonP_badRich.Branch("Rich_aerogel2_z", &kaonp_rich_tr3_z, "Rich_aerogel2_z/D");
    treeKaonP_badRich.Branch("Rich_aerogel2_path", &kaonp_rich_tr3_path, "Rich_aerogel2_path/D");
    treeKaonP_badRich.Branch("Rich_aerogel3_x", &kaonp_rich_tr4_x, "Rich_aerogel3_x/D");
    treeKaonP_badRich.Branch("Rich_aerogel3_y", &kaonp_rich_tr4_y, "Rich_aerogel3_y/D");
    treeKaonP_badRich.Branch("Rich_aerogel3_z", &kaonp_rich_tr4_z, "Rich_aerogel3_z/D");
    treeKaonP_badRich.Branch("Rich_aerogel3_path", &kaonp_rich_tr4_path, "Rich_aerogel3_path/D");

    // evnt != 321 & rich 321
    treeKaonP_goodRich.Branch("E", &kaonp_E, "E/D");
    treeKaonP_goodRich.Branch("px", &kaonp_px, "px/D");
    treeKaonP_goodRich.Branch("py", &kaonp_py, "py/D");
    treeKaonP_goodRich.Branch("pz", &kaonp_pz, "pz/D");
    treeKaonP_goodRich.Branch("Mom", &kaonp_Ph, "Mom/D");
    treeKaonP_goodRich.Branch("Beta", &kaonp_beta, "Beta/D");
    treeKaonP_goodRich.Branch("W", &kaonp_W, "W/D");
    treeKaonP_goodRich.Branch("vz", &kaonp_vz, "vz/D");
    treeKaonP_goodRich.Branch("Q2", &kaonp_Q2, "Q2/D");
    treeKaonP_goodRich.Branch("xF", &kaonp_xF, "xF/D");
    treeKaonP_goodRich.Branch("xB", &kaonp_xB, "xB/D");
    treeKaonP_goodRich.Branch("y", &kaonp_y, "y/D");
    treeKaonP_goodRich.Branch("z", &kaonp_z, "z/D");
    treeKaonP_goodRich.Branch("Pt", &kaonp_Pt, "Pt/D");
    treeKaonP_goodRich.Branch("PhT", &kaonp_PhT, "PhT/D");
    treeKaonP_goodRich.Branch("Phi_Lab", &kaonp_Phi, "Phi_Lab/D");
    treeKaonP_goodRich.Branch("Theta_Lab", &kaonp_Theta, "Theta_Lab/D");
    treeKaonP_goodRich.Branch("Pseudorapidity", &kaonp_eta, "Pseudorapidity/D");
    treeKaonP_goodRich.Branch("Phi_h", &kaonp_Phi_h, "Phi_h/D");
    treeKaonP_goodRich.Branch("helicity", &kaonp_helicity, "helicity/D");
    treeKaonP_goodRich.Branch("Mx", &kaonp_Mx, "Mx/D");
    treeKaonP_goodRich.Branch("PID_event", &kaonp_PID_event, "PID_event/D");
    treeKaonP_goodRich.Branch("Rich_Id", &kaonp_rich_id, "Rich_Id/D");
    treeKaonP_goodRich.Branch("Rich_PID", &kaonp_rich_pid, "Rich_PID/D");
    treeKaonP_goodRich.Branch("Rich_RL", &kaonp_rich_RL, "Rich_RL/D");
    treeKaonP_goodRich.Branch("Rich_RQ", &kaonp_rich_RQ, "Rich_RQ/D");
    treeKaonP_goodRich.Branch("Rich_nTot", &kaonp_rich_nTot, "Rich_nTot/D");
    treeKaonP_goodRich.Branch("Rich_ch", &kaonp_rich_ch, "Rich_ch/D");
    treeKaonP_goodRich.Branch("Rich_chi2", &kaonp_rich_chi2, "Rich_chi2/D");
    treeKaonP_goodRich.Branch("Polarization", &kaonp_Pol, "Polarization/D");
    treeKaonP_goodRich.Branch("epsilon", &kaonp_epsilon, "epsilon/D");
    treeKaonP_goodRich.Branch("chi2_pid", &kaonp_chi2pid, "chi2_pid/D");
    treeKaonP_goodRich.Branch("EventBuilder_PID", &kaonp_track_pid, "EventBuilder_PID/D");
    treeKaonP_goodRich.Branch("Rich_PMT_x", &kaonp_rich_tr1_x, "Rich_PMT_x/D");
    treeKaonP_goodRich.Branch("Rich_PMT_y", &kaonp_rich_tr1_y, "Rich_PMT_y/D");
    treeKaonP_goodRich.Branch("Rich_aerogel1_x", &kaonp_rich_tr2_x, "Rich_aerogel1_x/D");
    treeKaonP_goodRich.Branch("Rich_aerogel1_y", &kaonp_rich_tr2_y, "Rich_aerogel1_y/D");
    treeKaonP_goodRich.Branch("Rich_aerogel1_z", &kaonp_rich_tr2_z, "Rich_aerogel1_z/D");
    treeKaonP_goodRich.Branch("Rich_aerogel1_path", &kaonp_rich_tr2_path, "Rich_aerogel1_path/D");
    treeKaonP_goodRich.Branch("Rich_aerogel2_x", &kaonp_rich_tr3_x, "Rich_aerogel2_x/D");
    treeKaonP_goodRich.Branch("Rich_aerogel2_y", &kaonp_rich_tr3_y, "Rich_aerogel2_y/D");
    treeKaonP_goodRich.Branch("Rich_aerogel2_z", &kaonp_rich_tr3_z, "Rich_aerogel2_z/D");
    treeKaonP_goodRich.Branch("Rich_aerogel2_path", &kaonp_rich_tr3_path, "Rich_aerogel2_path/D");
    treeKaonP_goodRich.Branch("Rich_aerogel3_x", &kaonp_rich_tr4_x, "Rich_aerogel3_x/D");
    treeKaonP_goodRich.Branch("Rich_aerogel3_y", &kaonp_rich_tr4_y, "Rich_aerogel3_y/D");
    treeKaonP_goodRich.Branch("Rich_aerogel3_z", &kaonp_rich_tr4_z, "Rich_aerogel3_z/D");
    treeKaonP_goodRich.Branch("Rich_aerogel3_path", &kaonp_rich_tr4_path, "Rich_aerogel3_path/D");
    // KAON+ PLOT
    //dirKaonp->cd();
    double bin = 300;
    double chi_min = 0.1;
    double chi_max = 2000;
    auto make_bins = [](int bins, double min, double max) {
        return CreateLogBinning(bins, min, max);
      };
    const auto log_chi2 = make_bins(bin, chi_min, chi_max);
    // Mom
    TH1D kp_evnt_chi2 ("_evnt_chi2", "#chi^{2} EventBuilder PID | only EventBuilder | 1.2 < Mom < 8 GeV ; #chi^{2}; count", 300, -8, 8);
    TH1D kp_m ("_best_mass", "m extracted from #beta | 1.2 < Mom < 8 GeV ; m [GeV]; count", 300, 0, 1);
    TH1D kp_deltaB ("_delta_beta", "#beta_{meas} - #beta_{th} | 1.2 < Mom < 8 GeV ; #Delta_{#beta}; count", 300, -0.05, 0.05);
    TH1D kp_Mx ("_missing_mass", "missing mass | 1.2 < Mom < 8 GeV ; M_{x} [GeV]; count", 300, 0.4, 2);
    //TH1D kp_rich_chi2 ("_rich_chi2", "#chi^{2} RICH PID ; #chi^{2}; count", 300, -3, 3);
    TH2D kp_MomVsPhT ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  K+  | with EventBuilder + RICH ; P_{hT} [GeV]; Mom [GeV]", 400, 0, 1.2, 400, 1, 8);
    TH2D kp_MomVsXb ("_MomVsXb", "Correlation Mom vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; Mom [GeV]", 400, 0, 0.8, 400, 1, 8);
    TH2D kp_MomVsXf ("_MomVsXf", "Correlation Mom vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; Mom [GeV]", 400, 0, 0.5, 400, 1, 8);
    TH2D kp_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  K+  | with EventBuilder + RICH ; z; Mom [GeV]", 400, 0.2, 0.9, 400, 1, 8);
    TH2D kp_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  K+  | with EventBuilder + RICH ; y; Mom [GeV]", 400, 0.2, 0.75, 400, 1, 8);
    TH2D kp_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  K+  | with EventBuilder + RICH ; Eta; Mom [GeV]", 400, 1.5, 3.0, 400, 1, 8);
    TH2D kp_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  K+  | with EventBuilder + RICH ; Mom [GeV]; Theta [Rad]", 400, 1, 8, 400, 0.09, 0.4);
    TH2D kp_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Mom [GeV]", 400, -TMath::Pi(), TMath::Pi(), 400, 1, 8);
    TH2D kp_MomVsMx ("_MomVsMx", "correlation Mom vs M^{2}_{x} | K+ |; Mom [GeV]; M^{2}_{x} [GeV^{2}]", 400, 1, 8, 400, 2.45, 11.5);
    TH2D kp_MomVsBeta ("_MomVsBeta", "correlation Mom vs #beta | K+ |; Mom [GeV]; #beta", 400, 1, 8, 400, 0.88, 1.02);
    TH2D kp_MomVsBeta_ToF1 ("_MomVsBeta_ToF_goodChi2", "correlation Mom vs #beta | K+ | ToF #Chi^{2} < 3; Mom [GeV]; #beta", 400, 1, 8, 400, 0.88, 1.02);
    TH2D kp_MomVsBeta_ToF2 ("_MomVsBeta_ToF_badChi2", "correlation Mom vs #beta | K+ | ToF #Chi^{2} > 3; Mom [GeV]; #beta", 400, 1, 8, 400, 0.88, 1.02);
    TH2D kp_MomVsCh ("_MomVsCh", "correlation Mom vs #theta_{ch} | K+ |; Mom [GeV]; #theta_{ch} [rad]", 400, 1, 8, 400, 0.0, 0.35);
    TH2D kp_MomVsMass_Rich ("_MomVsMass_RICH", "correlation Mom vs Mass | K+ |; Mom [GeV]; Mass [GeV]", 400, 1, 8, 400, -0.2, 1.2);
    TH2D kp_MomVsCh_RICH1 ("_MomVsCh_RICH_goodChi2", "correlation Mom vs #theta_{ch} | K+ | RICH #Chi^{2} < 4; Mom [GeV]; #theta_{ch} [rad]", 400, 1, 8, 400, 0.0, 0.35);
    TH2D kp_MomVsCh_RICH2 ("_MomVsCh_RICH_badChi2", "correlation Mom vs #theta_ch | K+ | RICH #Chi^{2} > 4; Mom [GeV]; #theta_{ch} [rad]", 400, 1, 8, 400, 0.0, 0.35);
    // Q2
    TH2D kp_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; Q^{2} [GeV^{2}]", 400, 0, 0.8, 400, 1, 10);
    TH2D kp_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; Q^{2} [GeV^{2}]", 400, 0, 0.6, 400, 1, 10);
    TH2D kp_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  K+  | with EventBuilder + RICH ; Mom [GeV]; Q^{2} [GeV^{2}]", 400, 0.9, 6, 400, 1, 10);
    TH2D kp_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  K+  | with EventBuilder + RICH ; P_{hT} [GeV]; Q^{2} [GeV^{2}]", 400, 0, 1.2, 400, 1, 10);
    TH2D kp_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs Z  |  K+  | with EventBuilder + RICH ; z; Q^{2} [GeV^{2}]", 400, 0.2, 0.9, 400, 1, 10);
    TH2D kp_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  K+  | with EventBuilder + RICH ; y; Q^{2} [GeV^{2}]", 400, 0.2, 0.75, 400, 0.9, 10);
    TH2D kp_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  K+  | with EventBuilder + RICH ; Eta; Q^{2} [GeV^{2}]", 400, 1.5, 3.0, 400, 1, 10);
    TH2D kp_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", 400, -TMath::Pi(), TMath::Pi(), 400, 1, 10);
    // PhT
    TH2D kp_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs Z  |  K+  | with EventBuilder + RICH ; z; P_{hT} [GeV]", 400, 0.2, 1, 400, 0, 1.2);
    TH2D kp_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; P_{hT} [GeV]", 400, 0, 0.8, 400, 0, 1.2);
    TH2D kp_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  K+  | with EventBuilder + RICH ; Eta; P_{hT} [GeV]", 400, 1.5, 3.0, 400, 0, 1.2);
    TH2D kp_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; P_{hT} [GeV]", 400, -TMath::Pi(), TMath::Pi(), 400, 0, 1.2);
    // Z
    TH2D kp_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; z", 400, 0, 0.8, 400, 0.2, 0.9);
    TH2D kp_zVsXf ("_zVsXf", "Correlation Z vs x_{F}  |  K+  | with EventBuilder + RICH ; x_{F}; z", 400, 0, 0.6, 400, 0.2, 0.9);
    TH2D kp_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  K+  | with EventBuilder + RICH ; Eta; z", 400, 1.5, 3.0, 400, 0.2, 0.9);
    TH2D kp_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; z", 400, -TMath::Pi(), TMath::Pi(), 400, 0.2, 0.9);
    //
    TH2D kp_xBvsY ("_xBvsY", "Correlation y vs x_{B}  |  K+  | with EventBuilder + RICH ; x_{B}; y", 400, 0, 0.8, 400, 0.25, 0.75);
    // Angles
    TH2D kp_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  K+  | with EventBuilder + RICH ; #Phi_{h} [Rad]; Theta [Rad]", 400, -TMath::Pi(), TMath::Pi(), 400, 0.1, 0.4);
    TH2D kp_ThetaVsPhi_Lab ("_ThetaVsPhi_Lab", "Correlation #theta vs #Phi_{Lab} | K+ | with EventBuilder + RICH ; #Phi_{Lab} [Rad]; #theta [Rad]", 400, -TMath::Pi(), TMath::Pi(), 400, 0.05, 0.4);
    //
    TH2D empty_canvas ("____________________", "empty space;", 1, 0, 1, 1, 0, 1);
    TH1D empty_canvas_grich ("_________ good RICH ___________", "empty space;", 1, 0, 1);
    TH1D kp_evnt_pid_grich ("_evnt_PID_good_rich", "EventBuilder PID | K+ for RICH; PID", 300, 0, 2300);
    TH1D kp_evnt_chi2_grich ("_evnt_chi2_good_rich", "EventBuilder #chi^{2} | K+ for RICH; #chi^{2}", 300, -5, 5);
    TH1D kp_rich_chi2_grich ("_rich_chi2_good_rich", "RICH #chi^{2} | K+ for RICH; #chi^{2}", 300, 0, 14);
    TH1D kp_rich_Mchi2_grich ("_rich_Mchi2_good_rich", "RICH M#chi^{2} | K+ for RICH; M#chi^{2}", 300, 0, 14);
    TH1D kp_rich_RL_grich ("_rich_RL_good_rich", "RICH RL | K+ for RICH; RL", 300, 0, 8);
    TH1D kp_rich_RQ_grich ("_rich_RQ_good_rich", "RICH RQ | K+ for RICH; RQ", 300, 0, 0.7);
    TH1D kp_rich_nTot_grich ("_rich_nTot_good_rich", "RICH nTot | K+ for RICH; nTot", 80, 0, 40);
    TH1D kp_rich_ch_grich ("_rich_ch_good_rich", "RICH Cherenkov angle | K+ for RICH; #theta [Rad]", 300, 0.05, 0.35);
    TH2D kp_MomVsBeta_grich ("_MomVsBeta_good_rich", "correlation Mom vs #beta | K+ for RICH; Mom [GeV]; #beta", 300, 1, 8, 300, 0.88, 1.02);
    TH2D kp_MomVsCh_grich ("_MomVsCh_good_rich", "correlation Mom vs #theta_{ch} | K+ for RICH; Mom [GeV]; #theta_{ch} [rad]", 300, 1, 8, 300, 0.0, 0.35);
    TH2D kp_MomVsRL_grich ("_MomVsRL_good_rich", "correlation Mom vs RL | K+ for RICH; RL; Mom [GeV]", 300, 0, 8, 300, 1, 8);
    TH2D kp_ChvsChi2_grich ("_ChVsChi2_good_rich", "correlation #theta_{ch} vs #Chi^{2} | K+ for RICH; #Chi^{2}; #theta_{ch} [rad]", 300, log_chi2.data(), 300, 0, 0.35);
    TH2D kp_nTotvsChi2_grich ("_nTotVsChi2_good_rich", "correlation n_{tot} vs #Chi^{2} | K+ for RICH; #Chi^{2}; n_{tot}", 300, log_chi2.data(), 40, 0, 40);
    TH2D kp_MomvsChi2_grich ("_MomVsChi2_good_rich", "correlation Mom vs #Chi^{2} | K+ for RICH; #Chi^{2}; Mom [GeV]", 300, log_chi2.data(), 300, 1, 8);
    TH2D kp_RLvsChi2_grich ("_RLvsChi2_good_rich", "correlation RL vs #Chi^{2} | K+ for RICH; #Chi^{2}; RL", 300, log_chi2.data(), 300, 0, 8);
    TH2D kp_RLvsChi2_low_grich ("_RLvsChi2_low_good_rich", "correlation RL vs #Chi^{2} with #Chi^{2} < 80| K+ for RICH; #Chi^{2}; RL", 300, 0, 50, 300, 0, 8);
    TH2D kp_RLvsChi2_high_grich ("_RLvsChi2_high_good_rich", "correlation RL vs #Chi^{2} with #Chi^{2} > 50| K+ for RICH; #Chi^{2}; RL", 300, 50, 1200, 300, 0, 8);
    TH2D kp_RLvsRQ_grich ("_RLvsRQ_good_rich", "correlation RL vs RQ | K+ for RICH; RL; RQ", 300, 0, 8, 300, 0.0, 0.7);
    TH2D kp_RLvsCh_grich ("_RLvsCh_good_rich", "correlation RL vs #theta_{ch} | K+ for RICH; RL; #theta_{ch}", 300, 0, 8, 300, 0.0, 0.35);
    TH2D kp_RQvsCh_grich ("_RQvsCh_good_rich", "correlation RQ vs #theta_{ch} | K+ for RICH; RQ; #theta_{ch}", 300, 0.0, 0.7, 300, 0.0, 0.35);
    //TH2D kp_MomVsnTot_grich ("_MomVsnTot_good_rich", "correlation Mom vs #n_{Tot} | K+ for RICH; Mom [GeV]; #n_{Tot} [rad]", 300, 1, 8, 30, 0, 30);
    //
    TH1D empty_canvas_brich ("_________ bad RICH ___________", "empty space;", 1, 0, 1);
    TH1D kp_rich_pid_brich ("_rich_PID_bad_rich", "RICH PID | not K+ for RICH; PID", 300, 0, 2300);
    TH1D kp_evnt_chi2_brich ("_evnt_chi2_bad_rich", "EventBuilder #chi^{2} | not K+ for RICH; #chi^{2}", 300, -5, 5);
    TH1D kp_rich_Mchi2_brich ("_rich_Mchi2_bad_rich", "RICH M#chi^{2} | not K+ for RICH; M#chi^{2}", 300, 0, 14);
    TH1D kp_rich_chi2_brich ("_rich_chi2_bad_rich", "RICH #chi^{2} | not K+ for RICH; #chi^{2}", 300, 0, 14);
    TH1D kp_rich_RL_brich ("_rich_RL_bad_rich", "RICH RL | not K+ for RICH; RL", 300, 0, 8);
    TH1D kp_rich_RQ_brich ("_rich_RQ_bad_rich", "RICH RQ | not K+ for RICH; RQ", 300, 0, 0.7);
    TH1D kp_rich_nTot_brich ("_rich_nTot_bad_rich", "RICH nTot | not K+ for RICH; nTot", 80, 0, 40);
    TH1D kp_rich_ch_brich ("_rich_ch_bad_rich", "RICH Cherenkov angle | not K+ for RICH; #theta [Rad]", 300, 0.05, 0.35);
    TH2D kp_MomVsBeta_brich ("_MomVsBeta_bad_rich", "correlation Mom vs #beta | not K+ for RICH; Mom [GeV]; #beta", 300, 1, 8, 300, 0.88, 1.02);
    TH2D kp_MomVsCh_brich ("_MomVsCh_bad_rich", "correlation Mom vs #theta_{ch} | not K+ for RICH; Mom [GeV]; #theta_{ch} [rad]", 300, 1, 8, 300, 0.0, 0.35);
    TH2D kp_MomVsRL_brich ("_MomVsRL_bad_rich", "correlation Mom vs RL | not K+ for RICH; RL; Mom [GeV]", 200, 0, 8, 200, 1, 8);
    TH2D kp_RLvsChi2_brich ("_RLvsChi2_bad_rich", "correlation RL vs #Chi^{2} | not K+ for RICH; #Chi^{2}; RL", 300, log_chi2.data(), 300, 0, 8);
    TH2D kp_ChvsChi2_brich ("_ChVsChi2_bad_rich", "correlation #theta_{ch} vs #Chi^{2} | not K+ for RICH; #Chi^{2}; #theta_{ch} [rad]", 300, log_chi2.data(), 300, 0, 0.35);
    TH2D kp_nTotvsChi2_brich ("_nTotVsChi2_bad_rich", "correlation n_{tot} vs #Chi^{2} | not K+ for RICH; #Chi^{2}; n_{tot}", 300, log_chi2.data(), 40, 0, 40);
    TH2D kp_MomvsChi2_brich ("_MomVsChi2_bad_rich", "correlation Mom vs #Chi^{2} | not K+ for RICH; #Chi^{2}; Mom [GeV]", 300, log_chi2.data(), 300, 1, 8);
    TH2D kp_RLvsChi2_low_brich ("_RLvsChi2_low_bad_rich", "correlation RL vs #Chi^{2} with #Chi^{2} < 80| not K+ for RICH; #Chi^{2}; RL", 300, 0, 50, 300, 0, 8);
    TH2D kp_RLvsChi2_high_brich ("_RLvsChi2_high_bad_rich", "correlation RL vs #Chi^{2} with #Chi^{2} > 50| not K+ for RICH; #Chi^{2}; RL", 300, 50, 1200, 300, 0, 8);
    TH2D kp_RLvsRQ_brich ("_RLvsRQ_bad_rich", "correlation RL vs RQ | not K+ for RICH; RL; RQ", 300, 0, 8, 300, 0.0, 0.7);
    TH2D kp_RLvsCh_brich ("_RLvsCh_bad_rich", "correlation RL vs #theta_{ch} | not K+ for RICH; RL; #theta_{ch}", 300, 0, 8, 300, 0.0, 0.35);
    TH2D kp_RQvsCh_brich ("_RQvsCh_bad_rich", "correlation RQ vs #theta_{ch} | not K+ for RICH; RQ; #theta_{ch}", 300, 0.0, 0.7, 300, 0.0, 0.35);
    TH2D kp_Blob1VsCh_brich ("_Blob1VsCh_bad_rich", "correlation Blob1 vs #theta_{ch} | not K+ for RICH; RL; #theta_{ch} [rad]", 200, 0, 2.5, 200, 0.0, 0.35);
    TH2D kp_Blob2VsCh_brich ("_Blob2VsCh_bad_rich", "correlation Blob2 vs #theta_{ch} | not K+ for RICH; RL; #theta_{ch} [rad]", 200, 2.5, 8, 200, 0.0, 0.35);
    //
    TH1D empty_canvas_mrich ("_________ maybe RICH ___________", "empty space;", 1, 0, 1);
    TH1D kp_evnt_pid_mrich ("_evnt_PID_maybe_rich", "EventBuilder PID | K+ only for RICH; PID", 4600, -2300, 2300);
    TH1D kp_evnt_chi2_mrich ("_evnt_chi2_maybe_rich", "EventBuilder #chi^{2} | K+ only for RICH; #chi^{2}", 300, -5, 5);
    TH1D kp_rich_chi2_mrich ("_rich_chi2_maybe_rich", "RICH #chi^{2} | K+ only for RICH; #chi^{2}", 300, 0, 14);
    TH1D kp_rich_Mchi2_mrich ("_rich_Mchi2_maybe_rich", "RICH M#chi^{2} | K+ only for RICH; M#chi^{2}", 300, 0, 14);
    TH1D kp_rich_RL_mrich ("_rich_RL_maybe_rich", "RICH RL | K+ only for RICH; RL", 300, 0, 8);
    TH1D kp_rich_RQ_mrich ("_rich_RQ_maybe_rich", "RICH RQ | K+ only for RICH; RQ", 300, 0, 0.7);
    TH1D kp_rich_nTot_mrich ("_rich_nTot_maybe_rich", "RICH nTot | K+ only for RICH; nTot", 80, 0, 40);
    TH1D kp_rich_ch_mrich ("_rich_ch_maybe_rich", "RICH Cherenkov angle | K+ only for RICH; #theta [Rad]", 300, 0.05, 0.35);
    TH2D kp_MomVsBeta_mrich ("_MomVsBeta_maybe_rich", "correlation Mom vs #beta | K+ only for RICH; Mom [GeV]; #beta", 300, 1, 8, 300, 0.88, 1.02);
    TH2D kp_MomVsCh_mrich ("_MomVsCh_maybe_rich", "correlation Mom vs #theta_{ch} | K+ only for RICH; Mom [GeV]; #theta_{ch} [rad]", 300, 1, 8, 300, 0.0, 0.35);
    TH2D kp_MomVsnTot_mrich ("_MomVsnTot_maybe_rich", "correlation Mom vs #n_{Tot} | K+ only for RICH; Mom [GeV]; #n_{Tot} [rad]", 30, 1, 8, 30, 5, 35);
    TH2D kp_MomVsRL_mrich ("_MomVsRL_maybe_rich", "correlation Mom vs RL | K+ only for RICH; RL; Mom [GeV]", 300, 0, 8, 300, 1, 8);
    TH2D kp_ChvsChi2_mrich ("_ChVsChi2_maybe_rich", "correlation #theta_{ch} vs #Chi^{2} | K+ only for RICH; #Chi^{2}; #theta_{ch} [rad]", 300, log_chi2.data(), 300, 0, 0.35);
    TH2D kp_nTotvsChi2_mrich ("_nTotVsChi2_maybe_rich", "correlation n_{tot} vs #Chi^{2} | K+ only for RICH; #Chi^{2}; n_{tot}", 300, log_chi2.data(), 40, 0, 40);
    TH2D kp_RLvsChi2_mrich ("_RLvsChi2_maybe_rich", "correlation RL vs #Chi^{2} | K+ only for RICH; #Chi^{2}; RL", 300,log_chi2.data(), 300, 0, 8);
    TH2D kp_MomvsChi2_mrich ("_MomVsChi2_maybe_rich", "correlation Mom vs #Chi^{2} | K+ only for RICH; #Chi^{2}; Mom [GeV]", 300, log_chi2.data(), 300, 1, 8);
    TH2D kp_RLvsChi2_low_mrich ("_RLvsChi2_low_maybe_rich", "correlation RL vs #Chi^{2} with #Chi^{2} < 80| K+ only for RICH; #Chi^{2}; RL", 300, 0, 50, 300, 0, 8);
    TH2D kp_RLvsChi2_high_mrich ("_RLvsChi2_high_maybe_rich", "correlation RL vs #Chi^{2} with #Chi^{2} > 50| K+ only for RICH; #Chi^{2}; RL", 300, 50, 1200, 300, 0, 8);
    TH2D kp_RLvsRQ_mrich ("_RLvsRQ_maybe_rich", "correlation RL vs RQ | K+ only for RICH; RL; RQ", 300, 0, 8, 300, 0.0, 0.7);
    TH2D kp_RLvsCh_mrich ("_RLvsCh_maybe_rich", "correlation RL vs #theta_{ch} | K+ only for RICH; RL; #theta_{ch}", 300, 0, 8, 300, 0.0, 0.35);
    TH2D kp_RQvsCh_mrich ("_RQvsCh_maybe_rich", "correlation RQ vs #theta_{ch} | K+ only for RICH; RQ; #theta_{ch}", 300, 0.0, 0.7, 300, 0.0, 0.35);
    TH2D empty_canvas2 ("______________________", "empty space;", 1, 0, 1, 1, 0, 1);


    //
    TH1D kp_rich_chi2 ("_rich_chi2", "rich chi2", 300, 0, 1200);
    TH1D kp_rich_Mchi2 ("_rich_Mchi2", "rich Mchi2", 300, 0, 14);
    TH1D kp_rich_traj1_x ("_rich_traj1_x", "traj 1 - x; x [cm]", 50, -160, 5);
    TH1D kp_rich_traj1_y ("_rich_traj1_y", "traj 1 - y; y [cm]", 50, -40, 40);
    TH2D kp_rich_pmt_xy ("rich_pmt_xy", "PMT RICH 4th sector, xy plane | t+1 & K+; x [cm]; y [cm]", 500, -170, -50, 500, -70, 70);
    TH1D kp_rich_traj2_x ("_rich_traj2_x", "traj 2 - x; x [cm]", 50, -120, -80);
    TH1D kp_rich_traj2_y ("_rich_traj2_y", "traj 2 - y; y [cm]", 50, -1, 1);
    TH1D kp_rich_traj2_z ("_rich_traj2_z", "traj 2 - z; z [cm]", 50, 490, 540);
    TH1D kp_rich_traj2_path ("_rich_traj2_path", "traj 2 - path; path [cm]", 50, 500, 620);
    TH1D kp_rich_traj3_x ("_rich_traj3_x", "traj 3 - x; x [cm]", 50, -300, -100);
    TH1D kp_rich_traj3_y ("_rich_traj3_y", "traj 3 - y; y [cm]", 50, -100, 100);
    TH1D kp_rich_traj3_z ("_rich_traj3_z", "traj 3 - z; z [cm]", 50, 510, 570);
    TH1D kp_rich_traj3_path ("_rich_traj3_path", "traj 3 - path; path [cm]", 50, 510, 620);
    TH1D kp_rich_traj4_x ("_rich_traj4_x", "traj 4 - x; x [cm]", 50, -250, -150);
    TH1D kp_rich_traj4_y ("_rich_traj4_y", "traj 4 - y; y [cm]", 50, -110, 100);
    TH1D kp_rich_traj4_z ("_rich_traj4_z", "traj 4 - z; z [cm]", 50, 500, 540);
    TH1D kp_rich_traj4_path ("_rich_traj4_path", "traj 4 - path; path [cm]", 50, 510, 620);
    TH1D kp_rich_aerogel_x ("rich_aerogel_x", "traj aerogel 34; x [cm]", 50, -250, -100);
    TH1D kp_rich_aerogel_y ("rich_aerogel_y", "traj aerogel 34; y [cm]", 50, -110, 100);
    TH2D kp_rich_aerogel_xy ("rich_aerogel_xy", "Aerogel RICH 4th sector, xy plane | t-1 & K+; x [cm]; y [cm]", 400, -230, -70, 400, -100, 100);
    TH2D kp_rich_aerogel_l1_xy ("rich_aerogel_l1_xy", "Aerogel RICH 4th sector, layer 1, xy plane | t-1 & K+; x [cm]; y [cm]", 500, -170, -70, 500, -80, 80);
    TH2D kp_rich_aerogel_l2_xy ("rich_aerogel_l2_xy", "Aerogel RICH 4th sector, layer 2, xy plane | t-1 & K+; x [cm]; y [cm]", 500, -170, -70, 500, -80, 80);
    TH2D kp_rich_aerogel_xp ("rich_aerogel_xp", "Aerogel RICH 4th sector, x-path plane | t-1 & K+; x [cm]; path [cm]", 500, -170, -70, 500, 490, 640);
    TH2D kp_rich_aerogel_yp ("rich_aerogel_yp", "Aerogel RICH 4th sector, y-path plane | t-1 & K+; y [cm]; path [cm]", 500, -100, 70, 500, 490, 640);
    TH2D kp_rich_aerogel_xz ("rich_aerogel_xz", "Aerogel RICH 4th sector, xz plane | t-1 & K+; x [cm]; z [cm]", 100, -170, -70, 300, 100, 580);

    TH1D kp_Mx_tof ("_missing_mass_ToF", "missing mass | 1.2 < Mom < 8 GeV | ToF ; M_{x} [GeV]; count", 300, 0.4, 2);
    TH1D kp_Mx_rich ("_missing_mass_RICH", "missing mass | 1.2 < Mom < 8 GeV | RICH ; M_{x} [GeV]; count", 300, 0.4, 2);
    // ORA RIEMPI I GRAFICI
    // Kaon+
    Long64_t nEntries_kp = chainKaonP.GetEntries();
    for (Long64_t i = 0; i < nEntries_kp; i++) {
        // RICH cut
        //if (electron_ass != 321) continue;
        chainKaonP.GetEntry(i);
        //if(kaonp_rich_nTot <= 4) continue;
        //if(kaonp_Ph <= 3) continue;
        //if(kaonp_rich_chi2 > 50 || kaonp_Ph > 3) continue;
        // || kaonp_rich_Mchi2 == 0
        if (kaonp_track_pid == -11) continue;
        //if((kaonp_rich_RQ >= 0.1 && kaonp_rich_pid == 321 && kaonp_Ph >= 1.2 && kaonp_Ph < 8) || (kaonp_rich_pid != 321 && kaonp_Ph >= 1.2 && kaonp_Ph < 3 && std::abs(kaonp_chi2pid)<3 && KinematicPID_Vertex(kaonp_vz, torus, 321))){
        if (std::abs(kaonp_chi2pid) < 8 && KinematicPID_Vertex(kaonp_vz, torus, 321) && kaonp_rich_pid != 321){
            //treeKaonP_badRich.Fill();
            double m = 0.1625;
            double q = -0.4875;
            double y_limit = m*kaonp_rich_RL + q;
            double x_limit = (kaonp_rich_RQ - q)/m; // kaonp_rich_RQ < y_limit && kaonp_rich_RL > x_limit
            //if(kaonp_rich_RQ >= y_limit) continue;
            if(kaonp_rich_nTot == 0) continue;
            if(passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph)){
                kp_evnt_chi2_brich.Fill(kaonp_chi2pid);
                kp_rich_chi2_brich.Fill(kaonp_rich_chi2);
                kp_rich_Mchi2_brich.Fill(kaonp_rich_Mchi2);
                kp_rich_RL_brich.Fill(kaonp_rich_RL);
                kp_rich_RQ_brich.Fill(kaonp_rich_RQ);
                kp_rich_nTot_brich.Fill(kaonp_rich_nTot);
                kp_rich_ch_brich.Fill(kaonp_rich_ch);
                kp_rich_pid_brich.Fill(kaonp_rich_pid);
                kp_MomVsBeta_brich.Fill(kaonp_Ph, kaonp_beta);
                kp_MomVsCh_brich.Fill(kaonp_Ph, kaonp_rich_ch);
                kp_RLvsChi2_brich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
                kp_MomVsRL_brich.Fill(kaonp_rich_RL, kaonp_Ph);
                kp_ChvsChi2_brich.Fill(kaonp_rich_chi2, kaonp_rich_ch);
                kp_nTotvsChi2_brich.Fill(kaonp_rich_chi2, kaonp_rich_nTot);
                kp_MomvsChi2_brich.Fill(kaonp_rich_chi2, kaonp_Ph);
                if (kaonp_rich_chi2 < 50) kp_RLvsChi2_low_brich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
                if (kaonp_rich_chi2 >= 50) kp_RLvsChi2_high_brich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
                kp_RLvsRQ_brich.Fill(kaonp_rich_RL, kaonp_rich_RQ);
                kp_RLvsCh_brich.Fill(kaonp_rich_RL, kaonp_rich_ch);
                kp_RQvsCh_brich.Fill(kaonp_rich_RQ, kaonp_rich_ch);
                if (kaonp_rich_RL <= 2.5) kp_Blob1VsCh_brich.Fill(kaonp_rich_RL, kaonp_rich_ch);
                if (kaonp_rich_RL > 2.5) kp_Blob2VsCh_brich.Fill(kaonp_rich_RL, kaonp_rich_ch);
            }
            
        }
        if(std::abs(kaonp_chi2pid) < 5 && KinematicPID_Vertex(kaonp_vz, torus, 321) && kaonp_rich_pid == 321 && kaonp_track_pid == 321){
            if (kaonp_rich_RL > 6 || kaonp_rich_nTot <= 2 || kaonp_rich_RQ < 0.1) continue;
            //treeKaonP_goodRich.Fill();
            kp_evnt_chi2_grich.Fill(kaonp_chi2pid);
            kp_rich_chi2_grich.Fill(kaonp_rich_chi2);
            kp_rich_Mchi2_grich.Fill(kaonp_rich_Mchi2);
            kp_rich_RL_grich.Fill(kaonp_rich_RL);
            kp_rich_RQ_grich.Fill(kaonp_rich_RQ);
            kp_rich_nTot_grich.Fill(kaonp_rich_nTot);
            kp_rich_ch_grich.Fill(kaonp_rich_ch);
            kp_evnt_pid_grich.Fill(kaonp_track_pid);
            kp_MomVsBeta_grich.Fill(kaonp_Ph, kaonp_beta);
            kp_MomVsCh_grich.Fill(kaonp_Ph, kaonp_rich_ch);
            kp_RLvsChi2_grich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
            kp_MomVsRL_grich.Fill(kaonp_rich_RL, kaonp_Ph);
            kp_ChvsChi2_grich.Fill(kaonp_rich_chi2, kaonp_rich_ch);
            kp_nTotvsChi2_grich.Fill(kaonp_rich_chi2, kaonp_rich_nTot);
            kp_MomvsChi2_grich.Fill(kaonp_rich_chi2, kaonp_Ph);
            if (kaonp_rich_chi2 < 50) kp_RLvsChi2_low_grich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
            if (kaonp_rich_chi2 >= 50) kp_RLvsChi2_high_grich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
            kp_RLvsRQ_grich.Fill(kaonp_rich_RL, kaonp_rich_RQ);
            kp_RLvsCh_grich.Fill(kaonp_rich_RL, kaonp_rich_ch);
            kp_RQvsCh_grich.Fill(kaonp_rich_RQ, kaonp_rich_ch);
        }
        if(std::abs(kaonp_chi2pid) < 5 && KinematicPID_Vertex(kaonp_vz, torus, 321) && kaonp_rich_pid == 321 && kaonp_track_pid != 321){
            double m = 0.218;
            double q = -0.695;
            double y_limit = m*kaonp_rich_RL + q;
            double x_limit = (kaonp_rich_RQ - q)/m; // kaonp_rich_RQ < y_limit or kaonp_rich_RL > x_limit
            if (passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph)){
                //treeKaonP_goodRich.Fill();
                kp_evnt_chi2_mrich.Fill(kaonp_chi2pid);
                kp_rich_chi2_mrich.Fill(kaonp_rich_chi2);
                kp_rich_Mchi2_mrich.Fill(kaonp_rich_Mchi2);
                kp_rich_RL_mrich.Fill(kaonp_rich_RL);
                kp_rich_RQ_mrich.Fill(kaonp_rich_RQ);
                kp_rich_nTot_mrich.Fill(kaonp_rich_nTot);
                kp_rich_ch_mrich.Fill(kaonp_rich_ch);
                kp_evnt_pid_mrich.Fill(kaonp_track_pid);
                kp_MomVsBeta_mrich.Fill(kaonp_Ph, kaonp_beta);
                kp_MomVsCh_mrich.Fill(kaonp_Ph, kaonp_rich_ch);
                kp_MomVsnTot_mrich.Fill(kaonp_Ph, kaonp_rich_nTot);
                kp_RLvsChi2_mrich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
                kp_MomVsRL_mrich.Fill(kaonp_rich_RL, kaonp_Ph);
                kp_ChvsChi2_mrich.Fill(kaonp_rich_chi2, kaonp_rich_ch);
                kp_nTotvsChi2_mrich.Fill(kaonp_rich_chi2, kaonp_rich_nTot);
                kp_MomvsChi2_mrich.Fill(kaonp_rich_chi2, kaonp_Ph);
                if (kaonp_rich_chi2 < 50) kp_RLvsChi2_low_mrich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
                if (kaonp_rich_chi2 >= 50) kp_RLvsChi2_high_mrich.Fill(kaonp_rich_chi2, kaonp_rich_RL);
                kp_RLvsRQ_mrich.Fill(kaonp_rich_RL, kaonp_rich_RQ);
                kp_RLvsCh_mrich.Fill(kaonp_rich_RL, kaonp_rich_ch);
                kp_RQvsCh_mrich.Fill(kaonp_rich_RQ, kaonp_rich_ch);
            }
        }
        //
        if(std::abs(kaonp_chi2pid) < 5 && KinematicPID_Vertex(kaonp_vz, torus, 321)){
            if(kaonp_track_pid == 321) kp_Mx_tof.Fill(kaonp_Mx);
            if(kaonp_rich_pid == 321) kp_Mx_rich.Fill(kaonp_Mx);
            // passRichSelection(kaonp_rich_RL,kaonp_rich_RQ,kaonp_rich_nTot,kaonp_rich_pid,kaonp_track_pid,kaonp_rich_chi2,kaonp_Ph)
            if(kaonp_rich_pid == 321){
                double ph2 = kaonp_Ph*kaonp_Ph;
                double beta2 = kaonp_beta*kaonp_beta;
                kaonp_m = sqrt((ph2)*(1-beta2)/beta2);
                //
                kp_evnt_chi2.Fill(kaonp_chi2pid);
                kp_m.Fill(kaonp_bestMass);
                kp_MomVsMass_Rich.Fill(kaonp_Ph, kaonp_bestMass);
                double beta_th = kaonp_Ph/(sqrt(ph2+0.01948816));
                double delta_beta = kaonp_beta - beta_th;
                kp_deltaB.Fill(delta_beta);
                kp_rich_Mchi2.Fill(kaonp_rich_Mchi2);
                kp_Mx.Fill(kaonp_Mx);
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
                if(kaonp_chi2pid < 3) kp_MomVsBeta_ToF1.Fill(kaonp_Ph, kaonp_beta);
                if(kaonp_chi2pid >= 3) kp_MomVsBeta_ToF2.Fill(kaonp_Ph, kaonp_beta);
                kp_MomVsCh.Fill(kaonp_Ph, kaonp_rich_ch);
                if(kaonp_rich_chi2 < 4) kp_MomVsCh_RICH1.Fill(kaonp_Ph, kaonp_rich_ch);
                if(kaonp_rich_chi2 >= 4) kp_MomVsCh_RICH2.Fill(kaonp_Ph, kaonp_rich_ch);
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
                kp_rich_chi2.Fill(kaonp_rich_chi2);
                kp_rich_traj1_x.Fill(kaonp_rich_tr1_x);
                kp_rich_traj1_y.Fill(kaonp_rich_tr1_y);
                kp_rich_pmt_xy.Fill(kaonp_rich_tr1_x, kaonp_rich_tr1_y);
                kp_rich_traj2_x.Fill(kaonp_rich_tr2_x);
                kp_rich_traj2_y.Fill(kaonp_rich_tr2_y);
                kp_rich_traj2_z.Fill(kaonp_rich_tr2_z);
                kp_rich_traj2_path.Fill(kaonp_rich_tr2_path);
                kp_rich_traj3_x.Fill(kaonp_rich_tr3_x);
                kp_rich_traj3_y.Fill(kaonp_rich_tr3_y);
                kp_rich_traj3_z.Fill(kaonp_rich_tr3_z);
                kp_rich_traj3_path.Fill(kaonp_rich_tr3_path);
                kp_rich_traj4_x.Fill(kaonp_rich_tr4_x);
                kp_rich_traj4_y.Fill(kaonp_rich_tr4_y);
                kp_rich_traj4_z.Fill(kaonp_rich_tr4_z);
                kp_rich_traj4_path.Fill(kaonp_rich_tr4_path);
                kp_rich_aerogel_x.Fill(kaonp_rich_tr3_x);
                kp_rich_aerogel_x.Fill(kaonp_rich_tr4_x);
                kp_rich_aerogel_y.Fill(kaonp_rich_tr3_y);
                kp_rich_aerogel_y.Fill(kaonp_rich_tr4_y);
                kp_rich_aerogel_xy.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_y);
                kp_rich_aerogel_xy.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_y);
                kp_rich_aerogel_xy.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_y);
                //kp_rich_aerogel_xy.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_y);
                kp_rich_aerogel_l1_xy.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_y);
                kp_rich_aerogel_l2_xy.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_y);
                kp_rich_aerogel_xp.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_path);
                kp_rich_aerogel_xp.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_path);
                //kp_rich_aerogel_xp.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_path);
                kp_rich_aerogel_yp.Fill(kaonp_rich_tr2_y, kaonp_rich_tr2_path);
                kp_rich_aerogel_yp.Fill(kaonp_rich_tr3_y, kaonp_rich_tr3_path);
                kp_rich_aerogel_xz.Fill(kaonp_rich_tr2_x, kaonp_rich_tr2_z);
                kp_rich_aerogel_xz.Fill(kaonp_rich_tr3_x, kaonp_rich_tr3_z);
                //kp_rich_aerogel_xz.Fill(kaonp_rich_tr4_x, kaonp_rich_tr4_z);
            }
        }
    }

    
    //dirKaonp->cd();
    //outFile.cd("KaonPlus");
    // I prefer to have the tree before the th2d plot
    treeKaonP_badRich.Write();
    treeKaonP_goodRich.Write();
    treeKaonP.Write();

    kp_evnt_chi2.Write();
    kp_rich_chi2.Write();
    kp_rich_Mchi2.Write();
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

    std::vector<TH2D*> hists_kp = {
        &kp_MomVsPhT, &kp_MomVsXb, &kp_MomVsXf, &kp_MomVsZ, &kp_MomVsY, &kp_MomVsEta,
        &kp_MomVsTheta, &kp_MomVsPhi_h, &kp_MomVsMx, &kp_MomVsBeta, &kp_MomVsBeta_ToF1, &kp_MomVsBeta_ToF2, 
        &kp_MomVsCh, &kp_MomVsCh_RICH1, &kp_MomVsCh_RICH2,
        &kp_Q2VsXb, &kp_Q2VsXf, &kp_Q2VsMom, &kp_Q2VsPhT, &kp_Q2VsZ, &kp_Q2VsY, &kp_Q2VsEta, &kp_Q2VsPhi_h,
        &kp_PhTvsZ, &kp_PhTvsXb, &kp_PhTvsEta, &kp_PhTvsPhi_h,
        &kp_zVsXb, &kp_zVsXf, &kp_zVsEta, &kp_zVsPhi_h,
        &kp_xBvsY, &kp_ThetaVsPhi_h, &kp_ThetaVsPhi_Lab, &empty_canvas,
        &kp_MomVsMass_Rich,
        &kp_rich_pmt_xy, &kp_rich_aerogel_xy, &kp_rich_aerogel_l1_xy, &kp_rich_aerogel_l2_xy, &kp_rich_aerogel_xp, &kp_rich_aerogel_yp, &kp_rich_aerogel_xz
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

    // good rich
    empty_canvas_grich.Write();
    kp_evnt_pid_grich.Write();
    kp_evnt_chi2_grich.Write();
    kp_rich_chi2_grich.Write();
    kp_rich_Mchi2_grich.Write();
    kp_rich_RL_grich.Write();
    kp_rich_RQ_grich.Write();
    kp_rich_nTot_grich.Write();
    kp_rich_ch_grich.Write();
    kp_MomVsBeta_grich.Write();
    kp_MomVsCh_grich.Write();
    kp_MomVsRL_grich.Write();
    kp_MomvsChi2_grich.Write();
    kp_RLvsChi2_grich.Write();
    kp_RLvsChi2_low_grich.Write();
    kp_RLvsChi2_high_grich.Write();
    kp_ChvsChi2_grich.Write();
    kp_nTotvsChi2_grich.Write();
    kp_RLvsRQ_grich.Write();
    kp_RLvsCh_grich.Write();
    kp_RQvsCh_grich.Write();
    // maybe rich
    empty_canvas_mrich.Write();
    kp_evnt_pid_mrich.Write();
    kp_evnt_chi2_mrich.Write();
    kp_rich_chi2_mrich.Write();
    kp_rich_Mchi2_mrich.Write();
    kp_rich_RL_mrich.Write();
    kp_rich_RQ_mrich.Write();
    kp_rich_nTot_mrich.Write();
    kp_rich_ch_mrich.Write();
    kp_MomVsBeta_mrich.Write();
    kp_MomVsCh_mrich.Write();
    kp_MomVsnTot_mrich.Write();
    kp_MomVsRL_mrich.Write();
    kp_MomvsChi2_mrich.Write();
    kp_RLvsChi2_mrich.Write();
    kp_RLvsChi2_low_mrich.Write();
    kp_RLvsChi2_high_mrich.Write();
    kp_ChvsChi2_mrich.Write();
    kp_nTotvsChi2_mrich.Write();
    kp_RLvsRQ_mrich.Write();
    kp_RLvsCh_mrich.Write();
    kp_RQvsCh_mrich.Write();
    // bad rich
    empty_canvas_brich.Write();
    kp_rich_pid_brich.Write();
    kp_evnt_chi2_brich.Write();
    kp_rich_chi2_brich.Write();
    kp_rich_Mchi2_brich.Write();
    kp_rich_RL_brich.Write();
    kp_rich_RQ_brich.Write();
    kp_rich_nTot_brich.Write();
    kp_rich_ch_brich.Write();
    kp_MomVsBeta_brich.Write();
    kp_MomVsCh_brich.Write();
    kp_MomVsRL_brich.Write();
    kp_MomvsChi2_brich.Write();
    kp_RLvsChi2_brich.Write();
    kp_RLvsChi2_low_brich.Write();
    kp_RLvsChi2_high_brich.Write();
    kp_ChvsChi2_brich.Write();
    kp_nTotvsChi2_brich.Write();
    kp_RLvsRQ_brich.Write();
    kp_RLvsCh_brich.Write();
    kp_RQvsCh_brich.Write();
    kp_Blob1VsCh_brich.Write();
    kp_Blob2VsCh_brich.Write();
    //
    empty_canvas2.Write();

    std::set<std::string> invertOrder = {"c_rich_chi2", "c_rich_RL"};
    std::set<std::string> invertOrder2 = {"c_rich_Mchi2", "c_rich_ch"};
    // Funzione helper per disegnare due istogrammi su stessa canvas
    auto NormalizeHist = [](TH1D& h) {
        double integral = h.Integral();
        if (integral > 0) h.Scale(1.0 / integral);
    };
    auto DrawComparison = [&](TH1D& h_good, TH1D& h_bad, TH1D& h_maybe, const char* cname) {
        TCanvas* c = new TCanvas(cname, cname, 800, 600);

        // Creo copie così non modifico gli originali
        TH1D h_good_norm = h_good;
        TH1D h_bad_norm  = h_bad;
        TH1D h_maybe_norm  = h_maybe;
        h_good_norm.SetStats(0); 
        h_bad_norm.SetStats(0);
        h_maybe_norm.SetStats(0);

        NormalizeHist(h_good_norm), NormalizeHist(h_bad_norm), NormalizeHist(h_maybe_norm);
        h_good_norm.SetLineColor(kRed+1), h_bad_norm.SetLineColor(kAzure+1), h_maybe_norm.SetLineColor(kGreen+1);

        if (invertOrder.count(cname)) {
            h_good_norm.Draw("HIST");
            h_bad_norm.Draw("HIST SAME");
            h_maybe_norm.Draw("HIST SAME");
        } else if (invertOrder2.count(cname)){
            h_maybe_norm.Draw("HIST");
            h_good_norm.Draw("HIST SAME");
            h_bad_norm.Draw("HIST SAME");
            h_maybe_norm.Draw("HIST SAME");
        } else {
            h_bad_norm.Draw("HIST");
            h_good_norm.Draw("HIST SAME");
            h_maybe_norm.Draw("HIST SAME");
        }

        // Aggiungo legenda
        TLegend* leg = new TLegend(0.72,0.75,0.88,0.88);
        leg->AddEntry(&h_good_norm, "RICH K+ Evnt K+", "l");
        leg->AddEntry(&h_maybe_norm, "RICH K+ Evnt not K+", "l");
        leg->AddEntry(&h_bad_norm, "RICH not K+ Evnt K+", "l");
        leg->Draw();

        gPad->Update();
        c->Write();
    };

    DrawComparison(kp_evnt_chi2_grich, kp_evnt_chi2_brich, kp_evnt_chi2_mrich, "c_evnt_chi2");
    DrawComparison(kp_rich_chi2_grich, kp_rich_chi2_brich, kp_rich_chi2_mrich, "c_rich_chi2");
    DrawComparison(kp_rich_Mchi2_grich, kp_rich_Mchi2_brich, kp_rich_Mchi2_mrich, "c_rich_Mchi2");
    DrawComparison(kp_rich_RL_grich, kp_rich_RL_brich, kp_rich_RL_mrich, "c_rich_RL");
    DrawComparison(kp_rich_RQ_grich, kp_rich_RQ_brich, kp_rich_RQ_mrich, "c_rich_RQ");
    DrawComparison(kp_rich_nTot_grich, kp_rich_nTot_brich, kp_rich_nTot_mrich, "c_rich_nTot");
    DrawComparison(kp_rich_ch_grich, kp_rich_ch_brich, kp_rich_ch_mrich, "c_rich_ch");
    
    //outFile.Write();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
