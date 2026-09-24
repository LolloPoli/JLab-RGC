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
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
//#include "HipoChain.h"
#include "ccdb_reader.h"
#include "rcdb_reader.h"
#include "rich.h"
#include "particle.h"
#include "particle_detector.h"
#include "QADB.h"

using namespace clas12;
using namespace std;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

float BeamPolarization(Int_t run, Bool_t v) {
  if      (run>=16137 && run<=16148) return v ? 0.630 : 0.0; 
  else if (run>=16156 && run<=16178) return v ? -0.585 : 0.0; 
  else if (run>=16211 && run<=16228) return v ? 0.690 : 0.0; 
  else if (run>=16231 && run<=16260) return v ? -0.580 : 0.0; 
  else if (run>=16318 && run<=16333) return v ? 0.610 : 0.0; 
  else if (run>=16335 && run<=16357) return v ? -0.520 : 0.0; 
  else if (run>=16658 && run<=16675) return v ? 0.520 : 0.0; 
  else if (run>=16676 && run<=16720) return v ? 0.574 : 0.0; 
  else if (run>=16723 && run<=16766) return v ? -0.512 : 0.0; 
  else if (run>=16767 && run<=16772) return v ? 0.582 : 0.0; 
  else {
    fprintf(stderr,"WARNING: polarization unknown for run %d\n",run);
    return 0.0;
  }
}

struct RGC_FALL {
    float run;
    // float total_charge;
    // float positive_beam_charge;
    // float negative_beam_charge;
    float polariz;
    // float polariz_unc;
};
std::map<int, double> polariz_map;
void LoadPolarization(const std::string& filename) {
    std::ifstream infile(filename);
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string field;
        int run;
        double penultimate;
        std::string last;
        // Colonna 1: run
        std::getline(ss, field, ',');
        run = std::stoi(field);
        // Colonne 2-4: salta
        for (int i = 0; i < 3; ++i) {
            std::getline(ss, field, ',');
        }
        // Colonna 5: penultima (quella che ti serve)
        std::getline(ss, field, ',');
        penultimate = std::stod(field);
        // Colonna 6: ultima (opzionale, solo per completare parsing)
        std::getline(ss, last, ',');
        polariz_map[run] = penultimate;
    }
}

double BeamPolarization_file(Int_t run){
  if (polariz_map.empty()) {
    //std::cout << "Caricamento tabella polarizzazione..." << std::endl;
    LoadPolarization("spring23_NH3_goodRuns.txt"); // fall and summer 
  }
  return polariz_map[run];
}

// Electron PID
// B1
bool ElectronPID_ForwardDetector(int status){
    return (status > -4000 && status <= 2000);
}
// B2
bool ElectronPID_nPhe(int nphe){
    return nphe > 2;
}
// B3
bool ElectronPID_PCAL(double pcal){
    return pcal > 0.07;
}
// B4
bool ElectronPID_CalSFcut(int sector, int runnum, double p, double cal_energy){
    // calorimeter sampling fraction cut
    // getSector(), Nrun, electron mom, electron dep E in every cal, ECAL, PCAL, ECIN? or idk the itter o pre shower
  double scale = 3.5; // how many std away from mean to cut on
  // Common calculation for mean and std
  double mean = 0.0;
  double std = 0.0;
  std::vector<std::vector<double>> e_cal_sampl_mu;
  std::vector<std::vector<double>> e_cal_sampl_sigma;

  if ((runnum == 11) || ((runnum >= 5032 && runnum <= 5666) || (runnum >= 6616 && runnum <= 6783))) { // RGA
    e_cal_sampl_mu = {
      {0.2531, 0.2550, 0.2514, 0.2494, 0.2528, 0.2521},
      {-0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703},
      {4.939, 5.350, 5.102, 6.440, 6.149, 4.957}
    };

    e_cal_sampl_sigma = {
      {0.002726, 0.004157, 0.00522, 0.005398, 0.008453, 0.006553},
      {1.062, 0.859, 0.5564, 0.6576, 0.3242, 0.4423},
      {-4.089, -3.318, -2.078, -2.565, -0.8223, -1.274}
    };
  } else if (runnum >= 6120 && runnum <= 6604) { // RGB
    e_cal_sampl_mu = {
      {0.2520, 0.2520, 0.2479, 0.2444, 0.2463, 0.2478},
      {-0.8615, -0.8524, -0.6848, -0.5521, -0.5775, -0.7327},
      {5.596, 6.522, 5.752, 5.278, 6.430, 5.795}
    };

    e_cal_sampl_sigma = {
      {-0.02963, -0.1058, -0.05087, -0.04524, -0.02951, -0.01769},
      {20.4, 129.3, 0.6191, 0.6817, 20.84, 8.44},
      {-41.44, -101.6, -2.673, -2.606, -42.67, -21.73}
    };
  } else if ((runnum >= 11323 && runnum <= 21571) || (runnum >= 11093 && runnum <= 11300)) {
    // RGB winter 2020 // (also using this for RGB fall 2019, but it should be updated! TODO)
    e_cal_sampl_mu = {
      {0.2433, 0.2421, 0.2415, 0.2486, 0.2419, 0.2447},
      {-0.8052, -1.0495, -1.1747, -0.5170, -0.6840, -0.9022},
      {5.2750, 4.4886, 4.4935, 5.9044, 5.6716, 4.9288}
    };
    e_cal_sampl_sigma = {
      {0.0120, 0.0164, 0.0120, 0.0108, 0.0147, 0.0077},
      {0.1794, 0.1519, 0.1379, 0.1838, 0.0494, 0.3509},
      {-0.0695, 0.1553, 0.3300, 0.4330, 1.1032, -0.7996}
    };
  }

  // Calculation of mean and std
  mean = e_cal_sampl_mu[0][sector] + (e_cal_sampl_mu[1][sector] / 1000) * (p - e_cal_sampl_mu[2][sector]) * (p - e_cal_sampl_mu[2][sector]);
  std = e_cal_sampl_sigma[0][sector] + e_cal_sampl_sigma[1][sector] / (10 * (p - e_cal_sampl_sigma[2][sector]));

  // Return result
  return ((cal_energy / p) > (mean - scale * std)) && ((cal_energy / p) < (mean + scale * std));
}
// B5
bool ElectronPID_DiagonalCut(clas12::region_part_ptr pa, TLorentzVector *p){
  bool response = true;
  double einner = pa->cal(ECIN)->getEnergy();
  //double einner = pa->cal(ECIN)->getEnergy()+pa->cal(ECOUT)->getEnergy();
  if(p->P() > 4.5){
    double xxx = einner/p->P();
    double yyy = pa->cal(PCAL)->getEnergy()/p->P();
    if(xxx + yyy < (0.2) )response = false;
    //cout <<Form("Diagonal cut check: %lf + %lf = %lf > 0.2) ? %d",xxx,yyy,xxx+yyy,response);
    //cin.get();
  }
  return response;
}
//B6
bool ElectronPID_VertexCut(double vz, double vx, double vy){
    return (vz > -7.5 && vz < 2.5 && abs(vx) <= 2 && abs(vy) <= 2); // summer -9, 1, fall and spring use -7.5, 2.5
}

// EC lv,lu,lw cut
bool ElectronEC_cut(double lu, double lv,  double lw){
  if(lu > 19 && lu < 400 && lv > 9 && lw > 9) return true;
  else return false;
  
}
// strip removal
bool ElectronECin_strip(double lv, int sector){
  if(sector == 1){
    if(lv < 72 || lv > 94) return true;
    else return false;
  }
  else return true;
}
bool ElectronECout_strip(double lu, double lw, int sector){
  if(sector == 2){
    if(lw < 68 || lw > 84) return true;
    else return false;
  }
  else if (sector == 5){
    if(lu < 200 || lu > 220) return true;
    else return false;
  }
  else return true;
}
bool ElectronPCAL_cut(double lv, double lw, int sector){
  if(sector == 1){
    return (lv > 22.5 && lw > 22.5);
  }
  else{
    return (lv > 13.5 && lw > 13.5);
  }
}

bool ElectronPID_DC(double edge1, double edge2, double edge3, int torus){
  if(torus == -1){  // INBENDING ELECTRON 
    if(edge1 > 4 && edge2 > 5 && edge3 > 8) return true; // from Derek Holmberg
    else return false;
  }
  else if(torus == 1){  // OUTBENDING ELECTRON
    if(edge1 > 3 && edge2 > 3 && edge3 > 10) return true;
    else return false;
  }
  else return false;
}


// HADRON PID 
// C2
bool HadronPID_Reson(double mom, double m2x){
  if(mom > 1.1){
    if(m2x > 2.45) return true;
    else return false;
  }  
  else return false;
 
}

bool HadronPID_Q2(double q2){
  return q2 >= 0.95;
}

bool KinematicPID_xF(double xF){
  return xF > 0;
}

bool KinematicPID_W(double w){
  return w > 1.95;
}

bool HadronPID_DC(double edge1, double edge2, double edge3, int torus){
  if(torus == -1){  // INBENDING HADRON 
    if(edge1 > 3 && edge2 > 3 && edge3 > 7) return true;
    else return false;
  }
  else if(torus == 1){  // OUTBENDING HADRON  
    if(edge1 > 3 && edge2 > 3 && edge3 > 9) return true;
    else return false;
  }
  else return false;
}


bool KinematicPID_Vertex(double vz, int torus, int pid){
  bool signal = false;
  if(torus == -1){  // INBENDING
    if(pid > 20){
      if(-10 < vz && vz < 1) signal = true; // maybe better -10, 2
    } 
    if(pid < 20){
      if(-9 < vz && vz < 1) signal = true; // summer -9, 1, fall and spring use -7.5, 2.5
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

bool HadronPID_Chi2Pid(double chi2){
  return (std::abs(chi2) <= 3);
}

bool HadronPID_z(double z){
  return z > -1;
}


void rgc(const char* fileList){
    // the file list has a first line with the output file, and the follow lines are the input file
    clas12root::HipoChain chain;
    std::ifstream file(fileList);
    if (!file) {
        std::cerr << "Error: impossible to open the file" << fileList << std::endl;
        return;
    }
    
    std::string outputFile;
    std::vector<std::string> inputFiles;
    std::string line;
    bool firstLine = true;
    while (std::getline(file, line)) {
      if (firstLine) {
          outputFile = line;  // Prima riga = file di output
          firstLine = false;
      } else {
          inputFiles.push_back(line);  // Altre righe = file di input
      }
    }
    file.close();
    if (inputFiles.empty()) {
      std::cerr << "Error: no input inside " << fileList << std::endl;
      return;
    }
    for (int i = 0; i < inputFiles.size(); i++){
      chain.Add(inputFiles[i].c_str());
    }
    // Create the TTree
    TFile outFile(outputFile.c_str(), "RECREATE");
    TTree Kaon_tree("Kaon+", "kaon");
    double event = 0;
    double event_mc = 0;
    // ELECTRON
    double electron_px, electron_py, electron_pz, electron_mom, electron_eng;
    double electron_px_mc, electron_py_mc, electron_pz_mc, electron_mom_mc, electron_eng_mc;
    double electron_theta, electron_phi, electron_theta_mc, electron_phi_mc;
    int electron_Nphe, electron_status, electron_sector;
    double electron_PCAL, electron_ECAL, electron_ECIN, electron_CAL_Tot, electron_vz, electron_vx;
    double electron_edge1, electron_edge2, electron_edge3, electron_ass;
    double electron_ECin_lu, electron_ECin_lv, electron_ECin_lw;
    double electron_PCAL_lu, electron_PCAL_lv, electron_PCAL_lw;
    double electron_ECout_lu, electron_ECout_lv, electron_ECout_lw;
    double CVT_edge12, CVT_edge1, CVT_edge3, CVT_edge5, CVT_edge7;
    // elettrone
    double electron_Phi, electron_E, electron_W, electron_Q2;
    // gamma
    double gamma_px, gamma_py, gamma_pz, gamma_nu;
    // KAON
    double kaon_px, kaon_py, kaon_pz, kaon_mom, kaon_eng;
    double kaon_px_mc, kaon_py_mc, kaon_pz_mc, kaon_mom_mc, kaon_eng_mc;
    double kaon_epsilon, kaon_gamma, kaon_epsilon_mc, kaon_gamma_mc;
    double kaon_Q2, kaon_xB, kaon_xF, kaon_y, kaon_Mx, kaon_z, kaon_Pt, kaon_s, kaon_W;
    double kaon_Q2_mc, kaon_xB_mc, kaon_xF_mc, kaon_y_mc, kaon_Mx_mc, kaon_z_mc, kaon_Pt_mc, kaon_s_mc, kaon_W_mc;
    double kaon_theta, kaon_phi_lab, kaon_phi_h, kaon_eta, kaon_Pt_ratio;
    double kaon_theta_mc, kaon_phi_lab_mc, kaon_phi_h_mc, kaon_eta_mc, kaon_Pt_ratio_mc;
    double kaon_parent, kaon_parent_mc, kaon_parent_idx, kaon_vz;
    double kaon_rich_tr1_x, kaon_rich_tr1_y, kaon_rich_tr1;
    double kaon_rich_tr2_x, kaon_rich_tr2_y, kaon_rich_tr2;
    double kaon_rich_tr3_x, kaon_rich_tr3_y, kaon_rich_tr3;
    double kaon_rich_tr4_x, kaon_rich_tr4_y, kaon_rich_tr4;
    double kaon_best_mass, kaon_rich_Id, kaon_rich_PID, kaon_rich_ntot, kaon_rich_chi2, kaon_rich_RQ, kaon_rich_RL;
    double kaon_rich_ch, kaon_rich_Mchi2;
    double kaon_chi2pid, kaon_polariz;
    double clas_pid, kaon_beta;
    double beam_helicity;
    double kaon_tt;
    double photon_en;

    
    // VECTOR CMS
    TVector3 beta_cms_kaon;
    TVector3 beta_cms_kaon_mc;
    // KAON+
    // REC
    //Kaon_tree.Branch("event", &event);
    Kaon_tree.Branch("el_px", &electron_px), Kaon_tree.Branch("el_py", &electron_py), Kaon_tree.Branch("el_pz", &electron_pz);
    Kaon_tree.Branch("el_mom", &electron_mom), Kaon_tree.Branch("el_W", &electron_W);
    Kaon_tree.Branch("el_vz", &electron_vz), Kaon_tree.Branch("el_vx", &electron_vx);
    Kaon_tree.Branch("el_theta", &electron_theta), Kaon_tree.Branch("el_phi", &electron_phi);
    Kaon_tree.Branch("gamma_E", &photon_en);
    Kaon_tree.Branch("s", &kaon_s), Kaon_tree.Branch("clas_chi2", &kaon_chi2pid, "clas_chi2/D");
    Kaon_tree.Branch("clas_pid", &clas_pid), Kaon_tree.Branch("beta", &kaon_beta);
    Kaon_tree.Branch("Q2", &kaon_Q2), Kaon_tree.Branch("xB", &kaon_xB), Kaon_tree.Branch("xF", &kaon_xF), Kaon_tree.Branch("y", &kaon_y), Kaon_tree.Branch("z", &kaon_z), Kaon_tree.Branch("-t", &kaon_tt);
    Kaon_tree.Branch("kaon_px", &kaon_px, "kaon_px/D"), Kaon_tree.Branch("kaon_py", &kaon_py, "kaon_py/D"), Kaon_tree.Branch("kaon_pz", &kaon_pz, "kaon_pz/D");
    Kaon_tree.Branch("kaon_mom", &kaon_mom), Kaon_tree.Branch("kaon_Pt", &kaon_Pt), Kaon_tree.Branch("Pt_over_zQ", &kaon_Pt_ratio);
    Kaon_tree.Branch("kaon_theta", &kaon_theta), Kaon_tree.Branch("kaon_phi_lab", &kaon_phi_lab), Kaon_tree.Branch("kaon_phi_h", &kaon_phi_h);
    Kaon_tree.Branch("eta", &kaon_eta), Kaon_tree.Branch("kaon_vz", &kaon_vz) , Kaon_tree.Branch("kaon_parent", &kaon_parent);
    Kaon_tree.Branch("gamma", &kaon_gamma), Kaon_tree.Branch("epsilon", &kaon_epsilon);
    Kaon_tree.Branch("Mx", &kaon_Mx), Kaon_tree.Branch("W", &kaon_W);
    Kaon_tree.Branch("Polarization", &kaon_polariz, "Polarization/D"), Kaon_tree.Branch("Helicity", &beam_helicity, "Helicity/D");
    //CAL
    /*
    Kaon_tree.Branch("el_PCAL_e", &electron_PCAL, "el_PCAL_e/D"), Kaon_tree.Branch("el_status", &electron_status, "el_status/I"), Kaon_tree.Branch("el_sector", &electron_sector, "el_sector/I");
    Kaon_tree.Branch("el_PCAL_lu", &electron_PCAL_lu, "el_PCAL_lu/D"), Kaon_tree.Branch("el_PCAL_lv", &electron_PCAL_lv, "el_PCAL_lv/D"), Kaon_tree.Branch("el_PCAL_lw", &electron_PCAL_lw, "el_PCAL_lw/D");
    Kaon_tree.Branch("el_ECin_lu", &electron_ECin_lu, "el_ECin_lu/D"), Kaon_tree.Branch("el_ECin_lv", &electron_ECin_lv, "el_ECin_lv/D"), Kaon_tree.Branch("el_ECin_lw", &electron_ECin_lw, "el_ECin_lw/D");
    Kaon_tree.Branch("el_ECout_lu", &electron_ECout_lu, "el_ECout_lu/D"), Kaon_tree.Branch("el_ECout_lv", &electron_ECout_lv, "el_ECout_lv/D"), Kaon_tree.Branch("el_ECout_lw", &electron_ECout_lw, "el_ECout_lw/D");
    */
    //
    //Kaon_tree.Branch("CVT_edge1", &CVT_edge1, "CVT_edge1/D"); 
    //Kaon_tree.Branch("CVT_edge3", &CVT_edge3, "CVT_edge3/D"), Kaon_tree.Branch("CVT_edge5", &CVT_edge5, "CVT_edge5/D");
    //Kaon_tree.Branch("CVT_edge7", &CVT_edge7, "CVT_edge7/D"), Kaon_tree.Branch("CVT_edge12", &CVT_edge12, "CVT_edge12/D");
    //
    Kaon_tree.Branch("Rich_Id", &kaon_rich_Id, "Rich_Id/D");
    Kaon_tree.Branch("Rich_mass", &kaon_best_mass, "Rich_mass/D");
    Kaon_tree.Branch("Rich_PID", &kaon_rich_PID, "Rich_PID/D");
    Kaon_tree.Branch("Rich_RL", &kaon_rich_RL, "Rich_RL/D");
    Kaon_tree.Branch("Rich_RQ", &kaon_rich_RQ, "Rich_RQ/D");
    Kaon_tree.Branch("Rich_nTot", &kaon_rich_ntot, "Rich_nTot/D");
    Kaon_tree.Branch("Rich_ch", &kaon_rich_ch, "Rich_ch/D");
    Kaon_tree.Branch("Rich_chi2", &kaon_rich_chi2, "Rich_chi2/D");
    Kaon_tree.Branch("Rich_Mchi2", &kaon_rich_Mchi2, "Rich_Mchi2/D");
    Kaon_tree.Branch("Rich_PMT_edge", &kaon_rich_tr1);
    Kaon_tree.Branch("Rich_PMT_x", &kaon_rich_tr1_x, "Rich_PMT_x/D"), Kaon_tree.Branch("Rich_PMT_y", &kaon_rich_tr1_y, "Rich_PMT_y/D");
    Kaon_tree.Branch("Rich_aerogel1_edge", &kaon_rich_tr2);
    Kaon_tree.Branch("Rich_aerogel1_x", &kaon_rich_tr2_x, "Rich_aerogel1_x/D"), Kaon_tree.Branch("Rich_aerogel1_y", &kaon_rich_tr2_y, "Rich_aerogel1_y/D");
    Kaon_tree.Branch("Rich_aerogel2_edge", &kaon_rich_tr3);
    Kaon_tree.Branch("Rich_aerogel2_x", &kaon_rich_tr3_x, "Rich_aerogel2_x/D"), Kaon_tree.Branch("Rich_aerogel2_y", &kaon_rich_tr3_y, "Rich_aerogel2_y/D");
    Kaon_tree.Branch("Rich_aerogel3_edge", &kaon_rich_tr4);
    Kaon_tree.Branch("Rich_aerogel3_x", &kaon_rich_tr4_x, "Rich_aerogel3_x/D"), Kaon_tree.Branch("Rich_aerogel3_y", &kaon_rich_tr4_y, "Rich_aerogel3_y/D");

    auto db2 = TDatabasePDG::Instance();
    // Usefull Lorentz vector
    TLorentzVector beam(0, 0, 10.6, 10.6); // Fascio con energia 10.6 GeV
    TLorentzVector target(0, 0, 0, db2->GetParticle(2212)->Mass());  // Target
    TLorentzVector el_4vec(0, 0, 0, db2->GetParticle(11)->Mass());   // Scattered electron
    TLorentzVector el_4vec_mc(0, 0, 0, db2->GetParticle(11)->Mass());
    TLorentzVector pr(0, 0, 0, db2->GetParticle(2212)->Mass());      // Proton
    TLorentzVector gm(0, 0, 0, db2->GetParticle(22)->Mass());        // Photon
    TLorentzVector pip_4v(0,0,0,db2->GetParticle(211)->Mass());      // Pion+
    TLorentzVector pip_4v_mc(0,0,0,db2->GetParticle(211)->Mass());
    TLorentzVector pim(0,0,0,db2->GetParticle(-211)->Mass());        // Pion-
    TLorentzVector kp_4v(0,0,0,db2->GetParticle(-321)->Mass());       // Kaon+
    TLorentzVector kp_4v_mc(0,0,0,db2->GetParticle(321)->Mass());
    TLorentzVector km_4v(0,0,0,db2->GetParticle(-321)->Mass());         // Kaon-
    TLorentzVector positron(0, 0, 0, db2->GetParticle(-11)->Mass());   // Positron
    TLorentzVector Lab = beam + target;
    double ProtonMass = db2->GetParticle(2212)->Mass();

    int nbin = 200;
    //TH2D kp_PxVsPy ("_PxVsPy", "Correlation P_{x} vs P_{y}  |  K+  | with EventBuilder + RICH ; P_{x} [GeV]; P_{y} [GeV]", nbin, -2, 2, nbin, -2, 2);
    //TH2D el_PxVsPy ("el_PxVsPy", "Correlation P_{x} vs P_{y}  |  e-  | with EventBuilder + RICH ; P_{x} [GeV]; P_{y} [GeV]", nbin, -2, 2, nbin, -2, 2);
    //TH2D el_ECin_lvlw ("el_ECin_lvlw", "ECin lv vs lw  ; lw; lv", nbin, 0, 450, nbin, 0, 450);
    //TH2D el_ECout_lulw ("el_ECout_lulw", "ECout lu vs lw  ; lu; lw", nbin, 0, 450, nbin, 0, 450);
    //TH1D elpos_selection ("elpos_selection", "Electron and Positron selection; M_{e^{-}e^{+}} [GeV]", 100, 0, 3);
    //TH1D elpos_vz ("elpos_vz", "Electron and Positron vertex position; v_{z}^{e^{-}} - v_{z}^{e^{+}} [cm]", 100, -20, 20);
    //TH1D elpos_theta_ee ("elpos_theta_ee", "Electron-Positron opening angle; #theta_{e^{-}e^{+}} [deg]", 100, 0, 60);
    //TH2D elpos_MvsTheta ("elpos_MvsTheta", "Electron-Positron invariant mass vs opening angle; #theta_{e^{-}e^{+}} [deg]; M_{e^{-}e^{+}} [GeV]", 100, 0, 10, 100, 0, 0.3);
    //TH2D elpos_Mvsvz ("elpos_Mvsvz", "Electron-Positron invariant mass vs vertex position; v_{z}^{e^{-}} - v_{z}^{e^{+}} [cm]; M_{e^{-}e^{+}} [GeV]", 100, -20, 20, 100, 0, 3);

    auto config_c12=chain.GetC12Reader();
    chain.SetReaderTags({0});
    config_c12->addAtLeastPid(11,1);
    //config_c12->addAtLeastPid(321,1);
    auto& c12=chain.C12ref();
    clas12databases::SetRCDBRemoteConnection();
    clas12databases db;
    c12->connectDataBases(&db);

    if (config_c12->rcdb()) {
        cout << "rcdb found" << endl;
        auto& rcdbData = config_c12->rcdb()->current();
    }
    // Calibration Constants DataBase, e electron_sf contiene informazioni delle frazioni energetiche degli elettroni nei calorimentri
    // copia e incollata così, è una tabella di calibrazione

    if (config_c12->ccdb()) {
        cout << "ccdb found" << endl;
        auto& ccdbElSF = config_c12->ccdb()->requestTableDoubles("/calibration/eb/electron_sf");
    }

    //TH2D kp_rich_pmt_xy ("rich_pmt_xy", "PMT RICH 4th sector, xy plane | t+1 & K+; x [cm]; y [cm]", 100, -170, 170, 100, -70, 70);
    //TH2D kp_rich_aerogel_xy ("rich_aerogel_xy", "Aerogel RICH 4th sector, xy plane | t+1 & K+; x [cm]; y [cm]", 100, -170, 220, 100, -85, 85);



    QA::QADB qa("latest");
    qa.CheckForDefect("TotalOutlier");
    qa.CheckForDefect("TerminalOutlier");
    qa.CheckForDefect("MarginalOutlier");
    qa.CheckForDefect("SectorLoss");
    qa.CheckForDefect("LowLiveTime");
    qa.CheckForDefect("Misc");
    qa.CheckForDefect("ChargeHigh");
    qa.CheckForDefect("ChargeNegative");
    qa.CheckForDefect("ChargeUnknown");
    qa.CheckForDefect("PossiblyNoBeam");

    auto setupQA = [](QA::QADB& q) {
    q.CheckForDefect("TotalOutlier");
    q.CheckForDefect("TerminalOutlier");
    q.CheckForDefect("MarginalOutlier");
    q.CheckForDefect("SectorLoss");
    q.CheckForDefect("LowLiveTime");
    q.CheckForDefect("Misc");
    q.CheckForDefect("ChargeHigh");
    q.CheckForDefect("ChargeNegative");
    q.CheckForDefect("ChargeUnknown");
    q.CheckForDefect("PossiblyNoBeam");

    // RGC run con Misc lecito — empty target / He per dilution factor
    std::vector<int> allow_misc = {
        16194, 16089, 16185, 16308, 16184, 16307, 16309, // RGC Su22 He/ET
        16872, 16975,                                      // RGC Fa22 He/ET
        17763, 17764, 17765, 17766, 17767, 17768,          // RGC Sp23 He/ET
        17179, 17180, 17181, 17182, 17183, 17188, 17189,   // RICH off/partially down
        17252
        };
        for (auto run : allow_misc) q.AllowMiscBit(run);
    };

    setupQA(qa);

    int count = 0;
    int countK = 0;
    int eventCount = 0;
    while (chain.Next()){
        eventCount++;
        auto electrons = c12->getByID(11);  
        auto kaonp = c12->getByID(321);
        auto virtual_gamma = c12->getByID(22);
        vector<clas12::region_part_ptr> electrons_vec;
        vector<clas12::region_part_ptr> kaon_vec;
        vector<clas12::region_part_ptr> kaon_min_vec;
        vector<clas12::region_part_ptr> positrons_vec;
        vector<clas12::region_part_ptr> pion_vec;
        for(auto &p:c12->getDetParticles()){
            int pidd = p->getPid();
            if(p->getPid() == 11) electrons_vec.push_back(p);
            if(p->getPid() == 321) kaon_vec.push_back(p);
            if(p->getPid() == -321) kaon_min_vec.push_back(p);
            else if(p->getPid() == -11) positrons_vec.push_back(p);
            else if(p->getPid() == 211) pion_vec.push_back(p);
            else if(p->rich() && p->rich()->getBest_PID()==-321) kaon_min_vec.push_back(p);
            else if(p->rich() && p->rich()->getBest_PID()==321) kaon_vec.push_back(p);
            else if(p->rich() && p->rich()->getBest_PID()==11) electrons_vec.push_back(p);
            else if(p->rich() && p->rich()->getBest_PID()==-11) positrons_vec.push_back(p);
            else if(p->rich() && p->rich()->getBest_PID()==211) pion_vec.push_back(p);
        }
        auto N_run = c12->runconfig()->getRun(); 
        auto N_event = c12->runconfig()->getEvent();
        auto b_helicity = c12->event()->getHelicity();
        beam_helicity = b_helicity * qa.CorrectHelicitySign(N_run, N_event);
        double target_spinstate = BeamPolarization_file(N_run);
        auto torus = -1;
        // QADB
        bool qa_passed = qa.Pass(N_run, N_event);
        if (N_run > 16600 && N_run < 16700) qa_passed = false; // Hall C bleedthrough
        if (N_run > 17768 && N_run <= 17811) qa_passed = false; // RGC Sp23 outbending
        if (N_run == 17331 || N_run == 16987 ||
            N_run == 17079 || N_run == 17190 ||
            N_run == 17639) qa_passed = false;                   // low live time
        if (N_run == 16850 || N_run == 16851 || N_run == 16852 ||
            N_run == 16855 || N_run == 16879) qa_passed = false; // luminosity scans
        if (qa_passed){
          qa.AccumulateCharge();
          qa.AccumulateChargeHL();  
        }
        // REC ELECTRONS
        for (auto& e : electrons_vec) {
            event++;
            SetLorentzVector(el_4vec, e);          // update el_4vec with px,py,pz.
            TLorentzVector Elec = el_4vec;
            TLorentzVector q = beam - el_4vec;
            electron_px = e->getPx(), electron_py = e->getPy(), electron_pz = e->getPz();
            electron_mom = el_4vec.P(), electron_eng = el_4vec.E();
            electron_theta = el_4vec.Theta(), electron_phi = el_4vec.Phi();
            gamma_nu = beam.E() - el_4vec.E();
            double electron_chi2 = e->getChi2Pid();
            double electron_y = target.Dot(q) / target.Dot(beam);
            electron_Q2 = -q.M2();
            electron_W = pow(pow(pr.M(),2)+2*pr.M()*gamma_nu - electron_Q2, 0.5);
            electron_status = e->getStatus();
            electron_sector = e->getSector();
            electron_vz = e->par()->getVz();
            electron_vx = e->par()->getVx();
            double electron_vy = e->par()->getVy();
            electron_PCAL = e->cal(PCAL)->getEnergy();    
            electron_ECAL = e->cal(ECAL)->getEnergy();
            electron_ECIN = e->cal(ECIN)->getEnergy();
            electron_CAL_Tot = electron_ECAL + electron_PCAL + electron_ECIN;
            electron_edge1 = e->traj(6,6)->getEdge();
            electron_edge2 = e->traj(6,18)->getEdge();
            electron_edge3 = e->traj(6,36)->getEdge();
            //
            electron_PCAL_lu = e->cal(PCAL)->getLu();
            electron_PCAL_lv = e->cal(PCAL)->getLv(); 
            electron_PCAL_lw = e->cal(PCAL)->getLw();
            electron_ECin_lu = e->cal(ECIN)->getLu();
            electron_ECin_lv = e->cal(ECIN)->getLv();
            electron_ECin_lw = e->cal(ECIN)->getLw();
            electron_ECout_lu = e->cal(ECOUT)->getLu();
            electron_ECout_lv = e->cal(ECOUT)->getLv();
            electron_ECout_lw = e->cal(ECOUT)->getLw();

            // ElectronPID_ForwardDetector yes?
            // ElectronPID_CalSFcut(int sector, int runnum, double p, double cal_energy) ADD
            // ElectronPID_DiagonalCut(e, &Elec) necessary for RGA, RGC?
            // ElectronEC_cut(electron_ECout_lu ,electron_ECout_lv, electron_ECout_lw) && toglie metà statistica...
            if(ElectronPID_PCAL(electron_PCAL) && KinematicPID_W(electron_W) && std::fabs(electron_chi2) <= 3 &&
            ElectronPID_VertexCut(electron_vz, electron_vx, electron_vy) && ElectronPID_DC(electron_edge1, electron_edge2, electron_edge3, torus) &&
            ElectronPCAL_cut(electron_PCAL_lv, electron_PCAL_lw, electron_sector)  && electron_y < 0.8 && electron_Q2 > 0.95 &&
            ElectronECin_strip(electron_ECin_lv, electron_sector) && ElectronECout_strip(electron_ECout_lu, electron_ECout_lw, electron_sector) && qa_passed){
                /*
                bool pos_cut = false;
                for(auto& pos : positrons_vec){
                    double pos_px = pos->getPx(), pos_py = pos->getPy(), pos_pz = pos->getPz();
                    double pos_mom = pow(pos_px*pos_px + pos_py*pos_py + pos_pz*pos_pz, 0.5);
                    SetLorentzVector(positron, pos);
                    double inv_mass = (el_4vec + positron).M();
                    double pos_theta = positron.Theta(), pos_phi = positron.Phi();
                    double theta_ee = (el_4vec.Vect().Angle(positron.Vect())) * 180 / TMath::Pi();
                    double pos_vz = pos->par()->getVz();
                    //
                    //elpos_selection.Fill(inv_mass);
                    //elpos_vz.Fill(electron_vz - pos_vz);
                    //elpos_theta_ee.Fill(theta_ee);
                    //elpos_MvsTheta.Fill(theta_ee, inv_mass);
                    //elpos_Mvsvz.Fill(electron_vz - pos_vz, inv_mass);
                    if (inv_mass <= 0.1 && theta_ee <= 3){
                      count++;
                      pos_cut = true;
                      continue;
                    }
                }
                */
                // REC KAON+
                for(auto& kp : kaon_vec){
                  //if (pos_cut) continue;
                    SetLorentzVector(kp_4v, kp);
                    kaon_y = target.Dot(q) / target.Dot(beam);
                    TLorentzVector Mx = (Lab - el_4vec - kp_4v);
                    kaon_Mx = Mx.M(), kaon_mom = kp_4v.P(), kaon_eng = kp_4v.E();
                    kaon_z = target.Dot(kp_4v) / target.Dot(q);
                    kaon_s = (beam + target).M2();
                    beta_cms_kaon = Lab.BoostVector(); // beta for cms
                    TLorentzVector kaon_cms_4v = kp_4v;
                    kaon_cms_4v.Boost(-beta_cms_kaon);
                    kaon_xF = (2 * kaon_cms_4v.Pz()) / sqrt(kaon_s);
                    //kaon_polariz = BeamPolarization(N_run, true);
                    kaon_polariz = BeamPolarization_file(N_run);
                    kaon_vz = kp->par()->getVz();
                    clas_pid = kp->getPid();
                    kaon_beta = kp->getBeta();
                    kaon_px = kp->getPx(), kaon_py = kp->getPy(), kaon_pz = kp->getPz();
                    kaon_Q2 = -q.M2();
                    kaon_xB = kaon_Q2 / (2 * target.Dot(q));      
                    kaon_Pt = kp_4v.Perp(q.Vect()); 
                    kaon_theta = kp_4v.Theta(), kaon_phi_lab = kp_4v.Phi(), kaon_eta = kp_4v.PseudoRapidity();
                    kaon_chi2pid = kp->getChi2Pid();
                    kaon_tt = (q - kp_4v).M2();
                    photon_en = q.E(); // cut for gamma for best photon efficiency is 0.35 GeV
                    //
                    double kaon_edge1 = kp->traj(6,6)->getEdge();
                    double kaon_edge2 = kp->traj(6,18)->getEdge();
                    double kaon_edge3 = kp->traj(6,36)->getEdge();
                    // CVT - non riesco ad accedere
                    CVT_edge1 = kp->traj(5, 1)->getEdge();
                    CVT_edge3 = kp->traj(5, 3)->getEdge();
                    CVT_edge5 = kp->traj(5, 5)->getEdge();
                    CVT_edge7 = kp->traj(5, 7)->getEdge();
                    CVT_edge12 = kp->traj(5, 12)->getEdge();
                    // RICH
                    kaon_rich_tr1 = kp->traj(18,1)->getEdge();
                    kaon_rich_tr1_x = kp->traj(18,1)->getX();
                    kaon_rich_tr1_y = kp->traj(18,1)->getY();
                    kaon_rich_tr2 = kp->traj(18,2)->getEdge();
                    kaon_rich_tr2_x = kp->traj(18,2)->getX();
                    kaon_rich_tr2_y = kp->traj(18,2)->getY();
                    kaon_rich_tr3 = kp->traj(18,3)->getEdge();
                    kaon_rich_tr3_x = kp->traj(18,3)->getX();
                    kaon_rich_tr3_y = kp->traj(18,3)->getY();
                    kaon_rich_tr4 = kp->traj(18,4)->getEdge();
                    kaon_rich_tr4_x = kp->traj(18,4)->getX();
                    kaon_rich_tr4_y = kp->traj(18,4)->getY();

                    // Trento Convention - Scattering Plane Axis calculation
                    TVector3 zAxis = q.Vect().Unit();           // virtual photon direction
                    TVector3 l_vect = beam.Vect();              // spatial component of the beam lepton
                    TVector3 sl_vect = el_4vec.Vect();               // spatial component of the scattered lepton           
                    TVector3 yAxis = (l_vect.Cross(sl_vect)).Unit();
                    TVector3 xAxis = yAxis.Cross(zAxis);        // x-Axis of the scattering plane
                    TVector3 Ph_T = kp_4v.Vect() - (kp_4v.Vect() * zAxis) * zAxis; 
                    TVector3 Ph_T_hat = Ph_T.Unit();
                    double kaon_Ph_x = Ph_T.Dot(xAxis);
                    double kaon_Ph_y = Ph_T.Dot(yAxis);
                    kaon_phi_h = TMath::ATan2(kaon_Ph_y, -kaon_Ph_x);
                    //
                    kaon_best_mass = kp->rich()->getBest_mass();
                    kaon_rich_Id = kp->rich()->getId();
                    kaon_rich_PID = kp->rich()->getBest_PID();
                    kaon_rich_ch = kp->rich()->getBest_ch();
                    kaon_rich_chi2 = kp->rich()->getBest_c2(); // Best_c2
                    kaon_rich_Mchi2 = kp->rich()->getMchi2();
                    kaon_rich_RQ = kp->rich()->getRQ();
                    kaon_rich_RL = kp->rich()->getBest_RL();
                    kaon_rich_ntot = kp->rich()->getBest_ntot();

                    kaon_W = sqrt((kaon_eng + electron_eng)*(kaon_eng + electron_eng) - (kaon_mom + electron_mom)*(kaon_mom + electron_mom));
                    kaon_Pt_ratio = kaon_Pt / (kaon_z * sqrt(kaon_Q2));
                    kaon_gamma = (2*kaon_xB*ProtonMass)/(sqrt(kaon_Q2));
                    kaon_epsilon = (1- kaon_y - 0.25 * pow(kaon_gamma, 2) * pow(kaon_y, 2))/(1 - kaon_y + 0.5 * pow(kaon_y,2) + 0.25*pow(kaon_gamma,2)*pow(kaon_y,2));
                    //if(kaon_y <= 0.8 && kaon_mom >= 1 && kaon_z >= 0.2 && kaon_xF > 0){
                    if(HadronPID_Q2(kaon_Q2) && HadronPID_DC(kaon_edge1, kaon_edge2, kaon_edge3, -1) && kaon_y <= 0.8 && kaon_mom >= 1 &&
                    HadronPID_Chi2Pid(kaon_chi2pid) && KinematicPID_Vertex(kaon_vz, torus, 321) && kaon_z >= 0.2 && kaon_xF > 0){
                      //kp_rich_pmt_xy.Fill(kaon_rich_tr1_x, kaon_rich_tr1_y);
                      //kp_rich_aerogel_xy.Fill(kaon_rich_tr2_x, kaon_rich_tr2_y);
                      //kp_rich_aerogel_xy.Fill(kaon_rich_tr3_x, kaon_rich_tr3_y);
                      //kp_rich_aerogel_xy.Fill(kaon_rich_tr4_x, kaon_rich_tr4_y);
                      //kp_PxVsPy.Fill(kaon_px, kaon_py);
                      Kaon_tree.Fill();
                    }
                }
            }
        }
    }
    
    double Q = qa.GetAccumulatedCharge();
    double Q_p = qa.GetAccumulatedChargeHL(+1);
    double Q_m = qa.GetAccumulatedChargeHL(-1);
    cout << "total charge: " << Q << endl;
    cout << "positive charge: " << Q_p << endl;
    cout << "negative charge: " << Q_m << endl;
    /*
    double mb = std::min(C_pp + C_pm, C_mp + C_mm);
    double mt = std::min(C_pp + C_mp, C_pm + C_mm);

    double w_pp = (mb * mt) / ((C_pp + C_pm) * (C_pp + C_mp));
    double w_pm = (mb * mt) / ((C_pp + C_pm) * (C_pm + C_mm));
    double w_mp = (mb * mt) / ((C_mp + C_mm) * (C_pp + C_mp));
    double w_mm = (mb * mt) / ((C_mp + C_mm) * (C_pm + C_mm));

    cout << "wpp: " << w_pp << endl;
    cout << "wpm: " << w_pm << endl;
    cout << "wmp: " << w_mp << endl;  
    cout << "wmm: " << w_mm << endl;
    */



    outFile.Write();
    outFile.Close();
    
    cout << "" << endl;
    cout << "events saved: " << Kaon_tree.GetEntries() << endl;
    cout << "ROOT output file: " << outputFile << endl;
    cout << "" << endl;
}


