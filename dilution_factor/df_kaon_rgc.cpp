#include <numeric>
// RGC NH3 K+ dilution factor - draft, requires clas12root, ROOT, QADB on ifarm.
// ROOT usage: root -l -b -q 'df_kaon_rgc.cpp("summer22","summer2022_NH3_goodRuns.txt",true)'
// Third argument true corresponds to --pf-cell (convert cell PF to bath PF).
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <glob.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "clas12reader.h"
#include "QADB.h"
#include "rcdb_reader.h"
#include "ccdb_reader.h"
using namespace std;
using namespace clas12;

bool ElectronPID_VertexCut(double vz, double vx, double vy, const std::string& campaign){
    if (campaign == "summer22")
        return (vz > -9.0 && vz < 1.0 &&
                std::abs(vx) <= 2 && std::abs(vy) <= 2);
    return (vz > -7.5 && vz < 2.5 &&
            std::abs(vx) <= 2 && std::abs(vy) <= 2);
}

bool ElectronPID_PCAL(double pcal){
    return pcal > 0.07;
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
      if(-20 < vz && vz < 10) signal = true; // summer 22 -10 to +2.5 for all, even fall22
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



struct Period {string key, directory; double pf, dpf; vector<int> carbon, ch2, mt, foil;};
Period period(const string& s) {
 const string base="/cache/clas12/rg-c/production/";
 if(s=="summer22") return {s,base+"summer22/pass1/10.5gev/",.5364,.0003,{16290,16291,16292,16293,16296,16297,16700,16701,16702,16704},{16298,16299,16300,16302,16303},{16309},{16194}};
 if(s=="fall22" || s=="fall22_pre" || s=="fall22_post") return {s,base+"fall22/pass1/",s=="fall22_pre"?.5133:(s=="fall22_post"?.5135:.5048),s=="fall22_post"?.0001:.0003,{17138,17139,17140,17170},{17132,17133,17134,17135},{16975},{16872}};
 if(s=="spring23") return {s,base+"spring23/pass1/",.5152,.0003,{17566,17567,17569,17744,17745,17746},{17625,17626,17627,17635,17637,17638},{17766,17767,17768},{17763,17764,17765}};
 throw runtime_error("period must be summer22, fall22, fall22_pre, fall22_post, spring23");
}
map<int,double> loadNH3(const string& file) {
 ifstream in(file); if(!in) throw runtime_error("cannot open "+file);
 map<int,double> runs; string line;
 while(getline(in,line)) {if(line.empty()||line[0]=='#')continue; replace(line.begin(),line.end(),',',' '); istringstream ss(line);int r;double q;
  if(ss>>r>>q && q>0) runs[r]=q;}
 return runs;
}
vector<string> filesForRun(const string& dir, const string& target, int run) {
 char num[32]; snprintf(num,sizeof(num),"%06d",run);
 string pat=dir+target+"/dst/train/sidisdvcs/*"+num+"*.hipo";
 glob_t g{}; vector<string> out;
 if(glob(pat.c_str(),0,nullptr,&g)==0) for(size_t i=0;i<g.gl_pathc;++i)out.emplace_back(g.gl_pathv[i]);
 globfree(&g); return out;
}
int bin2d(double x,double q) {
 const vector<pair<double,double>> qs={{1,3},{3,5},{5,7},{7,11}};
 const vector<vector<double>> xs={{0,.12,.18,.24,.30,.8},{0,.24,.30,.36,.44,.8},{0,.44,.8},{0,.55,.8}};
 int idx=0; for(size_t i=0;i<qs.size();++i){
  if(q>=qs[i].first && q<(i==qs.size()-1?qs[i].second+1e-10:qs[i].second)){
   for(size_t j=0;j+1<xs[i].size();++j)if(x>=xs[i][j] && x<xs[i][j+1])return idx+(int)j;
   return -1;
  }idx+=(int)xs[i].size()-1;
 }return -1;
}
// z-Pt SIDIS bins (1-based), ordered by Pt slice then z slice.
const vector<pair<double,double>> PT_EDGES={{0,.25},{.25,.5},{.5,.8},{.8,1.4}};
const vector<vector<pair<double,double>>> Z_EDGES={
 {{.2,.3},{.3,.4},{.4,.5},{.5,.6},{.6,.7},{.7,.8},{.8,1.0}},
 {{.2,.3},{.3,.4},{.4,.5},{.5,.6},{.6,.7},{.7,.8},{.8,1.0}},
 {{.2,.3},{.3,.4},{.4,.5},{.5,.6},{.6,.7},{.7,.8},{.8,1.0}},
 {{.2,.5},{.5,.6},{.6,.7},{.7,1.0}}
};
constexpr int NZPT=25;
int binZPt(double z,double pt) {
 int offset=0;
 for(size_t ip=0;ip<PT_EDGES.size();++ip){
  // Half-open Pt intervals prevent double counting at Pt boundaries.
  if(pt>=PT_EDGES[ip].first && pt<PT_EDGES[ip].second){
   for(size_t iz=0;iz<Z_EDGES[ip].size();++iz)
    if(z>=Z_EDGES[ip][iz].first && z<Z_EDGES[ip][iz].second)
     return offset+(int)iz;
   return -1;
  }
  offset+=(int)Z_EDGES[ip].size();
 }
 return -1;
}
constexpr int NX=20, N2=14;
struct Sample {vector<double> x=vector<double>(NX,0), two=vector<double>(N2,0), zpt=vector<double>(NZPT,0); double charge=0; long long events=0,selected=0;};
void configureQA(QA::QADB& qa) {
 for(auto s:{"TotalOutlier","TerminalOutlier","MarginalOutlier","SectorLoss","LowLiveTime","Misc","ChargeHigh","ChargeNegative","ChargeUnknown","PossiblyNoBeam"})qa.CheckForDefect(s);
 for(int r:{16194,16089,16185,16308,16184,16307,16309,16872,16975,16979,17763,17764,17765,17766,17767,17768,17179,17180,17181,17182,17183,17188,17189,17252})qa.AllowMiscBit(r);
}
bool acceptedQA(int r) {
 if(r>16600&&r<16700)return false;
 if(r>17768&&r<=17811)return false;
 for(int x:{17331,16987,17079,17190,17639,16850,16851,16852,16855,16879})if(r==x)return false;
 return true;
}
Sample process(const Period& p,const string& target,const vector<int>& runs) {
 Sample result;
 for(int run:runs) {
  const double beforeRunX=accumulate(result.x.begin(),result.x.end(),0.0);
  auto paths=filesForRun(p.directory,target,run);
  if(paths.empty()){cerr<<"WARNING missing HIPO: "<<p.key<<" "<<target<<" "<<run<<"\n";continue;}
  clas12root::HipoChain chain;for(const auto& path:paths)chain.Add(path.c_str());
  auto conf=chain.GetC12Reader();chain.SetReaderTags({0});conf->addAtLeastPid(11,1);
  auto& c12=chain.C12ref();
  clas12databases::SetRCDBRemoteConnection(); clas12databases db;c12->connectDataBases(&db);
  QA::QADB qa("latest");configureQA(qa);
  TLorentzVector beam(0,0,10.6,10.6), proton(0,0,0,.9382720813), el(0,0,0,.00051099895), k(0,0,0,.493677);
  while(chain.Next()) {
   ++result.events;
   int r=c12->runconfig()->getRun(), ev=c12->runconfig()->getEvent();
   if(r!=run)continue;
   bool good=acceptedQA(r)&&qa.Pass(r,ev);
   if(!good)continue;
   qa.AccumulateCharge(); // QADB deduplicates charge by QA interval; NOT per event.
   vector<clas12::region_part_ptr> electrons,kaons;
   for(auto& part:c12->getDetParticles()) {
    int pid=part->getPid();
    if(pid==11 || (part->rich()&&part->rich()->getBest_PID()==11))electrons.push_back(part);
    if(pid==321 || (part->rich()&&part->rich()->getBest_PID()==321))kaons.push_back(part);
   }
   for(auto& e:electrons) {
    if(!e->cal(PCAL)||!e->cal(ECIN)||!e->cal(ECOUT)||!e->traj(6,6)||!e->traj(6,18)||!e->traj(6,36))continue;
    el.SetXYZM(e->getPx(),e->getPy(),e->getPz(),.00051099895);
    TLorentzVector photon=beam-el;
    double q2=-photon.M2(), nu=photon.E(), y=nu/beam.E();
    double w2=proton.M2()+2*proton.M()*nu-q2;
    if(w2<=0 || !ElectronPID_PCAL(e->cal(PCAL)->getEnergy()) || !KinematicPID_W(sqrt(w2)) || fabs(e->getChi2Pid())>3 ||
       !ElectronPID_VertexCut(e->par()->getVz(),e->par()->getVx(),e->par()->getVy(),p.key) ||
       !ElectronPID_DC(e->traj(6,6)->getEdge(),e->traj(6,18)->getEdge(),e->traj(6,36)->getEdge(),-1) ||
       !ElectronPCAL_cut(e->cal(PCAL)->getLv(),e->cal(PCAL)->getLw(),e->getSector()) ||
       y>=.8 || q2<=.95 ||
       !ElectronECin_strip(e->cal(ECIN)->getLv(),e->getSector()) ||
       !ElectronECout_strip(e->cal(ECOUT)->getLu(),e->cal(ECOUT)->getLw(),e->getSector()))continue;
    double xb=q2/(2*proton.M()*nu);
    if(!(xb>=0&&xb<.8))continue;
    for(auto& kp:kaons) {
     if(!kp->traj(6,6)||!kp->traj(6,18)||!kp->traj(6,36))continue;
     k.SetXYZM(kp->getPx(),kp->getPy(),kp->getPz(),.493677);
     double z=k.E()/nu;
     TLorentzVector gammaP=photon+proton, ks=k, qs=photon;
     if(gammaP.M2()<=0)continue;
     TVector3 boost=gammaP.BoostVector();ks.Boost(-boost);qs.Boost(-boost);
     double xf=2*ks.Vect().Dot(qs.Vect().Unit())/gammaP.M();
     if(!HadronPID_Q2(q2) || !HadronPID_DC(kp->traj(6,6)->getEdge(),kp->traj(6,18)->getEdge(),kp->traj(6,36)->getEdge(),-1) ||
        k.P()<1 || !HadronPID_Chi2Pid(kp->getChi2Pid()) ||
        !KinematicPID_Vertex(kp->par()->getVz(),-1,321) || z<.2 || xf<=0)continue;
     // Pt transverse to virtual-photon direction in the target rest frame.
     double pt=k.Vect().Perp(photon.Vect());
     int ix=int(xb/.04), iq=bin2d(xb,q2), izpt=binZPt(z,pt);
     if(ix>=0&&ix<NX)result.x[ix]++;
     if(iq>=0)result.two[iq]++;
     if(izpt>=0)result.zpt[izpt]++;
     ++result.selected;
    }
   }
  }
  const double runCharge=qa.GetAccumulatedCharge();
  const double runCountsX=accumulate(result.x.begin(),result.x.end(),0.0);
  // Difference from pre-run totals gives per-run selected K+ candidates in xB bins.
  const double runN=runCountsX-beforeRunX;
  result.charge+=runCharge;
  cout<<"[DF_RUN] "<<p.key<<" target="<<target<<" run="<<run
      <<" files="<<paths.size()<<" Q="<<runCharge<<" N_xB="<<runN
      <<" yield_xB="<<(runCharge>0?runN/runCharge:0.0)<<"\n";
 }
 return result;
}
// Equation (9) of RGC note; uncertainty from independent Poisson counts, charges fixed.
constexpr double A=-.888233,B=.170484,C=-.282251,G=.505587;
struct Value{double df=0,err=0;bool valid=false;};
Value compute(const array<double,5>& counts,const array<double,5>& charges,double pf) {
 for(double q:charges)if(q<=0)return {};
 array<double,5> n{};for(int i=0;i<5;++i)n[i]=counts[i]/charges[i];
 if(n[0]<=0)return {};
 double t=n[2]+A*n[1]+B*n[3]+C*n[4];
 double df=pf*t/(G*n[0]);
 double v=0;
 // derivatives wrt raw Poisson counts
 for(int i=0;i<5;++i){double d=0;
  if(i==0)d=-df/counts[0];
  if(i==1)d=pf*A/(G*n[0]*charges[1]);
  if(i==2)d=pf/(G*n[0]*charges[2]);
  if(i==3)d=pf*B/(G*n[0]*charges[3]);
  if(i==4)d=pf*C/(G*n[0]*charges[4]);
  v+=d*d*counts[i];
 }
 return {df,sqrt(v),true};
}
int df_kaon_rgc(const char* campaign="summer22",
                const char* nh3_file="summer2022_NH3_goodRuns.txt",
                bool pfcell=false)try{
 Period p=period(campaign);auto nh3=loadNH3(nh3_file);
 if(p.key=="fall22_pre" || p.key=="fall22_post"){
  const bool pre=(p.key=="fall22_pre");
  for(auto it=nh3.begin();it!=nh3.end();){
   const bool keep=pre?(it->first<17185):(it->first>=17185);
   if(!keep)it=nh3.erase(it);else ++it;
  }
  if(nh3.empty())throw runtime_error("no NH3 runs in selected Fall22 solenoid subset");
  cerr<<"WARNING: DIAGNOSTIC ONLY: auxiliary target runs are unchanged and NOT split/matched by solenoid polarity. "
      <<"Do not interpret pre/post differences as a validated physics comparison.\n";
 }

 double pf=p.pf*(pfcell?5.0/5.86:1.0);
 cerr<<"PF for Eq.9="<<pf<<" ("<<(pfcell?"input PFcell converted":"input PFbath assumed")<<")\n";
 if(p.key=="fall22")cerr<<"NOTE Fall22 solenoid flips at run 17185; this draft uses combined PF=0.5048. Split extraction requires configuration-matched auxiliary samples.\n";
 vector<int> nr;for(auto& [r,q]:nh3)nr.push_back(r);
 const vector<vector<int>> groups={nr,p.carbon,p.ch2,p.mt,p.foil};
 const vector<string> names={"NH3","C","CH2","ET","ET"};
 array<Sample,5> data;
 for(int i=0;i<5;++i)data[i]=process(p,names[i],groups[i]);
 array<double,5> charge{};for(int i=0;i<5;++i)charge[i]=data[i].charge;
 for(int i=0;i<5;++i)if(charge[i]<=0)throw runtime_error("zero QA charge for "+names[i]+"; check files/QADB");

 cout<<"\n[DF_CHECK] ===== "<<p.key<<" PF="<<pf<<" G="<<G<<" =====\n";
 const array<string,5> labels={"NH3","C","CH2","MT","FOIL"};
 for(int i=0;i<5;++i){
  double total=accumulate(data[i].x.begin(),data[i].x.end(),0.0);
  cout<<"[DF_TARGET] "<<labels[i]<<" Q="<<charge[i]
      <<" N_xB="<<total<<" n_xB="<<total/charge[i]
      <<" files_events="<<data[i].events<<" selected_K="<<data[i].selected<<"\n";
 }
 cout<<"[DF_CHECK] binning bin N_NH3 N_C N_CH2 N_MT N_F"
     <<" n_NH3 n_C n_CH2 n_MT n_F term_CH2 term_C term_MT term_F T R DF\n";

 string stem="DF_Kplus_"+p.key;
 TFile root((stem+".root").c_str(),"RECREATE");
 TH1D h1("DF_xB","K+ dilution factor;x_{B};DF",NX,0,.8);
 TH1D h2("DF_xQ2_bin","K+ dilution factor;custom xB-Q2 bin (1-based);DF",N2,.5,N2+.5);
 TH1D h3("DF_zPt_bin","K+ dilution factor;custom z-P_{T} bin (1-based);DF",NZPT,.5,NZPT+.5);
 ofstream csv(stem+".csv");csv<<"period,binning,bin,xB_min,xB_max,Q2_min,Q2_max,DF,stat_error,N_NH3,N_C,N_CH2,N_MT,N_FOIL,z_min,z_max,Pt_min,Pt_max\n";
 auto write=[&](const string& type,int index,double xmin,double xmax,double qmin,double qmax,TH1D& h,int mode,double zmin=-1,double zmax=-1,double ptmin=-1,double ptmax=-1){
  array<double,5> counts{};for(int i=0;i<5;++i)counts[i]=mode==2?data[i].two[index]:(mode==3?data[i].zpt[index]:data[i].x[index]);
  Value v=compute(counts,charge,pf);
  array<double,5> n{};
  for(int i=0;i<5;++i)n[i]=counts[i]/charge[i];
  const double termCH2=n[2],termC=A*n[1],termMT=B*n[3],termF=C*n[4];
  const double T=termCH2+termC+termMT+termF;
  const double R=n[0]>0?T/n[0]:0.0;
  cout<<"[DF_BIN] "<<type<<" "<<index+1;
  for(double c:counts)cout<<" "<<c;
  for(double v:n)cout<<" "<<v;
  cout<<" "<<termCH2<<" "<<termC<<" "<<termMT<<" "<<termF
      <<" "<<T<<" "<<R<<" "<<(v.valid?v.df:0.0)<<"\n";

  if(v.valid){h.SetBinContent(index+1,v.df);h.SetBinError(index+1,v.err);}
  csv<<p.key<<","<<type<<","<<index+1<<","<<xmin<<","<<xmax<<","<<qmin<<","<<qmax<<",";
  if(v.valid)csv<<v.df<<","<<v.err;else csv<<"nan,nan";
  for(double c:counts)csv<<","<<c;
  csv<<","<<zmin<<","<<zmax<<","<<ptmin<<","<<ptmax<<"\n";
 };
 for(int i=0;i<NX;++i)write("xB",i,.04*i,.04*(i+1),.95,1e9,h1,1);
 const vector<pair<double,double>> q={{1,3},{3,5},{5,7},{7,11}};
 const vector<vector<double>> x={{0,.12,.18,.24,.30,.8},{0,.24,.30,.36,.44,.8},{0,.44,.8},{0,.55,.8}};
 int j=0;for(size_t iq=0;iq<q.size();++iq)for(size_t ix=0;ix+1<x[iq].size();++ix,++j)write("xB_Q2",j,x[iq][ix],x[iq][ix+1],q[iq].first,q[iq].second,h2,2);
 int izpt=0;
 for(size_t ip=0;ip<PT_EDGES.size();++ip)
  for(const auto& zr:Z_EDGES[ip]){
   write("z_Pt",izpt,-1,-1,-1,-1,h3,3,zr.first,zr.second,PT_EDGES[ip].first,PT_EDGES[ip].second);
   ++izpt;
  }
 h1.Write();h2.Write();h3.Write();root.Close();
 cerr<<"Wrote "<<stem<<".root and "<<stem<<".csv\n";
 return 0;
}catch(const exception& e){cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
