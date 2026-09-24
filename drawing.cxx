#define MAXDATA 10000
#include <iostream>
#include <vector>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TError.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TSystem.h>

#include "definition.h"
#include "utility.h"
#include "computing.h"

bool saveRootFile=true; //set true to save the root file "displayMonitor.root"

static TH1D *hTime;
static vector <TH1D*> hHitCoinc;
static vector <TH1D*> hHitStart;
static vector <TH1D*> hHitEnd;
static vector <TH1D*> hHitDur;
static vector <TH2D*> hHitCorrDur;
static vector <TH2D*> hNotCalibVsCh;
static vector <TH2D*> hCalibVsCh;

static TH1D *hInterpolateGasRing=0;
static TH1D *hInterpolateNotGasRing=0;
static TH1D *hInterpolateIsProtonAndGasRing=0;
static TH1D *hInterpolateIsKaonAndGasRing=0;
static TH1D *hInterpolateIsPionAndNotGasRing=0;

static vector <TH1D*> vGasRadPMT;
static vector <TH1D*> vAeroRadPMT;
static vector <TH1D*> vGasnRadPMT;
static vector <TH1D*> vAeronRadPMT;

static TH1D* hDistanceExpectedPionAero = 0;
static TH1D* hDistanceExpectedKaonAero = 0;
static TH1D* hDistanceExpectedProtonAero = 0;
static TH1D* hDistanceExpectedPionGas = 0;
static TH1D* hDistanceExpectedKaonGas = 0;
static TH1D* hDistanceExpectedProtonGas = 0;

static int countEvent;

static TH2D *hHitCorr;
static TH2D *hMap;
static TH2D *htMap;
static TH2D *hTimeVsRadius; //CHECK THE REASON THIS CAUSES A SEGMENTATION FAULT.
static TH2D *hMapNC;
static TH2D *hnMap;
static vector<TH2D*> vPiMap;
static vector<TH2D*> vKMap;
static vector<TH2D*> vPrMap;
static TH1D *hRadius;
static vector<TH1D*> vRadiusMM;
static vector<TH1D*> vnRadiusMM;
static vector<TH1D*> vnRadius;
static vector<TH1D*> vPiRadius;
static vector<TH1D*> vKRadius;
static vector<TH1D*> vPrRadius;
static vector<TH1D*> vSpnPiRadius;
static vector<TH1D*> vSpnKRadius;
static vector<TH1D*> vSpnPrRadius;

static vector<TH1D*> vRadiusRMS;
static vector<TH1D*> vTimeRMS;
static vector<TH1D*> vRadiusRSD;
static vector<TH1D*> vTimeRSD;

static TH2D *hNewCenterGas;
static TH2D *hNewCenterAero;
static TH1D *hNewCenterGasX;
static TH1D *hNewCenterGasY;
static TH1D *hNewCenterAeroX;
static TH1D *hNewCenterAeroY;
static TH2D *hGasdXVsAngle;
static TH2D *hGasdYVsAngle;
static TH2D *hAerodXVsAngle;
static TH2D *hAerodYVsAngle;

static TH1D *hPiRad;
static TH1D *hKRad;
static TH1D *hPrRad;

static TH1D *htrig;
static TH1D *hx474;
static TH1D *hx519;
static TH1D *hx537;
static TH1D *hBeamCh;

static TH1D *hEdge;
static TH1D *hHitReco;

static TH2D *hUpGEM;
static TH2D *hDnGEM;
static TH2D *hBeam;
static TH2D *hBeamTheta;

static vector<TH1D*> vrsdRadius;
static vector<TH1D*> vrsdTime;

static vector<TH1D*> vspRadius;
static vector<TH1D*> vspTime;
static vector<TH1D*> vspPhoton;
static vector<TH1D*> vspnRadius;
static vector<TH1D*> vspnTime;
static vector<TH1D*> vspnPhoton;
static vector<TH1D*> vcutRadius;
static vector<TH1D*> vcutTime;
static vector<TH1D*> vcutPhoton;

static vector<TH1D*> vspSigPhoGas;
static vector<TH1D*> vspSigPhoAero;
static vector<TH1D*> vspnSigPhoGas;
static vector<TH1D*> vspnSigPhoAero;
static vector<TH1D*> vcutSigPhoGas;
static vector<TH1D*> vcutSigPhoAero;
static TH1D *hspSigVsPhoGas;
static TH1D *hspnSigVsPhoGas;
static TH1D *hcutSigVsPhoGas;
static TH1D *hspSigVsPhoAero;
static TH1D *hspnSigVsPhoAero;
static TH1D *hcutSigVsPhoAero;

static TH2D *hCalibVsTime;
static TH2D *hTimeVsCh;
static TH2D *hTimeWalkVsTime;
static TH2D *hStartVsDur;
static TH1D *hSingPhotRad;
static TH1D *hTimeRespectTrigger;

static vector<vector<TH1D*>> vNotTime;
static vector<vector<TH1D*>> vCalTime;

static TH1D *hSimPhotonLambdaGas;
static TH1D *hSimPhotonLambdaAero;

static int minPhoGas=2;
static int maxPhoGas=50;
static int minPhoAero=2;
static int maxPhoAero=15;

static TH1D *hA;//Numero di fotoni aerogel nel caso si veda ring sul gas.
static TH1D *hG;//Numero di fotoni gas nel caso si veda ring sull'aerogel.

static TH2D* hBaryCenter;
static TH2D* hBaryVsCenterGas;

const char *labelCh[3]={"Pi","K","Pr"};

double tmpSel = 187;

double piMass=0.13957;
double kMass=0.493677;
//double pMass=0.938272;
double pMass=0.000510999;
double eChPionAero = 0;
double eChKaonAero = 0;
double eChProtonAero = 0;
double eChPionGas = 0;
double eChKaonGas = 0;
double eChProtonGas = 0;


//---------------------------------------------------
void fillHistoMon(THeader *run){
  //---------------------------------------------------
  TCanvas *c = new TCanvas("c","c",1600,900);
  c->Draw();
  c->Divide(2);
  //c->Print(Form("lowRadEvent_%d.pdf[",run->runNum));

  gErrorIgnoreLevel=kWarning;
  TString fName=Form("%s/processed_data/integrated_dRICH_GEM_data/run_%04d_integrated.root",run->suite.c_str(),run->runNum);
  TFile *fIn = new TFile (fName,"READ");
  TTree *t = (TTree*) fIn->Get("dRICH");

  int evt, nedge, pol[MAXDATA], pmt[MAXDATA], spPhoton[10], spnPhoton[10], cutPhoton[10], fiber[MAXDATA], ch[MAXDATA], time[MAXDATA], otime[MAXDATA];
  double y[MAXDATA], x[MAXDATA], r[MAXDATA], rmm[MAXDATA], nx[MAXDATA], ny[MAXDATA], nr[MAXDATA], nrmm[MAXDATA], nt[MAXDATA], nttw[MAXDATA], dur[MAXDATA], rsdRadius[MAXDATA], rsdTime[MAXDATA], spRadius[10], spRadiusmm[10], spTime[10], spnRadius[10], spnTime[10], cutRadius[10], cutTime[10], xNCin, yNCin, xNCout, yNCout, photonWavelength[MAXDATA], rmsRadius[10], rmsPhoton[10], rmsTime[10];
  float gx0, gy0, gx1, gy1, gxa, gya, gxtheta, gytheta;
  bool trigSig[MAXDATA], goodHit[MAXDATA], goodPhoton[MAXDATA], coincPhoton[MAXDATA],externalPhoton[MAXDATA], goodSP[10], goodSPN[10], goodCUT[10];
  double trigtime, x474time, x519time, x537time;
  t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("trigtime",&trigtime);
  t->SetBranchAddress("x474time",&x474time);
  t->SetBranchAddress("x519time",&x519time);
  t->SetBranchAddress("x537time",&x537time);
  t->SetBranchAddress("nedge",&nedge);
  t->SetBranchAddress("pol",&pol);
  t->SetBranchAddress("fiber",&fiber);
  t->SetBranchAddress("ch",&ch);
  t->SetBranchAddress("pmt",&pmt);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("r",&r);
  t->SetBranchAddress("rmm",&rmm);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("otime",&otime);
  t->SetBranchAddress("nx",&nx);
  t->SetBranchAddress("ny",&ny);
  t->SetBranchAddress("nr",&nr);
  t->SetBranchAddress("nrmm",&nrmm);
  t->SetBranchAddress("nt",&nt);
  t->SetBranchAddress("nttw",&nttw);
  t->SetBranchAddress("dur",&dur);
  t->SetBranchAddress("trigSig",&trigSig);
  t->SetBranchAddress("goodHit",&goodHit);
  t->SetBranchAddress("goodPhoton",&goodPhoton);
  t->SetBranchAddress("coincPhoton",&coincPhoton);
  if(run->sensor=="SIMULATION")	t->SetBranchAddress("photonWavelength",&photonWavelength);
  t->SetBranchAddress("externalPhoton",&externalPhoton);
  t->SetBranchAddress("rmsRadius",&rmsRadius);
  t->SetBranchAddress("rmsTime",&rmsTime);
  t->SetBranchAddress("rsdRadius",&rsdRadius);
  t->SetBranchAddress("rsdTime",&rsdTime);
  t->SetBranchAddress("gx0",&gx0);
  t->SetBranchAddress("gy0",&gy0);
  t->SetBranchAddress("gx1",&gx1);
  t->SetBranchAddress("gy1",&gy1);
  t->SetBranchAddress("gxa",&gxa);
  t->SetBranchAddress("gya",&gya);
  t->SetBranchAddress("gxtheta",&gxtheta);
  t->SetBranchAddress("gytheta",&gytheta);
  t->SetBranchAddress("spRadius",&spRadius);
  t->SetBranchAddress("spRadiusmm",&spRadiusmm);
  t->SetBranchAddress("spTime",&spTime);
  t->SetBranchAddress("spPhoton",&spPhoton);
  t->SetBranchAddress("goodSP",&goodSP);
  t->SetBranchAddress("spnRadius",&spnRadius);
  t->SetBranchAddress("spnTime",&spnTime);
  t->SetBranchAddress("spnPhoton",&spnPhoton);
  t->SetBranchAddress("goodSPN",&goodSPN);
  t->SetBranchAddress("cutRadius",&cutRadius);
  t->SetBranchAddress("cutTime",&cutTime);
  t->SetBranchAddress("cutPhoton",&cutPhoton);
  t->SetBranchAddress("goodCUT",&goodCUT);
  t->SetBranchAddress("xNCin",&xNCin);
  t->SetBranchAddress("yNCin",&yNCin);
  t->SetBranchAddress("xNCout",&xNCout);
  t->SetBranchAddress("yNCout",&yNCout);





  double n[3] = {1+run->px474*(nN2-1), 1+run->px519*(nCO2-1), 1+run->px537*(nCO2-1)};
  double thsPi[3] = {mPi/sqrt(n[0]*n[0]-1), mPi/sqrt(n[1]*n[1]-1), mPi/sqrt(n[2]*n[2]-1)};
  double thsK[3] = {mK/sqrt(n[0]*n[0]-1), mK/sqrt(n[1]*n[1]-1), mK/sqrt(n[2]*n[2]-1)};
  double thsPr[3] = {mPr/sqrt(n[0]*n[0]-1), mPr/sqrt(n[1]*n[1]-1), mPr/sqrt(n[2]*n[2]-1)};

  bool x474Pi=false, x474K=false, x474Pr=false;
  if(run->energyGeV > thsPi[0])x474Pi=true;
  if(run->energyGeV > thsK[0])x474K=true;
  if(run->energyGeV > thsPr[0])x474Pr=true;

  bool x519Pi=false, x519K=false, x519Pr=false;
  if(run->energyGeV > thsPi[1])x519Pi=true;
  if(run->energyGeV > thsK[1])x519K=true;
  if(run->energyGeV > thsPr[1])x519Pr=true;

  bool x537Pi=false, x537K=false, x537Pr=false;
  if(run->energyGeV > thsPi[2])x537Pi=true;
  if(run->energyGeV > thsK[2])x537K=true;
  if(run->energyGeV > thsPr[2])x537Pr=true;

  cout <<Form("Pi thresholds: %lf %lf %lf Status %d %d %d\n",thsPi[0],thsPi[1],thsPi[2],x474Pi,x519Pi,x537Pi);
  cout <<Form("K thresholds: %lf %lf %lf Status %d %d %d\n",thsK[0],thsK[1],thsK[2],x474K,x519K,x537K);
  cout <<Form("Pr thresholds: %lf %lf %lf Status %d %d %d\n",thsPr[0],thsPr[1],thsPr[2],x474Pr,x519Pr,x537Pr);
  //cin.get();

  //Expected Cherenkov angles
  double beamMomentum = (double)run->energyGeV;
  double k = pow(beamMomentum/piMass,2);
  double betaPi = (k-1)/k; 
  k = pow(beamMomentum/kMass,2);
  double betaK = (k-1)/k; 
  k = pow(beamMomentum/pMass,2);
  double betaP = (k-1)/k; 
  double nAero = 1.02;
  double nGas = 1.00085;
  eChPionAero = 1000*acos(1/nAero/betaPi);
  eChPionGas = 1000*acos(1/nGas/betaPi);
  eChKaonAero = 1000*acos(1/nAero/betaK);
  eChKaonGas = 1000*acos(1/nGas/betaK);
  eChProtonAero = 1000*acos(1/nAero/betaP);
  eChProtonGas = 1000*acos(1/nGas/betaP);

  countEvent=0;
  cout <<Form("Filling monitoring histograms for run %d\n",run->runNum);
  for(int i = 0; i < t->GetEntries(); i++){
    if(i%100==0)printProgress((double)i/t->GetEntries());
    t->GetEntry(i);
   // cout <<Form("%d %d %lf %lf %lf %lf\n",i, evt, rmsRadius[9],gxa,gya,sqrt(gxtheta*gxtheta+gytheta*gytheta));

    double xc = 0;
    double yc = 0;
    int npho = 0;
    for(int j=0; j<nedge; j++){
      if(!goodPhoton[j]) continue;
      xc += nx[j];
      yc += ny[j];
      npho++;
    }
    if(npho>0){
      xc /= npho;
      yc /= npho;
    }
    hBaryCenter->Fill(xc,yc);
    double dist = sqrt((xc-xNCin)*(xc-xNCin) + (yc-yNCin)*(yc-yNCin));
    hBaryVsCenterGas->Fill(xc-xNCin,yc-yNCin);

    if(spPhoton[9] > 0)continue;
    if(goodSPN[4])hA->Fill(spnPhoton[9]);
    if(goodSPN[9])hG->Fill(spnPhoton[4]);
    if(rmsRadius[9]>4)continue;
    //GEM CUTS
    //printf("GEM %7.2f [%7.2f] %7.2f [%7.2f] %7.3f [%7.3f] \n",abs(gxa),GEM_CUT_X,abs(gya),GEM_CUT_Y,sqrt(gxtheta*gxtheta+gytheta*gytheta),GEM_CUT_R);
    //cout <<"Run Event "
    if(spPhoton[1]+spRadiusmm[2]>1 && spPhoton[3]+spRadiusmm[0]>1)hGasdXVsAngle->Fill(gxtheta,((spRadiusmm[3]+spRadiusmm[0]-spRadiusmm[1]-spRadiusmm[2])/2));
    if(spPhoton[0]+spRadiusmm[1]>1 && spPhoton[2]+spRadiusmm[3]>1)hGasdYVsAngle->Fill(gytheta,((spRadiusmm[0]+spRadiusmm[1]-spRadiusmm[2]-spRadiusmm[3])/2));
    if(spPhoton[6]>1 && spPhoton[8]>1)hAerodXVsAngle->Fill(gxtheta,((spRadiusmm[6]-spRadiusmm[8])/2));
    if(spPhoton[5]>1 && spPhoton[7]>1)hAerodYVsAngle->Fill(gytheta,((spRadiusmm[5]-spRadiusmm[7])/2));
    if(APPLY_GEM_CUT==true){
      //if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
      //if(abs(gx0) > GEM_CUT_X || abs(gy0) > GEM_CUT_Y)continue;
      //if(abs(gx1) > GEM_CUT_X || abs(gy1) > GEM_CUT_Y)continue;
      if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;
    }
    hNewCenterGas->Fill(xNCin,yNCin);
    hNewCenterGasX->Fill(xNCin);
    hNewCenterGasY->Fill(yNCin);
    hNewCenterAero->Fill(xNCout,yNCout);;
    hNewCenterAeroX->Fill(xNCout);
    hNewCenterAeroY->Fill(yNCout);


    /*auto kRadius=sqrt(pow(spRadiusmm[5]-spRadiusmm[7],2)+pow(spRadiusmm[6]-spRadiusmm[8],2));
      auto kAngle=sqrt(gxtheta*gxtheta+gytheta*gytheta);
      if(spPhoton[5]>0 && spPhoton[7]>0 && spPhoton[6]>0 && spPhoton[8]>0)hAerodYVsAngle->Fill(kAngle,kRadius);*/

    //if(abs(gxa) > GEM_CUT_X || abs(gya) > GEM_CUT_Y)continue;
    //if(sqrt(gxtheta*gxtheta+gytheta*gytheta) > GEM_CUT_R)continue;

    countEvent++;

    bool isPion   = false;
    bool isKaon   = false;
    bool isProton = false;
    if(run->beamChLogic==2){
      if(x474Pi == true && x519Pi==true && x537Pi==true && x474time!= 0 && x519time!=0 && x537time!=0){
        hBeamCh->Fill(0.,1);
        isPion = true;
      }
      if(x474K == false && x519K==false && x537K==true && x474time== 0 && x519time==0 && x537time!=0){
        hBeamCh->Fill(1.,1);
        isKaon = true;
      }
      if(x474Pr == false && x519Pr==false && x537Pr==false && x474time== 0 && x519time==0 && x537time==0){
        hBeamCh->Fill(2.,1);
        isProton = true;
      }
    }else if(run->beamChLogic==3){
      if(x474Pi == true && x519Pi==true && x537Pi==true && x474time!= 0 && x519time!=0 && x537time!=0){
        hBeamCh->Fill(0.,1);
        isPion = true;
      }
      if(x474K == false && x519K==true && x537K==true && x474time== 0 && x519time!=0 && x537time!=0){
        hBeamCh->Fill(1.,1);
        isKaon = true;
      }
      if(x474Pr == false && x519Pr==false && x537Pr==true && x474time== 0 && x519time==0 && x537time!=0){
        hBeamCh->Fill(2.,1);
        isProton = true;
      }
    }else if(run->beamChLogic==4){
      if(x474Pi && x474time && (x519Pi && x519time)){
        hBeamCh->Fill(0.,1);
        isPion=true;
      }else{
        hBeamCh->Fill(2.,1);
        isProton=true;
      }
    }else if(run->beamChLogic==5){
      if(x474Pi && x474time && (x519Pi && x519time)){
        hBeamCh->Fill(0.,1);
        isPion=true;
      }				
      if(!x519K && !x519time && (x474K && x474time)){
        hBeamCh->Fill(1.,1);
        isKaon=true;
      }				
      if(!x474Pr && !x474time && (!x519Pr && !x519time)){
        hBeamCh->Fill(2.,1);
        isProton=true;
      }
    }else if(run->beamChLogic==6){
      if(x474Pi && x474time){
        hBeamCh->Fill(0.,1);
        isPion=true;
      }				
      if(!x474Pr && !x474time){
        hBeamCh->Fill(2.,1);
        isProton=true;
      }

    }else if(run->beamChLogic==7){
      if((x474Pi && x474time) || (x519Pi && x519time)){
        hBeamCh->Fill(0.,1);
        isPion=true;
      }else{
        hBeamCh->Fill(2.,1);
        isProton=true;
      }
    }else if(run->beamChLogic==8){
      if(x519Pi && x519time){
        hBeamCh->Fill(0.,1);
        isPion=true;
      }else{
        hBeamCh->Fill(2.,1);
        isProton=true;
      }
    }else if(run->beamChLogic==9){
      if(x474Pi && x474time && (x519Pi && x519time)){
        hBeamCh->Fill(0.,1);
        isPion=true;
      }else if (x474Pi && x474time && !(x519Pi && x519time)){
        hBeamCh->Fill(1.,1);
        isKaon=true;
      }else if (!(x474Pi && x474time) && !(x519Pi && x519time)){
        hBeamCh->Fill(2.,1);
        isProton=true;
      }
    }else if(run->beamChLogic==10){
      if((x474Pi && x474time) && (x519Pi && x519time)){
        hBeamCh->Fill(0.,1);
        isPion=true;
      }else if ((x474Pi && x474time) && !(x519Pi && x519time)){
        hBeamCh->Fill(1.,1);
        isKaon=true;
      }else if (!(x474Pi && x474time) && !(x519Pi && x519time)){
        hBeamCh->Fill(2.,1);
        isProton=true;
      }
    }else {
      if(debug)cout <<"No beam cherenkov available for this run\n";
    }
    //	if(isPion==true)printf("Is_pion!\n");
    //	if(isKaon==true)printf("Is_kaon!\n");
    //	if(isProton==true)printf("Is_proton!\n");

    hUpGEM->Fill(gx0,gy0);
    hDnGEM->Fill(gx1,gy1);
    hBeam->Fill(gxa,gya);
    hBeamTheta->Fill(gxtheta,gytheta);
    if(trigtime!=0)htrig->Fill(trigtime);
    if(x474time!=0)hx474->Fill(x474time-trigtime);
    if(x519time!=0)hx519->Fill(x519time-trigtime);
    if(x537time!=0)hx537->Fill(x537time-trigtime);
    if(goodSP[4]==true && spPhoton[4] > minPhoGas && spPhoton[4]<maxPhoGas) vspSigPhoGas[spPhoton[4]]->Fill(spRadius[4]);
    if(goodSP[9]==true && spPhoton[9] > minPhoAero && spPhoton[9]<maxPhoAero) vspSigPhoAero[spPhoton[9]]->Fill(spRadius[9]);
    if(goodSPN[4]==true && spnPhoton[4] > minPhoGas && spnPhoton[4]<maxPhoGas) { vspnSigPhoGas[spnPhoton[4]]->Fill(spnRadius[4]); printf(" +++ fill %6d %7.2f \n",spnPhoton[4],spnRadius[4]);}
    if(goodSPN[9]==true && spnPhoton[9] > minPhoAero && spnPhoton[9]<maxPhoAero) vspnSigPhoAero[spnPhoton[9]]->Fill(spnRadius[9]);
    if(goodCUT[4]==true && cutPhoton[4] > minPhoGas && cutPhoton[4]<maxPhoGas) vcutSigPhoGas[cutPhoton[4]]->Fill(cutRadius[4]);
    if(goodCUT[9]==true && cutPhoton[9] > minPhoAero && cutPhoton[9]<maxPhoAero) vcutSigPhoAero[cutPhoton[9]]->Fill(cutRadius[9]); 

    /* if(goodCUT[9] == true && cutRadius[9]<145){
    //printf(" Plotta eve %3d %4d \n",i,evt);
    for(int j = 0; j < nedge; j++){
    if(goodPhoton[j] && externalPhoton[j]) //printf("%3d %3d %3d %3d (%4d, %4d) %7.2f %7.2f %lf \n",evt,j,externalPhoton[j],pol[j],fiber[j],ch[j],nttw[j],dur[j],nr[j]);
    }
    //cin.get();
    }*/
    for(int j = 0; j < 10; j++){
      if(goodSP[j]==true && spPhoton[j] > 2){
        vspRadius[j]->Fill(spRadius[j]);
        //if(j==9 && abs(spRadius[j]-185)<2) cout <<"Watch run " <<run->runNum  <<" event: " <<evt <<endl;
        vspTime[j]->Fill(spTime[j]);
        vspPhoton[j]->Fill(spPhoton[j]);
      }
      if(goodSPN[j]==true && spnPhoton[j] > 0){
        vspnRadius[j]->Fill(spnRadius[j]);
        vspnTime[j]->Fill(spnTime[j]);
        vspnPhoton[j]->Fill(spnPhoton[j]);
      }
      if(goodCUT[j]==true && cutPhoton[j] > 2){
        vcutRadius[j]->Fill(cutRadius[j]);
        vcutTime[j]->Fill(cutTime[j]);
        vcutPhoton[j]->Fill(cutPhoton[j]);
      }
    }

    if(goodSPN[4]==true && spnPhoton[4]>2){
      if(isPion==true)vSpnPiRadius[1]->Fill(spnRadius[4]);
      if(isPion==true)vSpnPiRadius[0]->Fill(spnRadius[4]);
      if(isKaon==true)vSpnKRadius[1]->Fill(spnRadius[4]);
      if(isKaon==true)vSpnKRadius[0]->Fill(spnRadius[4]);
      if(isProton==true)vSpnPrRadius[1]->Fill(spnRadius[4]);
      if(isProton==true)vSpnPrRadius[0]->Fill(spnRadius[4]);
      hInterpolateGasRing->Fill(spnRadius[9]);
    }
    if(goodSPN[9]==true && spnPhoton[9]>2){
      if(isPion==true)vSpnPiRadius[2]->Fill(spnRadius[9]);
      if(isPion==true)vSpnPiRadius[0]->Fill(spnRadius[9]);
      if(isKaon==true)vSpnKRadius[2]->Fill(spnRadius[9]);
      if(isKaon==true)vSpnKRadius[0]->Fill(spnRadius[9]);
      if(isProton==true)vSpnPrRadius[2]->Fill(spnRadius[9]);
      if(isProton==true)vSpnPrRadius[0]->Fill(spnRadius[9]);
    }

    if(!goodSPN[4])hInterpolateNotGasRing->Fill(spnRadius[9]);
    if(isPion && !goodSPN[4])hInterpolateIsPionAndNotGasRing->Fill(spnRadius[9]);
    if(isKaon && goodSPN[4])hInterpolateIsKaonAndGasRing->Fill(spnRadius[9]);
    if(isProton && goodSPN[4])hInterpolateIsProtonAndGasRing->Fill(spnRadius[9]);


    //if(i<10)printf(" Plotta eve %3d %4d \n",i,evt);
    int nHit = 0;
    htMap->Reset();
    hTimeVsRadius->Reset();
    bool printThis=false;
    for(int j = 0; j < nedge; j++){
      //if(i<10)printf(" %3d %3d (%4d, %4d) %7.2f %7.2f %2d ",j,pol[j],fiber[j],ch[j],nttw[j],dur[j],pol[j]);
      if(dur[j]>CUT_MIN_DUR && trigSig[j]==false){
        if(pol[j]==0)hHitStart[0]->Fill(nttw[j]);
        if(pol[j]==1)hHitEnd[0]->Fill(nttw[j]);
        if(externalPhoton[j]==false){
          if(pol[j]==0)hHitStart[1]->Fill(nttw[j]);
          if(pol[j]==1)hHitEnd[1]->Fill(nttw[j]);
          //if(i<10)printf(" --> in\n");
        }else{
          if(pol[j]==0)hHitStart[2]->Fill(nttw[j]);
          if(pol[j]==1)hHitEnd[2]->Fill(nttw[j]);
          //if(i<10)printf(" --> out \n");
        }
      }else{
        //                            if(i<10)printf("bad \n");
      }
      if(goodHit[j]==true && trigSig[j]==false && pol[j]==0){
        hHitDur[0]->Fill(dur[j]);
        hHitCorrDur[0]->Fill(dur[j],nttw[j]);
        if(externalPhoton[j]==false){
          hHitDur[1]->Fill(dur[j]);
          hHitCorrDur[1]->Fill(dur[j],nttw[j]);
        }else{
          hHitDur[2]->Fill(dur[j]);
          hHitCorrDur[2]->Fill(dur[j],nttw[j]);
        }
        int ref=reference(pmt[j],ch[j]);
        if(goodPhoton[j]==true){
          hCalibVsCh[0]->Fill(ref,nttw[j]);
          hNotCalibVsCh[0]->Fill(ref,time[j]);
          if(externalPhoton[j]==false){
            hCalibVsCh[1]->Fill(ref,nttw[j]);
            hNotCalibVsCh[1]->Fill(ref,time[j]);
          }else{
            hCalibVsCh[2]->Fill(ref,nttw[j]);
            hNotCalibVsCh[2]->Fill(ref,time[j]);
          }
        }
      }
      if(run->sensor=="SIMULATION"){
        if(goodPhoton[j]==true){
          if(externalPhoton[j]==true)hSimPhotonLambdaAero->Fill(photonWavelength[j]);
          else hSimPhotonLambdaGas->Fill(photonWavelength[j]);
        }
      }
      if(goodPhoton[j]==true){
        hHitCoinc[0]->Fill(nttw[j]);
        if(externalPhoton[j]==false){
          hHitCoinc[1]->Fill(nttw[j]);
          vnRadius[1]->Fill(nr[j]);
          vRadiusMM[1]->Fill(rmm[j]);
          vnRadiusMM[1]->Fill(nrmm[j]);
          if(isPion==true){
            vPiRadius[1]->Fill(nr[j]);
            vPiMap[1]->Fill(x[j],y[j]);
          }
          if(isKaon==true){
            vKRadius[1]->Fill(nr[j]);
            vKMap[1]->Fill(x[j],y[j]);
          }
          if(isProton==true){
            vPrRadius[1]->Fill(nr[j]);
            vPrMap[1]->Fill(x[j],y[j]);
          }
        }else{
          hHitCoinc[2]->Fill(nttw[j]);
          vnRadius[2]->Fill(nr[j]);
          vRadiusMM[2]->Fill(rmm[j]);
          vnRadiusMM[2]->Fill(nrmm[j]);
          if(isPion==true){
            vPiRadius[2]->Fill(nr[j]);
            vPiMap[2]->Fill(x[j],y[j]);
          }
          if(isKaon==true){
            vKRadius[2]->Fill(nr[j]);
            vKMap[2]->Fill(x[j],y[j]);
          }
          if(isProton==true){
            vPrRadius[2]->Fill(nr[j]);
            vPrMap[2]->Fill(x[j],y[j]);
          }
        }
        hMap->Fill(x[j],y[j]);
        if(spnRadius[9]>tmpSel && goodSPN[9]==true && spnPhoton[9] > 2 && externalPhoton[j]){ 
          vRadiusRSD[0]->Fill(rsdRadius[9]);
          vTimeRSD[0]->Fill(rsdTime[9]);
          printThis=true;
        }
        if(spnRadius[9]<tmpSel && goodSPN[9]==true && spnPhoton[9] > 2 && externalPhoton[j]){ 
          printThis=true;
          hTimeVsRadius->Fill(nrmm[j],nt[j]);
          hTimeVsRadius->Fill(mRadTomm(spnRadius[9],run->firstPath,run->firstMirrorPosition),spnTime[9]);
          //cout <<i <<" " <<nt[j] <<" " <<nrmm[j] <<endl;
          htMap->SetTitle(Form("Entry %d",i));
          htMap->Fill(nx[j],ny[j]);
          vRadiusRSD[1]->Fill(rsdRadius[9]);
          vTimeRSD[1]->Fill(rsdTime[9]);
        }
        hnMap->Fill(nx[j],ny[j]);
        hRadius->Fill(r[j]);
        vnRadius[0]->Fill(nr[j]);
        vRadiusMM[0]->Fill(rmm[j]);
        vnRadiusMM[0]->Fill(nrmm[j]);
        if(externalPhoton[j]==false){
          vGasRadPMT[pmt[j]]->Fill(r[j]);
          vGasnRadPMT[pmt[j]]->Fill(nr[j]);
        }
        else{
          vAeroRadPMT[pmt[j]]->Fill(r[j]);
          vAeronRadPMT[pmt[j]]->Fill(nr[j]);
        }
        if(isPion==true){
          vPiRadius[0]->Fill(nr[j]);
          vPiMap[0]->Fill(x[j],y[j]);
        }
        if(isKaon==true){
          vKRadius[0]->Fill(nr[j]);
          vKMap[0]->Fill(x[j],y[j]);
        }
        if(isProton==true){
          vPrRadius[0]->Fill(nr[j]);
          vPrMap[0]->Fill(x[j],y[j]);
        }
      }
    }
    if(spnRadius[9]>tmpSel && goodSPN[9]==true && spnPhoton[9] > 2 && printThis==true){
      vRadiusRMS[0]->Fill(rmsRadius[9]);
      vTimeRMS[0]->Fill(rmsTime[9]);
    } 
    if(spnRadius[9]<tmpSel && goodSPN[9]==true && spnPhoton[9] > 2 && printThis==true){
      vRadiusRMS[1]->Fill(rmsRadius[9]);
      vTimeRMS[1]->Fill(rmsTime[9]);
      c->cd(1);
      htMap->Draw("colz");
      TEllipse *aRing = new TEllipse(0,0,mRadTomm(spnRadius[9],run->firstPath,run->firstMirrorPosition));
      aRing->SetLineColor(4);
      aRing->SetFillStyle(0);
      aRing->Draw("same");
      TEllipse *geoCutRing = new TEllipse(0,0,run->geoCut);
      geoCutRing->SetLineColor(2);
      geoCutRing->SetFillStyle(0);
      geoCutRing->Draw("same");

      c->cd(2);
      hTimeVsRadius->Draw("colz");;
      //c->Print(Form("lowRadEvent_%d.pdf",run->runNum));
      printThis=false;
    }
    if(goodSP[4]){}
    if(goodSP[9] && !goodSP[4]){
      cout <<Form("Event %d\nMaterial Pion_exRad Kaon_exRad Prot_exRad\nGas %lf %lf %lf, measured %lf\nAero %lf %lf %lf, measured %lf\n",i, eChPionGas, eChKaonGas, eChProtonGas,spRadius[4], eChPionAero, eChKaonAero, eChProtonAero,spRadius[9]);
      if(isnan(eChPionAero)==false) hDistanceExpectedPionAero->Fill(spRadius[9]-eChPionAero);
      if(isnan(eChKaonAero)==false) hDistanceExpectedKaonAero->Fill(spRadius[9]-eChKaonAero);
      if(isnan(eChProtonAero)==false) hDistanceExpectedProtonAero->Fill(spRadius[9]-eChProtonAero);
      if(isnan(eChPionGas)==false) hDistanceExpectedPionGas->Fill(spRadius[4]-eChPionGas);
      if(isnan(eChPionGas)==false) hDistanceExpectedKaonGas->Fill(spRadius[4]-eChKaonGas);
      if(isnan(eChPionGas)==false) hDistanceExpectedProtonGas->Fill(spRadius[4]-eChProtonGas);
    }
  }
  printEnd();
  fIn->Close();
  //c->Print(Form("lowRadEvent_%d.pdf]",run->runNum));
  c->Close();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------
void displayMonitor(THeader *run){
  //---------------------------------------------------

  hBaryCenter = new TH2D("hBaryCenter","Barycenter of photons;x [mm];y [mm]",200,-100,100,200,-100,100);
  hBaryVsCenterGas = new TH2D("hBaryVsCenterGas","Barycenter - Center Gas;#Deltax [mm];#Deltay [mm]",200,-10,10,200,-10,10);

  gErrorIgnoreLevel=kWarning;
  string out_pdf0 = Form("%s/output/plot/%s/displayMonitor.pdf[",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf  = Form("%s/output/plot/%s/displayMonitor.pdf",run->suite.c_str(),run->outputDir.c_str());
  string out_pdf1 = Form("%s/output/plot/%s/displayMonitor.pdf]",run->suite.c_str(),run->outputDir.c_str());
  string out_root = Form("%s/output/plot/%s/displayMonitor.root",run->suite.c_str(),run->outputDir.c_str());
  string out_pass = Form("%s/output/plot/%s/passport.pdf",run->suite.c_str(),run->outputDir.c_str());

  TList *save = new TList();

  TCanvas *c1 = new TCanvas("c1","Timing gas",1600,900);
  c1->Divide(4,2);
  c1->Draw();
  c1->Print(out_pdf0.c_str());
  c1->Print(Form("%s[",out_pass.c_str()));
  for (int irad=1; irad<3; irad++){

    double tmin = run->timeInMin;
    if(run->timeOuMin<tmin)tmin=run->timeOuMin;
    double tmax = run->timeInMax;
    if(run->timeOuMax>tmax)tmax=run->timeOuMax;

    double tmin_rad = run->timeInMin; 
    double tmax_rad = run->timeInMax; 
    if(irad==2){ 
      tmin_rad = run->timeOuMin;
      tmax_rad = run->timeOuMax;
    }
    double dmin = CUT_MIN_DUR;
    //printf("Time ref (%7.2f:%7.2f) gas (%7.2f:%7.2f)   aerogel (%7.2f:%7.2f)",tmin,tmax,
    //		run->timeInMin,run->timeInMax,run->timeOuMin,run->timeOuMax);

    c1->cd(1);
    //Hit start and end, showing coincidence window extremes
    THStack *hsHitTime = new THStack("hsHitTime","Hit start and End - Duration > 35;[ns]");
    hHitEnd[irad]->SetLineColor(2);
    hsHitTime->Add(hHitStart[irad]);
    hsHitTime->Add(hHitEnd[irad]);
    hsHitTime->Draw("nostack");
    hsHitTime->GetXaxis()->SetRangeUser(run->timeOuMin-0.1*run->timeOuMin,run->timeInMax+110);
    //hsHitTime->GetXaxis()->SetRangeUser(300,600);
    hsHitTime->Draw("nostack");
    TLine *l1 = new TLine(tmin_rad,0,tmin_rad,hHitStart[0]->GetBinContent(hHitStart[0]->GetMaximumBin()));
    l1->SetLineColor(3);
    l1->Draw("same");
    TLine *l2 = new TLine(tmax_rad,0,tmax_rad,hHitStart[0]->GetBinContent(hHitStart[0]->GetMaximumBin()));
    l2->SetLineColor(3);
    l2->Draw("same");
    c1->Update();

    c1->cd(5);
    //Zoom on coincidence peak, fitted. Showing coincidence window extremes.
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    hHitCoinc[irad]->Draw();
    TF1 *fcoinc = new TF1("fcoinc","gaus(0)",200,400);
    fcoinc->SetParameters(20000,355,2);
    hHitCoinc[irad]->Fit(fcoinc,"Q","",run->timeOuMin,run->timeInMax);
    l1->SetY2(hHitCoinc[irad]->GetBinContent(hHitCoinc[irad]->GetMaximumBin())/2);
    l2->SetY2(hHitCoinc[irad]->GetBinContent(hHitCoinc[irad]->GetMaximumBin())/2);
    l1->Draw("same");
    l2->Draw("same");
    c1->Update();

    c1->cd(2);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);
    hx474->SetLineColor(2);
    hx474->GetXaxis()->SetRangeUser(-50,50);
    hx519->SetLineColor(3);
    hx519->GetXaxis()->SetRangeUser(-50,50);
    hx537->SetLineColor(4);
    hx537->GetXaxis()->SetRangeUser(-50,50);
    THStack *hsBeamCherenkov= new THStack("hsBeamCherenkov","Beam cherenkov;time [ns];counts [#]");
    save->Add(hsBeamCherenkov);
    gPad->SetLogy();	
    hsBeamCherenkov->Add(hx474);
    hsBeamCherenkov->Add(hx519);
    hsBeamCherenkov->Add(hx537);
    hsBeamCherenkov->Draw("nostack");
    hsBeamCherenkov->GetXaxis()->SetRangeUser(-75,75);
    hsBeamCherenkov->Draw("nostack");
    gPad->BuildLegend(.3,.7,.7,.9);
	
    c1->cd(6);
    htrig->GetXaxis()->SetRangeUser(htrig->GetBinCenter(htrig->FindFirstBinAbove(1))-10,htrig->GetBinCenter(htrig->FindLastBinAbove(1))+10);
    htrig->Draw();

    c1->cd(3);
    //Trigger time plot. If I can, start and end.
    hHitDur[irad]->GetXaxis()->SetRangeUser(0,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(1))+5);
    hHitDur[irad]->Draw();
    TLine *l3 = new TLine(dmin,0,dmin,hHitDur[irad]->GetBinContent(hHitDur[0]->GetMaximumBin()));
    l3->SetLineColor(3);
    l3->Draw("same");

    c1->cd(4);
    gPad->SetLogz();
    hHitCorrDur[irad]->GetXaxis()->SetRangeUser(0,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(3))+5);
    //hHitCorrDur[irad]->GetYaxis()->SetRangeUser(run->timeOuMin-20,run->timeInMax+20);
    hHitCorrDur[irad]->GetYaxis()->SetRangeUser(run->timeOuMin-20,run->timeOuMax+25);
    hHitCorrDur[irad]->Draw("colz");
    TLine *l4 = new TLine(0,tmin_rad,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(3))+5,tmin_rad);
    l4->SetLineColor(3);
    l4->Draw("same");
    TLine *l5 = new TLine(0,tmax_rad,hHitDur[irad]->GetBinCenter(hHitDur[irad]->FindLastBinAbove(3))+5,tmax_rad);
    l5->SetLineColor(3);
    l5->Draw("same");
    TLine *l6 = new TLine(dmin,run->timeOuMin-20,dmin,run->timeInMax+20);
    l6->SetLineColor(3);
    l6->Draw("same");

    c1->cd(7);
    //Time not calibrated VS Channel
    hNotCalibVsCh[irad]->GetYaxis()->SetRangeUser(run->timeOuMin,run->timeInMax+5);
    //hNotCalibVsCh[irad]->GetXaxis()->SetRangeUser(0,300);
    hNotCalibVsCh[irad]->Draw("colz");
    TProfile *tp0 = hNotCalibVsCh[irad]->ProfileX();
    tp0->GetYaxis()->SetRangeUser(run->timeOuMin,run->timeInMax+5);
    tp0->SetLineColor(kWhite);
    //tp0->Draw("same");

    c1->cd(8);
    hCalibVsCh[irad]->GetYaxis()->SetRangeUser(run->timeOuMin,run->timeInMax+5);
    //hCalibVsCh[irad]->GetXaxis()->SetRangeUser(0,300);
    hCalibVsCh[irad]->Draw("colz");
    TProfile *tp1 = hCalibVsCh[irad]->ProfileX();
    tp1->GetYaxis()->SetRangeUser(run->timeOuMin,run->timeInMax+5);
    tp1->SetLineColor(kWhite);
    //tp1->Draw("same");

    c1->Update();
    c1->Print(out_pdf.c_str());

    TCanvas *cc1 = new TCanvas("cc1","cc1",1600,900);
    gStyle->SetOptStat(0);
    cc1->Draw();
    cc1->Divide(3,2);
    cc1->cd(1);
    hsHitTime->Draw("nostack");
    l1->Draw("same");
    l2->Draw("same");
    cc1->cd(2);
    hHitDur[irad]->Draw();
    l3->Draw("same");
    cc1->cd(3); 
    gPad->SetLogz();
    hHitCorrDur[irad]->Draw("colz");
    l4->Draw("same");
    l5->Draw("same");
    l6->Draw("same");
    cc1->cd(4);
    hHitCoinc[irad]->Draw();
    l1->Draw("same");
    l2->Draw("same");
    cc1->cd(5);
    hNotCalibVsCh[irad]->Draw("colz");
    cc1->cd(6);
    hCalibVsCh[irad]->Draw("colz");
    cc1->Update();
    cc1->Print(Form("%s",out_pass.c_str()));
    cc1->Close();

    TCanvas *cT  = new TCanvas(Form("cT%d",irad),"c0",1600,900);
    cT->Draw();
    hNotCalibVsCh[irad]->Draw("colz");
    cT->Update();
    cT->Print(Form("hNotCalivVsCh_%d.root",irad));
    hCalibVsCh[irad]->Draw("colz");
    cT->Update();
    cT->Print(Form("hCalivVsCh_%d.root",irad));
    cT->Close();
  }
  for (int j=1; j<=1024; j++){
    TProfile *tp_gas = hNotCalibVsCh[1]->ProfileX();
    TProfile *tp_aer = hNotCalibVsCh[2]->ProfileX();
    double off_gas = tp_gas->GetBinContent(j);
    double off_aer = tp_aer->GetBinContent(j);
    double err_gas = tp_gas->GetBinError(j);
    double err_aer = tp_aer->GetBinError(j);
    double zero = 0.0;
    /*printf(" jjj %3d %7.2f %7.2f %7.2f %7.2f ",j,off_gas,off_aer,err_gas,err_aer);
      if(off_gas!=0 && off_aer==0)printf(" prof %3d %7.2f \n",j,360-off_gas);
      if(off_gas==0 && off_aer!=0)printf(" prof %3d %7.2f \n",j,355-off_aer);
      if(off_gas!=0 && off_aer!=0){
      if(err_gas<=err_aer) printf(" prof %3d %7.2f \n",j,360-off_gas);
      if(err_gas>err_aer) printf(" prof %3d %7.2f \n",j,355-off_aer);
      }
      if(off_gas==0 && off_aer==0)printf(" prof %3d %7.2f \n",j,zero);*/
  }



  TCanvas *c2 = new TCanvas("c2","Radii canvas",1600,900);
  c2->Divide(3,3);
  c2->Draw();
  c2->cd(1);
  hBeamCh->SetMinimum(0);
  hBeamCh->Draw();
  //Cherenkov combination for pion 1
  c2->cd(2);
  //Cherenkov combination for pion 2
  c2->cd(3);
  //Cherenkov combination for pion+kaon

  c2->cd(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  //Fitted GAS radius for SP.
  vspRadius[4]->SetTitle("sp radius - gas");
  TH1D *cp00 = (TH1D*)vspRadius[4]->Clone("hspRadius_fitIn");
  TF1 *fspRadius_4=new TF1();
  applyFit(cp00,fspRadius_4,"fspRadius_4",false);
  cp00->GetXaxis()->SetRangeUser(10,60);
  cp00->Draw(); 
  c2->cd(5);
  //Fitted GAS radius for SPN.
  vspnRadius[4]->SetTitle("Corrected sp radius - gas");
  TH1D *cp01 = (TH1D*)vspnRadius[4]->Clone("hspnRadius_fitIn");
  TF1 *fspnRadius_4=new TF1();
  //applyFit(cp01,fspnRadius_4,"fspnRadius_4",false);
  fspnRadius_4=applyFit(cp01,fspnRadius_4,"fspnRadius_4",false,0);
  cp01->Draw(); 
  c2->cd(6);
  //Fitted GAS radius for CUT.
  vcutRadius[4]->SetTitle("RMS cut radius - gas");
  TH1D *cp02 = (TH1D*)vcutRadius[4]->Clone("hcutRadius_fitIn");
  TF1 *fcutRadius_4=new TF1();
  applyFit(cp02,fcutRadius_4,"fcutRadius_4",false);
  cp02->Draw(); 
  c2->cd(7);
  //Fitted AERO radius for SP.
  vspRadius[9]->SetTitle("Sp radius - aerogel");
  TH1D *cp10 = (TH1D*)vspRadius[9]->Clone("hspRadius_fitOut");
  TF1 *fspRadius_9 = new TF1();
  applyFit(cp10,fspRadius_9,"fspRadius_4",true);
  cp10->Draw();
  c2->cd(8);
  //Fitted AERO radius for SPN.
  vspnRadius[9]->SetTitle("Corrected sp radius - aerogel");
  TH1D *cp11 = (TH1D*)vspnRadius[9]->Clone("hspnRadius_fitOut");
  TF1 *fspnRadius_9 = new TF1();
  //applyFit(cp11,fspnRadius_9,"fspnRadius_9",true);
  fspnRadius_9=applyFit(cp11,fspnRadius_9,"fspnRadius_9",true,0);
  cp11->Draw();
  c2->cd(9);
  //Fitted AERO radius for CUT.
  vcutRadius[9]->SetTitle("RMS cut radius - aerogel");
  TH1D *cp12 = (TH1D*)vcutRadius[9]->Clone("hcutRadius_fitOut");
  TF1 *fcutRadius_9 = new TF1();
  applyFit(cp12,fcutRadius_9,"fcutRadius_4",true);
  cp12->Draw();
  c2->Update();
  c2->Print(out_pdf.c_str());
  c2->Clear();


  TCanvas *c3 = new TCanvas("c3","Time and photon for radius",1600,900);
  c3->Divide(3,4);
  c3->cd(1);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);
  vspTime[4]->SetTitle("Mean time of the event - gas - Single Particle");
  vspTime[4]->Draw();
  c3->cd(4);
  vspPhoton[4]->SetTitle("# photons for particle - gas - Single Particle");
  vspPhoton[4]->Draw();
  c3->cd(7);
  vspTime[9]->SetTitle("Mean time of the event - aerogel - Single Particle");
  vspTime[9]->Draw();
  c3->cd(10);
  vspPhoton[9]->SetTitle("# photons for particle - aerogel - Single Particle");
  vspPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspPhoton[9]->Draw();
  c3->cd(2);
  vspnTime[4]->SetTitle("Mean time of the event - gas - Corrected Single Particle");
  vspnTime[4]->Draw();
  c3->cd(5);
  vspnPhoton[4]->SetTitle("# photons for particle - gas - Corrected Single Particle");
  vspnPhoton[4]->Draw();
  c3->cd(8);
  vspnTime[9]->SetTitle("Mean time of the event - aerogel - Corrected Single Particle");
  vspnTime[9]->Draw();
  c3->cd(11);
  vspnPhoton[9]->SetTitle("# photons for particle - aerogel - Corrected Single Particle");
  vspnPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspnPhoton[9]->Draw();
  c3->cd(3);
  vcutTime[4]->SetTitle("Mean time of the event - gas - Cut on RMS");
  vcutTime[4]->Draw();
  c3->cd(6);
  vcutPhoton[4]->SetTitle("# photons for particle - gas - Cut on RMS");
  vcutPhoton[4]->Draw();
  c3->cd(9);
  vcutTime[9]->SetTitle("Mean time of the event - aerogel - Cut on RMS");
  vcutTime[9]->Draw();
  c3->cd(12);
  vcutPhoton[9]->SetTitle("# photons for particle - aerogel - Cut on RMS");
  vcutPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vcutPhoton[9]->Draw();
  c3->Update();
  c3->Print(out_pdf.c_str());

  TCanvas *cc2 = new TCanvas("cc2","cc2",1600,900);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  cc2->Divide(3,2);
  cc2->cd(1);
  cp01->SetTitle("Single particle Cherenkov angle - gas;#theta_{Ch} [mRad]; counts"); 
  cp01->Draw(); 
  cc2->cd(4);
  cp11->SetTitle("Single particle Cherenkov angle - aerogel;#theta_{Ch} [mRad]; counts"); 
  cp11->Draw(); 
  cc2->cd(2);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1);
  vspnTime[4]->SetTitle("Mean time of the event - gas");
  vspnTime[4]->Draw();
  cc2->cd(3);
  vspnPhoton[4]->SetTitle("# photons for particle - gas");
  vspnPhoton[4]->GetXaxis()->SetRangeUser(0,60);
  vspnPhoton[4]->Draw();
  cc2->cd(5);
  vspnTime[9]->SetTitle("Mean time of the event - aerogel");
  vspnTime[9]->Draw();
  cc2->cd(6);
  vspnPhoton[9]->SetTitle("# photons for particle - aerogel");
  vspnPhoton[9]->GetXaxis()->SetRangeUser(0,25);
  vspnPhoton[9]->Draw();
  cc2->Update();

  TCanvas *c4 = new TCanvas("c4","GEM canvas",1600,900);
  c4->Divide(2,2);
  c4->Draw();
  c4->cd(1);
  gPad->SetGrid();
  hUpGEM->GetXaxis()->SetRangeUser(-20,20);
  hUpGEM->GetYaxis()->SetRangeUser(-20,20);
  hUpGEM->Draw("colz");
  c4->cd(2);
  gPad->SetGrid();
  hDnGEM->GetXaxis()->SetRangeUser(-20,20);
  hDnGEM->GetYaxis()->SetRangeUser(-20,20);
  hDnGEM->Draw("colz");
  c4->cd(3);
  gPad->SetGrid();
  hBeam->GetXaxis()->SetRangeUser(-20,20);
  hBeam->GetYaxis()->SetRangeUser(-20,20);
  hBeam->Draw("colz");
  c4->cd(4);
  gPad->SetGrid();
  hBeamTheta->GetXaxis()->SetRangeUser(-0.005,0.005);
  hBeamTheta->GetYaxis()->SetRangeUser(-0.005,0.005);
  hBeamTheta->Draw("colz");
  c4->Update();
  c4->Print(out_pdf.c_str());

  gStyle->SetOptFit(1);
  TCanvas *cPhotBin = new TCanvas("cPhotBin","Photon binning",1600,900);
  cPhotBin->Draw();
  //cPhotBin->Print("plot_photon_binning.pdf[");
  //SIGMA vs Photon Number
  for(int i = 1; i < maxPhoGas; i++){
    if(vspSigPhoGas[i-1]->GetEntries()<5) continue; //Qui c'era i e non i-1 ed un == 0;
    TF1 *fspGas = getFun(vspSigPhoGas[i-1],false); //Qui c'era i-1;
    hspSigVsPhoGas->SetBinContent(i,fspGas->GetParameter(2));
    hspSigVsPhoGas->SetBinError(i,fspGas->GetParError(2));
    vspSigPhoGas[i-1]->Draw();
    //cPhotBin->Print("plot_photon_binning.pdf");
  }
  for(int i = 1; i < maxPhoGas; i++){
    if(vspnSigPhoGas[i-1]->GetEntries()<5) continue;
    TF1 *fspnGas = getFun(vspnSigPhoGas[i-1],false);
    //if(fspnGas->GetParError(2) > 0.1) continue;
    //hspnSigVsPhoGas->SetBinContent(i,fspnGas->GetParameter(2));
    //hspnSigVsPhoGas->SetBinError(i,fspnGas->GetParError(2));
    hspnSigVsPhoGas->SetBinContent(i,vspnSigPhoGas[i-1]->GetRMS());
    hspnSigVsPhoGas->SetBinError(i,0.02);
    vspnSigPhoGas[i-1]->Draw();
    save->Add(vspnSigPhoGas[i-1]);
    printf(" VSP %6d %4d %7.2f %7.2f [%7.2f %7.2f] --> [%7.2f %7.2f] \n",maxPhoGas,i,vspnSigPhoGas[i-1]->GetEntries(),vspnSigPhoGas[i-1]->GetBinCenter(vspnSigPhoGas[i-1]->GetMaximumBin()),
									vspnSigPhoGas[i-1]->GetMean(), vspnSigPhoGas[i-1]->GetRMS(), fspnGas->GetParameter(2),fspnGas->GetParError(2));
    //cPhotBin->Print("plot_photon_binning.pdf");
  }
  for(int i = 1; i < maxPhoGas; i++){
    if(vcutSigPhoGas[i-1]->GetEntries()<5) continue;
    TF1 *fcutGas = getFun(vcutSigPhoGas[i-1],false);
    hcutSigVsPhoGas->SetBinContent(i,fcutGas->GetParameter(2));
    hcutSigVsPhoGas->SetBinError(i,fcutGas->GetParError(2));
    vcutSigPhoGas[i-1]->Draw();
    //cPhotBin->Print("plot_photon_binning.pdf");
  }
  for(int i = 1; i < maxPhoAero; i++){
    if(vspSigPhoAero[i-1]->GetEntries()<5) continue;
    TF1 *fspAero = getFun(vspSigPhoAero[i-1],true);
    hspSigVsPhoAero->SetBinContent(i,fspAero->GetParameter(2));
    hspSigVsPhoAero->SetBinError(i,fspAero->GetParError(2));
    vspSigPhoAero[i-1]->Draw();
    //cPhotBin->Print("plot_photon_binning.pdf");
  }
  for(int i = 1; i < maxPhoAero; i++){
    if(vspnSigPhoAero[i-1]->GetEntries()<5) continue;
    TF1 *fspnAero = getFun(vspnSigPhoAero[i-1],true);
    if(fspnAero->GetParError(2) > 0.4) continue;
    hspnSigVsPhoAero->SetBinContent(i,fspnAero->GetParameter(2));
    hspnSigVsPhoAero->SetBinError(i,fspnAero->GetParError(2));
    vspnSigPhoAero[i-1]->Draw();
    save->Add(vspnSigPhoAero[i-1]);
    //cPhotBin->Print("plot_photon_binning.pdf");
    //if(i==6)cPhotBin->Print("plot.root");
  }
  for(int i = 1; i < maxPhoAero; i++){
    if(vcutSigPhoAero[i-1]->GetEntries()<5) continue;
    TF1 *fcutAero = getFun(vcutSigPhoAero[i-1],true);
    hcutSigVsPhoAero->SetBinContent(i,fcutAero->GetParameter(2));
    hcutSigVsPhoAero->SetBinError(i,fcutAero->GetParError(2));
    vcutSigPhoAero[i-1]->Draw();
    //cPhotBin->Print("plot_photon_binning.pdf");
  }
  //cPhotBin->Print("plot_photon_binning.pdf]");


  bool simFlag = false;
  if(run->sensor=="SIMULATION")simFlag=true;
  TCanvas *c5 = new TCanvas("c5","Single photon",1600,900);
  c5->Divide(4,2);
  c5->cd(1);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  //vnRadius[1]->Scale(countEvent);
  vnRadius[1]->Rebin(4);
  vnRadius[1]->Fit("gaus","Q","",20,50);
  vnRadius[1]->Draw();
  c5->cd(5);
  //vnRadius[2]->Scale(1/(double)countEvent);
  //if(run->runNum==186)vnRadius[2]->Scale(0.15/0.18);
  vnRadius[2]->Rebin(4); 
  vnRadius[2]->Fit("gaus","Q","",130,240);
  //vnRadius[2]->SetMaximum(0.08);
  vnRadius[2]->Draw(); 
  c5->cd(2);
  fitSigma(hspSigVsPhoGas,false,simFlag);
  //hspSigVsPhoGas->GetXaxis()->SetRangeUser(0,maxPhoGas+5);
  hspSigVsPhoGas->GetYaxis()->SetRangeUser(0,1.5);
  hspSigVsPhoGas->Draw("E");
  c5->cd(6);
  fitSigma(hspSigVsPhoAero,true,simFlag);
  //hspSVsPhoAero->GetXaxis()->SetRangeUser(0,maxPhoAero+3);
  hspSigVsPhoAero->GetYaxis()->SetRangeUser(0,4);
  hspSigVsPhoAero->Draw("E");
  c5->cd(3);
  fitSigma(hspnSigVsPhoGas,false,simFlag);
  //hspnSigVsPhoGas->GetXaxis()->SetRangeUser(0,maxPhoGas+5);
  hspnSigVsPhoGas->GetYaxis()->SetRangeUser(0,1.5);
  hspnSigVsPhoGas->Draw("E");
  //hspnSigVsPhoGas->Print("gas.root");
  save->Add(hspnSigVsPhoGas);
  c5->cd(7);
  fitSigma(hspnSigVsPhoAero,true,simFlag);
  //hspnSigVsPhoAero->GetXaxis()->SetRangeUser(0,maxPhoAero+3);
  hspnSigVsPhoAero->GetYaxis()->SetRangeUser(0,4);
  hspnSigVsPhoAero->Draw("E");
  //hspnSigVsPhoAero->Print("aero.root");
  save->Add(hspnSigVsPhoAero);
  c5->cd(4);
  fitSigma(hcutSigVsPhoGas,false,simFlag);
  //hcutSigVsPhoGas->GetXaxis()->SetRangeUser(0,maxPhoGas+5);
  hcutSigVsPhoGas->GetYaxis()->SetRangeUser(0,1.5);
  hcutSigVsPhoGas->Draw("E");
  c5->cd(8);
  fitSigma(hcutSigVsPhoAero,true,simFlag);
  //hcutSigVsPhoAero->GetXaxis()->SetRangeUser(0,maxPhoAero+3);
  hcutSigVsPhoAero->GetYaxis()->SetRangeUser(0,4);
  hcutSigVsPhoAero->Draw("E");
  c5->Update();
  c5->Print(out_pdf.c_str());

  TCanvas *cc5 = new TCanvas("cc5","cc5",1600,900);
  cc5->Draw();
  cc5->Divide(2,2);
  cc5->cd(1);
  vnRadius[1]->SetTitle("Single photon Cherenkov angle - gas;#theta_{Ch} [mRad];counts");
  vnRadius[1]->Draw();
  cc5->cd(2);
  hspnSigVsPhoGas->SetTitle("#sigma_{#theta} vs photons for particle - gas;photons/particle;#sigma_{#theta} [mRad]");
  //hspnSigVsPhoGas->SetMaximum(.5);
  hspnSigVsPhoGas->Draw("E");
  cc5->cd(3);
  vnRadius[2]->SetTitle("Single photon Cherenkov angle - aerogel;#theta_{Ch} [mRad];counts");
  vnRadius[2]->Draw();
  cc5->cd(4);
  hspnSigVsPhoAero->SetTitle("#sigma_{#theta} vs photons for particle - aerogel;photons/particle;#sigma_{#theta} [mRad]");
  hspnSigVsPhoAero->Draw("E");

  TCanvas *cc = new TCanvas("cc","cc",1600,900);
  cc->Draw();
  hcutSigVsPhoGas->Draw();
  //cc->Print("cc.root");

  TCanvas *c6 = new TCanvas("c6","Selected particle canvas",1600,900);
  c6->Divide(3,3);
  c6->cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  vPiRadius[0]->Draw();
  c6->cd(2);
  vKRadius[0]->Draw();
  c6->cd(3);
  vPrRadius[0]->Draw();
  c6->cd(4);
  vPiRadius[1]->Fit("gaus","Q","",20,50);
  vPiRadius[1]->Draw();
  c6->cd(5);
  vKRadius[1]->Fit("gaus","Q","",20,50);
  vKRadius[1]->Draw();
  c6->cd(6);
  vPrRadius[1]->Fit("gaus","Q","",20,50);
  vPrRadius[1]->Draw();
  c6->cd(7);
  //vPiRadius[2]->Scale(1./countEvent);
  vPiRadius[2]->Fit("gaus","Q","",130,250);
  vPiRadius[2]->Draw();
  c6->cd(8);
  vKRadius[2]->Fit("gaus","Q","",130,250);
  vKRadius[2]->Draw();
  c6->cd(9);
  vPrRadius[2]->Fit("gaus","Q","",130,250);
  vPrRadius[2]->Draw();
  c6->Update();
  if(run->sensor!="SIMULATION")
    c6->Print(out_pdf.c_str());

  TCanvas *c8 = new TCanvas("c8","Selected particle THStack",1600,900);
  c8->Draw();
  c8->Divide(2);
  c8->cd(1);
  gPad->SetLogy();
  THStack *hsCut1 = new THStack("hsCut1","Single particle radius - Gas;Radius [mRad];Counts [#]");
  vSpnKRadius[1]->SetLineColor(2);
  vSpnPrRadius[1]->SetLineColor(3);
  hsCut1->Add(vSpnPiRadius[1]);
  hsCut1->Add(vSpnKRadius[1]);
  hsCut1->Add(vSpnPrRadius[1]);
  hsCut1->Draw("nostack hist");
  hsCut1->GetXaxis()->SetRangeUser(25,50);
  hsCut1->Draw("nostack hist");
  TLegend *lhsCut1 = new TLegend(0.7,0.7,0.9,0.9);
  lhsCut1->AddEntry(vSpnPiRadius[1],"Pion","lp");
  lhsCut1->AddEntry(vSpnKRadius[1],"Kaon","lp");
  lhsCut1->AddEntry(vSpnPrRadius[1],"Proton","lp");
  lhsCut1->Draw("same");
  c8->cd(2);
  gPad->SetLogy();
  THStack *hsCut2 = new THStack("hsCut2","Single particle radius - Aerogel;Radius [mRad];Counts [#]");
  save->Add(hsCut2);
  save->Add(vSpnPiRadius[2]);
  save->Add(vSpnKRadius[2]);
  save->Add(vSpnPrRadius[2]);
  vSpnKRadius[2]->SetLineColor(2);
  vSpnPrRadius[2]->SetLineColor(3);
  hsCut2->Add(vSpnPiRadius[2]);
  hsCut2->Add(vSpnKRadius[2]);
  hsCut2->Add(vSpnPrRadius[2]);
  hsCut2->Draw("nostack hist");
  TLegend *lhsCut2 = new TLegend(0.1,0.7,0.1,0.9);
  lhsCut2->AddEntry(vSpnPiRadius[2],"Pion","lp");
  lhsCut2->AddEntry(vSpnKRadius[2],"Kaon","lp");
  lhsCut2->AddEntry(vSpnPrRadius[2],"Proton","lp");
  lhsCut2->Draw("same");
  c8->Update();
  if(run->sensor!="SIMULATION")
    c8->Print(out_pdf.c_str());




  TCanvas *cSim = new TCanvas("cSim","Simulation canvas",1600,900);
  if(run->sensor=="SIMULATION"){
    cSim->Draw();
    cSim->Divide(2);
    cSim->cd(1);
    hSimPhotonLambdaGas->Draw();
    save->Add(hSimPhotonLambdaGas);
    cSim->cd(2);
    hSimPhotonLambdaAero->Draw();
    save->Add(hSimPhotonLambdaAero);
    cSim->Update();
    cSim->Print(out_pdf.c_str());
  }


  TCanvas *c7 = new TCanvas("c7","Selected particle single radius",1600,900);
  c7->Draw();
  c7->Divide(3,3);
  c7->cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  vSpnPiRadius[0]->Draw();
  c7->cd(2);
  vSpnKRadius[0]->Draw();
  c7->cd(3);
  vSpnPrRadius[0]->Draw();
  c7->cd(4);
  TF1 *fSpnPiRadius1 = new TF1();
  applyFit(vSpnPiRadius[1],fSpnPiRadius1,"fSpnPiRadius1",false);
  vSpnPiRadius[1]->Draw();
  c7->cd(5);
  vSpnKRadius[1]->SetLineColor(4);;
  TF1 *fSpnKRadius1 = new TF1();
  applyFit(vSpnKRadius[1],fSpnKRadius1,"fSpnKRadius1",false);
  vSpnKRadius[1]->Draw();
  c7->cd(6);
  vSpnPrRadius[1]->SetLineColor(4);;
  TF1 *fSpnPrRadius1 = new TF1();
  applyFit(vSpnPrRadius[1],fSpnKRadius1,"fSpnKRadius1",false);
  vSpnPrRadius[1]->Draw();
  c7->cd(7);
  TF1 *fSpnPiRadius2 = new TF1();
  applyFit(vSpnPiRadius[2],fSpnPiRadius2,"fSpnPiRadius2",true);
  vSpnPiRadius[2]->Draw();
  c7->cd(8);
  vSpnKRadius[2]->SetLineColor(4);;
  TF1 *fSpnKRadius2 = new TF1();
  applyFit(vSpnKRadius[2],fSpnKRadius2,"fSpnKRadius2",true);
  vSpnKRadius[2]->Draw();
  c7->cd(9);
  vSpnPrRadius[2]->SetLineColor(4);;
  TF1 *fSpnPrRadius2 = new TF1();
  applyFit(vSpnPrRadius[2],fSpnPrRadius2,"fSpnPrRadius2",true);
  vSpnPrRadius[2]->Draw();
  c7->Update();
  if(run->sensor!="SIMULATION")
    c7->Print(out_pdf.c_str());


  TCanvas *c9 = new TCanvas("c9","Selected particle maps",1600,900);
  c9->Divide(3,3);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1); 
  c9->cd(1);
  vPiMap[0]->SetTitle("Pions - Full map");
  vPiMap[0]->Draw("colz");
  c9->cd(2);
  vKMap[0]->SetTitle("Kaons - Full map");
  vKMap[0]->Draw("colz");
  c9->cd(3);
  vPrMap[0]->SetTitle("Protons - Full map");
  vPrMap[0]->Draw("colz");

  c9->cd(4);
  vPiMap[1]->SetTitle("Pions - Gas");
  vPiMap[1]->Draw("colz");
  c9->cd(5);
  vKMap[1]->SetTitle("Kaons - Gas");
  vKMap[1]->Draw("colz");
  c9->cd(6);
  vPrMap[1]->SetTitle("Protons - Gas");
  vPrMap[1]->Draw("colz");

  c9->cd(7);
  vPiMap[2]->SetTitle("Pions - Aerogel");
  vPiMap[2]->Draw("colz");
  c9->cd(8);
  vKMap[2]->SetTitle("Kaons - Aerogel");
  vKMap[2]->Draw("colz");
  c9->cd(9);
  vPrMap[2]->SetTitle("Protons - Aerogel");
  vPrMap[2]->Draw("colz");
  c9->Update();
  if(run->sensor!="SIMULATION")
    c9->Print(out_pdf.c_str());


  TCanvas *c10 = new TCanvas("c10","Map canvas",1600,900);
  c10->Divide(2);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0); 
  c10->Draw();
  c10->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  hMap->SetMaximum(hnMap->GetMaximum());
  hMap->Draw("colz");
  /*    if(radius[j] > run->geoCut){
        inPath=run->firstPath;
        zMir=run->firstMirrorPosition;
        }else{
        inPath=run->secondPath;
        zMir=run->secondMirrorPosition;
        }*/
  TEllipse *gRing = new TEllipse(0,0,mRadTomm(fspnRadius_4->GetParameter(1),run->secondPath,run->secondMirrorPosition));
  gRing->SetLineColor(4);
  gRing->SetFillStyle(0);
  gRing->Draw("same");
  TEllipse *aRing = new TEllipse(0,0,mRadTomm(fspnRadius_9->GetParameter(1),run->firstPath,run->firstMirrorPosition));
  aRing->SetLineColor(4);
  aRing->SetFillStyle(0);
  aRing->Draw("same");
  TEllipse *geoCutRing = new TEllipse(0,0,run->geoCut);
  geoCutRing->SetLineColor(2);
  geoCutRing->SetFillStyle(0);
  geoCutRing->Draw("same");
  c10->cd(2);
  gPad->SetLogz();
  hnMap->Draw("colz");
  aRing->Draw("same");
  gRing->Draw("same");
  c10->Update();
  c10->Print(out_pdf.c_str());


  TCanvas *c11 = new TCanvas("c11","PMT radius",1600,900);
  c11->Divide(4,2);
  c11->Draw();
  c11->cd(1);
  THStack *hsGasPMT0 = new THStack("hsGasPMT0","PMT north - GAS");
  vGasRadPMT[0]->SetTitle("not corrected north");
  vGasRadPMT[0]->SetLineColor(3);
  vGasnRadPMT[0]->SetTitle("corrected north");
  vGasnRadPMT[0]->SetLineColor(4);
  vGasRadPMT[2]->SetTitle("not corrected south");
  vGasRadPMT[2]->SetLineColor(1);
  vGasnRadPMT[2]->SetTitle("corrected south");
  vGasnRadPMT[2]->SetLineColor(2);
  hsGasPMT0->Add(vGasRadPMT[0]);
  hsGasPMT0->Add(vGasnRadPMT[0]);
  hsGasPMT0->Add(vGasRadPMT[2]);
  hsGasPMT0->Add(vGasnRadPMT[2]);
  hsGasPMT0->Draw("nostack");
  gPad->BuildLegend(.15,.7,.3,.9);
  c11->cd(2);
  THStack *hsGasPMT2 = new THStack("hsGasPMT2","PMT south - GAS");
  vGasRadPMT[2]->SetLineColor(3);
  hsGasPMT2->Add(vGasRadPMT[2]);
  hsGasPMT2->Add(vGasnRadPMT[2]);
  hsGasPMT2->Draw("nostack");
  c11->cd(3);
  THStack *hsGasPMT1 = new THStack("hsGasPMT1","PMT est - GAS");
  vGasRadPMT[1]->SetLineColor(3);
  vGasnRadPMT[3]->SetLineColor(2);
  hsGasPMT1->Add(vGasRadPMT[1]);
  hsGasPMT1->Add(vGasnRadPMT[1]);
  hsGasPMT1->Add(vGasnRadPMT[3]);
  hsGasPMT1->Draw("nostack");
  c11->cd(4);
  THStack *hsGasPMT3 = new THStack("hsGasPMT3","PMT west - GAS");
  vGasRadPMT[3]->SetLineColor(3);
  hsGasPMT3->Add(vGasRadPMT[3]);
  hsGasPMT3->Add(vGasnRadPMT[3]);
  hsGasPMT3->Draw("nostack");
  c11->cd(5);
  THStack *hsAeroPMT0 = new THStack("hsAeroPMT0","PMT north - AERO");
  vAeroRadPMT[0]->SetLineColor(3);
  vAeronRadPMT[2]->SetLineColor(2);
  hsAeroPMT0->Add(vAeroRadPMT[0]);
  hsAeroPMT0->Add(vAeronRadPMT[0]);
  hsAeroPMT0->Add(vAeronRadPMT[2]);
  hsAeroPMT0->Draw("nostack");
  c11->cd(6);
  THStack *hsAeroPMT2 = new THStack("hsAeroPMT2","PMT south - AERO");
  vAeroRadPMT[2]->SetLineColor(3);
  hsAeroPMT2->Add(vAeroRadPMT[2]);
  hsAeroPMT2->Add(vAeronRadPMT[2]);
  hsAeroPMT2->Draw("nostack");
  c11->cd(7);
  THStack *hsAeroPMT1 = new THStack("hsAeroPMT1","PMT est - AERO");
  vAeroRadPMT[1]->SetLineColor(3);
  vAeronRadPMT[3]->SetLineColor(2);
  hsAeroPMT1->Add(vAeroRadPMT[1]);
  hsAeroPMT1->Add(vAeronRadPMT[1]);
  hsAeroPMT1->Add(vAeronRadPMT[3]);
  hsAeroPMT1->Draw("nostack");
  c11->cd(8);
  THStack *hsAeroPMT3 = new THStack("hsAeroPMT3","PMT west - AERO");
  vAeroRadPMT[3]->SetLineColor(3);
  hsAeroPMT3->Add(vAeroRadPMT[3]);
  hsAeroPMT3->Add(vAeronRadPMT[3]);
  hsAeroPMT3->Draw("nostack");
  c11->Update();
  c11->Print(out_pdf.c_str());


  TCanvas *c12 = new TCanvas("c12","Radius corrections",1600,900);
  c12->Divide(3,2);
  c12->cd(1);
  hGasdXVsAngle->Draw("colz");
  c12->cd(2);
  hGasdYVsAngle->Draw("colz");
  c12->cd(3);
  //hNewCenterGasX->Draw();
  //c12->cd(4);
  //hNewCenterGasY->Draw();
  hNewCenterGas->Draw("colz");

  c12->cd(4);
  hAerodXVsAngle->Draw("colz");
  c12->cd(5);
  hAerodYVsAngle->Draw("colz");
  c12->cd(6);
  //hNewCenterAeroX->Draw();;
  //c12->cd(8);
  //hNewCenterAeroY->Draw();;
  hNewCenterAero->Draw("colz");;
  c12->Update();
  c12->Print(out_pdf.c_str());

  TCanvas *c13 = new TCanvas("c13","Radius in mm",1600,900);
  c13->Divide(3,2);
  c13->cd(1);
  vRadiusMM[0]->SetTitle("Radius");
  vRadiusMM[0]->Draw();
  c13->cd(2);
  vRadiusMM[1]->GetXaxis()->SetRangeUser(5,60);
  vRadiusMM[1]->SetTitle("Radius - Gas");
  vRadiusMM[1]->Draw();
  c13->cd(3);
  vRadiusMM[2]->GetXaxis()->SetRangeUser(50,95);
  vRadiusMM[2]->SetTitle("Radius - Aerogel");
  vRadiusMM[2]->Draw();
  c13->cd(4);
  vnRadiusMM[0]->SetTitle("Corrected radius");
  vnRadiusMM[0]->Draw();
  c13->Update();
  c13->cd(5);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  vnRadiusMM[1]->GetXaxis()->SetRangeUser(25,60);
  vnRadiusMM[1]->SetTitle("Corrected radius - Gas");
  vnRadiusMM[1]->Fit("gaus","","",vnRadiusMM[1]->GetBinCenter(vnRadiusMM[1]->GetMaximumBin())-12,vnRadius[1]->GetBinCenter(vnRadiusMM[1]->GetMaximumBin())+12);
  vnRadiusMM[1]->Draw();
  c13->cd(6);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  vnRadiusMM[2]->GetXaxis()->SetRangeUser(50,95);
  vnRadiusMM[2]->SetTitle("Corrected radius - Aerogel");
  vnRadiusMM[2]->Fit("gaus","","",68,76);
  vnRadiusMM[2]->Draw();
  c13->Update();
  c13->Print(out_pdf.c_str());

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas *cc13 = new TCanvas("cc13","cc13",1600,900);
  cc13->Divide(2);
  cc13->cd(1);
  vnRadiusMM[1]->GetXaxis()->SetRangeUser(10,60);
  vnRadiusMM[1]->SetTitle("Single photon radius mm - Gas");
  vnRadiusMM[1]->Draw();
  cc13->cd(2);
  vnRadiusMM[2] ->SetTitle("Single photon radius mm - Aerogel");
  vnRadiusMM[2]->Draw();
  cc13->Update(); 

  c10->Print(Form("%s",out_pass.c_str()));
  c4->Print(Form("%s",out_pass.c_str()));
  cc13->Print(Form("%s",out_pass.c_str()));
  gStyle->SetOptStat(0);
  cc2->Print(Form("%s",out_pass.c_str()));
  cc5->Print(Form("%s",out_pass.c_str()));
  //c12->Print(Form("%s",out_pass.c_str()));
  if(run->sensor!="SIMULATION") c7->Print(Form("%s",out_pass.c_str()));
  vSpnKRadius[1]->SetLineColor(2);
  vSpnPrRadius[1]->SetLineColor(3);
  vSpnKRadius[2]->SetLineColor(2);
  vSpnPrRadius[2]->SetLineColor(3);
  if(run->sensor!="SIMULATION") c8->Print(Form("%s",out_pass.c_str()));
  if(run->sensor=="SIMULATION") cSim->Print(Form("%s",out_pass.c_str()));
  cc13->Print(Form("%s]",out_pass.c_str()));


  TCanvas *c14 = new TCanvas("c14","RMS and RSD",1600,900);
  c14->Divide(4,2);
  c14->cd(1);
  //gPad->SetLogy();
  vRadiusRMS[0]->SetTitle(Form("Radius RMS > %.0lf mRad",tmpSel));
  vRadiusRMS[0]->Draw();
  c14->cd(2);
  //gPad->SetLogy();
  vRadiusRSD[0]->SetTitle(Form("Radius RSD > %.0lf mRad",tmpSel));
  vRadiusRSD[0]->Draw();
  c14->cd(3);
  //gPad->SetLogy();
  vTimeRMS[0]->SetTitle(Form("Time RMS > %.0lf mRad",tmpSel));
  vTimeRMS[0]->Draw();
  c14->cd(4);
  //gPad->SetLogy();
  vTimeRSD[0]->SetTitle(Form("Time RSD > %.0lf mRad",tmpSel));
  vTimeRSD[0]->Draw();
  c14->cd(5);
  //gPad->SetLogy();
  vRadiusRMS[1]->SetTitle(Form("Radius RMS < %.0lf mRad",tmpSel));
  vRadiusRMS[1]->Draw();
  c14->cd(6);
  //gPad->SetLogy();
  vRadiusRSD[1]->SetTitle(Form("Radius RSD < %.0lf mRad",tmpSel));
  vRadiusRSD[1]->Draw();
  c14->cd(7);
  //gPad->SetLogy();
  vTimeRMS[1]->SetTitle(Form("Time RMS < %.0lf mRad",tmpSel));
  vTimeRMS[1]->Draw();
  c14->cd(8);
  //gPad->SetLogy();
  vTimeRSD[1]->SetTitle(Form("Time RSD < %.0lf mRad",tmpSel));
  vTimeRSD[1]->Draw();

  c14->Update();
  c14->Print(out_pdf.c_str());

  TCanvas *c15 = new TCanvas("c15","Measured vs expected Cherenkov angle",1600,900);
  c15->Draw();
  gStyle->SetOptStat(1);
  c15->Divide(3,2);
  c15->cd(1);
  hDistanceExpectedPionAero->Draw();
  c15->cd(2);
  hDistanceExpectedKaonAero->Draw();
  c15->cd(3);
  hDistanceExpectedProtonAero->Draw();
  c15->cd(4);
  hDistanceExpectedPionGas->Draw();
  c15->cd(5);
  hDistanceExpectedKaonGas->Draw();
  c15->cd(6);
  hDistanceExpectedProtonGas->Draw();
  c15->Print(out_pdf.c_str());

  int rbFactor=2;
  TCanvas *c16 = new TCanvas("c16","Canvas to interpolate",1600,900);
  gStyle->SetOptStat(0);
  c16->Draw();
  c16->Divide(2,2);
  c16->cd(1);
  THStack *hsInterpolatePion = new THStack("hsInterpolatePion","Pion analysis;r [mRad];counts");
  auto cloneSpnPiRadius=(TH1D*) vSpnPiRadius[2]->Clone("cloneSpnPiRadius");
  cloneSpnPiRadius->SetTitle("#pi tagged");
  cloneSpnPiRadius->SetLineColor(4);
  cloneSpnPiRadius->Rebin(rbFactor);
  hInterpolateGasRing->SetLineColor(1);
  hInterpolateGasRing->Rebin(rbFactor);
  hsInterpolatePion->Add(cloneSpnPiRadius);
  hsInterpolatePion->Add(hInterpolateGasRing);
  hsInterpolatePion->Draw("nostack hist");
  gPad->BuildLegend(.7,.7,.9,.9);
  c16->cd(2);
  THStack *hsInterpolateKaon = new THStack("hsInterpolateKaon","Kaon analysis;r [mRad]; counts");
  auto cloneSpnKRadius=(TH1D*) vSpnKRadius[2]->Clone("cloneSpnKRadius");
  cloneSpnKRadius->SetTitle("k tagged");
  cloneSpnKRadius->SetLineColor(2); 
  cloneSpnKRadius->Rebin(rbFactor);
  hInterpolateNotGasRing->SetLineColor(1);
  hInterpolateNotGasRing->Rebin(rbFactor);
  hsInterpolateKaon->Add(cloneSpnKRadius);
  hsInterpolateKaon->Add(hInterpolateNotGasRing);
  hsInterpolateKaon->Draw("nostack hist");
  gPad->BuildLegend(.7,.7,.9,.9);
  c16->cd(3);
  THStack *hsInterpolateProton = new THStack("hsInterpolateProton","Proton analysis;r [mRad]; counts");
  auto cloneSpnPrRadius=(TH1D*) vSpnPrRadius[2]->Clone("cloneSpnPrRadius");
  cloneSpnPrRadius->SetTitle("p tagged");
  cloneSpnPrRadius->SetLineColor(3); 
  cloneSpnPrRadius->Rebin(rbFactor);
  hsInterpolateProton->Add(cloneSpnPrRadius);
  hsInterpolateProton->Add(hInterpolateNotGasRing);
  hsInterpolateProton->Draw("nostack hist"); 
  gPad->BuildLegend(.7,.7,.9,.9);
  c16->cd(4);
  THStack *hsNegation = new THStack("hsNegation","Inefficiency;#theta [mRad]; counts");
  hInterpolateIsProtonAndGasRing->Rebin(rbFactor);
  hInterpolateIsKaonAndGasRing->Rebin(rbFactor);
  hInterpolateIsPionAndNotGasRing->Rebin(rbFactor);
  hInterpolateIsProtonAndGasRing->SetLineColor(3);
  hInterpolateIsKaonAndGasRing->SetLineColor(2);
  hInterpolateIsPionAndNotGasRing->SetLineColor(4);
  hsNegation->Add(hInterpolateIsProtonAndGasRing);
  hsNegation->Add(hInterpolateIsKaonAndGasRing);
  hsNegation->Add(hInterpolateIsPionAndNotGasRing);
  hsNegation->Draw("nostack hist");
  gPad->BuildLegend(.7,.7,.9,.9);
  c16->Update();  
  c16->Print(out_pdf.c_str());
  c16->Close();

  // Baricenter
  TCanvas *c17 = new TCanvas("c17","Barycenter canvas",1600,900);
  c17->Divide(2);
  c17->cd(1);
  gStyle->SetOptStat(1);
  hBaryCenter->Draw("colz");
  c17->cd(2);
  hBaryVsCenterGas->Draw("colz");
  c17->Update();
  c17->Print(out_pdf.c_str());

  save->Add(vnRadius[1]);
  save->Add(vnRadius[2]);
  save->Add(vspnRadius[4]);
  save->Add(vspnRadius[9]);
  save->Add(vspnPhoton[4]);
  save->Add(vspnPhoton[9]);
  save->Add(vRadiusMM[0]);
  save->Add(vRadiusMM[1]);
  save->Add(vRadiusMM[2]);
  save->Add(hG);
  save->Add(hA);
  save->Add(hBaryCenter);
  save->Add(hBaryVsCenterGas);

  TCanvas *c20 = new TCanvas("c20","Canvas to interpolate",1600,900);
  gStyle->SetOptStat(1);
  for(int i = 1; i < maxPhoGas; i++){
    vspnSigPhoGas[i-1]->Draw();
    c20->Print(out_pdf.c_str());
  }
  

  c1->Print(out_pdf1.c_str());
  if(saveRootFile==true){
    TFile *fOut = new TFile(out_root.c_str(),"RECREATE");
    save->Write();
    fOut->Close();
    //TFile *fOutTmp = new TFile("tmp.root","RECREATE");
    //hspnSigVsPhoGas->Write();
    //hspnSigVsPhoAero->Write();
    //cout <<"Nbins: "<<vnRadius[1]->GetNbinsX() <<endl;
    //fOutTmp->Close();
  }
  cout <<Form("You can open the monitoring plots typing: evince %s\n",out_pass.c_str());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void inizializePlot(THeader *run){
  if(run->sensor=="SIMULATION" && APPLY_QUANTUM_EFFICIENCY==false){
    minPhoGas=40;
    maxPhoGas=150;
    minPhoAero=5;
    maxPhoAero=40;
  }
  //Create /output directory
  string out_dir = Form("%s/output/plot/%s",run->suite.c_str(),run->outputDir.c_str());
  cout <<out_dir.c_str() <<endl;
  if(gSystem->AccessPathName(out_dir.c_str())){
    if(std::system(Form("mkdir -p %s",out_dir.c_str()))){
      cout <<"[ERROR] Output directory not created\n";
      exit(EXIT_FAILURE);
    }
  }

  hTime = new TH1D("hTime all","Time distribution for all hits",1000,0,1000);

  for (int i=0; i<3; i++){
    TH1D *hStart = new TH1D(Form("hHitStart_%d",i),Form("Hit start - %d",i),2000,0,2000);
    hHitStart.push_back(hStart);
    TH1D *hEnd = new TH1D(Form("hHitEnd_%d",i),Form("Hit end - %d",i),2000,0,2000);
    hHitEnd.push_back(hEnd);
    TH1D *hDur = new TH1D(Form("hHitDur_%d",i),Form("Hit duration - %d; Dur [ns]",i),101,-.50,100.5);
    hHitDur.push_back(hDur);
    TH2D *hCorrDur = new TH2D(Form("hHitCorrDur_%d",i),Form("Hit correlation; Dur [ns]; Start [ns] - %d",i),71,-.5,70.5,2000,0,2000); 
    hHitCorrDur.push_back(hCorrDur);
    TH1D *hCoinc = new TH1D(Form("hHitCoinc_%d",i),Form("Coincidence peak - %d",i),5*(int)(run->timeInMax-run->timeOuMin)+10,run->timeOuMin-5,run->timeInMax+5);
    hHitCoinc.push_back(hCoinc);
    TH2D *hCal = new TH2D(Form("hCalibVsCh_%d",i),Form("#Delta time calibrated - trigger time Vs Channel - %d",i),1026,-0.5,1025.5,run->timeInMax-run->timeOuMin+60,run->timeOuMin-30,run->timeInMax+30);
    TH2D *hNoCal = new TH2D(Form("hNotCalibVsCh_%d",i),Form("#Delta time not calibrated - trigger time Vs Channel - %d",i),1026,-0.5,1025.5,run->timeInMax-run->timeOuMin+60,run->timeOuMin-30,run->timeInMax+30);
    hCalibVsCh.push_back(hCal);
    hNotCalibVsCh.push_back(hNoCal);
  }

  hHitCorr = new TH2D("hHitCorr","Hit correlation; Start [ns]; End[ns]",1000,0,1000,1000,0,1000); 
  hHitReco = new TH1D("hHitReco","Number of reconstructed hit for event",151,-0.5,150.5);
  hEdge = new TH1D("hEdge","Number of edge for event",151,-0.5,150.5);

  for(int i = 0; i < 4; i++){
    TH1D *hGasRadPMT = new TH1D(Form("hGasRadPMT%d",i),"hGasRadPMT;r [mRad]",160,20,60);
    vGasRadPMT.push_back(hGasRadPMT);
    TH1D *hAeroRadPMT = new TH1D(Form("hAeroRadPMT%d",i),"hAeroRadPMT,r [mRad]",440,140,250);
    vAeroRadPMT.push_back(hAeroRadPMT);
    TH1D *hGasnRadPMT = new TH1D(Form("hGasnRadPMT%d",i),"hGasnRadPMT;r [mRad]",160,20,60);
    vGasnRadPMT.push_back(hGasnRadPMT);
    TH1D *hAeronRadPMT = new TH1D(Form("hAeronRadPMT%d",i),"hAeronRadPMT,r [mRad]",440,140,250);
    vAeronRadPMT.push_back(hAeronRadPMT);
  }


  if(run->sensor=="MPPC") hMap = new TH2D("hMap","Hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT"|| run->sensor=="SIMULATION") hMap = new TH2D("hMap",Form("Hit position MAPMT - run %d;x [mm];y [mm]",run->runNum),sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMap = new TH2D("hMap","Hit position Other;x [mm];y [mm]",180,-90,90,180,-90,90);

  htMap = (TH2D*) hMap->Clone("htmpMap");

  for(int i = 0; i < 3; i++){
    TH2D *hPiMap = (TH2D*) hMap->Clone(Form("hPiMap%d",i));
    TH2D *hKMap = (TH2D*) hMap->Clone(Form("hKMap%d",i));
    TH2D *hPrMap = (TH2D*) hMap->Clone(Form("hPrMap%d",i));
    vPiMap.push_back(hPiMap);
    vKMap.push_back(hKMap);
    vPrMap.push_back(hPrMap);
  }

  hTimeVsRadius = new TH2D("hTimeVsRadius","Time Vs radius;r[mm];t[ns]",50,40,90,100,420,440);

  for(int i = 0; i < 3; i++){
    TH1D *hRadiusRMS = new TH1D(Form("hRadiusRMS_%d",i),"Radius RMS;RMS [mRad]",100,0,20);
    vRadiusRMS.push_back(hRadiusRMS);
    TH1D *hTimeRMS = new TH1D(Form("hTimeRMS_%d",i),"Time RMS; t [ns]",100,0,2);
    vTimeRMS.push_back(hTimeRMS);
    TH1D *hRadiusRSD= new TH1D(Form("hRadiusRSD_%d",i),"Radius Residui; RSD[mRad]",100,0,20);
    vRadiusRSD.push_back(hRadiusRSD);
    TH1D *hTimeRSD= new TH1D(Form("hTimeRSD_%d",i),"Time Residui; t[ns]",100,0,2.5);
    vTimeRSD.push_back(hTimeRSD);
  }

  hDistanceExpectedPionAero = new TH1D("hDistanceExpectedPionAero","#Delta #theta_{Ch} = measured - expected for #pi - Aerogel;#Delta #theta_{Ch} [mRad];counts",250,-50,50);
  hDistanceExpectedKaonAero = new TH1D("hDistanceExpectedKaonAero","#Delta #theta_{Ch} = measured - expected for K - Aerogel;#Delta #theta_{Ch} [mRad];counts",250,-50,50);
  hDistanceExpectedProtonAero = new TH1D("hDistanceExpectedProtonAero","#Delta #theta_{Ch} = measured - expected for p - Aerogel;#Delta #theta_{Ch} [mRad];counts",400,-50,55);
  hDistanceExpectedPionGas = new TH1D("hDistanceExpectedPionGas","#Delta #theta_{Ch} = measured - expected for #pi - Gas;#Delta #theta_{Ch} [mRad];counts",250,-50,50);
  hDistanceExpectedKaonGas = new TH1D("hDistanceExpectedKaonGas","#Delta #theta_{Ch} = measured - expected for K - Gas;#Delta #theta_{Ch} [mRad];counts",250,-50,50);
  hDistanceExpectedProtonGas = new TH1D("hDistanceExpectedProtonGas","#Delta #theta_{Ch} = measured - expected for p - Gas;#Delta #theta_{Ch} [mRad];counts",250,-50,50);

  if(run->sensor=="MPPC") hMapNC = new TH2D("hMapNC","All hit position MPPC;x [mm];y [mm]",sizeof(xBinMPPC)/sizeof(*xBinMPPC)-1,xBinMPPC,sizeof(yBinMPPC)/sizeof(*yBinMPPC)-1,yBinMPPC);
  else if(run->sensor=="MAPMT") hMapNC = new TH2D("hMapNC","All hit position MAPMT;x [mm];y [mm]",sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);
  else hMapNC = new TH2D("hMapNC","All hit position Other;x [mm];y [mm]",180,-90,90,180,-90,90);

  //hnMap = new TH2D("hnMap","Corrected positions of hit;x [mm];y [mm]",180,-90,90,180,-90,90);
  hnMap = new TH2D("hnMap","Corrected position MAPMT;x [mm];y [mm]",sizeof(xBinMAPMT)/sizeof(*xBinMAPMT)-1,xBinMAPMT,sizeof(yBinMAPMT)/sizeof(*yBinMAPMT)-1,yBinMAPMT);

  //STD GEM HISTO
  hUpGEM = new TH2D("hUpGEM","Upstream GEM; x_0[mm];y_0[mm]",501,-100.2,100.2,501,-100.2,100.2);
  hDnGEM = new TH2D("hDnGEM","Downstream GEM; x_0[mm];y_0[mm]",501,-100,100,501,-100,100);
  hBeam = new TH2D("hBeam","Beam profile at aerogel; x_0[mm];y_0[mm]",501,-100,100,501,-100,100);
  hBeamTheta = new TH2D("hBeamTheta","Beam divergence; x_0[Rad];y_0[Rad]",501,-.02,.02,501,-.02,.02);


  hRadius = new TH1D("hRadius","Single photon radius - before corrections;r [mRad]",500,0,250);
  //	hnRadius = new TH1D("hnRadius","Single photon radius - after corrections;r [mRad]",400,0,200);
  TH1D *hnRadius0 = new TH1D("hnRadius0","Single photon radius - after corrections;r [mRad]",500,0,250);
  TH1D *hnRadius1 = new TH1D("hnRadius1","Single photon radius - Gas;r [mRad]",1200,0,60);
  TH1D *hnRadius2 = new TH1D("hnRadius2","Single photon radius - Aerogel;r [mRad]",1040,120,250);
  vnRadius.push_back(hnRadius0);
  vnRadius.push_back(hnRadius1);
  vnRadius.push_back(hnRadius2);


  for(int i = 0; i < 3; i++){
    TH1D *hRadiusMM = new TH1D(Form("hRadiusMM_%d",i),"Radius in mm;r[mm]",400,0,100);
    TH1D *hnRadiusMM = new TH1D(Form("hnRadiusMM_%d",i),"Corrected radius in mm;r[mm]",400,0,100);
    vRadiusMM.push_back(hRadiusMM);
    vnRadiusMM.push_back(hnRadiusMM);
  }

  TH1D *hPiRadius0 = new TH1D("hPiRadius0","Pi - Single photon radius - after corrections;r [mRad]",500,0,250);
  TH1D *hPiRadius1 = new TH1D("hPiRadius1","Pi - Single photon radius - Gas;r [mRad]",400,0,100);
  TH1D *hPiRadius2 = new TH1D("hPiRadius2","Pi - Single photon radius - Aerogel;r [mRad]",520,120,250);
  vPiRadius.push_back(hPiRadius0);
  vPiRadius.push_back(hPiRadius1);
  vPiRadius.push_back(hPiRadius2);
  TH1D *hKRadius0 = new TH1D("hKRadius0","K - Single photon radius - after corrections;r [mRad]",500,0,250);
  TH1D *hKRadius1 = new TH1D("hKRadius1","K - Single photon radius - Gas;r [mRad]",400,0,100);
  TH1D *hKRadius2 = new TH1D("hKRadius2","K - Single photon radius - Aerogel;r [mRad]",520,120,250);
  vKRadius.push_back(hKRadius0);
  vKRadius.push_back(hKRadius1);
  vKRadius.push_back(hKRadius2);
  TH1D *hPrRadius0 = new TH1D("hPrRadius0","Pr - Single photon radius - after corrections;r [mRad]",500,0,250);
  TH1D *hPrRadius1 = new TH1D("hPrRadius1","Pr - Single photon radius - Gas;r [mRad]",400,0,100);
  TH1D *hPrRadius2 = new TH1D("hPrRadius2","Pr - Single photon radius - Aerogel;r [mRad]",520,120,250);
  vPrRadius.push_back(hPrRadius0);
  vPrRadius.push_back(hPrRadius1);
  vPrRadius.push_back(hPrRadius2);


  TH1D *hSpnPiRadius0 = new TH1D("hSpnPiRadius0","Pi - Single particle radius - after corrections;r [mRad]",500,0,250);
  TH1D *hSpnPiRadius1 = new TH1D("hSpnPiRadius1","Pi - Single particle radius - Gas;r [mRad]",400,0,100);
  TH1D *hSpnPiRadius2 = new TH1D("hSpnPiRadius2","Pi - Single particle radius - Aerogel;r [mRad]",520,120,250);
  vSpnPiRadius.push_back(hSpnPiRadius0);
  vSpnPiRadius.push_back(hSpnPiRadius1);
  vSpnPiRadius.push_back(hSpnPiRadius2);
  TH1D *hSpnKRadius0 = new TH1D("hSpnKRadius0","K - Single particle radius - after corrections;r [mRad]",500,0,250);
  TH1D *hSpnKRadius1 = new TH1D("hSpnKRadius1","K - Single particle radius - Gas;r [mRad]",400,0,100);
  TH1D *hSpnKRadius2 = new TH1D("hSpnKRadius2","K - Single particle radius - Aerogel;r [mRad]",520,120,250);
  vSpnKRadius.push_back(hSpnKRadius0);
  vSpnKRadius.push_back(hSpnKRadius1);
  vSpnKRadius.push_back(hSpnKRadius2);
  TH1D *hSpnPrRadius0 = new TH1D("hSpnPrRadius0","Pr - Single particle radius - after corrections;r [mRad]",500,0,250);
  TH1D *hSpnPrRadius1 = new TH1D("hSpnPrRadius1","Pr - Single particle radius - Gas;r [mRad]",400,0,100);
  TH1D *hSpnPrRadius2 = new TH1D("hSpnPrRadius2","Pr - Single particle radius - Aerogel;r [mRad]",520,120,250);
  vSpnPrRadius.push_back(hSpnPrRadius0);
  vSpnPrRadius.push_back(hSpnPrRadius1);
  vSpnPrRadius.push_back(hSpnPrRadius2);

  for(int i = 0; i < 10; i++){
    TH1D *hspRadius = new TH1D(Form("hspRadius_%d",i),Form("Single particle radius - %d - before corrections;radius [mRad]",i),2000,0,250);
    vspRadius.push_back(hspRadius);
    TH1D *hspTime = new TH1D(Form("hspTime_%d",i),Form("Single particle radius - %d - before corrections;time [ns]",i),5*(int)(run->timeInMax-run->timeOuMin),run->timeOuMin,run->timeInMax);
    vspTime.push_back(hspTime);
    TH1D *hspPhoton = new TH1D(Form("hspPhoton_%d",i),Form("Single particle radius - %d - before corrections;photon [#]",i),150,0,150);
    vspPhoton.push_back(hspPhoton);

    TH1D *hspnRadius = new TH1D(Form("hspnRadius_%d",i),Form("Single particle radius - %d - after corrections;radius [mRad]",i),2000,0,250);
    vspnRadius.push_back(hspnRadius);
    TH1D *hspnTime = new TH1D(Form("hspnTime_%d",i),Form("Single particle radius - %d - after corrections;time [ns]",i),5*(int)(run->timeInMax-run->timeOuMin),run->timeOuMin,run->timeInMax);
    vspnTime.push_back(hspnTime);
    TH1D *hspnPhoton = new TH1D(Form("hspnPhoton_%d",i),Form("Single particle radius - %d - after corrections;photon [#]",i),150,0,150);
    vspnPhoton.push_back(hspnPhoton);

    TH1D *hcutRadius = new TH1D(Form("hcutRadius_%d",i),Form("Single particle radius - %d - after rms cuts;radius [mRad]",i),2000,0,250);
    vcutRadius.push_back(hcutRadius);
    TH1D *hcutTime = new TH1D(Form("hcutTime_%d",i),Form("Single particle radius - %d - after rms cuts;time [ns]",i),5*(int)(run->timeInMax-run->timeOuMin),run->timeOuMin,run->timeInMax);
    vcutTime.push_back(hcutTime);
    TH1D *hcutPhoton = new TH1D(Form("hcutPhoton_%d",i),Form("Single particle radius - %d - after rms cuts;photon [#]",i),150,0,150);
    vcutPhoton.push_back(hcutPhoton);
  }

  for(int i = 0; i < 2; i++){
    TH1D *hrsdRadius = new TH1D(Form("hrsdRadius_%d",i),Form("Radius residui - %d;rsd_r [mRad];counts [#]",i),120,0,12);
    vrsdRadius.push_back(hrsdRadius);
    TH1D *hrsdTime = new TH1D(Form("hrsdTime_%d",i),Form("Time residui - %d;rsd_t [ns];counts [#]",i),120,0,12);
    vrsdTime.push_back(hrsdTime);
  }

  for(int i = 0; i < maxPhoGas; i++){
    int nBin=2000;
    TH1D *hspSigPhoGas= new TH1D(Form("hspSigPhoGas_%02d",i),Form("hspSigPhoGas_%02d",i),nBin,0,250);
    vspSigPhoGas.push_back(hspSigPhoGas);
    TH1D *hspnSigPhoGas= new TH1D(Form("hspnSigPhoGas_%02d",i),Form("hspnSigPhoGas_%02d",i),nBin,0,250);
    vspnSigPhoGas.push_back(hspnSigPhoGas);
    TH1D *hcutSigPhoGas= new TH1D(Form("hcutSigPhoGas_%02d",i),Form("hcutSigPhoGas_%02d",i),nBin,0,250);
    vcutSigPhoGas.push_back(hcutSigPhoGas);
  }
  for(int i = 0; i < maxPhoAero; i++){
    int nBin=1000;
    TH1D *hspSigPhoAero= new TH1D(Form("hspSigPhoAero_%02d",i),Form("hspSigPhoAero_%02d",i),nBin,0,250);
    vspSigPhoAero.push_back(hspSigPhoAero);
    TH1D *hspnSigPhoAero= new TH1D(Form("hspnSigPhoAero_%02d",i),Form("hspnSigPhoAero_%02d",i),nBin,0,250);
    vspnSigPhoAero.push_back(hspnSigPhoAero);
    TH1D *hcutSigPhoAero= new TH1D(Form("hcutSigPhoAero_%02d",i),Form("hcutSigPhoAero_%02d",i),nBin,0,250);
    vcutSigPhoAero.push_back(hcutSigPhoAero);
  }
  hspSigVsPhoGas = new TH1D("hspSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hspnSigVsPhoGas= new TH1D("hsnpSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hcutSigVsPhoGas = new TH1D("hcutSigVsPhoGas","#sigma_{r} vs photon per particle - Gas;Photon/particle [#];#sigma_{r} [mRad]",maxPhoGas+1,-0.5,maxPhoGas+0.5);
  hspSigVsPhoAero = new TH1D("hspSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);
  hspnSigVsPhoAero= new TH1D("hsnpSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);
  hcutSigVsPhoAero = new TH1D("hcutSigVsPhoAero","#sigma_{r} vs photon per particle - Aerogel;Photon/particle [#];#sigma_{r} [mRad]",maxPhoAero+1,-0.5,maxPhoAero+0.5);



  //hCoinc = new TH1D("hCoinc","Coincidence peak",5*(int)(run->timeOuMax-run->timeInMin)+10,run->timeInMin-5,run->timeOuMax+5);

  ////// MONITORING PLOT DEFINITION
  //hCalibVsTime = new TH2D("hCalibVsTime","Calibrated time vs acquired time",(int)(run->timeOuMax-run->timeInMin)+11,run->timeInMin-5.5,run->timeOuMax+5.5,(int)(run->timeOuMax-run->timeInMin)+11,run->timeInMin-5.5,run->timeOuMax+5.5);
  //hTimeWalkVsTime = new TH2D("hTimeWalkVsTime","TimeWalk corrected time vs acquired time",(int)(run->timeOuMax-run->timeInMin)+11,run->timeInMin-5.5,run->timeOuMax+5.5,5*(int)(run->timeOuMax-run->timeInMin)+50,run->timeInMin-5,run->timeOuMax+5);
  hCalibVsTime = new TH2D("hCalibVsTime","Calibrated time vs acquired time",51,349.5,400.1,(int)(run->timeInMax-run->timeOuMin)+11,run->timeOuMin-5.5,run->timeInMax+5.5);
  hTimeVsCh = new TH2D("hTimeVsCh","#Delta time not calibrated - trigger time Vs Channel",1026,-0.5,1025.5,70,330,400);
  hTimeWalkVsTime = new TH2D("hTimeWalkVsTime","TimeWalk corrected time vs acquired time",51,349.5,401.5,5*(int)(run->timeInMax-run->timeOuMin)+50,run->timeOuMin-5,run->timeInMax+5);
  //hDur = new TH1D("hDur","Hit duration",101,-.50,100.5);
  hStartVsDur = new TH2D("hStartVsDur","DuratioVsStart; Start [ns]; Duration[ns]",(int)(run->timeInMax-run->timeOuMin)+11,run->timeOuMin-5.5,run->timeInMax+5.5,101,-.5,100.5); 
  hSingPhotRad = new TH1D("hSingPhotRad","Single photon radius; Radius [mRad]",400,0,200);

  for(int i = 0; i < 32; i++){
    vector<TH1D*> tmpNot;
    vector<TH1D*> tmpCal;
    for(int j = 0; j < 1025; j++){
      TH1D* hNot = new TH1D(Form("hNot_%d_%d",i+4,j),"Time not calibrated",1001,-.5,1000.5);
      tmpNot.push_back(hNot);
      TH1D* hCal = new TH1D(Form("hCal_%d_%d",i+4,j),"Time calibrated",1001,-.5,1000.5);
      tmpCal.push_back(hCal);
    }
    vNotTime.push_back(tmpNot);
    vCalTime.push_back(tmpCal);
  } 
  hx474=new TH1D("hx474","Cherenkov x474",3000,-1000,2000);
  hx519=new TH1D("hx519","Cherenkov x519",3000,-1000,2000);
  hx537=new TH1D("hx537","Cherenkov x537",3000,-1000,2000);
  htrig=new TH1D("htrig","Trigger;time [ns]",3000,-1000,2000);


  hPiRad = new TH1D("hPiRad","#pi radius - gas; radius [mRad]",200,0,100);
  hKRad = new TH1D("hKRad","K radius - gas; radius [mRad]",200,0,100);
  hPrRad = new TH1D("hPrRad","p radius - gas; radius [mRad]",200,0,100);

  hBeamCh = new TH1D("hBeamCh","Beam Cherenkov",3,-.5,2.5);

  hNewCenterGas = new TH2D("hNewCenterGas","Map of new center of the event - Gas;x [mm];y [mm]",100,-5,5,100,-5,5);
  hNewCenterAero = new TH2D("hNewCenterAero","Map of new center of the event - Aerogel;x [mm];y [mm]",100,-5,5,100,-5,5);
  hNewCenterGasX= new TH1D("hNewCenterGasX","New center X - Gas;x [mm];counts",100,-5,5);
  hNewCenterGasY= new TH1D("hNewCenterGasY","New center Y - Gas;y [mm];counts",100,-5,5);
  hNewCenterAeroX= new TH1D("hNewCenterAeroX","New center X - Aero;x [mm];counts",100,-5,5);
  hNewCenterAeroY= new TH1D("hNewCenterAeroY","New center Y - Aero;y [mm];counts",100,-5,5);
  hGasdXVsAngle = new TH2D("hGasdXVsAngle","Gas - X mean semidifference Vs Angle; #theta_{x} [Rad];(#delta x)/3 [mm]",100,-0.005,0.005,100,-8,8);
  hGasdYVsAngle = new TH2D("hGasdYVsAngle","Gas - Y mean semidifference Vs Angle; #theta_{y} [Rad];(#delta y)/2 [mm]",100,-0.005,0.005,100,-8,8);
  hAerodXVsAngle = new TH2D("hAerodXVsAngle","Aero - X mean semidifference Vs Angle; #theta{x} [Rad];(#delta x)/2 [mm]",100,-0.005,0.005,100,-8,8);
  hAerodYVsAngle = new TH2D("hAerodYVsAngle","Aero - Y mean semidifference Vs Angle; #theta_{y} [Rad];(#delta y)/2 [mm]",100,-0.005,0.005,100,-8,8);
  //hAerodYVsAngle = new TH2D("hAerodYVsAngle","Aero - Y mean semidifference Vs Angle; #theta_{y} [Rad];(#delta y)/2 [mm]",50,-0,0.012,40,0,16);

  hSimPhotonLambdaGas = new TH1D("hSimPhotonLambdaGas","Simulation photon wavelength distribution - Gas;#lambda [nm];Hits [#]",600,150,750);
  hSimPhotonLambdaAero = new TH1D("hSimPhotonLambdaAero","Simulation photon wavelength distribution - Aero;#lambda [nm];Hits [#]",600,150,750);

  hA = new TH1D("hA","Number of photons - aerogel;photons/particle;counts",16,-.5,15.5);
  hG = new TH1D("hG","Number of photons - gas;photons/particle;counts",41,-.5,40.5);

  double AAA=2;
  double hInterpolateMin=120, hInterpolateMax=250;
  hInterpolateGasRing = new TH1D("hInterpolateGasRing","Gas ring exists;r [mRad]; counts",AAA*(hInterpolateMax-hInterpolateMin),hInterpolateMin,hInterpolateMax);
  hInterpolateNotGasRing =  new TH1D("hInterpolateNotGasRing","No gas ring;r [mRad]; counts",AAA*(hInterpolateMax-hInterpolateMin),hInterpolateMin,hInterpolateMax);
  hInterpolateIsProtonAndGasRing =  new TH1D("hInterpolateIsProtonAndGasRing","p tag w/ gas;r [mRad]; counts",AAA*(hInterpolateMax-hInterpolateMin),hInterpolateMin,hInterpolateMax);
  hInterpolateIsKaonAndGasRing =  new TH1D("hInterpolateIsKaonAndGasRing","k tag w/ gas;r [mRad]; counts",AAA*(hInterpolateMax-hInterpolateMin),hInterpolateMin,hInterpolateMax);
  hInterpolateIsPionAndNotGasRing = new TH1D("hInterpolateIsPionAndNotGasRing","#pi tag w/o gas;r [mRad]; counts",AAA*(hInterpolateMax-hInterpolateMin),hInterpolateMin,hInterpolateMax);

}




void displayBase(THeader *run){}
void displaySP(THeader *run){}
void displaySPN(THeader *run){}
void displayCUT(THeader *run){}
void displayRSD(THeader *run){}
void displayPhotonAnalysis(THeader *run){}
void fillHisto(THeader *run){}

