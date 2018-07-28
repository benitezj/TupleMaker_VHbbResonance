#ifndef BASETUPLEMAKER_H
#define BASETUPLEMAKER_H

// ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TLorentzVector.h>

// Event Loop includes
#include "EventLoop/Algorithm.h"
#include "EventLoop/StatusCode.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "PileupReweighting/PileupReweightingTool.h"
#include "PileupReweighting/TPileupReweighting.h"


#ifndef __MAKECINT__
#include "MuonEfficiencyCorrections/MuonTriggerScaleFactors.h"
#include "CxAODTools/CommonProperties.h"
#include "CxAODTools_DB/DBProperties.h"
#include "CxAODTools/TriggerTool.h"
#endif // not __MAKECINT__

class TFile;
class TTree;

#include "CxAODTools/ConfigStore.h"

#include "TupleMaker_VHbbResonance/BaseTuple.h"

class BaseTupleMaker : public EL::Algorithm
{
public:
  BaseTupleMaker ();//Default Do not use.
  BaseTupleMaker (std::string configPath);
  virtual ~BaseTupleMaker();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  
  // this is needed to distribute the algorithm to the workers
  ClassDef(BaseTupleMaker, 1);


protected:

  ConfigStore* config = 0;

  //for debugging
  bool m_debug;
  int m_printTruth;
  int m_maxEvents;
  int m_randRun;

  std::vector<std::string> m_systNames;
  std::string m_currentSyst;
  
  BaseTuple * m_tuple = 0; //!
 
  xAOD::TEvent * m_event; //!
  TTree * m_tree; //!
  std::string m_filename; // output filename 

  CP::TPileupReweighting* m_tprw; //!
  CP::PileupReweightingTool *m_pileupreweighting; //!
  std::string m_ilumicalcFile ;
  std::string m_mcPileUpFile ;

  std::string m_EventInfoIn;
  const xAOD::EventInfo * m_eventInfo;//!
  EL::StatusCode processBaseEventInfo();

  virtual EL::StatusCode processEvent();

  void incrementCounter(const char * name, float w=1.);
  float getCounter(const char * name);

  void PrintEventFlags(){
    cout<<"PrintEventFlags: isCleanEvent="<<m_tuple->eve_isCleanEvent<<"  passTrig="<<m_tuple->eve_passTrig<<" passTrigMatch="<<m_tuple->eve_passTrigMatch<<endl;
  }
  void printMCTruth();

  double deltaR(float meta, float mphi, float peta, float pphi) {
    double deta= fabs(meta - peta), dphi= fabs(mphi - pphi);
    if (dphi > TMath::Pi() ) dphi = 2*(TMath::Pi()) - dphi;
    return sqrt((dphi*dphi)+(deta*deta));
  }
  double mT(TLorentzVector l1, TLorentzVector l2) {
    return sqrt( pow(l1.Pt() + l2.Pt(), 2) -
		 pow(l1.Px() + l2.Px(), 2) -
		 pow(l1.Py() + l2.Py(), 2) );
  }

  double smaller(double a, double b){
    if(a<b) return a;
    if(b<a) return b;
    return a;
  }


  //Trigger
  std::vector<std::string> m_triggerPaths;

  //Trigger SF
#ifndef __MAKECINT_
  TriggerTool * m_triggerTool;// !
#endif // not __MAKECINT_
  
  CP::MuonTriggerScaleFactors *  m_trig_sfmuon; 


  ///Particle masses
  float m_HMass=125000.;
  float m_ZMass=91200.;
  float m_WMass=80400.;
  

private:
 

  //sytematic variations
  std::vector< TTree* > m_systTuples;
     
  std::vector<unsigned int> m_runsToAnalyze;//These two should be filled in sync
  std::vector<ULong64_t> m_eventsToAnalyze;

  std::vector<string> counterNames;
  std::vector<float> counters;
  void printCounters();

  //event counts 
  TH1D* HEventCounts; //!

  //Pileup histos
  TH1F* HPileUp; //!
  
  ///for pile-up weights
  std::vector<std::string> pileUpMCFiles;
  std::vector<std::string> pileUpDataFiles;
  TH1F* PileUpRatios[BaseTuple::MAXMCPUWEIGHTS]; //!
  EL::StatusCode initPileUpWeights(TFile * f);
  EL::StatusCode fillPileUpWeights();

  ///Generator info
  std::string m_MCTruthIn;
  EL::StatusCode processMCTruth();
  void fillFakeComposite(int * p, int l1, int l2, int pdgid);
  void fillTruthParticle(int index, TLorentzVector P, int status, int pdg);
  void fillMCProcesses();
  
  TH1F* HGen_Zmm_m; //!
  TH1F* HGen_Zmm_pt; //!
  TH1F* HGen_Zmm_eta; //!
  TH1F* HGen_Zmm_m1_pt; //!
  TH1F* HGen_Zmm_m1_eta; //!
  TH1F* HGen_Zmm_m2_pt; //!
  TH1F* HGen_Zmm_m2_eta; //!

  TH1F* HGen_Zee_m; //!
  TH1F* HGen_Zee_pt; //!
  TH1F* HGen_Zee_eta; //!
  TH1F* HGen_Zee_e1_pt; //!
  TH1F* HGen_Zee_e1_eta; //!
  TH1F* HGen_Zee_e2_pt; //!
  TH1F* HGen_Zee_e2_eta; //!

  TH1F* HGen_Ztt_m; //!
  TH1F* HGen_Ztt_pt; //!
  TH1F* HGen_Ztt_eta; //!
  TH1F* HGen_Ztt_t1_pt; //!
  TH1F* HGen_Ztt_t1_eta; //!
  TH1F* HGen_Ztt_t2_pt; //!
  TH1F* HGen_Ztt_t2_eta; //!

  TH1F* HGen_Wmv_m; //!
  TH1F* HGen_Wmv_pt; //!
  TH1F* HGen_Wmv_eta; //!
  TH1F* HGen_Wmv_m_pt; //!
  TH1F* HGen_Wmv_m_eta; //!
  TH1F* HGen_Wmv_v_pt; //!
  TH1F* HGen_Wmv_v_eta; //!

  TH1F* HGen_Wev_m; //!
  TH1F* HGen_Wev_pt; //!
  TH1F* HGen_Wev_eta; //!
  TH1F* HGen_Wev_e_pt; //!
  TH1F* HGen_Wev_e_eta; //!
  TH1F* HGen_Wev_v_pt; //!
  TH1F* HGen_Wev_v_eta; //!

  TH1F* HGen_Wtv_m; //!
  TH1F* HGen_Wtv_pt; //!
  TH1F* HGen_Wtv_eta; //!
  TH1F* HGen_Wtv_t_pt; //!
  TH1F* HGen_Wtv_t_eta; //!
  TH1F* HGen_Wtv_v_pt; //!
  TH1F* HGen_Wtv_v_eta; //!

  TH1F* HGen_H0bb_m; //!
  TH1F* HGen_H0bb_pt; //!
  TH1F* HGen_H0bb_eta; //!
  TH1F* HGen_H0bb_b1_pt; //!
  TH1F* HGen_H0bb_b1_eta; //!
  TH1F* HGen_H0bb_b2_pt; //!
  TH1F* HGen_H0bb_b2_eta; //!

  TH1F* HGen_ZH_m; //!
  TH1F* HGen_ZH_pt; //!
  TH1F* HGen_ZH_eta; //!

  TH1F* HGen_WH_m; //!
  TH1F* HGen_WH_pt; //!
  TH1F* HGen_WH_eta; //!


  ///VH Resonance analysis
  TH1F* HGen_RToVH_m; //!
  TH1F* HGen_RToVH_pt; //!
  TH1F* HGen_RToVH_eta; //!

  TH1F* HGen_RToVH_V_m; //!
  TH1F* HGen_RToVH_V_pt; //!
  TH1F* HGen_RToVH_V_eta; //!

  TH1F* HGen_RToVH_H_m; //!
  TH1F* HGen_RToVH_H_pt; //!
  TH1F* HGen_RToVH_H_eta; //!

  TH1F* HGen_RToVH_H_b1_m; //!
  TH1F* HGen_RToVH_H_b1_pt; //!
  TH1F* HGen_RToVH_H_b1_eta; //!

  TH1F* HGen_RToVH_H_b2_m; //!
  TH1F* HGen_RToVH_H_b2_pt; //!
  TH1F* HGen_RToVH_H_b2_eta; //!

  TH1F* HGen_RToVH_H_p; //!
  TH1F* HGen_RToVH_H_dR; //!
  TH2F* HGen_RToVH_H_dRVsp; //!

  EL::StatusCode bookGenHistograms(TFile *f);
  void fillGenHistograms();


};

#endif
