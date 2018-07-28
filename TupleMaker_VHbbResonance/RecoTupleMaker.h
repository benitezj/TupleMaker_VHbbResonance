#ifndef RECOTUPLEMAKER_H
#define RECOTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/BaseTupleMaker.h"
#include "TupleMaker_VHbbResonance/RecoTuple.h"

#include "xAODMuon/Muon.h"
#include "xAODEgamma/Electron.h"
#include "xAODJet/Jet.h"
#include "xAODMissingET/MissingETComponentMap.h"
#include "xAODMissingET/MissingETContainer.h"

#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingET.h"

#ifndef __MAKECINT__
#include "CxAODTools/BTaggingTool.h"
#endif // not __MAKECINT__


class RecoTupleMaker : public BaseTupleMaker
{
public:
  RecoTupleMaker ();//Default Do not use.
  RecoTupleMaker (std::string configPath);
  ~RecoTupleMaker ();

  virtual EL::StatusCode initialize ();
 
  // this is needed to distribute the algorithm to the workers
  ClassDef(RecoTupleMaker, 1);


protected:

  RecoTuple * m_tuple = 0; //!

#ifndef __MAKECINT_
  BTaggingTool *m_bTagTool;  // !
#endif // not __MAKECINT_

  EL::StatusCode processRecoEventInfo();
  virtual EL::StatusCode processEvent();


  bool passMuJetForOR(int mu, int jet);
  bool passMuonOR(int idx);
  bool passElectronOR(int idx);
  bool passJetOR(int idx);
  bool passLooseMuon(int idx,int iso=1);
  bool passLooseElectron(int idx,int iso=1);
  int countLooseLeptons();

  bool passCentralJet(int i){
    if( fabs(m_tuple->jet_eta[i]) > 2.5) return false;
    if( m_tuple->jet_pt[i] < 20000 ) return false;
    if( fabs(m_tuple->jet_eta[i]) < 2.4  && m_tuple->jet_pt[i] < 60000 && m_tuple->jet_jvt[i] < 0.59 ) return false;
    //if( fabs(m_tuple->jet_eta[i]) < 2.4  && m_tuple->jet_pt[i] < 50000 && m_tuple->jet_jvt[i] < 0.64 ) return false;
    return true;
  } 
  bool passFWJet(int i){		       
    if( fabs(m_tuple->jet_eta[i]) <= 2.5 ) return false;
    if( fabs(m_tuple->jet_eta[i]) > 4.5 ) return false;
    if( m_tuple->jet_pt[i] < 30000 ) return false;
    return true;
  }

  bool passFatJet(int i){
    if( fabs(m_tuple->fatjet_eta[i]) > 2.0 ) return false;
    if( m_tuple->fatjet_pt[i] < 200000 ) return false;
    return true;
  }

  bool passTrackJet(int i){
    //used in Xbb Note
    if( fabs(m_tuple->trkjet_eta[i]) > 2.5 ) return false;//NOTE: Xbb code uses abs(), but gives same result
    if( m_tuple->trkjet_pt[i] < 10000 ) return false;
    if( m_tuple->trkjet_numConstituents[i] < 2 ) return false;
    return true;
  } 

  bool passBTag(float val,int algo){
    //// mv2c20 %70, cut from  JetSubStructureUtils/data/config_13TeV_Htagging_MC15_Prerecommendations_20150812.dat
    //if(algo==2) return static_cast<int>(val  > -0.3098) == 1 ? true : false ;   
    if(algo==2) return (val  > -0.3098) ;  // gives same result as above
    return false;
  }
  

  bool passBadJetVeto();//use Akt0.4EM jets to check if there is bad one then remove event

  
  void PrintElectron(int i){
    std::cout<<"PrintElectron: pT="<<m_tuple->ele_pt[i]<<" eta="<<m_tuple->ele_eta[i]<<" phi="<<m_tuple->ele_phi[i]<<" PV="<<m_tuple->ele_isFromPV[i] <<" shtrk="<<m_tuple->ele_sharedtrk[i]<<" ptiso="<<m_tuple->ele_ptiso[i]<<" etiso="<<m_tuple->ele_etiso[i]<<std::endl;
  }
  void PrintMuon(int i){
    std::cout<<"PrintMuon: pT="<<m_tuple->muo_pt[i]<<" eta="<<m_tuple->muo_eta[i]<<" phi="<<m_tuple->muo_phi[i]<<" PV="<<m_tuple->muo_isFromPV[i]<<" ptiso="<<m_tuple->muo_ptiso[i]<<" etiso="<<m_tuple->muo_etiso[i]<<std::endl;
  }
  void PrintJet(int i){
    std::cout<<"PrintJet: pT="<<m_tuple->jet_pt[i]<<" eta="<<m_tuple->jet_eta[i]<<" phi="<<m_tuple->jet_phi[i]<<" JVT="<<m_tuple->jet_jvt[i]<<" goodJet="<<m_tuple->jet_good[i]<<" SumPtTrkPt500PV="<<m_tuple->jet_SumPtTrkPt500PV[i]<<" NumTrkPt500PV="<<m_tuple->jet_NumTrkPt500PV[i]<<" btag="<<m_tuple->jet_btag[i]<<std::endl;
  }
  void PrintTrkJet(int i){
    std::cout<<"PrintTrkJet: pT="<<m_tuple->trkjet_pt[i]<<" eta="<<m_tuple->trkjet_eta[i]<<" phi="<<m_tuple->trkjet_phi[i]<<" numConst="<<m_tuple->trkjet_numConstituents[i]<<" btag="<<m_tuple->trkjet_btag[i]<<std::endl;
  }
  void PrintFatJet(int i){
    std::cout<<"PrintFatJet: pT="<<m_tuple->fatjet_pt[i]<<" eta="<<m_tuple->fatjet_eta[i]<<" phi="<<m_tuple->fatjet_phi[i]<<" M="<<m_tuple->fatjet_m[i]<<" Mcor="<<m_tuple->fatjet_m_cor[i]<<" nTrkJ="<<m_tuple->fatjet_ntrkjet[i];
    if( m_tuple->fatjet_ltrkjet[i]>=0){std::cout<<"(pT1="<<m_tuple->trkjet_pt[m_tuple->fatjet_ltrkjet[i]]<<")";}
    if( m_tuple->fatjet_sltrkjet[i]>=0){std::cout<<"(pT2="<<m_tuple->trkjet_pt[m_tuple->fatjet_sltrkjet[i]]<<")";}  
    std::cout<<std::endl;
  }
  void PrintMM(int i){
    std::cout<<"PrintMM: pT="<<m_tuple->mm_pt[i]<<" eta="<<m_tuple->mm_eta[i]<<" phi="<<m_tuple->mm_phi[i]<<" M="<<m_tuple->mm_m[i]<<" dR="<<m_tuple->mm_dR[i]<<std::endl;
    PrintMuon(m_tuple->mm_leg1[i]);
    PrintMuon(m_tuple->mm_leg2[i]);
  }
  void PrintEE(int i){
    std::cout<<"PrintEE: pT="<<m_tuple->ee_pt[i]<<" eta="<<m_tuple->ee_eta[i]<<" phi="<<m_tuple->ee_phi[i]<<" M="<<m_tuple->ee_m[i]<<" dR="<<m_tuple->ee_dR[i]<<std::endl;
    PrintElectron(m_tuple->ee_leg1[i]);
    PrintElectron(m_tuple->ee_leg2[i]);
  }

  void PrintEventParticles();
  void PrintElectrons();
  void PrintMuons();
  void PrintJets();
  void PrintFatJets();
  void PrintTrkJets();
  void PrintMET();
  

  const xAOD::ElectronContainer * m_Electrons = 0;
  const xAOD::MuonContainer * m_Muons = 0;
  const xAOD::JetContainer * m_Jets = 0;
  const xAOD::JetContainer * m_FatJets = 0;
  const xAOD::JetContainer * m_TrackJets = 0;
  const xAOD::MissingETContainer* m_MissingET = 0;


private:

  std::string m_electronsNameIn;
  std::string m_muonsNameIn;
  std::string m_tausNameIn;
  std::string m_jetsNameIn;
  std::string m_fatJetsNameIn;
  std::string m_trkJetsNameIn;
  std::string m_missingETIn;
  
  std::vector<std::string> m_electronsSysts;
  std::vector<std::string> m_muonsSysts;
  std::vector<std::string> m_tausSysts;  
  std::vector<std::string> m_jetsSysts;
  std::vector<std::string> m_fatJetsSysts;
  std::vector<std::string> m_trkJetsSysts;
  std::vector<std::string> m_missingETSysts;


  EL::StatusCode processElectrons();
  EL::StatusCode processMuons();
  // EL::StatusCode processTaus();
  EL::StatusCode processJets();
  EL::StatusCode processFatJets();
  EL::StatusCode processTrackJets();
  EL::StatusCode processMET();
  EL::StatusCode processEE();  
  EL::StatusCode processMM();  
  EL::StatusCode processEM();  
  EL::StatusCode processJJ();  
  EL::StatusCode processEMET();  
  EL::StatusCode processMUMET();  
  
  EL::StatusCode processParticles();
  //bool passRecoEventCleaning();


};

#endif
