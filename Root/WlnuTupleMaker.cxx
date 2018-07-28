#include "TupleMaker_VHbbResonance/WlnuTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(WlnuTupleMaker)

WlnuTupleMaker :: WlnuTupleMaker (){}

WlnuTupleMaker :: WlnuTupleMaker (std::string configPath) :
  RecoTupleMaker(configPath),
  m_tuple(0),
  leptonIsoOption(0)
{
  
  leptonIsoOption = config->get<int>("tuple.leptonIsoOption");

}

EL::StatusCode WlnuTupleMaker :: initialize ()
{
  if(!m_tuple){
    //tuple must be created in a higher level
    cout<<"WlnuTupleMaker :: initialize() : Error, m_tuple should have been defined at this point"<<endl; 
    return EL::StatusCode::FAILURE;
  }
  RecoTupleMaker::m_tuple = m_tuple;
  
  return RecoTupleMaker::initialize();
}


EL::StatusCode WlnuTupleMaker :: processEvent(){
  
  //Nothing to do here

  return RecoTupleMaker :: processEvent();
}

EL::StatusCode WlnuTupleMaker :: fillJets ()
{

  /////Central jets
  m_tuple->v_njet=0;
  m_tuple->v_jet=-1;
  float maxpt=0.;
  for(int i=0;i<m_tuple->njet;i++){
    if(!passCentralJet(i)) continue;
    
    TLorentzVector Jp4;
    Jp4.SetPtEtaPhiM(m_tuple->jet_pt[i],m_tuple->jet_eta[i],m_tuple->jet_phi[i],m_tuple->jet_m[i]);

    //check overlap with muons
    if(Jp4.DeltaR(l1p4)<0.4) continue;

    m_tuple->v_njet++;
 
    if(m_tuple->jet_pt[i] > maxpt){
      m_tuple->v_jet = i ;
      maxpt = m_tuple->jet_pt[i];
    }
  }


  /////Forward jets
  m_tuple->v_njetfw=0;
  m_tuple->v_jetfw=-1;
  float maxptfw=0.;
  for(int i=0;i<m_tuple->njet;i++){
    if(!passFWJet(i)) continue;
    
    TLorentzVector Jp4;
    Jp4.SetPtEtaPhiM(m_tuple->jet_pt[i],m_tuple->jet_eta[i],m_tuple->jet_phi[i],m_tuple->jet_m[i]);

    //check overlap with muons
    if(Jp4.DeltaR(l1p4)<0.4) continue;

    m_tuple->v_njetfw++;
 
    if(m_tuple->jet_pt[i] > maxptfw){
      m_tuple->v_jetfw = i ;
      maxptfw = m_tuple->jet_pt[i];
    }
  }


  /////Fat jets
  m_tuple->v_nfatjet=0;
  m_tuple->v_fatjet=-1;
  float maxptfat=0.;
  for(int i=0;i<m_tuple->nfatjet;i++){
    
    TLorentzVector Jp4;
    Jp4.SetPtEtaPhiM(m_tuple->jet_pt[i],m_tuple->jet_eta[i],m_tuple->jet_phi[i],m_tuple->jet_m[i]);

    if(Jp4.DeltaR(l1p4)<1.0) continue;
    m_tuple->v_nfatjet++;
 
    if(m_tuple->fatjet_pt[i] > maxptfat){
      m_tuple->v_fatjet = i ;
      maxptfat = m_tuple->fatjet_pt[i];
    }
  }

  return EL::StatusCode::SUCCESS;
}


