#include "TupleMaker_VHbbResonance/ZllTupleMaker.h"


// this is needed to distribute the algorithm to the workers
ClassImp(ZllTupleMaker)


ZllTupleMaker :: ZllTupleMaker ()
{
  ///Do not use this constructor
}

ZllTupleMaker :: ZllTupleMaker (std::string configPath) :
  RecoTupleMaker(configPath),
  m_tuple(0),
  leptonIsoOption(0)
{
  leptonIsoOption = config->get<int>("tuple.leptonIsoOption");
}


EL::StatusCode ZllTupleMaker :: initialize ()
{
  
  if(!m_tuple){
    //tuple must be created in a fully defined TupleMaker
    cout<<"ZllTupleMaker :: initialize() : Error, m_tuple should have been defined at this point"<<endl; 
    return EL::StatusCode::FAILURE;
  }
  RecoTupleMaker::m_tuple = m_tuple;
  
  return RecoTupleMaker::initialize();
}


EL::StatusCode ZllTupleMaker :: processEvent(){
  
  //Nothing to do here

  return RecoTupleMaker :: processEvent();
}


EL::StatusCode ZllTupleMaker :: fillJets ()
{
  //Note: l1p4 and l2p4 must be filled before

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
    if(Jp4.DeltaR(l2p4)<0.4) continue;

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

    if(Jp4.DeltaR(l1p4)<0.4) continue;
    if(Jp4.DeltaR(l2p4)<0.4) continue;

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
    Jp4.SetPtEtaPhiM(m_tuple->fatjet_pt[i],m_tuple->fatjet_eta[i],m_tuple->fatjet_phi[i],m_tuple->fatjet_m[i]);

    if(Jp4.DeltaR(l1p4)<1.0) continue;
    if(Jp4.DeltaR(l2p4)<1.0) continue;

    m_tuple->v_nfatjet++;
 
    if(m_tuple->fatjet_pt[i] > maxptfat){
      m_tuple->v_fatjet = i ;
      maxptfat = m_tuple->fatjet_pt[i];
    }
  }



  /////trk jets
  m_tuple->v_ntrkjet=0;
  m_tuple->v_trkjet=-1;
  float maxpttrk=0.;
  for(int i=0;i<m_tuple->ntrkjet;i++){
    
    TLorentzVector Jp4;
    Jp4.SetPtEtaPhiM(m_tuple->trkjet_pt[i],m_tuple->trkjet_eta[i],m_tuple->trkjet_phi[i],m_tuple->trkjet_m[i]);

    if(Jp4.DeltaR(l1p4)<.3) continue;
    if(Jp4.DeltaR(l2p4)<.3) continue;

    m_tuple->v_ntrkjet++;
 
    if(m_tuple->trkjet_pt[i] > maxpttrk){
      m_tuple->v_trkjet = i ;
      maxpttrk = m_tuple->trkjet_pt[i];
    }
  }

  return EL::StatusCode::SUCCESS;
}


