#include "TupleMaker_VHbbResonance/ZHmmbbTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(ZHmmbbTupleMaker)


ZHmmbbTupleMaker :: ZHmmbbTupleMaker ()
{
  ///Do not use this constructor
}

ZHmmbbTupleMaker :: ZHmmbbTupleMaker (std::string configPath) :
  ZHllbbTupleMaker(configPath),
  m_tuple(0)
{
}

EL::StatusCode ZHmmbbTupleMaker :: initialize ()
{
  //this is the top level, ntuple must be defined here
  m_tuple = new ZHmmbbTuple();
  ZHllbbTupleMaker::m_tuple = m_tuple;
    
  return ZHllbbTupleMaker::initialize(); 
}



EL::StatusCode ZHmmbbTupleMaker :: processMMBB(){
  //Note: Muon array must be filled prior
  m_tuple->nvh=0;

  int ivh=0;
  for(int z = 0; z < m_tuple->nmm ; z++){
    for(int h = 0; h < m_tuple->njj ; h++){
      if( ivh < VHTuple::MAXVH) {

	TLorentzVector zp4;
	zp4.SetPtEtaPhiM(m_tuple->mm_pt[z],m_tuple->mm_eta[z],m_tuple->mm_phi[z],m_tuple->mm_m[z]);
	TLorentzVector hp4;
	hp4.SetPtEtaPhiM(m_tuple->jj_pt[h],m_tuple->jj_eta[h],m_tuple->jj_phi[h],m_tuple->jj_m[h]);
	TLorentzVector vhp4 = zp4 + hp4;

	m_tuple->vh_charge[ivh] = m_tuple->mm_charge[z] + m_tuple->jj_charge[h];
	m_tuple->vh_E[ivh] = vhp4.E();
	m_tuple->vh_p[ivh] = vhp4.P();
	m_tuple->vh_pt[ivh] = vhp4.Pt(); 
	m_tuple->vh_phi[ivh] = vhp4.Phi();
	m_tuple->vh_eta[ivh] = vhp4.Eta();
	m_tuple->vh_dR[ivh] = zp4.DeltaR(hp4);
          
	m_tuple->vh_v[ivh] = z;
	m_tuple->vh_h[ivh] = h;

	m_tuple->vh_dPhiHMET[ivh]= m_tuple->jj_phi[h] - m_tuple->met_phi;

	TLorentzVector hp4Corr;
	hp4Corr.SetPtEtaPhiM(m_tuple->jj_pt[h] * 125000 / m_tuple->jj_m[h], m_tuple->jj_eta[h],m_tuple->jj_phi[h], 125000);

	m_tuple->vh_m[ivh] = (zp4 + hp4Corr).M();
	m_tuple->vh_m1[ivh] = vhp4.M();          

	++ivh;
      }else{
	Info( "processMMBB()" , " Trying to build too many mm-jj pairs, max is %i" , VHTuple::MAXVH);
      }
    }
  }
  m_tuple->nvh=ivh;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ZHmmbbTupleMaker :: processEvent(){

  if( ZHllbbTupleMaker :: processEvent()  == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;
  

  //>2 Loose Leptons Veto
  if(countLooseLeptons()>2) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_LooseLeptonsVeto");

  
  //init lists
  std::vector<int> list;
  std::vector<int> listsel;
  for(int z = 0; z < m_tuple->nmm ; z++){
    list.push_back(z);
  }
  if(list.size()==0) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoDiLepton");


  //lepton id
  for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
    if( !passLooseMuon( m_tuple->mm_leg1[*it], leptonIsoOption) ) continue;
    if( !passLooseMuon( m_tuple->mm_leg2[*it], leptonIsoOption) ) continue;

    listsel.push_back(*it);
  }
  if(listsel.size()==0) return EL::StatusCode::FAILURE;
  list = listsel;
  listsel.clear();
  incrementCounter("eventCounter_2LeptonId");


  if(leptonIsoOption == 100){
    //Lepton isolation with fix
    for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
      if( (m_tuple->muo_ptiso[m_tuple->mm_leg1[*it]]-m_tuple->muo_trkpt[m_tuple->mm_leg2[*it]]*(m_tuple->mm_dR[*it]<0.2))/m_tuple->muo_pt[m_tuple->mm_leg1[*it]] > 0.1 ) continue;
      if( (m_tuple->muo_ptiso[m_tuple->mm_leg2[*it]]-m_tuple->muo_trkpt[m_tuple->mm_leg1[*it]]*(m_tuple->mm_dR[*it]<0.2))/m_tuple->muo_pt[m_tuple->mm_leg2[*it]] > 0.1 ) continue;
      listsel.push_back(*it);
    }
    if(listsel.size()==0) return EL::StatusCode::FAILURE;
    list=listsel;
    listsel.clear();
    incrementCounter("eventCounter_LeptonIsolation");
  }


  int l1=m_tuple->mm_leg1[list[0]];
  int l2=m_tuple->mm_leg2[list[0]];

  l1p4.SetPtEtaPhiM(m_tuple->muo_pt[l1],
  		    m_tuple->muo_eta[l1],
  		    m_tuple->muo_phi[l1],
  		    m_tuple->muo_m[l1]);
  l2p4.SetPtEtaPhiM(m_tuple->muo_pt[l2],
  		    m_tuple->muo_eta[l2],
  		    m_tuple->muo_phi[l2],
  		    m_tuple->muo_m[l2]);


  //Find the b-jets and fill extra jets
  if( selectJets() == EL::StatusCode::FAILURE ) 
    return EL::StatusCode::FAILURE;


  ///Create the llbb candidates
  if( processMMBB()  == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;
  

  //Find the selected VH : Note leptons and Jets must be pT ordered
  m_tuple->vh=-1;
  for(int i=0;i<m_tuple->nvh;i++){
    if(m_tuple->mm_leg1[m_tuple->vh_v[i]] == l1 &&
       m_tuple->mm_leg2[m_tuple->vh_v[i]] == l2 &&
       m_tuple->jj_leg1[m_tuple->vh_h[i]] == jet1 &&
       m_tuple->jj_leg2[m_tuple->vh_h[i]] == jet2 
       ) 
      m_tuple->vh=i;
  }
  if( m_tuple->vh==-1){ 
    cout<<"processEvent: something is wrong with the VH logic"<<endl;
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
}


