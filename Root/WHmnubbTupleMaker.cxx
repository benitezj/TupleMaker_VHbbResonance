#include "TupleMaker_VHbbResonance/WHmnubbTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(WHmnubbTupleMaker)


WHmnubbTupleMaker :: WHmnubbTupleMaker ()
{
  ///Do not use this constructor
}

WHmnubbTupleMaker :: WHmnubbTupleMaker (std::string configPath) :
  WHlnubbTupleMaker(configPath),
  m_tuple(0)
{
}

EL::StatusCode WHmnubbTupleMaker :: initialize ()
{
  //this is the top level, ntuple must be defined here
  m_tuple = new WHmnubbTuple();
  WHlnubbTupleMaker::m_tuple = m_tuple;
    
  return WHlnubbTupleMaker::initialize(); 
}



EL::StatusCode WHmnubbTupleMaker :: processMuMETBB(){
  //Note: Muon array must be filled prior
  m_tuple->nvh=0;

  int ivh=0;
  for(int w = 0; w < m_tuple->nmuMet ; w++){
    for(int h = 0; h < m_tuple->njj ; h++){
      if( ivh < VHTuple::MAXVH) {

	TLorentzVector wp4;
	wp4.SetPtEtaPhiM(m_tuple->muMet_pt[w],m_tuple->muMet_eta[w],m_tuple->muMet_phi[w],m_tuple->muMet_m[w]);

	TLorentzVector hp4;
	hp4.SetPtEtaPhiM(m_tuple->jj_pt[h],m_tuple->jj_eta[h],m_tuple->jj_phi[h],m_tuple->jj_m[h]);

	TLorentzVector vhp4 = wp4 + hp4;

	m_tuple->vh_charge[ivh] = m_tuple->muMet_charge[w] + m_tuple->jj_charge[h];
	m_tuple->vh_E[ivh] = vhp4.E();
	m_tuple->vh_p[ivh] = vhp4.P();
	m_tuple->vh_pt[ivh] = vhp4.Pt(); 
	m_tuple->vh_phi[ivh] = vhp4.Phi();
	m_tuple->vh_eta[ivh] = vhp4.Eta();
	m_tuple->vh_dR[ivh] = wp4.DeltaR(hp4);
	
	m_tuple->vh_v[ivh] = w;
	m_tuple->vh_h[ivh] = h;

	m_tuple->vh_dPhiHMET[ivh]= m_tuple->jj_phi[h] - m_tuple->met_phi;


        TLorentzVector hp4Corr;
        hp4Corr.SetPtEtaPhiM(m_tuple->jj_pt[h] * 125000 / m_tuple->jj_m[h], m_tuple->jj_eta[h],m_tuple->jj_phi[h], 125000);

        m_tuple->vh_m[ivh] = (wp4 + hp4Corr).M();
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

EL::StatusCode WHmnubbTupleMaker :: processEvent(){
  
  if( WHlnubbTupleMaker :: processEvent()  == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;
  

  //>1 Loose Leptons Veto
  if(countLooseLeptons()>1) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_LooseLeptonsVeto");

  //Tight Muon
  int l1=-1;
  for (int i = 0; i < m_tuple->nmuo ; i++){
    if( !passLooseMuon(i,leptonIsoOption) ) continue;
    l1=i;
  }
  if(l1==-1) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_1TightLepton");


  //Set the selected leptons
  lp4.SetPtEtaPhiM(m_tuple->muo_pt[l1],
		   m_tuple->muo_eta[l1],
		   m_tuple->muo_phi[l1],
		   m_tuple->muo_m[l1]);
  

  //Find the b-jets and fill extra jets
  if( selectJets() == EL::StatusCode::FAILURE ) 
    return EL::StatusCode::FAILURE;


  ///Create the lnubb candidates
  if( processMuMETBB()  == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;


  //Find the selected VH : Note leptons and Jets must be pT ordered
  m_tuple->vh=-1;
  for(int i=0;i<m_tuple->nvh;i++){
    if(m_tuple->muMet_mu[m_tuple->vh_v[i]] == l1 &&
       m_tuple->jj_leg1[m_tuple->vh_h[i]] == jet1 &&
       m_tuple->jj_leg2[m_tuple->vh_h[i]] == jet2 
       ) 
      m_tuple->vh=i;
  }
  if( m_tuple->vh==-1) {
    cout<<"processEvent: something is wrong with the VH logic"<<endl; 
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
}




  // //init list of candidates
  // std::vector<int> list;
  // std::vector<int> listsel;  
  // for(int i = 0; i < m_tuple->nmuo ; i++){
  //   list.push_back(i);
  // }
  // if(list.size()==0) return EL::StatusCode::FAILURE;
  // incrementCounter("eventCounter_NoLeptons");

  // //Lepton id 
  // for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
  //   if( m_tuple->muo_idVHl[*it] == false ) continue;
  //   listsel.push_back(*it);
  // }
  // if(listsel.size()==0) return EL::StatusCode::FAILURE;
  // list=listsel;
  // listsel.clear();
  // incrementCounter("eventCounter_LeptonId");
  
  // ///Lepton pT cut
  // for(std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
  //   if( m_tuple->muo_pt[*it] < 25000 ) continue;
  //   listsel.push_back(*it);
  // }
  // if(listsel.size()==0) return EL::StatusCode::FAILURE;
  // list=listsel;
  // listsel.clear();
  // incrementCounter("eventCounter_LeptonPt");

  // //Lepton eta cut
  // for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
  //   if( fabs(m_tuple->muo_eta[*it]) > 2.5  ) continue;
  //   listsel.push_back(*it);
  // }
  // if(listsel.size()==0) return EL::StatusCode::FAILURE;
  // list=listsel;
  // listsel.clear();
  // incrementCounter("eventCounter_LeptonEta");


  // //Lepton track isolation
  // for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
  //   if( m_tuple->muo_ptiso[*it]/m_tuple->muo_pt[*it] > 0.04  ) continue;
  //   listsel.push_back(*it);
  // }
  // if(listsel.size()==0) return EL::StatusCode::FAILURE;
  // list=listsel;
  // listsel.clear();
  // incrementCounter("eventCounter_LeptonPtIso");
  
  // //Lepton calo isolation
  // for (std::vector<int>::iterator it = list.begin() ; it != list.end(); ++it){
  //   if( m_tuple->muo_etiso[*it]/m_tuple->muo_pt[*it] > 0.04 ) continue;
  //   listsel.push_back(*it);
  // }
  // if(listsel.size()==0) return EL::StatusCode::FAILURE;
  // list=listsel;
  // listsel.clear();
  // incrementCounter("eventCounter_LeptonCaloIso");

  // //Remove events with more than 2 leptons
  // if(list.size()>1) return EL::StatusCode::FAILURE;
  // incrementCounter("eventCounter_ManyLeptons");
