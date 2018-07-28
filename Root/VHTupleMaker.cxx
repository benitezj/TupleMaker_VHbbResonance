#include "TupleMaker_VHbbResonance/VHTupleMaker.h"

// this is needed to distribute the algorithm to the workers
ClassImp(VHTupleMaker)


VHTupleMaker :: VHTupleMaker ()
{
  ///Do not use this constructor
}

VHTupleMaker :: VHTupleMaker (std::string configPath) :
  RecoTupleMaker(configPath),
  m_tuple(0),
  leptonIsoOption(0),
  jet1(-1),
  jet2(-1),
  bTagOption(1),
  fatJetPtCut(0.)
{

  leptonIsoOption = config->get<int>("tuple.leptonIsoOption");

  //resolved analysis
  bTagOption  = config->get<int>("tuple.bTag");

  //merged analysis
  fatJetPtCut  = config->get<float>("tuple.fatJetPtCut");

}

EL::StatusCode VHTupleMaker :: initialize ()
{
  if( m_debug) std::cout<<" VHTupleMaker :: initialize "<<std::endl;

  if(!m_tuple){
    //tuple must be created in a fully defined TupleMaker
    cout<<"VHTupleMaker :: initialize() : Error, m_tuple should have been defined at this point"<<endl; 
    return EL::StatusCode::FAILURE;
  }
  RecoTupleMaker::m_tuple = m_tuple;
  
  return RecoTupleMaker::initialize();
}


EL::StatusCode VHTupleMaker :: processEvent(){

  return RecoTupleMaker :: processEvent();
}


void VHTupleMaker :: fillFWJets ()
{
  m_tuple->vh_njetfw=0;
  m_tuple->vh_jetfw=-1;
  float maxptfw=0.;
  for(int i=0;i<m_tuple->njet;i++){
    if(!passFWJet(i)) continue;
    m_tuple->vh_njetfw++;
    if(m_tuple->jet_pt[i] > maxptfw){
      m_tuple->vh_jetfw = i ;
      maxptfw = m_tuple->jet_pt[i];
    }
  }
}



void VHTupleMaker :: fillCentralJetsVeto()
{
  m_tuple->vh_njetveto=0;
  for(int i=0;i<m_tuple->njet;i++){
    if(!passCentralJet(i)) continue;
    if(this->jetOverlap(i)) continue;
    m_tuple->vh_njetveto++;
  }
}


std::vector<int>  VHTupleMaker :: fillCentralJets ()
{
  std::vector<int> list;
  m_tuple->vh_njet=0;
  m_tuple->vh_jet1=-1;
  float maxpt=0.;
  for(int i=0;i<m_tuple->njet;i++){
    if(!passCentralJetVH(i)) continue;
    if(this->jetOverlap(i)) continue;
    m_tuple->vh_njet++;
    list.push_back(i);

    if(m_tuple->jet_pt[i] > maxpt){
      m_tuple->vh_jet1 = i ;
      maxpt = m_tuple->jet_pt[i];
    }

  }

  // if(m_tuple->vh_njet == 3) //check the pT ordering
  //   cout<<"j1="<<m_tuple->jet_pt[m_tuple->vh_jet1]
  // 	<<", j2="<<m_tuple->jet_pt[m_tuple->vh_jet2]
  // 	<<", j3="<<m_tuple->jet_pt[m_tuple->vh_jet3]<<endl;
  

  return list;
}


void VHTupleMaker :: fillBJets ()
{
  m_tuple->vh_nbjet=0;
  m_tuple->vh_bjet1=-1;
  m_tuple->vh_bjet2=-1;
  for(int i=0;i<m_tuple->njet;i++){
    if(!passCentralJetVH(i)) continue;
    if(this->jetOverlap(i)) continue;
    if(m_tuple->jet_btag[i] == 0) continue; 
    m_tuple->vh_nbjet++;
    if(m_tuple->vh_bjet1<0) m_tuple->vh_bjet1  = i ;
    else m_tuple->vh_bjet2 = i ;
  }

  //must be pT ordered for later
  if(m_tuple->vh_bjet2>0 && 
     m_tuple->jet_pt[m_tuple->vh_bjet2] > m_tuple->jet_pt[ m_tuple->vh_bjet1] )
    {
      int tmp= m_tuple->vh_bjet1;
      m_tuple->vh_bjet1 = m_tuple->vh_bjet2 ;
      m_tuple->vh_bjet2 = tmp;
    }
  
}

EL::StatusCode  VHTupleMaker :: selectJets ()
{
  
  ////////////////////////////////////
  ////Define the Forward Jets for veto
  ///////////////////////////////////
  fillFWJets();  
  
  
  // //Veto events with forward jets
  // if( m_tuple->vh_njetfw > 0) return EL::StatusCode::FAILURE;
  // incrementCounter("eventCounter_FWJetVeto");
  // This cut is not applied for 2-lep channel

  ////////////////////////
  ///Define the Central Jets for Veto
  ////////////////////////
  fillCentralJetsVeto();

  // //remove events with > 3 central Jets
  // if( m_tuple->vh_njetveto > 3 ) return EL::StatusCode::FAILURE;
  // incrementCounter("eventCounter_More3CentralJets");
  // This cut is not applied for 2-lep channel
  

  /////////////////////////
  ///Define the Central Jets
  // these have a different pT requirement
  ////////////////////////
  std::vector<int> centralJetList=fillCentralJets();

  //remove events with < 1 central Jets //Just to check single Jet efficiency
  if( m_tuple->vh_njet < 1 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_AtLeast1CentralJet");

  //remove events with < 2 central Jets
  if( m_tuple->vh_njet < 2 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_AtLeast2CentralJets");

  //pT cut on first Jet (jets are pT ordered)
  if( m_tuple->jet_pt[m_tuple->vh_jet1] < 45000 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoHighPtJet");

  
  /////////////////////
  ///Apply b-tagging 
  ////////////////////
  fillBJets();


  if( bTagOption==0 ){ 
    jet1 = centralJetList[0];
    jet2 = centralJetList[1];
  }
  

  if( bTagOption==1 ){    
    //remove events with no b-Jets
    if( m_tuple->vh_nbjet == 0 ) return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_0BJets");

    //remove events with > 2 b-Jets
    if( m_tuple->vh_nbjet > 2 ) return EL::StatusCode::FAILURE;
    incrementCounter("eventCounter_More2BJets");

    ////////////////////////
    //identify the two signal jets 
    //////////////////////
    jet1=-1;
    jet2=-1;
    if( m_tuple->vh_njet  == 2 ){
      jet1 = centralJetList[0];
      jet2 = centralJetList[1];
    }else if( m_tuple->vh_nbjet == 2 ){    
      jet1 = m_tuple->vh_bjet1;
      jet2 = m_tuple->vh_bjet2; 
    } else if( m_tuple->vh_njet  > 2 && m_tuple->vh_nbjet == 1 ){
      ///loop over central jets and find the highest pT which is not jet1
      ///Note this list already has central jet requirement
      jet1 = m_tuple->vh_bjet1;
      float maxpt=0;
      for (std::vector<int>::iterator it = centralJetList.begin() ; it != centralJetList.end(); ++it){
	if((*it)!=jet1 && m_tuple->jet_pt[*it] > maxpt){
	  jet2 = *it ;
	  maxpt = m_tuple->jet_pt[*it];
	}
      }        
    }

  }



  if(jet1==-1 || jet2==-1){
    cout<<"selectJets() : Something is wrong with the b-jets logic"<<endl;
    exit(0);
  }
  
  //must be pT ordered for later
  if(m_tuple->jet_pt[jet2] > m_tuple->jet_pt[jet1] ){
    int tmp=jet1;
    jet1 = jet2;
    jet2 = tmp;
  }


  return EL::StatusCode::SUCCESS;
}



////////////////////////////////////////////////////////////////
/////////Fat jet analysis
////////////////////////////////////////////////////////
EL::StatusCode  VHTupleMaker :: selectHiggs ()
{
 
  ////////////////////////////////////
  ////Define the Forward Jets for veto
  ////////////////////////////////////
  fillFWJets();  
  
  /////////////////////////
  ///select fat jets with basic selections
  ////////////////////////
  std::vector<int> fatJetList = fillFatJets();

  //remove events with no fat jets
  if( m_tuple->vh_nfatjet < 1 ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_NoFatJet");

  //select the Higgs candidate
  higgsJet  =  m_tuple->vh_fatjet1; //select the leading jet as the Higgs candidate

  //pT cut 
  if( m_tuple->fatjet_pt[higgsJet] < fatJetPtCut ) return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_FatJetPt");

  //bad jet veto
  if( ! passBadJetVeto())  return EL::StatusCode::FAILURE;
  incrementCounter("eventCounter_BadJetVeto");
  
  ///////////////////////////////////
  /////////NO MORE CUTS BELOW
  ////////////////////////////////
  
  ///fill info about additional b-jets (higgsJet must have been filled)
  fillAddBJets();
 

  ////Fill MC information
  if( m_tuple->eve_isMC  ){
    
    ////////////////////////////////////////////////
    ///  fill truth info for the Higgs candiate
    //////////////////////////////////////////////
    m_tuple->vh_h_truthdR[higgsJet] = 100.;
    if(m_tuple->truth_RToVH_H>=0) 
      m_tuple->vh_h_truthdR[higgsJet] = deltaR(m_tuple->fatjet_eta[higgsJet],
					       m_tuple->fatjet_phi[higgsJet],
					       m_tuple->truth_eta[m_tuple->truth_RToVH_H],
					       m_tuple->truth_phi[m_tuple->truth_RToVH_H]);
  

    ///b-tagging weight, Note this works for fat jets with only ONE track jet as well
    std::vector<const xAOD::Jet*> trackJetsForBTagSF;
    if(m_tuple->fatjet_ltrkjet[higgsJet] != -1)
      trackJetsForBTagSF.push_back(m_tuple->trkjet_address[m_tuple->fatjet_ltrkjet[higgsJet]]);
    if(m_tuple->fatjet_sltrkjet[higgsJet] != -1)
      trackJetsForBTagSF.push_back(m_tuple->trkjet_address[m_tuple->fatjet_sltrkjet[higgsJet]]);

 
    if( m_debug){
      PrintTrkJets();//print full list of track jets in the event
      std::cout<<" Track jets for BTag SF : "<<std::endl;
      for(unsigned int i=0;i<trackJetsForBTagSF.size();i++){
	cout<<"  pt="<<trackJetsForBTagSF[i]->pt()<<endl;
      }
    }

    //Need to add track jets outside fat jet if using additional btag veto
    //given a fat jet with index "higgsJet" identify the outer track jets
    for(int idx=0;idx<m_tuple->ntrkjet;idx++){
      if( ! passTrackJet(idx) ) continue;
      if( this->trkJetOverlap(idx) ) continue;
      if( trkJetLinkedHiggs(idx) ) continue;
      trackJetsForBTagSF.push_back(m_tuple->trkjet_address[idx]);
    }

    if( m_debug){
      std::cout<<" Track jets for BTag SF with additional b-tag veto: "<<std::endl;
      for(unsigned int i=0;i<trackJetsForBTagSF.size();i++){
	cout<<"  pt="<<trackJetsForBTagSF[i]->pt()<<endl;
      }
    }


    std::map<std::string, float> btagEffSFs = m_bTagTool->computeEventWeight(trackJetsForBTagSF);
    m_tuple->eve_btag_w = btagEffSFs["Nominal"];   
    m_tuple->eve_btag_w_FT_EFF_Eigen_B_0__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_B_0__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_B_1__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_B_1__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_B_2__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_B_2__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_0__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_0__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_1__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_1__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_2__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_2__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_3__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_C_3__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_0__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_0__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_1__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_1__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_2__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_2__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_3__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_3__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_4__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_Eigen_Light_4__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_extrapolation_from_charm__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_extrapolation_from_charm__1up = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_extrapolation__1down = m_tuple->eve_btag_w;
    m_tuple->eve_btag_w_FT_EFF_extrapolation__1up = m_tuple->eve_btag_w;

    if( btagEffSFs["FT_EFF_Eigen_B_0__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_B_0__1down = btagEffSFs["FT_EFF_Eigen_B_0__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_B_0__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_B_0__1up = btagEffSFs["FT_EFF_Eigen_B_0__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_B_1__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_B_1__1down = btagEffSFs["FT_EFF_Eigen_B_1__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_B_1__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_B_1__1up = btagEffSFs["FT_EFF_Eigen_B_1__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_B_2__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_B_2__1down = btagEffSFs["FT_EFF_Eigen_B_2__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_B_2__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_B_2__1up = btagEffSFs["FT_EFF_Eigen_B_2__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_0__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_0__1down = btagEffSFs["FT_EFF_Eigen_C_0__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_0__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_0__1up = btagEffSFs["FT_EFF_Eigen_C_0__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_1__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_1__1down = btagEffSFs["FT_EFF_Eigen_C_1__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_1__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_1__1up = btagEffSFs["FT_EFF_Eigen_C_1__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_2__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_2__1down = btagEffSFs["FT_EFF_Eigen_C_2__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_2__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_2__1up = btagEffSFs["FT_EFF_Eigen_C_2__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_3__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_3__1down = btagEffSFs["FT_EFF_Eigen_C_3__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_C_3__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_C_3__1up = btagEffSFs["FT_EFF_Eigen_C_3__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_0__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_0__1down = btagEffSFs["FT_EFF_Eigen_Light_0__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_0__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_0__1up = btagEffSFs["FT_EFF_Eigen_Light_0__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_1__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_1__1down = btagEffSFs["FT_EFF_Eigen_Light_1__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_1__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_1__1up = btagEffSFs["FT_EFF_Eigen_Light_1__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_2__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_2__1down = btagEffSFs["FT_EFF_Eigen_Light_2__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_2__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_2__1up = btagEffSFs["FT_EFF_Eigen_Light_2__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_3__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_3__1down = btagEffSFs["FT_EFF_Eigen_Light_3__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_3__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_3__1up = btagEffSFs["FT_EFF_Eigen_Light_3__1up"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_4__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_4__1down = btagEffSFs["FT_EFF_Eigen_Light_4__1down"] ;
    if( btagEffSFs["FT_EFF_Eigen_Light_4__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_Eigen_Light_4__1up = btagEffSFs["FT_EFF_Eigen_Light_4__1up"] ;
    if( btagEffSFs["FT_EFF_extrapolation_from_charm__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_extrapolation_from_charm__1down = btagEffSFs["FT_EFF_extrapolation_from_charm__1down"] ;
    if( btagEffSFs["FT_EFF_extrapolation_from_charm__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_extrapolation_from_charm__1up = btagEffSFs["FT_EFF_extrapolation_from_charm__1up"] ;
    if( btagEffSFs["FT_EFF_extrapolation__1down"] > 0. )  m_tuple->eve_btag_w_FT_EFF_extrapolation__1down = btagEffSFs["FT_EFF_extrapolation__1down"] ;
    if( btagEffSFs["FT_EFF_extrapolation__1up"] > 0. )  m_tuple->eve_btag_w_FT_EFF_extrapolation__1up = btagEffSFs["FT_EFF_extrapolation__1up"] ;

    
    // //only fill weight variations in "Nominal" tuple
    // if (m_currentSyst.compare("Nominal") ==0 && m_bTagTool->doWeightVar()) {
    //   for (auto effSF : btagEffSFs) {
    // 	Info("computeBTagSFWeight", "Relative weight for :%s: = %f", effSF.first.c_str(), effSF.second/btagEffSFs["Nominal"]);
    //   }
    // }
    // cout<<"m_tuple->eve_btag_w_FT_EFF_Eigen_B_0__1down: = "<<m_tuple->eve_btag_w_FT_EFF_Eigen_B_0__1down/m_tuple->eve_btag_w<<endl;

  }

 
  return EL::StatusCode::SUCCESS;
}


std::vector<int>  VHTupleMaker :: fillFatJets ()
{
  std::vector<int> list;
  m_tuple->vh_nfatjet=0;
  m_tuple->vh_fatjet1=-1;//leading fat jet after basic selections
  float maxpt=0.;
  for(int i=0;i<m_tuple->nfatjet;i++){
    if(!passFatJet(i)) continue;
    if(this->fatJetOverlap(i)) continue;//overlap with leptons
    
    m_tuple->vh_nfatjet++;
    list.push_back(i);

    if(m_tuple->fatjet_pt[i] > maxpt){
      m_tuple->vh_fatjet1 = i ;
      maxpt = m_tuple->fatjet_pt[i];
    }

  }


  return list;
}


void VHTupleMaker :: fillAddBJets()
{

  //additional b-tagged track jets
  m_tuple->vh_naddtrkbjet = 0;
  for(int i=0;i<m_tuple->ntrkjet;i++){
    if( ! passTrackJet(i) ) continue;
    if( trkJetLinkedHiggs(i) ) continue;
    if( m_tuple->trkjet_btag[i] == 0 ) continue; 
    m_tuple->vh_naddtrkbjet++;
  }

}


bool VHTupleMaker :: jetOverlapHiggs (int jetindex){
  if(deltaR(m_tuple->jet_eta[jetindex],m_tuple->jet_phi[jetindex],
	    m_tuple->fatjet_eta[higgsJet],m_tuple->fatjet_phi[higgsJet]) < 1.4) return 1;
  return 0; 
}

bool VHTupleMaker :: trkJetLinkedHiggs(int jetindex){

  std::vector<const xAOD::Jet*> linkedTrkJets;
  if(! m_tuple->fatjet_address[higgsJet]->getAssociatedObjects("GhostAntiKt2TrackJet", linkedTrkJets) )
    Error("trkJetLinkedHiggs", "Failed to retrieve GhostAntiKt2TrackJet");
  
  for(unsigned int t=0;t<linkedTrkJets.size();t++)
    if( m_tuple->trkjet_address[jetindex] == linkedTrkJets[t] ) 
      return 1;
  
  return 0;
}
