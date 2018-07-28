// packages includes
#include "TupleMaker_VHbbResonance/RecoTupleMaker.h"
#include <TRegexp.h>

ClassImp(RecoTupleMaker)

RecoTupleMaker :: RecoTupleMaker () {}
RecoTupleMaker :: ~RecoTupleMaker () {}

RecoTupleMaker :: RecoTupleMaker (std::string configPath) :
  BaseTupleMaker(configPath),
  m_tuple(0),
  m_bTagTool(nullptr),
  m_electronsNameIn("none"),
  m_muonsNameIn("none"),
  m_jetsNameIn("none"),
  m_fatJetsNameIn("none"),
  m_missingETIn("none")
{

}


EL::StatusCode RecoTupleMaker :: initialize ()
{
  if( m_debug) std::cout<<" RecoTupleMaker :: initialize "<<std::endl;

  if(!m_tuple){//tuple could have been defined upstream
    m_tuple = new RecoTuple();   
  }  
  BaseTupleMaker::m_tuple = m_tuple;
  

  //ConfigStore* config = ConfigStore::createStore(configPath);

  m_electronsNameIn   = config->get<std::string>("tuple.ElectronsIn");
  m_muonsNameIn       = config->get<std::string>("tuple.MuonsIn");
  m_jetsNameIn        = config->get<std::string>("tuple.JetsIn");
  m_missingETIn       = config->get<std::string>("tuple.missingETIn");
  m_fatJetsNameIn     = config->get<std::string>("tuple.FatJetsIn");
  m_trkJetsNameIn     = config->get<std::string>("tuple.TrackJetsIn");  

    
  //Fill systematic variations
  config->getif< std::vector< std::string > >("tuple.ElectronsSysts",m_electronsSysts);
  for(unsigned int i=0;i<m_electronsSysts.size();i++)
    m_systNames.push_back(m_electronsSysts[i]);

  config->getif< std::vector< std::string > >("tuple.MuonsSysts", m_muonsSysts);
  for(unsigned int i=0;i<m_muonsSysts.size();i++)
    m_systNames.push_back(m_muonsSysts[i]);

  config->getif< std::vector< std::string > >("tuple.JetsSysts", m_jetsSysts);
  for(unsigned int i=0;i<m_jetsSysts.size();i++)
    m_systNames.push_back(m_jetsSysts[i]);

  config->getif< std::vector< std::string > >("tuple.FatJetsSysts", m_fatJetsSysts);
  for(unsigned int i=0;i<m_fatJetsSysts.size();i++)
    m_systNames.push_back(m_fatJetsSysts[i]);

  config->getif< std::vector< std::string > >("tuple.TrackJetsSysts", m_trkJetsSysts);
  for(unsigned int i=0;i<m_trkJetsSysts.size();i++)
    m_systNames.push_back(m_trkJetsSysts[i]);

  config->getif< std::vector< std::string > >("tuple.missingETSysts", m_missingETSysts);
  for(unsigned int i=0;i<m_missingETSysts.size();i++)
    m_systNames.push_back(m_missingETSysts[i]);




  //b-tagging
  //---------
  std::vector<std::string> bTagToolConfigs;
  config->getif<std::vector<std::string> >("bTagToolConfigs", bTagToolConfigs);
  if (bTagToolConfigs.size() >= 4) {
    if( m_debug) std::cout<<" RecoTupleMaker :: initialize : bTagToolConfigs.size() >= 4 "<<std::endl;

    BTaggingTool::Config_t args{
      { "TaggerName", bTagToolConfigs[0]},
	{ "OperatingPoint", bTagToolConfigs[1]},
	  { "JetAuthor", bTagToolConfigs[2]},
	    { "Scheme", bTagToolConfigs[3]},
	      { "rScheme", "Medium"}
    };
    //if (!m_bTagTool) 
    m_bTagTool = new BTaggingTool();
    //EL_CHECK("initializeTools()",);
    m_bTagTool->initialize(args);
    m_bTagTool->setWeightVar(true);
  } else {
    Warning("initializeTools()",
	    "Could not initialize BTaggingTool due to invalid bTagToolConfigs in config!");
    return EL::StatusCode::FAILURE;
  }
	
  if( m_debug) std::cout<<" RecoTupleMaker :: initialize : Done  "<<std::endl;
  return BaseTupleMaker::initialize();
}

EL::StatusCode RecoTupleMaker :: processEvent() 
{

  //switch the containers
  TString electrons(m_electronsNameIn.c_str());
  for(unsigned int i=0;i<m_electronsSysts.size();i++)
    if(m_currentSyst.compare(m_electronsSysts[i]) == 0)
      electrons(TRegexp("Nominal")) = m_currentSyst.c_str() ;

  TString muons(m_muonsNameIn.c_str());
  for(unsigned int i=0;i<m_muonsSysts.size();i++)
    if(m_currentSyst.compare(m_muonsSysts[i]) == 0)
      muons(TRegexp("Nominal")) = m_currentSyst.c_str() ;

  TString jets(m_jetsNameIn.c_str());
  for(unsigned int i=0;i<m_jetsSysts.size();i++)
    if(m_currentSyst.compare(m_jetsSysts[i]) == 0)
      jets(TRegexp("Nominal")) = m_currentSyst.c_str() ;

  TString fatjets(m_fatJetsNameIn.c_str());
  for(unsigned int i=0;i<m_fatJetsSysts.size();i++)
    if(m_currentSyst.compare(m_fatJetsSysts[i]) == 0)
      fatjets(TRegexp("Nominal")) = m_currentSyst.c_str() ;

  TString trkjets(m_trkJetsNameIn.c_str());
  for(unsigned int i=0;i<m_trkJetsSysts.size();i++)
    if(m_currentSyst.compare(m_trkJetsSysts[i]) == 0)
      trkjets(TRegexp("Nominal")) = m_currentSyst.c_str() ;

  TString met(m_missingETIn.c_str());
  for(unsigned int i=0;i<m_missingETSysts.size();i++)
    if(m_currentSyst.compare(m_missingETSysts[i]) == 0)
      met(TRegexp("Nominal")) = m_currentSyst.c_str() ;


  //std::cout<<"Config: "<<electrons<<" "<<muons <<" "<<jets <<" "<<fatjets <<" "<<trkjets <<" "<<met<<std::endl;


  ////////////////////
  if( !m_event->retrieve( m_Electrons, electrons.Data() ).isSuccess() ) {
    Error("execute", "Failed to retrieve %s contained. Exiting.", electrons.Data());
    return EL::StatusCode::FAILURE;
  }
  
  if( !m_event->retrieve( m_Muons, muons.Data() ).isSuccess() ) {
    Error("execute", "Failed to retrieve %s contained. Exiting.", muons.Data());
    return EL::StatusCode::FAILURE;
  }

  if( !m_event->retrieve( m_Jets, jets.Data() ).isSuccess() ) {
    Error("execute", "Failed to retrieve %s contained. Exiting.", jets.Data());
    return EL::StatusCode::FAILURE;
  }
 
  if( !m_event->retrieve( m_FatJets, fatjets.Data() ).isSuccess() ) {
    Error("execute", "Failed to retrieve %s contained. Exiting.", fatjets.Data());
    return EL::StatusCode::FAILURE;
  }

  if( !m_event->retrieve( m_TrackJets, trkjets.Data() ).isSuccess() ) {
    Error("execute", "Failed to retrieve %s contained. Exiting.", trkjets.Data());
    return EL::StatusCode::FAILURE;
  }
  
  if (!m_event->retrieve(m_MissingET, met.Data() ).isSuccess()) {
    Error("execute()", "Failed to retrieve missing ET collection  %s. Exiting.", met.Data()); 
    return EL::StatusCode::FAILURE;
  }                                                          


  /////////////////////////////////
  if( BaseTupleMaker::processEvent() == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;

  if( processRecoEventInfo() == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;

  if( processParticles() == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;


  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RecoTupleMaker :: processParticles ()
{
  
  if(processElectrons() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  if(processMuons() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  if(processJets() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  if(processTrackJets() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;//order with fat jets matters
  if(processFatJets() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  if(processMET() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  if(processEE() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  if(processMM() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  if(processEM() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  //if(processJJ() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  //if(processMUMET() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;
  //if(processEMET() == EL::StatusCode::FAILURE) return EL::StatusCode::FAILURE;

 
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode RecoTupleMaker :: processRecoEventInfo() 
{
 
  //m_tuple->nvtx = Props::NVtx2Trks.get( m_eventInfo);//disabled during move to AB2.3.41 (to be compatible with DB00-06
  m_tuple->vtx_z[0] = Props::ZPV.get( m_eventInfo);

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode RecoTupleMaker :: processElectrons() 
{

  m_tuple->nele=0;

  // loop over electrons
  int i = 0;
  for( xAOD::ElectronContainer::const_iterator iter = m_Electrons->begin();
       iter!= m_Electrons->end(); ++iter) {
    if( i < RecoTuple::MAXELECTRONS) {
      
      if( (*iter)->pt() < 7000 ) continue;

      m_tuple->ele_charge[i] = 0;
      m_tuple->ele_m[i] = (*iter)->m();
      m_tuple->ele_pt[i] = (*iter)->pt();
      m_tuple->ele_phi[i] = (*iter)->phi();
      m_tuple->ele_eta[i] = (*iter)->eta();
      m_tuple->ele_E[i] = (*iter)->e();
      m_tuple->ele_cluseta[i] = DBProps::ElClusterEta.get(*iter); 
      m_tuple->ele_trkpt[i] = Props::innerTrackPt.get(*iter);

      m_tuple->ele_z0[i] = Props::z0.get(*iter);
      m_tuple->ele_z0sintheta[i] = Props::z0sinTheta.get(*iter); 
      m_tuple->ele_d0[i] = Props::d0.get(*iter);
      m_tuple->ele_d0sig[i] = Props::d0sigBL.get(*iter);
      m_tuple->ele_isFromPV[i] = DBProps::isFromPV.get(*iter);
      m_tuple->ele_sharedtrk[i] = Props::hasSharedTrack.get(*iter);

      m_tuple->ele_ptiso[i] = (*iter)->auxdata<float>("ptvarcone20");
      m_tuple->ele_etiso[i] = (*iter)->auxdata<float>("topoetcone20");
      m_tuple->ele_isoWP1[i] = m_tuple->ele_ptiso[i] / m_tuple->ele_pt[i] < 0.1;
      m_tuple->ele_isoWP2[i] = Props::isLooseTrackOnlyIso.get(*iter);
      m_tuple->ele_isoWP3[i] = Props::isLooseIso.get(*iter);
      m_tuple->ele_isoWP4[i] = Props::isTightIso.get(*iter);
      m_tuple->ele_isoWP5[i] = Props::isGradientIso.get(*iter);
      m_tuple->ele_isoWP6[i] = Props::isGradientLooseIso.get(*iter);
      
      m_tuple->ele_idVHl[i] = DBProps::isVVCxAODMakerElectron.get(*iter)==1; 
      m_tuple->ele_isLooseIDElectron[i] = DBProps::isLooseIDElectron.get(*iter);
      m_tuple->ele_isMediumIDElectron[i] = DBProps::isMediumIDElectron.get(*iter);
      m_tuple->ele_isTightIDElectron[i] = DBProps::isTightIDElectron.get(*iter);

      m_tuple->ele_matchHLT_e24_lhmedium_L1EM18VH[i] = Props::matchHLT_e24_lhmedium_L1EM18VH.get(*iter);
      m_tuple->ele_matchHLT_e24_lhmedium_L1EM20VH[i] = Props::matchHLT_e24_lhmedium_L1EM20VH.get(*iter);
      m_tuple->ele_matchHLT_e26_lhtight_iloose[i] =    Props::matchHLT_e26_lhtight_iloose.get(*iter);
      m_tuple->ele_matchHLT_e60_lhmedium[i] =          Props::matchHLT_e60_lhmedium.get(*iter);
      m_tuple->ele_matchHLT_e120_lhloose[i] =          Props::matchHLT_e120_lhloose.get(*iter);

      ++i;
    }else{
      Info( "processElectrons()" , " Trying to build too many electrons, max is %i" , RecoTuple::MAXELECTRONS);
    }
  }
  m_tuple->nele=i;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RecoTupleMaker :: processMuons() 
{


  m_tuple->nmuo=0;
  
  // loop over muons
  int i = 0;
  for( xAOD::MuonContainer::const_iterator iter = m_Muons->begin();
       iter!= m_Muons->end(); ++iter) {
    if( i < RecoTuple::MAXMUONS) {
      
      if( (*iter)->pt() < 7000 ) continue;
      

      m_tuple->muo_charge[i] = 0;//(*iter)->charge();
      m_tuple->muo_m[i] = (*iter)->m();
      m_tuple->muo_pt[i] = (*iter)->pt();
      m_tuple->muo_phi[i] = (*iter)->phi();
      m_tuple->muo_eta[i] = (*iter)->eta();
      m_tuple->muo_E[i] = (*iter)->e();
      m_tuple->muo_trkpt[i] = Props::innerTrackPt.get(*iter);      

      m_tuple->muo_z0[i] = Props::z0.get(*iter); 
      m_tuple->muo_z0sintheta[i] = Props::z0sinTheta.get(*iter); 
      m_tuple->muo_d0[i] = Props::d0.get(*iter); 
      m_tuple->muo_d0sig[i] = Props::d0sigBL.get(*iter); 
      m_tuple->muo_isFromPV[i] = DBProps::isFromPV.get(*iter);

      //(*iter)->isolation(m_tuple->muo_ptiso[i], xAOD::Iso::ptcone20);
      //(*iter)->isolation(m_tuple->muo_etiso[i], xAOD::Iso::topoetcone30);
      m_tuple->muo_ptiso[i] = (*iter)->auxdata<float>("ptvarcone30");
      m_tuple->muo_etiso[i] = (*iter)->auxdata<float>("topoetcone20");
      m_tuple->muo_isoWP1[i] =  m_tuple->muo_ptiso[i] / m_tuple->muo_pt[i] < 0.1;
      m_tuple->muo_isoWP2[i] = Props::isLooseTrackOnlyIso.get(*iter);
      m_tuple->muo_isoWP3[i] = Props::isLooseIso.get(*iter);
      m_tuple->muo_isoWP4[i] = Props::isTightIso.get(*iter);
      m_tuple->muo_isoWP5[i] = Props::isGradientIso.get(*iter);
      m_tuple->muo_isoWP6[i] = Props::isGradientLooseIso.get(*iter);

      m_tuple->muo_idVHl[i] = DBProps::isVVCxAODMakerMuon.get(*iter)==1; 
      m_tuple->muo_isLooseIDMuon[i] = DBProps::isLooseIDMuon.get(*iter);
      m_tuple->muo_isMediumIDMuon[i] = DBProps::isMediumIDMuon.get(*iter);
      m_tuple->muo_isTightIDMuon[i] = DBProps::isTightIDMuon.get(*iter);

      //HLT_mu26_imedium
      m_tuple->muo_matchHLT_mu20_iloose_L1MU15[i] = Props::matchHLT_mu20_iloose_L1MU15.get(*iter); 
      m_tuple->muo_matchHLT_mu26_imedium[i] = Props::matchHLT_mu26_imedium.get(*iter); 
      m_tuple->muo_matchHLT_mu50[i] = Props::matchHLT_mu50.get(*iter); 

      ++i;
    }else{
      Info( "processMuons()" , " Trying to build too many muons, max is %i" , RecoTuple::MAXMUONS);
    }
  }
  m_tuple->nmuo=i;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RecoTupleMaker :: processJets() 
{

  m_tuple->njet=0;

  // loop over jets
  int i = 0;
  for( xAOD::JetContainer::const_iterator iter = m_Jets->begin();
       iter!= m_Jets->end(); ++iter) {
    if( i < RecoTuple::MAXJETS) {

      if( (*iter)->pt() < 15000 ) continue;


      m_tuple->jet_charge[i] = 0;//currently not set
      m_tuple->jet_m[i] = (*iter)->m();
      m_tuple->jet_pt[i] = (*iter)->pt();
      m_tuple->jet_phi[i] = (*iter)->phi();
      m_tuple->jet_eta[i] = (*iter)->eta();
      m_tuple->jet_E[i] = (*iter)->e();
      m_tuple->jet_jvf[i] = Props::jvf0.get(*iter);
      m_tuple->jet_jvt[i] = Props::Jvt.get(*iter);
      m_tuple->jet_mv2c00[i] = Props::MV2c00.get(*iter);
      m_tuple->jet_mv2c10[i] = Props::MV2c10.get(*iter);
      m_tuple->jet_mv2c20[i] = Props::MV2c20.get(*iter);
      m_tuple->jet_sv1ip3d[i] = Props::SV1_IP3D.get(*iter);
      m_tuple->jet_mvb[i] = Props::MVb.get(*iter);
      m_tuple->jet_mv1[i] = Props::MV1.get(*iter);
      m_tuple->jet_good[i] = Props::goodJet.get(*iter);
      m_tuple->jet_ntrk[i] = Props::NumTrkPt1000PV.get(*iter);
      m_tuple->jet_NumTrkPt500PV[i] = Props::NumTrkPt500PV.get(*iter);
      m_tuple->jet_SumPtTrkPt500PV[i] = Props::SumPtTrkPt500PV.get(*iter);

      m_tuple->jet_btag[i] = passBTag(m_tuple->jet_mv2c20[i],2) ? 1 : 0;

      ++i;
    }else{
      Info( "processJets()" , " Trying to build too many jets, max is %i" , RecoTuple::MAXJETS);
    }
  }

  m_tuple->njet=i;
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode RecoTupleMaker :: processFatJets() 
{
  m_tuple->nfatjet=0;


  // loop over jets
  int i = 0;
  for( xAOD::JetContainer::const_iterator iter = m_FatJets->begin(); iter!= m_FatJets->end(); ++iter) {
    if( i < RecoTuple::MAXJETS) {

      if( (*iter)->pt() < 100000 ) continue;

      m_tuple->fatjet_address[i]=(*iter);//address needed by derived classes to check links

      /////////Nominal, no muon correction
      m_tuple->fatjet_m[i] = (*iter)->m();
      m_tuple->fatjet_pt[i] = (*iter)->pt();
      m_tuple->fatjet_phi[i] = (*iter)->phi();
      m_tuple->fatjet_eta[i] = (*iter)->eta();

      //"OneMu" correction not saved in last 00-06 production 
      //For now will use Xbb correction which only applies to jets with >=2 track jets
      m_tuple->fatjet_m_cor[i]   =  (*iter)->auxdata< float >("OneMu_m");
      m_tuple->fatjet_pt_cor[i]  =  (*iter)->auxdata< float >("OneMu_pt");
      m_tuple->fatjet_phi_cor[i] =  (*iter)->auxdata< float >("OneMu_phi");
      m_tuple->fatjet_eta_cor[i] =  (*iter)->auxdata< float >("OneMu_eta");

      // ///Muon correction based on Xbb tagger, does not work for jets with ony one track jet
      // m_tuple->fatjet_m_cor[i]   =  (*iter)->auxdata< float >("XbbCorrectedJetP4_m");//corrected_jet.M();
      // m_tuple->fatjet_pt_cor[i]  =  (*iter)->auxdata< float >("XbbCorrectedJetP4_pt");//corrected_jet.Pt();
      // m_tuple->fatjet_phi_cor[i] =  (*iter)->auxdata< float >("XbbCorrectedJetP4_phi");//corrected_jet.Phi();
      // m_tuple->fatjet_eta_cor[i] =  (*iter)->auxdata< float >("XbbCorrectedJetP4_eta");//corrected_jet.Eta();

      //btagging based on Xbb tagger
      m_tuple->fatjet_XbbL2b[i]   = Props::isXbbLoose2b.get(*iter); 
      m_tuple->fatjet_XbbL2bmH[i] = Props::XbbLoose2bmH.get(*iter); 
      m_tuple->fatjet_XbbnTag[i]  = Props::XbbnTag.get(*iter); 

      ///get the list of linked track jets
      std::vector<const xAOD::Jet*> linkedTrkJets;
      if(!(*iter)->getAssociatedObjects("GhostAntiKt2TrackJet", linkedTrkJets)){
      	Error("processFatJets", "Failed to retrieve GhostAntiKt2TrackJet");
      }
      //for(unsigned int t=0;t<linkedTrkJets.size();t++)
      ///std::cout<<i<<" : "<<t<<" "<<linkedTrkJets[t]->pt()<<std::endl; -->Crashes (missing links?)

      // loop over track jets, use the linked ones only.
      int ntjet=0;
      int ntjetb=0;
      int ltrkjet=-1;

      //this relies on track jet block being filled before
      for(int idx=0;idx<m_tuple->ntrkjet;idx++){
	//check quality
	if( !passTrackJet(idx) ) continue;

	//check link
	bool linked=0;
	for(unsigned int t=0;t<linkedTrkJets.size();t++)
	  if( linkedTrkJets[t] == m_tuple->trkjet_address[idx] ) 
	    linked = 1;
	if(!linked) continue;

	///count track jets 
	ntjet++;

	//check btag's
	if( m_tuple->trkjet_btag[idx] == 1)
	  ntjetb++;

	
	//find leading track jet
	if(ltrkjet ==-1 
	   || m_tuple->trkjet_pt[idx] > m_tuple->trkjet_pt[ltrkjet]
	   )
	  ltrkjet = idx; 

      }


      //find subleading trk jet
      int sltrkjet=-1;
      for(int idx=0;idx<m_tuple->ntrkjet;idx++){
	//check quality
	if( !passTrackJet(idx) ) continue;

	//check link
	bool linked=0;
	for(unsigned int t=0; t<linkedTrkJets.size(); t++)
	  if( linkedTrkJets[t] == m_tuple->trkjet_address[idx] ) 
	    linked = 1;
	if(!linked) continue;

	//do not consider the leading jet
	if( idx == ltrkjet ) continue;
	
	//find leading track jet
	if( sltrkjet ==-1 
	    || m_tuple->trkjet_pt[idx] > m_tuple->trkjet_pt[sltrkjet]
	    )
	  sltrkjet = idx; 
	
      }


      m_tuple->fatjet_ntrkjet[i] = ntjet; 
      m_tuple->fatjet_ntrkjetb[i] = ntjetb; 
      m_tuple->fatjet_ltrkjet[i] = ltrkjet;
      m_tuple->fatjet_sltrkjet[i] = sltrkjet;


      ++i;
    }else{
      Info( "processFatjets()" , " Trying to build too many jets, max is %i" , RecoTuple::MAXJETS);
    }
  }

  m_tuple->nfatjet=i;
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode RecoTupleMaker :: processTrackJets() 
{
  m_tuple->ntrkjet=0;


  // loop over jets
  int i = 0;
  for( xAOD::JetContainer::const_iterator iter = m_TrackJets->begin(); iter!= m_TrackJets->end(); ++iter) {
    if( i < RecoTuple::MAXJETS) {
      m_tuple->trkjet_address[i] = (*iter);
	
      m_tuple->trkjet_m[i] = (*iter)->m();
      m_tuple->trkjet_pt[i] = (*iter)->pt();
      m_tuple->trkjet_phi[i] = (*iter)->phi();
      m_tuple->trkjet_eta[i] = (*iter)->eta();
      m_tuple->trkjet_numConstituents[i] = Props::TrkJetnumConstituents.get(*iter);

      m_tuple->trkjet_sv1ip3d[i]=Props::SV1_IP3D.get(*iter);
      m_tuple->trkjet_mv2c00[i]=Props::MV2c00.get(*iter);
      m_tuple->trkjet_mv2c10[i]=Props::MV2c10.get(*iter);
      m_tuple->trkjet_mv2c20[i]=Props::MV2c20.get(*iter);

      m_tuple->trkjet_btag[i] = passBTag(m_tuple->trkjet_mv2c20[i],2) ? 1 : 0;

      ++i;
    }else{
      Info( "processTrackjets()" , " Trying to build too many jets, max is %i" , RecoTuple::MAXJETS);
    }
  }

  m_tuple->ntrkjet=i;
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode RecoTupleMaker :: processMET() 
{
  
  m_tuple->met_sumet = (*m_MissingET).at(0)->sumet();
  m_tuple->met_pt = (*m_MissingET).at(0)->met();
  m_tuple->met_eta = 0.;
  m_tuple->met_phi = (*m_MissingET).at(0)->phi();
  m_tuple->met_px = (*m_MissingET).at(0)->mpx();
  m_tuple->met_py = (*m_MissingET).at(0)->mpy();
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode RecoTupleMaker :: processEE() 
{
  //Note: Muon array must be filled prior
  m_tuple->nee=0;
  for(int e1 = 0; e1 < m_tuple->nele ; e1++){
    for(int e2 = 0; e2 < m_tuple->nele ; e2++){
      if( m_tuple->ele_pt[e1] <= m_tuple->ele_pt[e2] ) continue; //pT order, remove duplicates
      if( !passLooseElectron(e1,2) ) continue;
      if( !passLooseElectron(e2,2) ) continue;
      if( !passElectronOR(e1) ) continue;
      if( !passElectronOR(e2) ) continue;

      if( m_tuple->nee < RecoTuple::MAXEE) {

	//build the di-electron
	TLorentzVector e1p4;
	e1p4.SetPtEtaPhiM(m_tuple->ele_pt[e1],m_tuple->ele_eta[e1],m_tuple->ele_phi[e1],m_tuple->ele_m[e1]);
	TLorentzVector e2p4;
	e2p4.SetPtEtaPhiM(m_tuple->ele_pt[e2],m_tuple->ele_eta[e2],m_tuple->ele_phi[e2],m_tuple->ele_m[e2]);
	TLorentzVector eep4=e1p4 + e2p4;

	m_tuple->ee_charge[m_tuple->nee] = m_tuple->ele_charge[e1] + m_tuple->ele_charge[e2];
	m_tuple->ee_m[m_tuple->nee] = eep4.M();	  
	m_tuple->ee_E[m_tuple->nee] = eep4.E();
	m_tuple->ee_p[m_tuple->nee] = eep4.P();
	m_tuple->ee_pt[m_tuple->nee] = eep4.Pt(); 
	m_tuple->ee_phi[m_tuple->nee] = eep4.Phi();
	m_tuple->ee_eta[m_tuple->nee] = eep4.Eta();
	m_tuple->ee_dR[m_tuple->nee] = e1p4.DeltaR(e2p4);

	m_tuple->ee_leg1[m_tuple->nee] = e1;
	m_tuple->ee_leg2[m_tuple->nee] = e2;

	++m_tuple->nee;
      }else{
	Info( "processEE()" , " Trying to build too many ele pairs, max is %i" , RecoTuple::MAXEE);
      }

    }
  }

  return EL::StatusCode::SUCCESS;
}




EL::StatusCode RecoTupleMaker :: processMM() 
{
  //Note: Muon array must be filled prior
  m_tuple->nmm=0;
  for(int m1 = 0; m1 < m_tuple->nmuo ; m1++){
    for(int m2 = 0; m2 < m_tuple->nmuo ; m2++){
      if( m_tuple->muo_pt[m1] <= m_tuple->muo_pt[m2] ) continue; //pT order, remove duplicates
      if( !passLooseMuon(m1,2) ) continue;
      if( !passLooseMuon(m2,2) ) continue;
      if( !passMuonOR(m1) ) continue;
      if( !passMuonOR(m2) ) continue;

      if( m_tuple->nmm < RecoTuple::MAXMM) {

	//build the dimuon
	TLorentzVector m1p4;
	m1p4.SetPtEtaPhiM(m_tuple->muo_pt[m1],m_tuple->muo_eta[m1],m_tuple->muo_phi[m1],m_tuple->muo_m[m1]);
	TLorentzVector m2p4;
	m2p4.SetPtEtaPhiM(m_tuple->muo_pt[m2],m_tuple->muo_eta[m2],m_tuple->muo_phi[m2],m_tuple->muo_m[m2]);
	TLorentzVector mmp4=m1p4 + m2p4;

	m_tuple->mm_charge[m_tuple->nmm] = m_tuple->muo_charge[m1] + m_tuple->muo_charge[m2];
	m_tuple->mm_m[m_tuple->nmm] = mmp4.M();	  
	m_tuple->mm_E[m_tuple->nmm] = mmp4.E();
	m_tuple->mm_p[m_tuple->nmm] = mmp4.P();
	m_tuple->mm_pt[m_tuple->nmm] = mmp4.Pt(); 
	m_tuple->mm_phi[m_tuple->nmm] = mmp4.Phi();
	m_tuple->mm_eta[m_tuple->nmm] = mmp4.Eta();
	m_tuple->mm_dR[m_tuple->nmm] = m1p4.DeltaR(m2p4);
	  
	m_tuple->mm_leg1[m_tuple->nmm] = m1;
	m_tuple->mm_leg2[m_tuple->nmm] = m2;

	++m_tuple->nmm;
      }else{
	Info( "processMM()" , " Trying to build too many muon pairs, max is %i" , RecoTuple::MAXMM);
      }

    }
  }

  return EL::StatusCode::SUCCESS;
}




EL::StatusCode RecoTupleMaker :: processEM() 
{
  m_tuple->nem=0;
  int i=0;
  for(int l1 = 0; l1 < m_tuple->nele ; l1++){
    for(int l2 = 0; l2 < m_tuple->nmuo ; l2++){
      if( !passLooseElectron(l1,2) ) continue;
      if( !passLooseMuon(l2,2) ) continue;
      if( !passElectronOR(l1) ) continue;
      if( !passMuonOR(l2) ) continue;

      if( i < RecoTuple::MAXEM) {

	  TLorentzVector l1p4;
	  l1p4.SetPtEtaPhiM(m_tuple->ele_pt[l1],m_tuple->ele_eta[l1],m_tuple->ele_phi[l1],m_tuple->ele_m[l1]);
	  TLorentzVector l2p4;
	  l2p4.SetPtEtaPhiM(m_tuple->muo_pt[l2],m_tuple->muo_eta[l2],m_tuple->muo_phi[l2],m_tuple->muo_m[l2]);
	  TLorentzVector emp4=l1p4 + l2p4;

	  m_tuple->em_charge[i] = m_tuple->ele_charge[l1] + m_tuple->muo_charge[l2];
	  m_tuple->em_m[i] = emp4.M();	  
	  m_tuple->em_E[i] = emp4.E();
	  m_tuple->em_p[i] = emp4.P();
	  m_tuple->em_pt[i] = emp4.Pt(); 
	  m_tuple->em_phi[i] = emp4.Phi();
	  m_tuple->em_eta[i] = emp4.Eta();
	  m_tuple->em_dR[i] = l1p4.DeltaR(l2p4);

	  m_tuple->em_leg1[i] = l1;
	  m_tuple->em_leg2[i] = l2;

	  ++i;
	}else{
	  Info( "processEM()" , " Trying to build too many e-mu pairs, max is %i" , RecoTuple::MAXEM);
	}
    }
  }
  m_tuple->nem=i;
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode RecoTupleMaker :: processJJ() 
{

  m_tuple->njj=0;
  // loop over jets1
  int ijj = 0;
  for(int j1 = 0; j1 < m_tuple->njet ; j1++){
    for(int j2 = 0; j2 < m_tuple->njet ; j2++){
      if( ! passCentralJet(j1)) continue;
      if( ! passCentralJet(j2)) continue; 
      if( ! (m_tuple->jet_pt[j1] > m_tuple->jet_pt[j2]) ) continue; 

      if( ijj < RecoTuple::MAXJJ ) {
	  
	//build the dijet
	TLorentzVector j1p4;
	j1p4.SetPtEtaPhiM(m_tuple->jet_pt[j1],m_tuple->jet_eta[j1],m_tuple->jet_phi[j1],m_tuple->jet_m[j1]);
	TLorentzVector j2p4;
	j2p4.SetPtEtaPhiM(m_tuple->jet_pt[j2],m_tuple->jet_eta[j2],m_tuple->jet_phi[j2],m_tuple->jet_m[j2]);
	TLorentzVector jjp4=j1p4 + j2p4;

	m_tuple->jj_charge[ijj] = m_tuple->jet_charge[j1] + m_tuple->jet_charge[j2];
	m_tuple->jj_m[ijj] = jjp4.M();	  
	m_tuple->jj_E[ijj] = jjp4.E();
	m_tuple->jj_p[ijj] = jjp4.P();
	m_tuple->jj_pt[ijj] = jjp4.Pt(); 
	m_tuple->jj_phi[ijj] = jjp4.Phi();
	m_tuple->jj_eta[ijj] = jjp4.Eta();
	m_tuple->jj_dR[ijj] = j1p4.DeltaR(j2p4);

	m_tuple->jj_leg1[ijj] = j1;
	m_tuple->jj_leg2[ijj] = j2;

	++ijj;
      }else{
	Info( "processJJ()" , " Trying to build too many jet pairs, max is %i" , RecoTuple::MAXJJ);
	j1=m_tuple->njet;//break out of the loops
	j2=m_tuple->njet;
      }
    }
  }
  m_tuple->njj=ijj;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RecoTupleMaker :: processEMET() 
{
  //Note: Muon array must be filled prior
  m_tuple->neMet=0;

  int ieMet=0;
  for(int i = 0; i < m_tuple->nele ; i++){
    if( !passLooseElectron(i,2) ) continue;
    if( !passElectronOR(i) ) continue;

    if( ieMet < RecoTuple::MAXEMET) {
      
      //build the dimuon
      TLorentzVector ep4;
      ep4.SetPtEtaPhiM(m_tuple->ele_pt[i],m_tuple->ele_eta[i],m_tuple->ele_phi[i],m_tuple->ele_m[i]);
      TLorentzVector metp4;
      metp4.SetPtEtaPhiM(m_tuple->met_pt,m_tuple->met_eta,m_tuple->met_phi,m_tuple->met_pt);
      
      TLorentzVector eMetp4=ep4 + metp4;
      
      m_tuple->eMet_charge[ieMet] = m_tuple->ele_charge[i];
      m_tuple->eMet_m[ieMet] = eMetp4.M();	  
      m_tuple->eMet_E[ieMet] = eMetp4.E();
      m_tuple->eMet_p[ieMet] = eMetp4.P();
      m_tuple->eMet_pt[ieMet] = eMetp4.Pt(); 
      m_tuple->eMet_phi[ieMet] = eMetp4.Phi();
      m_tuple->eMet_eta[ieMet] = eMetp4.Eta();
      m_tuple->eMet_mT[ieMet] = mT(ep4,metp4);
      m_tuple->eMet_ele[ieMet] = i;
	    
      ++ieMet;
    }else{
      Info( "processEMET()" , " Trying to build too many e-met pairs, max is %i" , RecoTuple::MAXEMET);
    }
  }
  m_tuple->neMet=ieMet;
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode RecoTupleMaker :: processMUMET() 
{
  //Note: Muon array must be filled prior
  m_tuple->nmuMet=0;

  int imuMet=0;
  for(int i = 0; i < m_tuple->nmuo ; i++){
    if( !passLooseMuon(i,2) ) continue;
    if( !passMuonOR(i) ) continue;


    if( imuMet < RecoTuple::MAXMUMET) {
      
      //build the dimuon
      TLorentzVector mup4;
      mup4.SetPtEtaPhiM(m_tuple->muo_pt[i],m_tuple->muo_eta[i],m_tuple->muo_phi[i],m_tuple->muo_m[i]);
      TLorentzVector metp4;
      metp4.SetPtEtaPhiM(m_tuple->met_pt,m_tuple->met_eta,m_tuple->met_phi,m_tuple->met_pt);
      TLorentzVector muMetp4=mup4 + metp4;
      
      m_tuple->muMet_charge[imuMet] = m_tuple->muo_charge[i];
      m_tuple->muMet_m[imuMet] = muMetp4.M();	  
      m_tuple->muMet_E[imuMet] = muMetp4.E();
      m_tuple->muMet_p[imuMet] = muMetp4.P();
      m_tuple->muMet_pt[imuMet] = muMetp4.Pt(); 
      m_tuple->muMet_phi[imuMet] = muMetp4.Phi();
      m_tuple->muMet_eta[imuMet] = muMetp4.Eta();
      m_tuple->muMet_mT[imuMet] = mT(mup4,metp4);
      m_tuple->muMet_mu[imuMet] = i;
      
      ++imuMet;
    }else{
      Info( "processMUMET()" , " Trying to build too many mu-met pairs, max is %i" , RecoTuple::MAXMUMET);
    }
  }
  m_tuple->nmuMet=imuMet;
  return EL::StatusCode::SUCCESS;
}



/////////////////////Tools
void  RecoTupleMaker::PrintElectrons(){
  cout<<"PrintElectrons:"<<endl;
  for(int i = 0; i < m_tuple->nele ; i++)  {cout<<" ";PrintElectron(i);}
}
void  RecoTupleMaker::PrintMuons(){
  cout<<"PrintMuons:"<<endl;
  for(int i = 0; i < m_tuple->nmuo ; i++)   {cout<<" ";PrintMuon(i);}
}
void  RecoTupleMaker::PrintJets(){
  cout<<"PrintJets:"<<endl;
  for(int i = 0; i < m_tuple->njet ; i++) {cout<<" "; PrintJet(i);}
}
void  RecoTupleMaker::PrintFatJets(){
  cout<<"PrintFatJets:"<<endl;
  for(int i = 0; i < m_tuple->nfatjet ; i++)  {cout<<" ";PrintFatJet(i);}
}
void  RecoTupleMaker::PrintTrkJets(){
  cout<<"PrintTrkJets:"<<endl;
  for(int i = 0; i < m_tuple->ntrkjet ; i++)  {cout<<" ";PrintTrkJet(i);}
}
void  RecoTupleMaker::PrintMET(){
  cout<<"MET: sumet: "<<m_tuple->met_pt<<"  phi: "<<m_tuple->met_phi<<endl;
}


void  RecoTupleMaker::PrintEventParticles(){
  PrintElectrons();
  PrintMuons();
  PrintJets();
  PrintFatJets();
  PrintTrkJets();
}


bool RecoTupleMaker :: passElectronOR(int idx){

  ///this is the e-mu OR
  if( m_tuple->ele_sharedtrk[idx] == 1 ) return false;

  //this is the e-jet OR
  for(int i = 0; i < m_tuple->njet ; i++){
    if( ! passCentralJet(i) ) continue;
    if( ! passJetOR(i) ) continue;


    ///if e is within the cone of the above selected jet remove the electron
    if(m_tuple->ele_pt[idx]<=0.) continue;//cannot divide by 0
    if(deltaR(m_tuple->ele_eta[idx],m_tuple->ele_phi[idx],
	      m_tuple->jet_eta[i],m_tuple->jet_phi[i]) 
       < smaller(0.4,0.04+10000/m_tuple->ele_pt[idx])  ) return false;
    
  }

  return true;
}

bool RecoTupleMaker :: passMuJetForOR(int mu, int jet){
  
  if( m_tuple->jet_SumPtTrkPt500PV[jet] <= 0. || m_tuple->muo_pt[mu] <=0. ) return false; //cannot devide by 0.
  
  if ( m_tuple->jet_NumTrkPt500PV[jet] < 3 ) return true;
  
  if ( m_tuple->muo_pt[mu] / m_tuple->jet_SumPtTrkPt500PV[jet] > 0.7 
       && m_tuple->jet_pt[jet] / m_tuple->muo_pt[mu] < 2 ) return true;
  
  return false;
}

bool RecoTupleMaker :: passMuonOR(int idx){

  for(int i = 0; i < m_tuple->njet ; i++){
    if( ! passCentralJet(i) ) continue;
    if( ! passJetOR(i) ) continue;

    ///if muon is within the cone of the above selected jet remove the muon
    if(m_tuple->muo_pt[idx]<=0.) continue; //cannot divide by 0
    if(deltaR(m_tuple->muo_eta[idx],m_tuple->muo_phi[idx],
	      m_tuple->jet_eta[i],m_tuple->jet_phi[i]) 
       < smaller(0.4,0.04+10000/m_tuple->muo_pt[idx])  ) return false;
    
  }

  return true;
}

bool RecoTupleMaker :: passJetOR(int idx){

  ///remove muons
  for(int i = 0; i < m_tuple->nmuo ; i++){
    if( ! passLooseMuon(i,2) ) continue;
    if( ! passMuJetForOR(i,idx) ) continue; 
    
    if(deltaR(m_tuple->muo_eta[i],m_tuple->muo_phi[i],
	      m_tuple->jet_eta[idx],m_tuple->jet_phi[idx]) 
       < 0.2  ) return false;
    
  }
  

  //remove electrons
  for(int i = 0; i < m_tuple->nele ; i++){
    if( !passLooseElectron(i,2) ) continue;
    
    if(deltaR(m_tuple->ele_eta[i],m_tuple->ele_phi[i],
	      m_tuple->jet_eta[idx],m_tuple->jet_phi[idx]) 
       < 0.2  ) return false;
  }
  

  return true;
}


bool RecoTupleMaker :: passLooseMuon(int idx, int iso){
  if( m_tuple->muo_isLooseIDMuon[idx] !=1  ) return false;
  if( m_tuple->muo_isFromPV[idx] != 1 ) return false;
  if( m_tuple->muo_pt[idx] < 10000 ) return false;
  if( fabs(m_tuple->muo_eta[idx]) > 2.7  ) return false;
  //if( ! passMuonOR(idx) ) return false;

  if( iso==1 && m_tuple->muo_isoWP1[idx] != 1 ) return false;  
  if( iso==2 && m_tuple->muo_isoWP2[idx] != 1 ) return false;  
  if( iso==3 && m_tuple->muo_isoWP3[idx] != 1 ) return false;  
  if( iso==4 && m_tuple->muo_isoWP4[idx] != 1 ) return false;  
  if( iso==5 && m_tuple->muo_isoWP5[idx] != 1 ) return false;  
  if( iso==6 && m_tuple->muo_isoWP6[idx] != 1 ) return false;  

  return true;
}


bool RecoTupleMaker :: passLooseElectron(int idx, int iso){
  if( m_tuple->ele_isLooseIDElectron[idx] !=1 ) return false;
  if( m_tuple->ele_isFromPV[idx] != 1 ) return false;
  if( m_tuple->ele_pt[idx] < 10000 ) return false;
  if( fabs(m_tuple->ele_cluseta[idx]) > 2.47 ) return false;
  //if( ! passElectronOR(idx) ) return false;

  if( iso==1 && m_tuple->ele_isoWP1[idx] != 1 ) return false;  
  if( iso==2 && m_tuple->ele_isoWP2[idx] != 1 ) return false;  
  if( iso==3 && m_tuple->ele_isoWP3[idx] != 1 ) return false;  
  if( iso==4 && m_tuple->ele_isoWP4[idx] != 1 ) return false;  
  if( iso==5 && m_tuple->ele_isoWP5[idx] != 1 ) return false;  
  if( iso==6 && m_tuple->ele_isoWP6[idx] != 1 ) return false;  

  return true;
}


int RecoTupleMaker :: countLooseLeptons(){
  int Nmu=0;
  for (int i = 0; i < m_tuple->nmuo ; i++){
    if( !passLooseMuon(i,2) ) continue;
    if( !passMuonOR(i) ) continue;
    if( m_tuple->muo_pt[i] < 25000) continue;
    //PrintMuon(i);
    Nmu++;
  }

  int Ne=0;
  for (int i = 0; i < m_tuple->nele ; i++){
    if( !passLooseElectron(i,2) ) continue;
    if( !passElectronOR(i) ) continue;
    if( m_tuple->ele_pt[i] < 25000) continue;
    //PrintElectron(i);
    Ne++;
  }
  return Nmu+Ne;
}

bool RecoTupleMaker :: passBadJetVeto(){
  
  ///https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HowToCleanJets2015
  for(int i = 0; i < m_tuple->njet ; i++){
    if( ! passCentralJet(i) ) continue;
    if( ! passJetOR(i)  ) continue;
    if( m_tuple->jet_good[i] == 0 ) return false;
  }
  
  return true;
}
