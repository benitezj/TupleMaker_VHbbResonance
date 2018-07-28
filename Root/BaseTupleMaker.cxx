#include "TupleMaker_VHbbResonance/BaseTupleMaker.h"
#include <fstream>
#include <sstream>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EventLoop/Job.h"
#include "EventLoop/OutputStream.h"
#include "EventLoop/Worker.h"
#include "PathResolver/PathResolver.h"

#include "xAODRootAccess/Init.h"
#include "xAODTruth/TruthParticleContainer.h"



ClassImp(BaseTupleMaker)

BaseTupleMaker :: BaseTupleMaker (){}


BaseTupleMaker :: BaseTupleMaker (std::string configPath) :
  config(ConfigStore::createStore(configPath)),
  m_debug(false),
  m_printTruth(0),
  m_maxEvents(-1),
  m_randRun(-1),
  m_currentSyst("none"),
  m_tuple(0),
  m_event(0),
  m_tree(0),
  m_filename("none"),
  m_tprw(0),
  m_pileupreweighting(0),
  m_ilumicalcFile(""),
  m_mcPileUpFile(""),
  m_EventInfoIn("none"),
  m_eventInfo(0),
  m_triggerTool(0),
  m_MCTruthIn("none")
{

  m_debug       = config->get<bool>("debug");
  m_filename    = config->get<std::string>("tuple.Label");
  if( m_debug) cout<<"tuple.Label: "<<m_filename<<endl;
  m_EventInfoIn = config->get<std::string>("tuple.EventInfo");
  if( m_debug) cout<<"tuple.EventInfo: "<<m_EventInfoIn<<endl;
  m_MCTruthIn   = config->get<std::string>("tuple.MCTruth");
  if( m_debug) cout<<"tuple.MCTrut: "<<m_MCTruthIn<<endl;
  m_printTruth  = config->get<int>("tuple.printTruth");
  if( m_debug) cout<<"tuple.printTruth: "<<m_printTruth<<endl;
  m_maxEvents   = config->get<int>("tuple.maxEvents");
  if( m_debug) cout<<"tuple.maxEvent: "<<m_maxEvents<<endl;

  ///Trigger paths
  m_triggerPaths = config->get< std::vector< std::string > >("tuple.triggerPaths");
  for(unsigned int i=0; i< m_triggerPaths.size(); i++)
    cout<<"Trigger Path: "<<m_triggerPaths[i]<<endl;

  ///////PRW configuration
  std::vector<std::string> m_pileUpMCNames = config->get<  std::vector< std::string > >("tuple.pileupMCFiles");
  for(unsigned int i=0;i< m_pileUpMCNames.size(); i++)
    pileUpMCFiles.push_back(std::string(gSystem->Getenv("ROOTCOREBIN"))+"/data/TupleMaker_VHbbResonance/"+m_pileUpMCNames[i]);  
  std::vector<std::string> m_pileUpDataNames = config->get< std::vector< std::string > >("tuple.pileupDataFiles");
  for(unsigned int i=0;i< m_pileUpDataNames.size(); i++)
    pileUpDataFiles.push_back(std::string(gSystem->Getenv("ROOTCOREBIN"))+"/data/TupleMaker_VHbbResonance/"+m_pileUpDataNames[i]);

  
  ////read in the list of events requested
  std::string eventListFile;
  config->getif< std::string  >("tuple.eventListFile",eventListFile);
  cout<<"EventListFile: "<<eventListFile<<endl;
  if(eventListFile.compare("")!=0){    
    std::ifstream infile(eventListFile.c_str());
    if(! infile.is_open()){
      cout<<"eventlist file was not found, but was requested in config."<<endl;
      exit(0);
    }

    std::string line;
    int run_id;
    ULong64_t event_id; 
    while (std::getline(infile, line)){
      std::istringstream iss(line);
      iss >> run_id >> event_id;
      if(0<run_id&&run_id<10000000){
	m_runsToAnalyze.push_back(run_id);
	m_eventsToAnalyze.push_back(event_id);
	if( m_debug) cout<<"EventList: "<<run_id<<" "<<event_id<<endl;
      }
    }
  }

  

  cout<<"BaseTupleMaker was created"<<endl;
}


BaseTupleMaker :: ~BaseTupleMaker (){
}


EL::StatusCode BaseTupleMaker :: setupJob (EL::Job& job)
{
 
  job.useXAOD();
  xAOD::Init( "BaseTupleMaker" ).ignore();
  
  EL::OutputStream out(m_filename);
  job.outputAdd(out);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode BaseTupleMaker :: histInitialize ()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode BaseTupleMaker :: initialize ()
{
  m_event = wk()->xaodEvent();
  Info("initialize()","Number of events = %lli", m_event->getEntries() );


  TFile * f = wk()->getOutputFile(m_filename);
  if(!f) { 
    Error( "initialize()" , "Failed to get output file");
    return EL::StatusCode::FAILURE;
  }


  ////Define the ntuple content, tuple could have been defined upstream  
  if(! m_tuple) 
    m_tuple = new BaseTuple();

  if(!m_tuple) { 
    Error( "initialize()" , "Failed to make m_tree or m_tuple");
    return EL::StatusCode::FAILURE;
  }
  
  ///create tree for each systematic variation
  m_systNames.push_back("Nominal");//Nominal tuple will always be made
  for(unsigned int i =0 ; i<m_systNames.size(); i++){
    m_tree = new TTree(TString("tuple_")+m_systNames[i].c_str(),m_systNames[i].c_str());
    m_tree->SetDirectory(f);
    m_tuple->DefineBranches(m_tree);
    m_systTuples.push_back(m_tree);
    cout<<"New Tree Created :"<<m_tree->GetName()<<endl;
  }

  Info("initialize()","Branches defined");

  ///pile-up tool 
  m_pileupreweighting = new CP::PileupReweightingTool("PileupReweightingTool");
  if( m_debug) m_pileupreweighting->msg().setLevel( MSG::DEBUG );

  m_tprw=new CP::TPileupReweighting("m_tprw");
  config->getif<std::string>("ilumicalcFile",m_ilumicalcFile);
  if(m_ilumicalcFile.compare("")!=0 ) {
    cout<<"BaseTupleMaker :: initialize m_ilumicalcFile. "<<m_ilumicalcFile.c_str()<<endl;

    std::vector<std::string> lc;
    lc.push_back(m_ilumicalcFile);
    m_pileupreweighting->setProperty("LumiCalcFiles",lc);
       
    std::vector<std::string> prw;
    config->getif<std::vector<std::string> >("prw",prw);
    if(prw.size() > 0 ) 
      m_pileupreweighting->setProperty("ConfigFiles",prw);

    m_pileupreweighting->setProperty("DefaultChannel",410000);
    m_pileupreweighting->initialize(); 


    //this is the sub-tool
    for(unsigned int j=0;j<prw.size();j++) 
      m_tprw->AddConfigFile(PathResolverFindCalibFile(prw[j].c_str()));
    m_tprw->AddLumiCalcFile(PathResolverFindCalibFile(m_ilumicalcFile.c_str()));
    m_tprw->SetDefaultChannel(410000);
    m_tprw->Initialize();
  }
  

  ///////Init MC Pile-up weights
  if( initPileUpWeights(f) == EL::StatusCode::FAILURE ){
    Error( "initialize()" , "Failed to initialize the pile-up weights");
    return EL::StatusCode::FAILURE;
  }


  ////Save the history of the sum of event weights
  cout<<"Clonning MetaData_EventCount from: "<<wk()->inputFile()->GetName()<<endl;
  TH1D* HEventCountsIn = dynamic_cast<TH1D*>(wk()->inputFile()->Get("MetaData_EventCount"));
  if (!HEventCountsIn) {
    Error("initialize()", (TString("MetaData_EventCount not found in ")+wk()->inputFile()->GetName()).Data());
    //return EL::StatusCode::FAILURE;

    ////When running on xAOD directly this histogram is not available, will create an empty one
    Error("initialize()", "Will make a empty MetaData_EventCount in the ouput");
    HEventCounts = new TH1D("MetaData_EventCount","MetaData_EventCount",10,0,10);
  }else {
    HEventCounts = (TH1D*)(HEventCountsIn->Clone("MetaData_EventCount"));
    //HEventCounts->Reset(); //do not reset to start adding first file here
  }
  HEventCounts->SetDirectory(f);
  



  ///Trigger SF tools
  m_triggerTool = new TriggerTool(*config,"lvqq");//currently use only single lepton trigger
  //EL_CHECK("BaseTupleMaker::init", );
  m_triggerTool->initialize();

  m_trig_sfmuon = new CP::MuonTriggerScaleFactors("mTrigSFClass");
  m_trig_sfmuon->setProperty("MuonQuality","Medium");
  m_trig_sfmuon->initialize();





  if(bookGenHistograms(f) == EL::StatusCode::FAILURE ){
    Error( "initialize()" , "Failed to book Generator histograms");
    return EL::StatusCode::FAILURE;
  }

  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode BaseTupleMaker::fileExecute(){
   
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode BaseTupleMaker::changeInput(bool firstFile){
  cout<<"Processing: "<<wk()->inputFile()->GetName()<<endl;
  
  //Add the sum of the event weights
  if(!firstFile){//HEventCounts does not exist for first file
    TH1D* HEventCountsIn = dynamic_cast<TH1D*>(wk()->inputFile()->Get("MetaData_EventCount"));
    if(!HEventCountsIn){
      Error("initialize()", (TString("MetaData_EventCount not found in ")+wk()->inputFile()->GetName()).Data());
      //return EL::StatusCode::FAILURE;
    }else {
      HEventCounts->Add(HEventCountsIn);
    }
  }
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode BaseTupleMaker :: execute ()
{
  //cout<<"In BaseTupleMaker :: execute "<<endl;
  
  ///skip events after requested
  if(int(getCounter("eventCounter_inputCounter")) >= m_maxEvents && m_maxEvents!=-1)
    return EL::StatusCode::SUCCESS;


  //retrieve the event info object 
  if( ! m_event->retrieve( m_eventInfo, m_EventInfoIn.c_str()).isSuccess() ) {
    Error("execute()","Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  
  if( int(getCounter("eventCounter_inputCounter")) % 100 == 0) 
    cout<<"Event number = "<<getCounter("eventCounter_inputCounter")
	<<" ,  selected events = "<<getCounter("eventCounter_outputCounter")
	<<" ,  runid = "<< m_eventInfo->runNumber() 
	<<" ,  eventid = "<< m_eventInfo->eventNumber() 
	<<endl;
  
  incrementCounter("eventCounter_inputCounter");
  

  /////here skip events not in the requested list
  if(m_eventsToAnalyze.size()>0){
    bool ev_found=false;
    for(unsigned int ev_counter=0;  ev_counter < m_eventsToAnalyze.size();  ev_counter++)
      if(m_eventInfo->runNumber() == m_runsToAnalyze[ev_counter] && 
	 m_eventInfo->eventNumber() == m_eventsToAnalyze[ev_counter] )
	ev_found=true;
    if(!ev_found) return EL::StatusCode::SUCCESS;
    incrementCounter("eventCounter_EventList");
  }  

  ///Debug info
  if( m_debug)
    cout<<"Analyzing Event : "<<m_eventInfo->runNumber()<<" "<<m_eventInfo->eventNumber()<<endl;

  
  ///loop over systematics
  for(unsigned int i =0 ; i<m_systTuples.size(); i++){

    // m_currentSyst is used by derived classes to grab the proper particle collections from the event
    m_currentSyst = m_systNames[i]; 

    //Only do systematics for MC
    if( static_cast<int>(m_eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) == 0 &&  m_currentSyst.compare("Nominal") != 0 ) 
      continue; 

    m_tree = m_systTuples[i]; 

    if( this->processEvent() != EL::StatusCode::FAILURE ){ 
      m_tree->Fill(); 
      incrementCounter("eventCounter_outputCounter");
    }

  }

  
  return EL::StatusCode::SUCCESS;
}


///////////////////////////

EL::StatusCode BaseTupleMaker :: processEvent() 
{

  //cout<<"In BaseTupleMaker :: processEvent "<<endl;

  if( processBaseEventInfo() == EL::StatusCode::FAILURE) 
    return EL::StatusCode::FAILURE;
  
  // if( processBaseEventCleaning() == EL::StatusCode::FAILURE) 
  // here just set a flag to be applied at plot level

  m_tuple->ntruth = 0;
  if( m_tuple->eve_isMC 
      && m_currentSyst.compare("Nominal") == 0  //only fill truth for Nominal tuple
      && processMCTruth() == EL::StatusCode::FAILURE ) 
    return EL::StatusCode::FAILURE;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode BaseTupleMaker :: processBaseEventInfo() 
{
  //cout<<"In BaseTupleMaker :: processBaseEventInfo "<<endl;
  m_tuple->eve_isMC = static_cast<int>(m_eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION));
  m_tuple->eve_num = m_eventInfo->eventNumber(); 
  m_tuple->eve_run = m_eventInfo->runNumber();
  m_tuple->eve_lb = (int) ( m_eventInfo->auxdata< unsigned int >( "lumiBlock" ) ) ;//why is there no lumiblock method?
  
  //generate a randoom run needed later
  if(m_tuple->eve_isMC) m_randRun = m_pileupreweighting->getRandomRunNumber(*m_eventInfo,false);
  else m_randRun = 0;

  //event cleaning only for data
  if(m_tuple->eve_isMC!=1) m_tuple->eve_isCleanEvent = Props::isCleanEvent.get(m_eventInfo);
  else m_tuple->eve_isCleanEvent = 1;
			     
  ///mu value needs correction in data
  m_tuple->eve_mu = m_eventInfo->averageInteractionsPerCrossing();
  if( ! m_tuple->eve_isMC &&  m_ilumicalcFile.compare("")!=0 ){    
    if( m_debug) std::cout<<" BaseTupleMaker :: processBaseEventInfo : "
			  <<" uncorr="<<m_eventInfo->averageInteractionsPerCrossing() 
			  <<" corr="<<m_pileupreweighting->getCorrectedMu(*m_eventInfo)
			  <<" corrTPRW="<<m_tprw->GetLumiBlockMu(m_tuple->eve_run,m_tuple->eve_lb)<<std::endl;
    m_tuple->eve_mu = m_tprw->GetLumiBlockMu(m_tuple->eve_run,m_tuple->eve_lb);
  }


  ///MC quantities
  m_tuple->eve_mc_w = 1.0;
  m_tuple->eve_mc_num = 0;
  m_tuple->eve_mc_chan = 0;
  if(m_tuple->eve_isMC){
    if( m_debug) std::cout<<"BaseTupleMaker :: processBaseEventInfo : MCEventWeight,mcEventNumber,mcChannelNumber"<<std::endl;
    if( Props::MCEventWeight.exists(m_eventInfo) ) m_tuple->eve_mc_w = Props::MCEventWeight.get(m_eventInfo); //under if() to be able to read xAOD's 
    m_tuple->eve_mc_num = m_eventInfo->mcEventNumber();
    m_tuple->eve_mc_chan = m_eventInfo->mcChannelNumber();
  }

  incrementCounter("inputSumW",m_tuple->eve_mc_w);
 
  //pile-up weights, for data will just be set to 1.
  if( fillPileUpWeights() == EL::StatusCode::FAILURE )
    return EL::StatusCode::FAILURE;  

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode  BaseTupleMaker :: fillPileUpWeights(){
  if( m_debug) std::cout<<"BaseTupleMaker :: fillPileUpWeights "<<std::endl;

  ///pile-up distribution
  HPileUp->Fill(m_tuple->eve_mu);

  //init the weights
  for(int i=0;i<m_tuple->nmcpuw;i++)
    m_tuple->eve_mc_puw[i] = 1.0;

  if(m_tuple->eve_isMC){
    if( m_debug) std::cout<<"BaseTupleMaker :: fillPileUpWeights for MC"<<std::endl;
    for(int i=0;i<m_tuple->nmcpuw;i++){

      if(!PileUpRatios[i]) return EL::StatusCode::FAILURE;
      int mubin= 1 + int(m_tuple->eve_mu);
      if(mubin<1 || mubin>PileUpRatios[i]->GetNbinsX())
	return EL::StatusCode::FAILURE;

      m_tuple->eve_mc_puw[i] = PileUpRatios[i]->GetBinContent(mubin);
    }    
  }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode BaseTupleMaker :: processMCTruth() {
  if( m_debug) std::cout<<" BaseTupleMaker :: processMCTruth"<<std::endl;

  if(! m_tuple->eve_isMC){
    cout<<" Calling processMCTruth() for Data. "<<endl;
    return EL::StatusCode::FAILURE;
  }
  
  if(m_printTruth > 0 ) printMCTruth();
  
  const xAOD::TruthParticleContainer* Truth = 0;
  if( ! m_event->retrieve( Truth , m_MCTruthIn.c_str() ).isSuccess() ){
    Error("execute()", "Failed to retrieve  TruthParticle collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  for (unsigned int i = 0 ; i < Truth->size(); ++i) {
    if( m_tuple->ntruth == BaseTuple::MAXTRUTH){
      //cout<<"processMCTruth(): number of truth particles is at max "<<BaseTuple::MAXTRUTH<<endl;
      break;
    }
    const xAOD::TruthParticle* part = Truth->at(i);
    
    if(part->status()==1 || part->status()==2 ) continue;

    fillTruthParticle(m_tuple->ntruth,part->p4(),part->status(),part->pdgId());
    m_tuple->ntruth++;    
  }  


  fillMCProcesses();
  fillGenHistograms();
 

  return EL::StatusCode::SUCCESS;
}


void BaseTupleMaker :: fillTruthParticle(int index, TLorentzVector P, int status, int pdg){
  m_tuple->truth_m[index]= P.M();
  m_tuple->truth_p[index]= P.P();
  m_tuple->truth_pt[index]= P.Pt();
  if(m_tuple->truth_pt[index] > 1 ) m_tuple->truth_eta[index]= P.Eta();
    else m_tuple->truth_eta[index] = -9999;
  m_tuple->truth_phi[index]= P.Phi();
  m_tuple->truth_status[index]= status;
  m_tuple->truth_pdg[index]= pdg;
}

void BaseTupleMaker :: fillFakeComposite(int * p, int l1, int l2, int pdgid){
    if( m_tuple->ntruth == BaseTuple::MAXTRUTH){
      cout<<"fillMCProcesses(): number of truth particles is at max. cannot create fake composite"
	  <<BaseTuple::MAXTRUTH<<endl;
      return;
    }
    TLorentzVector P1;
    P1.SetPtEtaPhiM(m_tuple->truth_pt[l1],m_tuple->truth_eta[l1],
		    m_tuple->truth_phi[l1],m_tuple->truth_m[l1]);
    TLorentzVector P2;
    P2.SetPtEtaPhiM(m_tuple->truth_pt[l2],m_tuple->truth_eta[l2],
		    m_tuple->truth_phi[l2],m_tuple->truth_m[l2]);

    fillTruthParticle(m_tuple->ntruth,P1+P2,-1,pdgid);
    *p = m_tuple->ntruth;
    m_tuple->ntruth++;
}


void BaseTupleMaker :: fillMCProcesses(){
  if( m_debug) std::cout<<" BaseTupleMaker :: fillMCProcesses"<<std::endl;
  if(! m_tuple->eve_isMC){
    cout<<" Calling fillMCProcesses() for Data. "<<endl;
    return;
  }

  
  ////////////////////////
  //Find Z-->mm
  ///////////////////////
  m_tuple->truth_Zmm=-1;
  m_tuple->truth_Zmm_m1=-1;
  m_tuple->truth_Zmm_m2=-1;

  ///Find the muons from the Z decay, order by pT
  for (int l1 = 0 ; l1 < m_tuple->ntruth; ++l1)
    for (int l2 = 0 ; l2 < m_tuple->ntruth; ++l2)
      if(m_tuple->truth_status[l1]==3 && m_tuple->truth_pdg[l1] == -13 //mu+
	 && m_tuple->truth_status[l2]==3  && m_tuple->truth_pdg[l2] == 13 //m-
	 ){
	if(m_tuple->truth_pt[l1] > m_tuple->truth_pt[l2]){
	  m_tuple->truth_Zmm_m1=l1;
	  m_tuple->truth_Zmm_m2=l2;
	}else{
	  m_tuple->truth_Zmm_m1=l2;
	  m_tuple->truth_Zmm_m2=l1;
	}	
      }


  for (int i = 0 ; i < m_tuple->ntruth; ++i) 
    if(m_tuple->truth_status[i]==62 && m_tuple->truth_pdg[i]==23)
      m_tuple->truth_Zmm=i;

  //in case Z was not found make a fake one from the muons
  if(m_tuple->truth_Zmm==-1 && m_tuple->truth_Zmm_m1>-1 && m_tuple->truth_Zmm_m2>-1)
    fillFakeComposite(&(m_tuple->truth_Zmm),m_tuple->truth_Zmm_m1,m_tuple->truth_Zmm_m2,23);


  ///////////////////////
  //Find Z-->ee
  ///////////////////////
  m_tuple->truth_Zee=-1;
  m_tuple->truth_Zee_e1=-1;
  m_tuple->truth_Zee_e2=-1;
    
  ///electrons 
  for (int l1 = 0 ; l1 < m_tuple->ntruth; ++l1)
    for (int l2 = 0 ; l2 < m_tuple->ntruth; ++l2)
      if(m_tuple->truth_status[l1]==3 && m_tuple->truth_pdg[l1] == -11
	 && m_tuple->truth_status[l2]==3  && m_tuple->truth_pdg[l2] == 11
	 ){
	if(m_tuple->truth_pt[l1] > m_tuple->truth_pt[l2]){
	  m_tuple->truth_Zee_e1=l1;
	  m_tuple->truth_Zee_e2=l2;
	}else{
	  m_tuple->truth_Zee_e1=l2;
	  m_tuple->truth_Zee_e2=l1;
	}	
      }

  for (int i = 0 ; i < m_tuple->ntruth; ++i) 
    if(m_tuple->truth_status[i]==62 && m_tuple->truth_pdg[i]==23)
      m_tuple->truth_Zee=i;

  //in case Z was not found make a fake one from 
  if(m_tuple->truth_Zee==-1 && m_tuple->truth_Zee_e1>-1 && m_tuple->truth_Zee_e2>-1)
    fillFakeComposite(&(m_tuple->truth_Zee),m_tuple->truth_Zee_e1,m_tuple->truth_Zee_e2,23);



  ////////////////////////
  //Find Z-->tau tau
  ///////////////////////
  m_tuple->truth_Ztt=-1;
  m_tuple->truth_Ztt_t1=-1;
  m_tuple->truth_Ztt_t2=-1;

  ///First loop in case the Z is not stored (Sherpa samples)
  for (int l1 = 0 ; l1 < m_tuple->ntruth; ++l1)
    for (int l2 = 0 ; l2 < m_tuple->ntruth; ++l2)
      if(m_tuple->truth_status[l1]==3 && m_tuple->truth_pdg[l1] == -15 
	 && m_tuple->truth_status[l2]==3  && m_tuple->truth_pdg[l2] == 15 
	 ){
	if(m_tuple->truth_pt[l1] > m_tuple->truth_pt[l2]){
	  m_tuple->truth_Ztt_t1=l1;
	  m_tuple->truth_Ztt_t2=l2;
	}else{
	  m_tuple->truth_Ztt_t1=l2;
	  m_tuple->truth_Ztt_t2=l1;
	}	
      }


  for (int i = 0 ; i < m_tuple->ntruth; ++i) 
    if(m_tuple->truth_status[i]==62 && m_tuple->truth_pdg[i]==23)
      m_tuple->truth_Ztt=i;


  //in case Z was not found make a fake one from 
  if(m_tuple->truth_Ztt==-1 && m_tuple->truth_Ztt_t1>-1 && m_tuple->truth_Ztt_t2>-1)
    fillFakeComposite(&(m_tuple->truth_Ztt),m_tuple->truth_Ztt_t1,m_tuple->truth_Ztt_t2,23);


  ///////////////////////
  //Find W->mu nu
  ///////////////////////
  m_tuple->truth_Wmv=-1;
  m_tuple->truth_Wmv_m=-1;
  m_tuple->truth_Wmv_v=-1;

  ///First loop in case the W is not stored (Sherpa samples)
  for (int i = 0 ; i < m_tuple->ntruth; ++i){
    if(m_tuple->truth_status[i]==3 && abs(m_tuple->truth_pdg[i]) == 13)
      m_tuple->truth_Wmv_m=i;
    if(m_tuple->truth_status[i]==3 && abs(m_tuple->truth_pdg[i]) == 14)
      m_tuple->truth_Wmv_v=i;
    if(m_tuple->truth_status[i]==62 && abs(m_tuple->truth_pdg[i])== 24)
      m_tuple->truth_Wmv=i;
  }

  //in case W was not found make a fake one 
  if(m_tuple->truth_Wmv==-1 && m_tuple->truth_Wmv_m>-1 && m_tuple->truth_Wmv_v>-1)
    fillFakeComposite(&(m_tuple->truth_Wmv),m_tuple->truth_Wmv_m,m_tuple->truth_Wmv_v,24);


  ///////////////////////
  //Find W->e nu
  ///////////////////////
  m_tuple->truth_Wev=-1;
  m_tuple->truth_Wev_e=-1;
  m_tuple->truth_Wev_v=-1;

  for (int i = 0 ; i < m_tuple->ntruth; ++i){
    if(m_tuple->truth_status[i]==3 && abs(m_tuple->truth_pdg[i]) == 11)
      m_tuple->truth_Wev_e=i;
    if(m_tuple->truth_status[i]==3 && abs(m_tuple->truth_pdg[i]) == 12)
      m_tuple->truth_Wev_v=i;
    if(m_tuple->truth_status[i]==62 && abs(m_tuple->truth_pdg[i])== 24 )
      m_tuple->truth_Wev=i;
  }
 
  //in case W was not found make a fake one 
  if(m_tuple->truth_Wev==-1 && m_tuple->truth_Wev_e>-1 && m_tuple->truth_Wev_v>-1)
    fillFakeComposite(&(m_tuple->truth_Wev),m_tuple->truth_Wev_e,m_tuple->truth_Wev_v,24);



  ///////////////////////
  //Find W->tau nu
  ///////////////////////
  m_tuple->truth_Wtv=-1;
  m_tuple->truth_Wtv_t=-1;
  m_tuple->truth_Wtv_v=-1;

  ///First loop in case the W is not stored (Sherpa samples)
  for (int i = 0 ; i < m_tuple->ntruth; ++i){
    if(m_tuple->truth_status[i]==3 && abs(m_tuple->truth_pdg[i]) == 15)
      m_tuple->truth_Wtv_t=i;
    if(m_tuple->truth_status[i]==3 && abs(m_tuple->truth_pdg[i]) == 16)
      m_tuple->truth_Wtv_v=i;
    if(m_tuple->truth_status[i]==62 && abs(m_tuple->truth_pdg[i])== 24 )
      m_tuple->truth_Wtv=i;
  }
 
  //in case W was not found make a fake one 
  if(m_tuple->truth_Wtv==-1 && m_tuple->truth_Wtv_t>-1 && m_tuple->truth_Wtv_v>-1)
    fillFakeComposite(&(m_tuple->truth_Wtv),m_tuple->truth_Wtv_t,m_tuple->truth_Wtv_v,24);


  //////////////////
  ///define truth_V index for convenience
  ////////////////
  m_tuple->truth_V=-1;
  if(m_tuple->truth_Zee>-1)m_tuple->truth_V=m_tuple->truth_Zee;
  if(m_tuple->truth_Zmm>-1)m_tuple->truth_V=m_tuple->truth_Zmm;
  if(m_tuple->truth_Ztt>-1)m_tuple->truth_V=m_tuple->truth_Ztt;
  if(m_tuple->truth_Wev>-1)m_tuple->truth_V=m_tuple->truth_Wev;
  if(m_tuple->truth_Wmv>-1)m_tuple->truth_V=m_tuple->truth_Wmv;
  if(m_tuple->truth_Wtv>-1)m_tuple->truth_V=m_tuple->truth_Wtv;

  ////////////////////////
  //Find H->bb
  ///////////////////////
  m_tuple->truth_H0bb=-1;
  m_tuple->truth_H0bb_b1=-1;
  m_tuple->truth_H0bb_b2=-1;

  ///First loop in case the Z is not stored (Sherpa samples)
  for (int l1 = 0 ; l1 < m_tuple->ntruth; ++l1)
    for (int l2 = 0 ; l2 < m_tuple->ntruth; ++l2)
      if(m_tuple->truth_status[l1]==3 && m_tuple->truth_pdg[l1] == 5 //b
	 && m_tuple->truth_status[l2]==3  && m_tuple->truth_pdg[l2] == -5 //b-bar
	 ){
	if(m_tuple->truth_pt[l1] > m_tuple->truth_pt[l2]){
	  m_tuple->truth_H0bb_b1=l1;
	  m_tuple->truth_H0bb_b2=l2;
	}else{
	  m_tuple->truth_H0bb_b1=l2;
	  m_tuple->truth_H0bb_b2=l1;
	}	
      }

  for (int i = 0 ; i < m_tuple->ntruth; ++i) 
    if(m_tuple->truth_status[i]==62 && m_tuple->truth_pdg[i]==25)
      m_tuple->truth_H0bb=i;
  
  //in case H0 was not found make a fake one 
  if(m_tuple->truth_H0bb==-1 && m_tuple->truth_H0bb_b1>-1 && m_tuple->truth_H0bb_b2>-1)
    fillFakeComposite(&(m_tuple->truth_H0bb),m_tuple->truth_H0bb_b1,m_tuple->truth_H0bb_b2,25);



  /////////////////////////////////////////////////////////////////
  //VH composites: V and H must be in the truth list (does not work for Sherpa samples)
  /////////////////////////////////////////////////////////////////////
  m_tuple->truth_ZHmmbb=-1;
  m_tuple->truth_ZHeebb=-1;
  m_tuple->truth_WHmvbb=-1;
  m_tuple->truth_WHevbb=-1;
  if(m_tuple->truth_Zmm>-1 && m_tuple->truth_H0bb>-1)
    fillFakeComposite(&(m_tuple->truth_ZHmmbb),m_tuple->truth_Zmm,m_tuple->truth_H0bb,131355);

  if(m_tuple->truth_Zee>-1 && m_tuple->truth_H0bb>-1)
    fillFakeComposite(&(m_tuple->truth_ZHeebb),m_tuple->truth_Zee,m_tuple->truth_H0bb,111155);

  if(m_tuple->truth_Wmv>-1 && m_tuple->truth_H0bb>-1)
    fillFakeComposite(&(m_tuple->truth_WHmvbb),m_tuple->truth_Wmv,m_tuple->truth_H0bb,131455);

  if(m_tuple->truth_Wev>-1 && m_tuple->truth_H0bb>-1)
    fillFakeComposite(&(m_tuple->truth_WHevbb),m_tuple->truth_Wev,m_tuple->truth_H0bb,111255);


  ////////////////////////
  //Find R->VH (V=Z/W), H->bb
  ///////////////////////
  m_tuple->truth_RToVH=-1;
  m_tuple->truth_RToVH_V=-1;
  m_tuple->truth_RToVH_H=-1;
  m_tuple->truth_RToVH_H_b1=-1;
  m_tuple->truth_RToVH_H_b2=-1;
  
  for (int i = 0 ; i < m_tuple->ntruth; ++i){
    if(m_tuple->truth_status[i]==62 && abs(m_tuple->truth_pdg[i])==9000001)
      m_tuple->truth_RToVH=i;
 
    if(m_tuple->truth_status[i]==22 && ( abs(m_tuple->truth_pdg[i])==23 ||  abs(m_tuple->truth_pdg[i])==24) )
      m_tuple->truth_RToVH_V=i;
 
    if(m_tuple->truth_status[i]==22 && abs(m_tuple->truth_pdg[i])==25)
      m_tuple->truth_RToVH_H=i;
 
  }
  if(m_printTruth >=2)
    cout<<"RToVH: "<<m_tuple->truth_RToVH<<" "<<m_tuple->truth_RToVH_V<<" "<<m_tuple->truth_RToVH_H<<endl;

  //Resonance was not kept in derivations, will create it from the V+H
  if( m_tuple->truth_RToVH == -1 && m_tuple->truth_RToVH_V > -1 && m_tuple->truth_RToVH_H > -1)
    fillFakeComposite(&(m_tuple->truth_RToVH),m_tuple->truth_RToVH_V,m_tuple->truth_RToVH_H,9000001);

  //find the decay of the Higgs : 
  m_tuple->truth_RToVH_H_decay=-1;
  for (int l1 = 0 ; l1 < m_tuple->ntruth; ++l1)
    for (int l2 = 0 ; l2 < m_tuple->ntruth; ++l2){
      if( l1 != l2 ) {
	if(m_tuple->truth_status[l1]==23 && m_tuple->truth_pdg[l1] == 4 //c
	   && m_tuple->truth_status[l2]==23  && m_tuple->truth_pdg[l2] == -4 //c-bar
	   ) m_tuple->truth_RToVH_H_decay=4;
	if(m_tuple->truth_status[l1]==23 && m_tuple->truth_pdg[l1] == 5 //b
	   && m_tuple->truth_status[l2]==23  && m_tuple->truth_pdg[l2] == -5 //b-bar
	   ) m_tuple->truth_RToVH_H_decay=5;
	if(m_tuple->truth_status[l1]==23 && m_tuple->truth_pdg[l1] == 21 //g
	   && m_tuple->truth_status[l2]==23  && m_tuple->truth_pdg[l2] == 21 //g
	   ) m_tuple->truth_RToVH_H_decay=21;
      }
    }


  //find the b's and pT order 
  for (int l1 = 0 ; l1 < m_tuple->ntruth; ++l1)
    for (int l2 = 0 ; l2 < m_tuple->ntruth; ++l2)
      if(m_tuple->truth_status[l1]==23 && m_tuple->truth_pdg[l1] == 5 //b
	 && m_tuple->truth_status[l2]==23  && m_tuple->truth_pdg[l2] == -5 //b-bar
	 ){
	if(m_tuple->truth_pt[l1] > m_tuple->truth_pt[l2]){
	  m_tuple->truth_RToVH_H_b1=l1;
	  m_tuple->truth_RToVH_H_b2=l2;
	}else{
	  m_tuple->truth_RToVH_H_b1=l2;
	  m_tuple->truth_RToVH_H_b2=l1;
	}	
      }
  
  if(m_tuple->truth_RToVH_H_b1>-1 && m_tuple->truth_RToVH_H_b2>-1)
    m_tuple->truth_RToVH_H_dR = deltaR(m_tuple->truth_eta[m_tuple->truth_RToVH_H_b1],
				       m_tuple->truth_phi[m_tuple->truth_RToVH_H_b1],
				       m_tuple->truth_eta[m_tuple->truth_RToVH_H_b2],
				       m_tuple->truth_phi[m_tuple->truth_RToVH_H_b2]);




  //Find t-tbar (needed to combine with mtt slices)
  for (int i = 0 ; i < m_tuple->ntruth; ++i){
    if(m_tuple->truth_status[i]==3 && m_tuple->truth_pdg[i] == 6)
      m_tuple->truth_ttbar_t1=i;
    if(m_tuple->truth_status[i]==3 && m_tuple->truth_pdg[i] == -6)
      m_tuple->truth_ttbar_t2=i;
  }
  if(m_tuple->truth_ttbar_t1>0 && m_tuple->truth_ttbar_t2>0)
    fillFakeComposite(&(m_tuple->truth_ttbar),m_tuple->truth_ttbar_t1,m_tuple->truth_ttbar_t2,9000066);


  ////Leave this for later
  // m_tuple->truth_top=-1;
  // m_tuple->truth_top_mm_m1=-1;
  // m_tuple->truth_top_mm_m2=-1;
  // m_tuple->truth_top_ee_e1=-1;
  // m_tuple->truth_top_ee_e2=-1;
  // m_tuple->truth_top_em_e=-1;
  // m_tuple->truth_top_em_m=-1;

}


EL::StatusCode BaseTupleMaker :: bookGenHistograms(TFile *f ){
  if( m_debug) std::cout<<" BaseTupleMaker :: bookGenHistograms"<<std::endl;

  ///NOTE: should set Sumw2() otherwise cannot use bin errors later

  HGen_Zmm_m = new TH1F("HGen_Zmm_m","HGen_Zmm_m",1000,0,1000);//in GeV
  HGen_Zmm_pt= new TH1F("HGen_Zmm_pt","HGen_Zmm_pt",1000,0,1000);
  HGen_Zmm_eta = new TH1F("HGen_Zmm_eta","HGen_Zmm_eta",200,-10,10);
  HGen_Zmm_m1_pt = new TH1F("HGen_Zmm_m1_pt","HGen_Zmm_m1_pt",1000,0,1000);
  HGen_Zmm_m1_eta = new TH1F("HGen_Zmm_m1_eta","HGen_Zmm_m1_eta",200,-10,10);
  HGen_Zmm_m2_pt = new TH1F("HGen_Zmm_m2_pt","HGen_Zmm_m2_pt",1000,0,1000);
  HGen_Zmm_m2_eta = new TH1F("HGen_Zmm_m2_eta","HGen_Zmm_m2_eta",200,-10,10);
  HGen_Zmm_m->SetDirectory(f);
  HGen_Zmm_pt->SetDirectory(f);
  HGen_Zmm_eta->SetDirectory(f);
  HGen_Zmm_m1_pt->SetDirectory(f);
  HGen_Zmm_m1_eta->SetDirectory(f);
  HGen_Zmm_m2_pt->SetDirectory(f);
  HGen_Zmm_m2_eta->SetDirectory(f);

  HGen_Zee_m = new TH1F("HGen_Zee_m","HGen_Zee_m",1000,0,1000);//in GeV
  HGen_Zee_pt= new TH1F("HGen_Zee_pt","HGen_Zee_pt",1000,0,1000);
  HGen_Zee_eta = new TH1F("HGen_Zee_eta","HGen_Zee_eta",200,-10,10);
  HGen_Zee_e1_pt = new TH1F("HGen_Zee_e1_pt","HGen_Zee_e1_pt",1000,0,1000);
  HGen_Zee_e1_eta = new TH1F("HGen_Zee_e1_eta","HGen_Zee_e1_eta",200,-10,10);
  HGen_Zee_e2_pt = new TH1F("HGen_Zee_e2_pt","HGen_Zee_e2_pt",1000,0,1000);
  HGen_Zee_e2_eta = new TH1F("HGen_Zee_e2_eta","HGen_Zee_e2_eta",200,-10,10);
  HGen_Zee_m->SetDirectory(f);
  HGen_Zee_pt->SetDirectory(f);
  HGen_Zee_eta->SetDirectory(f);
  HGen_Zee_e1_pt->SetDirectory(f);
  HGen_Zee_e1_eta->SetDirectory(f);
  HGen_Zee_e2_pt->SetDirectory(f);
  HGen_Zee_e2_eta->SetDirectory(f);

  HGen_Ztt_m = new TH1F("HGen_Ztt_m","HGen_Ztt_m",1000,0,1000);//in GeV
  HGen_Ztt_pt= new TH1F("HGen_Ztt_pt","HGen_Ztt_pt",1000,0,1000);
  HGen_Ztt_eta = new TH1F("HGen_Ztt_eta","HGen_Ztt_eta",200,-10,10);
  HGen_Ztt_t1_pt = new TH1F("HGen_Ztt_t1_pt","HGen_Ztt_t1_pt",1000,0,1000);
  HGen_Ztt_t1_eta = new TH1F("HGen_Ztt_t1_eta","HGen_Ztt_t1_eta",200,-10,10);
  HGen_Ztt_t2_pt = new TH1F("HGen_Ztt_t2_pt","HGen_Ztt_t2_pt",1000,0,1000);
  HGen_Ztt_t2_eta = new TH1F("HGen_Ztt_t2_eta","HGen_Ztt_t2_eta",200,-10,10);
  HGen_Ztt_m->SetDirectory(f);
  HGen_Ztt_pt->SetDirectory(f);
  HGen_Ztt_eta->SetDirectory(f);
  HGen_Ztt_t1_pt->SetDirectory(f);
  HGen_Ztt_t1_eta->SetDirectory(f);
  HGen_Ztt_t2_pt->SetDirectory(f);
  HGen_Ztt_t2_eta->SetDirectory(f);

  HGen_Wmv_m = new TH1F("HGen_Wmv_m","HGen_Wmv_m",1000,0,1000);//in GeV
  HGen_Wmv_pt= new TH1F("HGen_Wmv_pt","HGen_Wmv_pt",1000,0,1000);
  HGen_Wmv_eta = new TH1F("HGen_Wmv_eta","HGen_Wmv_eta",200,-10,10);
  HGen_Wmv_m_pt = new TH1F("HGen_Wmv_m_pt","HGen_Wmv_m_pt",1000,0,1000);
  HGen_Wmv_m_eta = new TH1F("HGen_Wmv_m_eta","HGen_Wmv_m_eta",200,-10,10);
  HGen_Wmv_v_pt = new TH1F("HGen_Wmv_v_pt","HGen_Wmv_v_pt",1000,0,1000);
  HGen_Wmv_v_eta = new TH1F("HGen_Wmv_v_eta","HGen_Wmv_v_eta",200,-10,10);
  HGen_Wmv_m->SetDirectory(f);
  HGen_Wmv_pt->SetDirectory(f);
  HGen_Wmv_eta->SetDirectory(f);
  HGen_Wmv_m_pt->SetDirectory(f);
  HGen_Wmv_m_eta->SetDirectory(f);
  HGen_Wmv_v_pt->SetDirectory(f);
  HGen_Wmv_v_eta->SetDirectory(f);

  HGen_Wev_m = new TH1F("HGen_Wev_m","HGen_Wev_m",1000,0,1000);//in GeV
  HGen_Wev_pt= new TH1F("HGen_Wev_pt","HGen_Wev_pt",1000,0,1000);
  HGen_Wev_eta = new TH1F("HGen_Wev_eta","HGen_Wev_eta",200,-10,10);
  HGen_Wev_e_pt = new TH1F("HGen_Wev_e_pt","HGen_Wev_e_pt",1000,0,1000);
  HGen_Wev_e_eta = new TH1F("HGen_Wev_e_eta","HGen_Wev_e_eta",200,-10,10);
  HGen_Wev_v_pt = new TH1F("HGen_Wev_v_pt","HGen_Wev_v_pt",1000,0,1000);
  HGen_Wev_v_eta = new TH1F("HGen_Wev_v_eta","HGen_Wev_v_eta",200,-10,10);
  HGen_Wev_m->SetDirectory(f);
  HGen_Wev_pt->SetDirectory(f);
  HGen_Wev_eta->SetDirectory(f);
  HGen_Wev_e_pt->SetDirectory(f);
  HGen_Wev_e_eta->SetDirectory(f);
  HGen_Wev_v_pt->SetDirectory(f);
  HGen_Wev_v_eta->SetDirectory(f);

  HGen_Wtv_m = new TH1F("HGen_Wtv_m","HGen_Wtv_m",1000,0,1000);//in GeV
  HGen_Wtv_pt= new TH1F("HGen_Wtv_pt","HGen_Wtv_pt",1000,0,1000);
  HGen_Wtv_eta = new TH1F("HGen_Wtv_eta","HGen_Wtv_eta",200,-10,10);
  HGen_Wtv_t_pt = new TH1F("HGen_Wtv_t_pt","HGen_Wtv_t_pt",1000,0,1000);
  HGen_Wtv_t_eta = new TH1F("HGen_Wtv_t_eta","HGen_Wtv_t_eta",200,-10,10);
  HGen_Wtv_v_pt = new TH1F("HGen_Wtv_v_pt","HGen_Wtv_v_pt",1000,0,1000);
  HGen_Wtv_v_eta = new TH1F("HGen_Wtv_v_eta","HGen_Wtv_v_eta",200,-10,10);
  HGen_Wtv_m->SetDirectory(f);
  HGen_Wtv_pt->SetDirectory(f);
  HGen_Wtv_eta->SetDirectory(f);
  HGen_Wtv_t_pt->SetDirectory(f);
  HGen_Wtv_t_eta->SetDirectory(f);
  HGen_Wtv_v_pt->SetDirectory(f);
  HGen_Wtv_v_eta->SetDirectory(f);

  HGen_H0bb_m = new TH1F("HGen_H0bb_m","HGen_H0bb_m",1000,0,1000);//in GeV
  HGen_H0bb_pt= new TH1F("HGen_H0bb_pt","HGen_H0bb_pt",1000,0,1000);
  HGen_H0bb_eta = new TH1F("HGen_H0bb_eta","HGen_H0bb_eta",200,-10,10);
  HGen_H0bb_b1_pt = new TH1F("HGen_H0bb_b1_pt","HGen_H0bb_b1_pt",1000,0,1000);
  HGen_H0bb_b1_eta = new TH1F("HGen_H0bb_b1_eta","HGen_H0bb_b1_eta",200,-10,10);
  HGen_H0bb_b2_pt = new TH1F("HGen_H0bb_b2_pt","HGen_H0bb_b2_pt",1000,0,1000);
  HGen_H0bb_b2_eta = new TH1F("HGen_H0bb_b2_eta","HGen_H0bb_b2_eta",200,-10,10);
  HGen_H0bb_m->SetDirectory(f);
  HGen_H0bb_pt->SetDirectory(f);
  HGen_H0bb_eta->SetDirectory(f);
  HGen_H0bb_b1_pt->SetDirectory(f);
  HGen_H0bb_b1_eta->SetDirectory(f);
  HGen_H0bb_b2_pt->SetDirectory(f);
  HGen_H0bb_b2_eta->SetDirectory(f);


  ///These will be filled with either ZH-->mmbb or ZH-->eebb (assume only one is created per event)
  HGen_ZH_m = new TH1F("HGen_ZH_m","HGen_ZH_m",1000,0,10000);//in GeV
  HGen_ZH_pt= new TH1F("HGen_ZH_pt","HGen_ZH_pt",1000,0,1000);
  HGen_ZH_eta = new TH1F("HGen_ZH_eta","HGen_ZH_eta",200,-10,10);
  HGen_ZH_m->SetDirectory(f);
  HGen_ZH_pt->SetDirectory(f);
  HGen_ZH_eta->SetDirectory(f);

  HGen_WH_m = new TH1F("HGen_WH_m","HGen_WH_m",1000,0,10000);//in GeV
  HGen_WH_pt= new TH1F("HGen_WH_pt","HGen_WH_pt",1000,0,1000);
  HGen_WH_eta = new TH1F("HGen_WH_eta","HGen_WH_eta",200,-10,10);
  HGen_WH_m->SetDirectory(f);
  HGen_WH_pt->SetDirectory(f);
  HGen_WH_eta->SetDirectory(f);


  //VH Resonance analysis
  HGen_RToVH_m = new TH1F("HGen_RToVH_m","HGen_RToVH_m",1000,0,10000);//in 10 GeV
  HGen_RToVH_pt= new TH1F("HGen_RToVH_pt","HGen_RToVH_pt",1000,0,1000);
  HGen_RToVH_eta = new TH1F("HGen_RToVH_eta","HGen_RToVH_eta",200,-10,10);
  HGen_RToVH_m->SetDirectory(f);
  HGen_RToVH_pt->SetDirectory(f);
  HGen_RToVH_eta->SetDirectory(f);

  HGen_RToVH_V_m = new TH1F("HGen_RToVH_V_m","HGen_RToVH_V_m",1000,0,1000);//in GeV
  HGen_RToVH_V_pt= new TH1F("HGen_RToVH_V_pt","HGen_RToVH_V_pt",2000,0,2000);
  HGen_RToVH_V_eta = new TH1F("HGen_RToVH_V_eta","HGen_RToVH_V_eta",200,-10,10);
  HGen_RToVH_V_m->SetDirectory(f);
  HGen_RToVH_V_pt->SetDirectory(f);
  HGen_RToVH_V_eta->SetDirectory(f);

  HGen_RToVH_H_m = new TH1F("HGen_RToVH_H_m","HGen_RToVH_H_m",1000,0,1000);//in GeV
  HGen_RToVH_H_pt= new TH1F("HGen_RToVH_H_pt","HGen_RToVH_H_pt",2000,0,2000);
  HGen_RToVH_H_eta = new TH1F("HGen_RToVH_H_eta","HGen_RToVH_H_eta",200,-10,10);
  HGen_RToVH_H_m->SetDirectory(f);
  HGen_RToVH_H_pt->SetDirectory(f);
  HGen_RToVH_H_eta->SetDirectory(f);

  HGen_RToVH_H_b1_m = new TH1F("HGen_RToVH_H_b1_m","HGen_RToVH_H_b1_m",1000,0,1000);//in GeV
  HGen_RToVH_H_b1_pt= new TH1F("HGen_RToVH_H_b1_pt","HGen_RToVH_H_b1_pt",1000,0,1000);
  HGen_RToVH_H_b1_eta = new TH1F("HGen_RToVH_H_b1_eta","HGen_RToVH_H_b1_eta",200,-10,10);
  HGen_RToVH_H_b1_m->SetDirectory(f);
  HGen_RToVH_H_b1_pt->SetDirectory(f);
  HGen_RToVH_H_b1_eta->SetDirectory(f);

  HGen_RToVH_H_b2_m = new TH1F("HGen_RToVH_H_b2_m","HGen_RToVH_H_b2_m",1000,0,1000);//in GeV
  HGen_RToVH_H_b2_pt= new TH1F("HGen_RToVH_H_b2_pt","HGen_RToVH_H_b2_pt",1000,0,1000);
  HGen_RToVH_H_b2_eta = new TH1F("HGen_RToVH_H_b2_eta","HGen_RToVH_H_b2_eta",200,-10,10);
  HGen_RToVH_H_b2_m->SetDirectory(f);
  HGen_RToVH_H_b2_pt->SetDirectory(f);
  HGen_RToVH_H_b2_eta->SetDirectory(f);

  HGen_RToVH_H_p= new TH1F("HGen_RToVH_H_p","HGen_RToVH_H_p",1000,0,10000);
  HGen_RToVH_H_dR= new TH1F("HGen_RToVH_H_dR","HGen_RToVH_H_dR",100,0,10);
  HGen_RToVH_H_dRVsp= new TH2F("HGen_RToVH_H_dRVsp","HGen_RToVH_H_dRVsp",1000,0,10000,100,0,10);
  HGen_RToVH_H_p->SetDirectory(f);
  HGen_RToVH_H_dR->SetDirectory(f);
  HGen_RToVH_H_dRVsp->SetDirectory(f);



  return EL::StatusCode::SUCCESS;
}

void BaseTupleMaker :: fillGenHistograms(){
  if(! m_tuple->eve_isMC){
    cout<<" Calling fillGenHistograms() for Data. "<<endl;
    return;
  }

  if(m_tuple->truth_Zmm>-1) HGen_Zmm_m->Fill(m_tuple->truth_m[m_tuple->truth_Zmm]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zmm>-1) HGen_Zmm_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Zmm]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zmm>-1) HGen_Zmm_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Zmm],m_tuple->eve_mc_w);
  if(m_tuple->truth_Zmm_m1>-1) HGen_Zmm_m1_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Zmm_m1]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zmm_m1>-1) HGen_Zmm_m1_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Zmm_m1],m_tuple->eve_mc_w);
  if(m_tuple->truth_Zmm_m2>-1) HGen_Zmm_m2_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Zmm_m2]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zmm_m2>-1) HGen_Zmm_m2_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Zmm_m2],m_tuple->eve_mc_w);

  //for sherpa samples need to calculate the parent quantities from the decay products
  if(m_tuple->truth_Zmm == -1 && m_tuple->truth_Zmm_m1>-1 && m_tuple->truth_Zmm_m2>-1 ){
    TLorentzVector L1; 
    L1.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Zmm_m1],
		    m_tuple->truth_eta[m_tuple->truth_Zmm_m1],
		    m_tuple->truth_phi[m_tuple->truth_Zmm_m1],
		    m_tuple->truth_m[m_tuple->truth_Zmm_m1]);

    TLorentzVector L2; 
    L2.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Zmm_m2],
		    m_tuple->truth_eta[m_tuple->truth_Zmm_m2],
		    m_tuple->truth_phi[m_tuple->truth_Zmm_m2],
		    m_tuple->truth_m[m_tuple->truth_Zmm_m2]);
    
    HGen_Zmm_m->Fill((L1+L2).M()/1000,m_tuple->eve_mc_w);
    HGen_Zmm_pt->Fill((L1+L2).Pt()/1000,m_tuple->eve_mc_w);
    if((L1+L2).Pt()>1.)
      HGen_Zmm_eta->Fill((L1+L2).Eta(),m_tuple->eve_mc_w);
    else HGen_Zmm_eta->Fill(-999,m_tuple->eve_mc_w);
  }



  if(m_tuple->truth_Zee>-1) HGen_Zee_m->Fill(m_tuple->truth_m[m_tuple->truth_Zee]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zee>-1) HGen_Zee_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Zee]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zee>-1) HGen_Zee_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Zee],m_tuple->eve_mc_w);
  if(m_tuple->truth_Zee_e1>-1) HGen_Zee_e1_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Zee_e1]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zee_e1>-1) HGen_Zee_e1_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Zee_e1],m_tuple->eve_mc_w);
  if(m_tuple->truth_Zee_e2>-1) HGen_Zee_e2_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Zee_e2]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Zee_e2>-1) HGen_Zee_e2_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Zee_e2],m_tuple->eve_mc_w);

  if(m_tuple->truth_Zee == -1 && m_tuple->truth_Zee_e1>-1 && m_tuple->truth_Zee_e2>-1 ){
    TLorentzVector L1; 
    L1.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Zee_e1],
		    m_tuple->truth_eta[m_tuple->truth_Zee_e1],
		    m_tuple->truth_phi[m_tuple->truth_Zee_e1],
		    m_tuple->truth_m[m_tuple->truth_Zee_e1]);

    TLorentzVector L2; 
    L2.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Zee_e2],
		    m_tuple->truth_eta[m_tuple->truth_Zee_e2],
		    m_tuple->truth_phi[m_tuple->truth_Zee_e2],
		    m_tuple->truth_m[m_tuple->truth_Zee_e2]);
    
    HGen_Zee_m->Fill((L1+L2).M()/1000,m_tuple->eve_mc_w);
    HGen_Zee_pt->Fill((L1+L2).Pt()/1000,m_tuple->eve_mc_w);
    if((L1+L2).Pt()>1.)
      HGen_Zee_eta->Fill((L1+L2).Eta(),m_tuple->eve_mc_w);
    else HGen_Zee_eta->Fill(-999,m_tuple->eve_mc_w);
  }


  if(m_tuple->truth_Ztt>-1) HGen_Ztt_m->Fill(m_tuple->truth_m[m_tuple->truth_Ztt]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Ztt>-1) HGen_Ztt_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Ztt]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Ztt>-1) HGen_Ztt_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Ztt],m_tuple->eve_mc_w);
  if(m_tuple->truth_Ztt_t1>-1) HGen_Ztt_t1_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Ztt_t1]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Ztt_t1>-1) HGen_Ztt_t1_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Ztt_t1],m_tuple->eve_mc_w);
  if(m_tuple->truth_Ztt_t2>-1) HGen_Ztt_t2_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Ztt_t2]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Ztt_t2>-1) HGen_Ztt_t2_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Ztt_t2],m_tuple->eve_mc_w);

  //for sherpa samples need to calculate the parent quantities from the decay products
  if(m_tuple->truth_Ztt == -1 && m_tuple->truth_Ztt_t1>-1 && m_tuple->truth_Ztt_t2>-1 ){
    TLorentzVector L1; 
    L1.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Ztt_t1],
		    m_tuple->truth_eta[m_tuple->truth_Ztt_t1],
		    m_tuple->truth_phi[m_tuple->truth_Ztt_t1],
		    m_tuple->truth_m[m_tuple->truth_Ztt_t1]);

    TLorentzVector L2; 
    L2.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Ztt_t2],
		    m_tuple->truth_eta[m_tuple->truth_Ztt_t2],
		    m_tuple->truth_phi[m_tuple->truth_Ztt_t2],
		    m_tuple->truth_m[m_tuple->truth_Ztt_t2]);
    
    HGen_Ztt_m->Fill((L1+L2).M()/1000,m_tuple->eve_mc_w);
    HGen_Ztt_pt->Fill((L1+L2).Pt()/1000,m_tuple->eve_mc_w);
    if((L1+L2).Pt()>1.)
      HGen_Ztt_eta->Fill((L1+L2).Eta(),m_tuple->eve_mc_w);
    else HGen_Ztt_eta->Fill(-999,m_tuple->eve_mc_w);
  }


  if(m_tuple->truth_Wmv>-1) HGen_Wmv_m->Fill(m_tuple->truth_m[m_tuple->truth_Wmv]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wmv>-1) HGen_Wmv_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wmv]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wmv>-1) HGen_Wmv_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wmv],m_tuple->eve_mc_w);
  if(m_tuple->truth_Wmv_m>-1) HGen_Wmv_m_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wmv_m]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wmv_m>-1) HGen_Wmv_m_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wmv_m],m_tuple->eve_mc_w);
  if(m_tuple->truth_Wmv_v>-1) HGen_Wmv_v_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wmv_v]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wmv_v>-1) HGen_Wmv_v_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wmv_v],m_tuple->eve_mc_w);

  if(m_tuple->truth_Wmv == -1 && m_tuple->truth_Wmv_m>-1 && m_tuple->truth_Wmv_v>-1 ){
    TLorentzVector L1; 
    L1.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Wmv_m],
		    m_tuple->truth_eta[m_tuple->truth_Wmv_m],
		    m_tuple->truth_phi[m_tuple->truth_Wmv_m],
		    m_tuple->truth_m[m_tuple->truth_Wmv_m]);

    TLorentzVector L2; 
    L2.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Wmv_v],
		    m_tuple->truth_eta[m_tuple->truth_Wmv_v],
		    m_tuple->truth_phi[m_tuple->truth_Wmv_v],
		    m_tuple->truth_m[m_tuple->truth_Wmv_v]);
    
    HGen_Wmv_m->Fill((L1+L2).M()/1000,m_tuple->eve_mc_w);
    HGen_Wmv_pt->Fill((L1+L2).Pt()/1000,m_tuple->eve_mc_w);
    if((L1+L2).Pt()>1.)
      HGen_Wmv_eta->Fill((L1+L2).Eta(),m_tuple->eve_mc_w);
    else HGen_Wmv_eta->Fill(-999,m_tuple->eve_mc_w);
  }


  if(m_tuple->truth_Wev>-1) HGen_Wev_m->Fill(m_tuple->truth_m[m_tuple->truth_Wev]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wev>-1) HGen_Wev_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wev]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wev>-1) HGen_Wev_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wev],m_tuple->eve_mc_w);
  if(m_tuple->truth_Wev_e>-1) HGen_Wev_e_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wev_e]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wev_e>-1) HGen_Wev_e_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wev_e],m_tuple->eve_mc_w);
  if(m_tuple->truth_Wev_v>-1) HGen_Wev_v_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wev_v]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wev_v>-1) HGen_Wev_v_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wev_v],m_tuple->eve_mc_w);

  if(m_tuple->truth_Wev == -1 && m_tuple->truth_Wev_e>-1 && m_tuple->truth_Wev_v>-1 ){
    TLorentzVector L1; 
    L1.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Wev_e],
		    m_tuple->truth_eta[m_tuple->truth_Wev_e],
		    m_tuple->truth_phi[m_tuple->truth_Wev_e],
		    m_tuple->truth_m[m_tuple->truth_Wev_e]);

    TLorentzVector L2; 
    L2.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Wev_v],
		    m_tuple->truth_eta[m_tuple->truth_Wev_v],
		    m_tuple->truth_phi[m_tuple->truth_Wev_v],
		    m_tuple->truth_m[m_tuple->truth_Wev_v]);
    
    HGen_Wev_m->Fill((L1+L2).M()/1000,m_tuple->eve_mc_w);
    HGen_Wev_pt->Fill((L1+L2).Pt()/1000,m_tuple->eve_mc_w);
    if((L1+L2).Pt()>1.)
      HGen_Wev_eta->Fill((L1+L2).Eta(),m_tuple->eve_mc_w);
    else HGen_Wev_eta->Fill(-999,m_tuple->eve_mc_w);
  }


  if(m_tuple->truth_Wtv>-1) HGen_Wtv_m->Fill(m_tuple->truth_m[m_tuple->truth_Wtv]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wtv>-1) HGen_Wtv_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wtv]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wtv>-1) HGen_Wtv_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wtv],m_tuple->eve_mc_w);
  if(m_tuple->truth_Wtv_t>-1) HGen_Wtv_t_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wtv_t]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wtv_t>-1) HGen_Wtv_t_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wtv_t],m_tuple->eve_mc_w);
  if(m_tuple->truth_Wtv_v>-1) HGen_Wtv_v_pt->Fill(m_tuple->truth_pt[m_tuple->truth_Wtv_v]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_Wtv_v>-1) HGen_Wtv_v_eta->Fill(m_tuple->truth_eta[m_tuple->truth_Wtv_v],m_tuple->eve_mc_w);

  if(m_tuple->truth_Wtv == -1 && m_tuple->truth_Wtv_t>-1 && m_tuple->truth_Wtv_v>-1 ){
    TLorentzVector L1; 
    L1.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Wtv_t],
		    m_tuple->truth_eta[m_tuple->truth_Wtv_t],
		    m_tuple->truth_phi[m_tuple->truth_Wtv_t],
		    m_tuple->truth_m[m_tuple->truth_Wtv_t]);

    TLorentzVector L2; 
    L2.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_Wtv_v],
		    m_tuple->truth_eta[m_tuple->truth_Wtv_v],
		    m_tuple->truth_phi[m_tuple->truth_Wtv_v],
		    m_tuple->truth_m[m_tuple->truth_Wtv_v]);
    
    HGen_Wtv_m->Fill((L1+L2).M()/1000,m_tuple->eve_mc_w);
    HGen_Wtv_pt->Fill((L1+L2).Pt()/1000,m_tuple->eve_mc_w);
    if((L1+L2).Pt()>1.)
      HGen_Wtv_eta->Fill((L1+L2).Eta(),m_tuple->eve_mc_w);
    else HGen_Wtv_eta->Fill(-999,m_tuple->eve_mc_w);
  }



  if(m_tuple->truth_H0bb>-1) HGen_H0bb_m->Fill(m_tuple->truth_m[m_tuple->truth_H0bb]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_H0bb>-1) HGen_H0bb_pt->Fill(m_tuple->truth_pt[m_tuple->truth_H0bb]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_H0bb>-1) HGen_H0bb_eta->Fill(m_tuple->truth_eta[m_tuple->truth_H0bb],m_tuple->eve_mc_w);
  if(m_tuple->truth_H0bb_b1>-1) HGen_H0bb_b1_pt->Fill(m_tuple->truth_pt[m_tuple->truth_H0bb_b1]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_H0bb_b1>-1) HGen_H0bb_b1_eta->Fill(m_tuple->truth_eta[m_tuple->truth_H0bb_b1],m_tuple->eve_mc_w);
  if(m_tuple->truth_H0bb_b2>-1) HGen_H0bb_b2_pt->Fill(m_tuple->truth_pt[m_tuple->truth_H0bb_b2]/1000,m_tuple->eve_mc_w);
  if(m_tuple->truth_H0bb_b2>-1) HGen_H0bb_b2_eta->Fill(m_tuple->truth_eta[m_tuple->truth_H0bb_b2],m_tuple->eve_mc_w);

  if(m_tuple->truth_H0bb == -1 && m_tuple->truth_H0bb_b1>-1 && m_tuple->truth_H0bb_b2>-1 ){
    TLorentzVector L1; 
    L1.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_H0bb_b1],
		    m_tuple->truth_eta[m_tuple->truth_H0bb_b1],
		    m_tuple->truth_phi[m_tuple->truth_H0bb_b1],
		    m_tuple->truth_m[m_tuple->truth_H0bb_b1]);

    TLorentzVector L2; 
    L2.SetPtEtaPhiM(m_tuple->truth_pt[m_tuple->truth_H0bb_b2],
		    m_tuple->truth_eta[m_tuple->truth_H0bb_b2],
		    m_tuple->truth_phi[m_tuple->truth_H0bb_b2],
		    m_tuple->truth_m[m_tuple->truth_H0bb_b2]);
    
    HGen_H0bb_m->Fill((L1+L2).M()/1000,m_tuple->eve_mc_w);
    HGen_H0bb_pt->Fill((L1+L2).Pt()/1000,m_tuple->eve_mc_w);
    if((L1+L2).Pt()>1.)
      HGen_H0bb_eta->Fill((L1+L2).Eta(),m_tuple->eve_mc_w);
    else HGen_H0bb_eta->Fill(-999,m_tuple->eve_mc_w);
  }
  
  if(m_tuple->truth_ZHmmbb>-1){
    HGen_ZH_m->Fill(m_tuple->truth_m[m_tuple->truth_ZHmmbb]/1000,m_tuple->eve_mc_w);
    HGen_ZH_pt->Fill(m_tuple->truth_pt[m_tuple->truth_ZHmmbb]/1000,m_tuple->eve_mc_w);
    HGen_ZH_eta->Fill(m_tuple->truth_eta[m_tuple->truth_ZHmmbb],m_tuple->eve_mc_w);
  }else if(m_tuple->truth_ZHeebb>-1){
    HGen_ZH_m->Fill(m_tuple->truth_m[m_tuple->truth_ZHeebb]/1000,m_tuple->eve_mc_w);
    HGen_ZH_pt->Fill(m_tuple->truth_pt[m_tuple->truth_ZHeebb]/1000,m_tuple->eve_mc_w);
    HGen_ZH_eta->Fill(m_tuple->truth_eta[m_tuple->truth_ZHeebb],m_tuple->eve_mc_w);
  }

  if(m_tuple->truth_WHmvbb>-1){
    HGen_WH_m->Fill(m_tuple->truth_m[m_tuple->truth_WHmvbb]/1000,m_tuple->eve_mc_w);
    HGen_WH_pt->Fill(m_tuple->truth_pt[m_tuple->truth_WHmvbb]/1000,m_tuple->eve_mc_w);
    HGen_WH_eta->Fill(m_tuple->truth_eta[m_tuple->truth_WHmvbb],m_tuple->eve_mc_w);
  }else if(m_tuple->truth_WHevbb>-1){
    HGen_WH_m->Fill(m_tuple->truth_m[m_tuple->truth_WHevbb]/1000,m_tuple->eve_mc_w);
    HGen_WH_pt->Fill(m_tuple->truth_pt[m_tuple->truth_WHevbb]/1000,m_tuple->eve_mc_w);
    HGen_WH_eta->Fill(m_tuple->truth_eta[m_tuple->truth_WHevbb],m_tuple->eve_mc_w);
  }



  ///Resonance analysis
  if(m_tuple->truth_RToVH>-1){
    HGen_RToVH_m->Fill(m_tuple->truth_m[m_tuple->truth_RToVH]/1000,m_tuple->eve_mc_w);
    HGen_RToVH_pt->Fill(m_tuple->truth_pt[m_tuple->truth_RToVH]/1000,m_tuple->eve_mc_w);
    HGen_RToVH_eta->Fill(m_tuple->truth_eta[m_tuple->truth_RToVH],m_tuple->eve_mc_w);
    

    if(m_tuple->truth_RToVH_V>-1){
      HGen_RToVH_V_m->Fill(m_tuple->truth_m[m_tuple->truth_RToVH_V]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_V_pt->Fill(m_tuple->truth_pt[m_tuple->truth_RToVH_V]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_V_eta->Fill(m_tuple->truth_eta[m_tuple->truth_RToVH_V],m_tuple->eve_mc_w);
    }
    if(m_tuple->truth_RToVH_H>-1){
      HGen_RToVH_H_m->Fill(m_tuple->truth_m[m_tuple->truth_RToVH_H]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_H_pt->Fill(m_tuple->truth_pt[m_tuple->truth_RToVH_H]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_H_eta->Fill(m_tuple->truth_eta[m_tuple->truth_RToVH_H],m_tuple->eve_mc_w);
      HGen_RToVH_H_p->Fill(m_tuple->truth_p[m_tuple->truth_RToVH_H]/1000,m_tuple->eve_mc_w);
    }
    if(m_tuple->truth_RToVH_H_b1>-1){
      HGen_RToVH_H_b1_m->Fill(m_tuple->truth_m[m_tuple->truth_RToVH_H_b1]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_H_b1_pt->Fill(m_tuple->truth_pt[m_tuple->truth_RToVH_H_b1]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_H_b1_eta->Fill(m_tuple->truth_eta[m_tuple->truth_RToVH_H_b1],m_tuple->eve_mc_w);
    }
    if(m_tuple->truth_RToVH_H_b2>-1){
      HGen_RToVH_H_b2_m->Fill(m_tuple->truth_m[m_tuple->truth_RToVH_H_b2]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_H_b2_pt->Fill(m_tuple->truth_pt[m_tuple->truth_RToVH_H_b2]/1000,m_tuple->eve_mc_w);
      HGen_RToVH_H_b2_eta->Fill(m_tuple->truth_eta[m_tuple->truth_RToVH_H_b2],m_tuple->eve_mc_w);
    }

    if(m_tuple->truth_RToVH_H_b1>-1 && m_tuple->truth_RToVH_H_b2>-1){
      HGen_RToVH_H_dR->Fill(m_tuple->truth_RToVH_H_dR,m_tuple->eve_mc_w);
      HGen_RToVH_H_dRVsp->Fill(m_tuple->truth_p[m_tuple->truth_RToVH_H]/1000, m_tuple->truth_RToVH_H_dR, m_tuple->eve_mc_w);
    }
    
  }


}




void BaseTupleMaker :: printMCTruth(){

  const xAOD::TruthParticleContainer* Truth = 0;
  if( ! m_event->retrieve( Truth , m_MCTruthIn.c_str() ).isSuccess() ){
    Error("execute()", "Failed to retrieve  TruthParticle collection. Exiting." );
    return;
  }

  cout<<"MC Truth List : List Size = "<<Truth->size()<<endl;
  for (unsigned int i = 0 ; i < Truth->size(); ++i) {
    const xAOD::TruthParticle * part = Truth->at(i);
    cout<<"Status="<<part->status() 
	<<",  PDG="<<part->pdgId()
	<<",  M="<<part->p4().M()
	<<",  Pt="<<part->p4().Pt();
    if(part->p4().Pt()>1)cout<<",  Eta="<<part->p4().Eta();
    cout<<",  Phi="<<part->p4().Phi()
	<<endl;
  }  

}



///////////////////
EL::StatusCode BaseTupleMaker :: postExecute ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode BaseTupleMaker :: finalize ()
{  
  printCounters();
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode BaseTupleMaker :: histFinalize ()
{
  //TFile * f = wk()->getOutputFile(m_filename);
  //f->ls();
  return EL::StatusCode::SUCCESS;
}

//////////////////////////////////
void  BaseTupleMaker::incrementCounter(const char * name, float w){
  if(  m_currentSyst.compare("Nominal") !=0 && m_currentSyst.compare("none") !=0) return;//only use counter for nominal

  //check if counter is already registered:
  bool registered=0;
  unsigned int i=0;
  for(;i<counterNames.size();i++){
    if(counterNames[i]==name){
      counters[i]+=w;
      registered=1;
      break;
    }
  }
  if(!registered){
    counterNames.push_back(name);
    float * counter=new float(w);
    counters.push_back(*counter);
  }

}

float  BaseTupleMaker::getCounter(const char * name){  
  for(unsigned int i=0;i<counterNames.size();i++){
    if(counterNames[i]==name)
      return counters[i];
  }
  return 0.;
}

void BaseTupleMaker::printCounters(){
  for(unsigned int i=0;i<counterNames.size();i++){
    cout<<counterNames[i]<<" = "<<counters[i]<<endl;
  }
}


//////////////

EL::StatusCode BaseTupleMaker :: initPileUpWeights(TFile * f){

  ///////////////////////////////////
  ///Save the pile up distribution for each processed sample ( this is not used for the weights)
  //////////////////////////////////
  HPileUp=new TH1F("HPileUp","HPileUp",100,0,100);
  HPileUp->SetDirectory(f);

  //////////////////////////////////////////////////////
  //Create the ratio histograms to be used for the weights
  //There will be one weight for every combination of Data and MC distribution
  ////////////////////////////////////////////////////////
  //now fill weight the histograms
  m_tuple->nmcpuw=0;
  for(unsigned int i=0;i<pileUpMCFiles.size();i++){
    TFile FMC(pileUpMCFiles[i].c_str(),"READ");
    TH1F* HMC=(TH1F*)FMC.Get("pileup");
    if(!HMC){
      Error("initPileUpWeights()", "No MC PileUp Histo found");
      return EL::StatusCode::FAILURE;
    }

    //check the histogram:
    for(int b=1;b<=HMC->GetNbinsX();b++)
      if(HMC->GetBinContent(b) < 0. //can't be negative
	 || HMC->GetBinContent(b) > 1.0 //can't be >1 because histo integral is 1.
	 ){
	cout<<"Pileup bin content is bad in "<<pileUpMCFiles[i].c_str()<<endl;
      }

    for(unsigned int j=0;j<pileUpDataFiles.size();j++){
      TFile FData(pileUpDataFiles[j].c_str(),"READ");
      TH1F* HData=(TH1F*)FData.Get("pileup");
      if(!HData){
	Error("initPileUpWeights()", "No Data PileUp Histo found");
	return EL::StatusCode::FAILURE;
      }


      if(HData->GetNbinsX()!=HMC->GetNbinsX()){
	Error("initPileUpWeights()", "PileUp Bins don't match");
	return EL::StatusCode::FAILURE;
      }

      //check the histogram:
      for(int b=1;b<=HData->GetNbinsX();b++)
	if(HData->GetBinContent(b)<0. //can't be negative
	   || HData->GetBinContent(b) >1.0 //can't be >1 because histo integral is 1.
	   ){
	  Error("initPileUpWeights()", "Pileup bin contents are bad");
	  return EL::StatusCode::FAILURE;
	}


      //Create the ratios and link them to the output file
      PileUpRatios[m_tuple->nmcpuw]=new TH1F(std::string("PileUpRatios")+(long)(m_tuple->nmcpuw),"pileup",
					     HPileUp->GetNbinsX(),
					     HPileUp->GetXaxis()->GetXmin(),
					     HPileUp->GetXaxis()->GetXmax());
      PileUpRatios[m_tuple->nmcpuw]->SetDirectory(f);

      //Calculate the ratio Data/MC
      //Also calculate the MC integral after reweighting
      float mcIntegral=0.;
      for(int b=1;b<=HMC->GetNbinsX();b++){
	if(HMC->GetBinContent(b)<=0. || HData->GetBinContent(b)<=0.)
	  PileUpRatios[m_tuple->nmcpuw]->SetBinContent(b,0.);
	else PileUpRatios[m_tuple->nmcpuw]->SetBinContent(b,HData->GetBinContent(b)/HMC->GetBinContent(b));
	
	mcIntegral+= PileUpRatios[m_tuple->nmcpuw]->GetBinContent(b) * HMC->GetBinContent(b) ;
      }

      //rescale ratio to correct MC integral
      if(mcIntegral>0.) PileUpRatios[m_tuple->nmcpuw]->Scale(1./mcIntegral);
      else {
	Error("initPileUpWeights()", "Failed to rescale PileUpRatio");
	return EL::StatusCode::FAILURE;
      }

      cout<<"PileUpWeight: "<<m_tuple->nmcpuw<<"  "<<pileUpMCFiles[i].c_str()<<" "<<pileUpDataFiles[j].c_str()<<endl;
      
      m_tuple->nmcpuw++;
      FData.Close();
    }
    FMC.Close();
  }

  return EL::StatusCode::SUCCESS;
}
