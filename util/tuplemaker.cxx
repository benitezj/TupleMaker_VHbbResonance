#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/DiskListEOS.h"
#include "SampleHandler/SampleLocal.h"

#include <TRegexp.h>
#include <RootCoreUtils/StringUtil.h>

#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
//#include "EventLoop/ProofDriver.h"
//#include "EventLoop/LSFDriver.h"
//#include "EventLoopGrid/PrunDriver.h"

#include <stdlib.h> 

#include "CxAODTools/ConfigStore.h"

#include "TupleMaker_VHbbResonance/BaseTupleMaker.h"
#include "TupleMaker_VHbbResonance/RecoTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZmmTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZeeTupleMaker.h"
#include "TupleMaker_VHbbResonance/WmunuTupleMaker.h"
#include "TupleMaker_VHbbResonance/WenuTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHmmbbTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHeebbTupleMaker.h"
#include "TupleMaker_VHbbResonance/WHmnubbTupleMaker.h"
#include "TupleMaker_VHbbResonance/WHenubbTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHmmJTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHeeJTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHemJTupleMaker.h"


using namespace std;

int main(int argc, char* argv[]) {

  if(argc < 2 ){ 
    std::cout<<"You need to provide an option"<<std::endl;
    exit(0);
  }
  std::string action = argv[1];

  std::string inPath(getenv("TUPLEINPATH"));
  cout<<"TUPLEINPATH  "<<inPath.c_str()<<endl;

  std::string sampleName(getenv("TUPLESAMPLENAME"));
  cout<<"TUPLESAMPLENAME  "<<sampleName.c_str()<<endl;

  std::string outPath(getenv("TUPLEOUTPATH"));
  cout<<"TUPLEOUTPATH  "<<outPath.c_str()<<endl;

  std::string tupleConfig(getenv("TUPLECONFIG"));
  cout<<"TUPLECONFIG  "<<tupleConfig.c_str()<<endl;

  std::string channel(getenv("TUPLECHANNEL"));
  cout<<"TUPLECHANNEL  "<<channel.c_str()<<endl;

  int filesperjob = atoi(getenv("TUPLENFILES"));
  cout<<"TUPLENFILES  "<<filesperjob<<endl;


  // Construct the sample to run on:
  SH::SampleHandler sampleHandler;

  if(action == "local" ){
    //SH::DiskListLocal list(inPath + "/" + sampleName);
    //SH::scanSingleDir(sampleHandler,"tuple",list,"*outputLabel_*.root*");

    SH::SampleLocal* sample = new SH::SampleLocal("tuple");
    SH::DiskListLocal list(inPath + "/" + sampleName);
    int counter=0;
    while(list.next() && ( counter<filesperjob || filesperjob == -1 ) ){
      //if(RCU::match_expr(RCU::glob_to_regexp("*.root*"),list.fileName().c_str())){///does not compile after moving to 2.3.42
      if(list.fileName().find(".root")!=list.fileName().npos){
	TFile File((inPath + "/" + sampleName+"/"+list.fileName()).c_str());
	if(File.Get("CollectionTree") != NULL){
	  sample->add(list.path());
	  counter++;
	}
      }
    }
    if(counter>0){
      sampleHandler.add(sample);
    }else{
      cout<<"Could not find any root files: "<<inPath + "/" + sampleName<<endl;
      return 0;
    }
  }else{
    cout<<"Invalid action required."<<endl;
    return 0;
  }
  
  sampleHandler.setMetaString("nc_tree", "CollectionTree");
  sampleHandler.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler(sampleHandler);
  //job.options()->setDouble (EL::Job::optMaxEvents,-1);

  // create algorithm, set job options, and add our analysis to the job:
  BaseTupleMaker* tupleMaker =0;
  if(channel=="base")tupleMaker = new BaseTupleMaker(tupleConfig);
  if(channel=="reco")tupleMaker = new RecoTupleMaker(tupleConfig);
  if(channel=="Zmm")tupleMaker = new ZmmTupleMaker(tupleConfig);
  if(channel=="Zee")tupleMaker = new ZeeTupleMaker(tupleConfig);
  if(channel=="Wmunu")tupleMaker = new WmunuTupleMaker(tupleConfig);
  if(channel=="Wenu")tupleMaker = new WenuTupleMaker(tupleConfig);
  if(channel=="ZHmmbb")tupleMaker = new ZHmmbbTupleMaker(tupleConfig);
  if(channel=="ZHeebb")tupleMaker = new ZHeebbTupleMaker(tupleConfig);
  if(channel=="WHmnubb")tupleMaker = new WHmnubbTupleMaker(tupleConfig);
  if(channel=="WHenubb")tupleMaker = new WHenubbTupleMaker(tupleConfig);
  if(channel=="ZHmmJ")tupleMaker = new ZHmmJTupleMaker(tupleConfig);
  if(channel=="ZHeeJ")tupleMaker = new ZHeeJTupleMaker(tupleConfig);
  if(channel=="ZHemJ")tupleMaker = new ZHemJTupleMaker(tupleConfig);

  if(!tupleMaker){cout<<"Wrong channel"<<endl; exit(0);}
  job.algsAdd(tupleMaker);
  
  if(action=="local" ){//LOCAL/EOS in lxplus
    cout<<"Processing locally."<<endl;
    EL::DirectDriver driver;
    driver.submitOnly(job,outPath + "/" + sampleName);
  }

  return 0;
}

