#ifndef VHTUPLEMAKER_H
#define VHTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/RecoTupleMaker.h"
#include "TupleMaker_VHbbResonance/VHTuple.h"

class VHTupleMaker : public RecoTupleMaker
{
public:
  VHTupleMaker ();//Default Do not use.
  VHTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(VHTupleMaker, 1);

protected:
 
  VHTuple * m_tuple = 0 ; //!

  int leptonIsoOption;//for leptons

  //resolved analysis
  int jet1; //signal jets pT ordered 
  int jet2; //
  int bTagOption;

  virtual EL::StatusCode processEvent();

  //Overlap check with the leptons
  virtual bool jetOverlap(int jetindex) = 0; //to be defined in higher class
  virtual bool fatJetOverlap(int jetindex) = 0; //to be defined in higher class
  virtual bool trkJetOverlap(int jetindex) = 0; //to be defined in higher class


  //Resolved selection
  EL::StatusCode  selectJets () ;
  std::vector<int> fillCentralJets(); 
  void fillCentralJetsVeto(); 
  void fillFWJets();
  void fillBJets(); 


  //merged analysis selection
  float fatJetPtCut;  
  int higgsJet; 
  EL::StatusCode  selectHiggs() ;
  std::vector<int> fillFatJets(); 
  void fillAddBJets(); 
  bool jetOverlapHiggs(int jetindex);
  bool trkJetLinkedHiggs(int jetindex);


  ///
  void PrintVHFlags(){
    std::cout<<"PrintVHFlags: passVHMultiLepVeto="<<m_tuple->eve_passVHMultiLepVeto<<" nAddTrkBjet="<<m_tuple->vh_naddtrkbjet<<std::endl;
  }


private:


  //this has an additional 30GeV cut
  bool passCentralJetVH(int i){
    if(!passCentralJet(i)) return false;
    if( m_tuple->jet_pt[i] < 30000 ) return false;
    return true;
  } 


  
};

#endif
