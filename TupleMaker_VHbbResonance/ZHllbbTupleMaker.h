#ifndef ZHLLBBTUPLEMAKER_H
#define ZHLLBBTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/VHTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHllbbTuple.h"

class ZHllbbTupleMaker : public VHTupleMaker
{
public:
  ZHllbbTupleMaker ();//Default Do not use.
  ZHllbbTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(ZHllbbTupleMaker, 1);

protected:
 
  ZHllbbTuple * m_tuple = 0 ; //!
  
  TLorentzVector l1p4; //4-vector of selected lepton
  TLorentzVector l2p4;

  virtual EL::StatusCode processEvent();

  bool jetOverlap(int jetindex); 
  bool fatJetOverlap(int jetindex); 
  bool trkJetOverlap(int jetindex); 

    
private:


  
};

#endif
