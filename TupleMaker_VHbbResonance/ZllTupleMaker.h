#ifndef ZLLTUPLEMAKER_H
#define ZLLTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/RecoTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZllTuple.h"

class ZllTupleMaker : public RecoTupleMaker
{
public:
  ZllTupleMaker ();//Default Do not use.
  ZllTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(ZllTupleMaker, 1);

protected:
 
  ZllTuple * m_tuple = 0 ; //!
  
  TLorentzVector l1p4;
  TLorentzVector l2p4;

  int leptonIsoOption;

  virtual EL::StatusCode processEvent();

  EL::StatusCode fillJets ();//must be called after filling leptons
  
private:


  
};

#endif
