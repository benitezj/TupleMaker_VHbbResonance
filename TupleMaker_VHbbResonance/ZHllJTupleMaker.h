#ifndef ZHLLJTUPLEMAKER_H
#define ZHLLJTUPLEMAKER_H

#include "TupleMaker_VHbbResonance/VHTupleMaker.h"
#include "TupleMaker_VHbbResonance/ZHllJTuple.h"

class ZHllJTupleMaker : public VHTupleMaker
{
public:
  ZHllJTupleMaker ();//Default Do not use.
  ZHllJTupleMaker (std::string configPath);

  virtual EL::StatusCode initialize ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(ZHllJTupleMaker, 1);

protected:
 
  ZHllJTuple * m_tuple = 0 ; //!
  
  TLorentzVector l1p4; //4-vector of selected lepton
  TLorentzVector l2p4;

  virtual EL::StatusCode processEvent();

  bool jetOverlap(int jetindex); 
  bool fatJetOverlap(int jetindex); 
  bool trkJetOverlap(int jetindex); 

    
private:


  
};

#endif
