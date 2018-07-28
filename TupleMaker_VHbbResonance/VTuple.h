#ifndef VTUPLE
#define VTUPLE

#include "RecoTuple.h"

class VTuple : public RecoTuple {
 public:
  VTuple();
  ~VTuple();

  int v_njet=0;//# of central jets after overlap cleaning
  int v_jet=0;//index of leading central jet after overlap clening

  int v_njetfw=0;//# of forward jets after overlap cleaning
  int v_jetfw=0;//index of leading fw jet after overlap cleaning

  int v_nfatjet=0;//# of fat jets after overlap cleaning
  int v_fatjet=0;//index of leading fat jet after overlap cleaning

  int v_ntrkjet=0;//# of trk jets after overlap cleaning
  int v_trkjet=0;//index of leading trk jet after overlap cleaning


  virtual void DefineBranches(TTree * tr); 

 private:
  
};
#endif 
