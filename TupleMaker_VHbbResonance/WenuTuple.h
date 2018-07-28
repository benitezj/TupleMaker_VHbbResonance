#ifndef WENUTUPLE
#define WENUTUPLE

#include "WlnuTuple.h"

class WenuTuple : public WlnuTuple {
 public:
  WenuTuple();
  ~WenuTuple();

  virtual void DefineBranches(TTree * tr); 

  //void SetBranchAddresses(TTree * tr);
  
 private:

  //template<class T> void SetBranchAddressSafe(TTree * tr, const char * label, T * var);
  
};
#endif 
