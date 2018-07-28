#ifndef ZEETUPLE
#define ZEETUPLE

#include "ZllTuple.h"

class ZeeTuple : public ZllTuple {
 public:
  ZeeTuple();
  ~ZeeTuple();
  
  //This is the top level.
  //Any variables specific to Z-->mu mu here:

  virtual void DefineBranches(TTree * tr); 

  //void SetBranchAddresses(TTree * tr);
  
 private:

  //template<class T> void SetBranchAddressSafe(TTree * tr, const char * label, T * var);
  
};
#endif 
