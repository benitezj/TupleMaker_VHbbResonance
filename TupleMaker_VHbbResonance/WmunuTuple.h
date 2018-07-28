#ifndef WMUNUTUPLE
#define WMUNUTUPLE

#include "WlnuTuple.h"

class WmunuTuple : public WlnuTuple {
 public:
  WmunuTuple();
  ~WmunuTuple();


  virtual void DefineBranches(TTree * tr); 

  //void SetBranchAddresses(TTree * tr);
  
 private:

  //template<class T> void SetBranchAddressSafe(TTree * tr, const char * label, T * var);
  
};
#endif 
