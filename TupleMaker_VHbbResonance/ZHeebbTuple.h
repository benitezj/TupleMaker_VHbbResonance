#ifndef ZHEEBBTUPLE
#define ZHEEBBTUPLE

#include "ZHllbbTuple.h"

class ZHeebbTuple : public ZHllbbTuple {
 public:
  ZHeebbTuple();
  ~ZHeebbTuple();
  
  //This is the top level.
  //Any variables specific to Z-->ee  here:

  virtual void DefineBranches(TTree * tr); 

 private:

};
#endif
