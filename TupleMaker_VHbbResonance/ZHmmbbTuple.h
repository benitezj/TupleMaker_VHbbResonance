#ifndef ZHMMBBTUPLE
#define ZHMMBBTUPLE

#include "ZHllbbTuple.h"

class ZHmmbbTuple : public ZHllbbTuple {
 public:
  ZHmmbbTuple();
  ~ZHmmbbTuple();
  
  //This is the top level.
  //Any variables specific to Z-->mu mu here:

  virtual void DefineBranches(TTree * tr); 

 private:

};
#endif
