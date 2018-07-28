#include "TupleMaker_VHbbResonance/ZHeeJTuple.h"

ZHeeJTuple::ZHeeJTuple():
  ZHllJTuple()
{
}

ZHeeJTuple::~ZHeeJTuple() {;}

void ZHeeJTuple::DefineBranches(TTree * tr) {
  
  ZHllJTuple::DefineBranches(tr);
  
}



