#include "TupleMaker_VHbbResonance/ZHmmJTuple.h"

ZHmmJTuple::ZHmmJTuple():
  ZHllJTuple()
{
}

ZHmmJTuple::~ZHmmJTuple() {;}

void ZHmmJTuple::DefineBranches(TTree * tr) {
  
  ZHllJTuple::DefineBranches(tr);
  
}



