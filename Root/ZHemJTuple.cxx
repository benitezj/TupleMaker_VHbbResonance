#include "TupleMaker_VHbbResonance/ZHemJTuple.h"

ZHemJTuple::ZHemJTuple():
  ZHllJTuple()
{
}

ZHemJTuple::~ZHemJTuple() {;}

void ZHemJTuple::DefineBranches(TTree * tr) {
  
  ZHllJTuple::DefineBranches(tr);
  
}



