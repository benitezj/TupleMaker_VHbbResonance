#include "TupleMaker_VHbbResonance/ZllTuple.h"

ZllTuple::ZllTuple():
  VTuple()
{
}

ZllTuple::~ZllTuple() {;}

void ZllTuple::DefineBranches(TTree * tr) {
  VTuple::DefineBranches(tr);

  tr->Branch("zll_ll",&zll_ll,"zll_ll/I");
  tr->Branch("zll_l1",&zll_l1,"zll_l1/I");
  tr->Branch("zll_l2",&zll_l2,"zll_l2/I");

}
