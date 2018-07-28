#include "TupleMaker_VHbbResonance/ZmmTuple.h"

ZmmTuple::ZmmTuple():
  ZllTuple()
{
}

ZmmTuple::~ZmmTuple() {;}

void ZmmTuple::DefineBranches(TTree * tr) {
  
  ZllTuple::DefineBranches(tr);
  
}


// template <class T> void ZmmTuple::SetBranchAddressSafe(TTree * tr, const char * label, T * var) {
//   if( tr->GetBranch(label) ) {
//     tr->SetBranchAddress(label,var);
//   }
// }
// void ZmmTuple::SetBranchAddresses(TTree * tr) {
//   SetBranchAddressSafe(tr,"eve_num",&eve_num);
// }


