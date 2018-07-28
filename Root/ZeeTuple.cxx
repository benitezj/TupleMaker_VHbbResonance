#include "TupleMaker_VHbbResonance/ZeeTuple.h"

ZeeTuple::ZeeTuple():
  ZllTuple()
{
}

ZeeTuple::~ZeeTuple() {;}

void ZeeTuple::DefineBranches(TTree * tr) {
  ZllTuple::DefineBranches(tr);

}


// template <class T> void ZeeTuple::SetBranchAddressSafe(TTree * tr, const char * label, T * var) {
//   if( tr->GetBranch(label) ) {
//     tr->SetBranchAddress(label,var);
//   }
// }
// void ZeeTuple::SetBranchAddresses(TTree * tr) {
//   SetBranchAddressSafe(tr,"eve_num",&eve_num);
// }


