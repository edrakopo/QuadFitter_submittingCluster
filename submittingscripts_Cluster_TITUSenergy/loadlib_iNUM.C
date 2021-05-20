#include "TROOT.h"
#include "TSystem.h"

#include <iostream>


void loadlib_iNUM()
{
  gSystem->Load("/exports/csce/datastore/ph/groups/PPE/titus/ts-WChRecoSandBox/lib/libWCLAnalysis.so");
  //gSystem->Load("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing/PhotonMask_12in_test_iNUM_C.so");
}
