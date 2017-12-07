#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "lardataobj/RecoBase/OpFlash.h"

namespace recob{

  // Note the > instead < so that largest flashes come first
  struct OpFlashSortByPE_t {
    bool operator() (const recob::OpFlash i, const recob::OpFlash j){ return i.TotalPE() > j.TotalPE(); }
  } OpFlashSortByPE;

  struct OpFlashPtrSortByTime_t {
    bool operator() (const art::Ptr<recob::OpFlash> i, const art::Ptr<recob::OpFlash> j){ return i->Time() < j->Time(); }
  } OpFlashPtrSortByTime;

  // Note the > instead < so that largest flashes come first
  struct OpFlashPtrSortByPE_t {
    bool operator() (const art::Ptr<recob::OpFlash> i, const art::Ptr<recob::OpFlash> j){ return i->TotalPE() > j->TotalPE(); }
  } OpFlashPtrSortByPE;

}


