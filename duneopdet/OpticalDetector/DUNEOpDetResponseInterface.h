// \file DUNEOpDetResponseInterface.h
// \brief service containing information about the response of optical detectors for dune specific uses.
//
// \author: JStock (jason.stock@mines.sdsmt.edu)
// Based on OpDetResponseInterface by AHimmel (ahimmel@fnal.gov).

#ifndef DUNE_OPDET_RESPONSE_INTERFACE_H
#define DUNE_OPDET_RESPONSE_INTERFACE_H

#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "dune/OpticalDetector/DUNEOpDetResponse.h"

namespace opdet
{
  class DUNEOpDetResponseInterface : public opdet::OpDetResponseInterface {
    public:
      virtual bool detectedLite(int OpDet, int &newOpChannel, int& hardwareChannel) const
      {
        return doDetectedLite( OpDet, newOpChannel, hardwareChannel);
      }
    private:
      virtual bool doDetectedLite(int OpDet, int &newOpChannel, int& hardwareChannel) const = 0;
  };
}

DECLARE_ART_SERVICE_INTERFACE(opdet::DUNEOpDetResponseInterface, LEGACY)

#endif //DUNE_OPDET_RESPONSE_INTERFACE_H
