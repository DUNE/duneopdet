// Tyler LaBree
// Northern Illinois University

#include "PhotonLibrary.h"
#include "LightYieldMap.h"
#include "DepthPlot.h"
void MakePlot() {
  Int_t dims[3] = {50,93,137};
  PhotonLibrary* lib = new PhotonLibrary(
    "/dune/app/users/wemark/photon-library/lib_dunevd10kt_1x8x14_3_20221025.root"
    , dims);
  
  // Define depth ranges on which to average over to make light yield plots.
  Int_t centerSlice[2] = {68,70};
  Int_t centerAvg[2] = {60,78};
  Int_t fullAvg[2] = {6,131};

  /* Make a light yield plot perpendicular to the z-dir, averaged over the
   * range 60-78 voxels, using the photon library made above.
   */
  LightYieldMap* lymap = new LightYieldMap(lib, 2, centerAvg, "Average over ");
  lymap->Draw();

  // Draw average and minimum light yield across detector depth.
  DrawDepthPlot(lib, 6, 131);
}
