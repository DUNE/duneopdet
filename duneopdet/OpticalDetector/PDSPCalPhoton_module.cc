////////////////////////////////////////////////////////////////////////
// Class:       PDSPCalPhoton
// Plugin Type: producer (Unknown Unknown)
// File:        PDSPCalPhoton_module.cc
//
// Generated at Tue Aug  9 15:32:19 2022 by Tingjun Yang using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TH1D.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TF1.h"

#include <vector>
#include <memory>
#include <complex>
#include <fftw3.h>

double elecrespfun(double *x, double *par){
  double y = 0;
  if (x[0]<par[3]){
    y = 0;
  }
  else if (x[0]<par[4]){
    y = par[0]*(1-exp(-(x[0]-par[3])/par[1]));
  }
  else{
    y = par[0]*(1-exp(-(par[4]-par[3])/par[1]))*exp(-(x[0]-par[4])/par[2]);
  }
  return y;
}

namespace pdsp {
  class PDSPCalPhoton;
}

class pdsp::PDSPCalPhoton : public art::EDProducer {
public:
  explicit PDSPCalPhoton(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPCalPhoton(PDSPCalPhoton const&) = delete;
  PDSPCalPhoton(PDSPCalPhoton&&) = delete;
  PDSPCalPhoton& operator=(PDSPCalPhoton const&) = delete;
  PDSPCalPhoton& operator=(PDSPCalPhoton&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fRawOpWaveformLabel;
  std::string   fResponseFile;
  std::string   fFiltFunc;
  std::vector<double> fPECalib;
  int ntbin;
  TH1D *wienerfilter;
  TH1D *elecres;
  std::vector<double> vwiener;
  std::vector<std::complex<double>> vinvres;
  TF1 *filtfun;
};

double CalcMedian(std::vector<short int> scores){

  size_t size = scores.size();

  if (size == 0){
    return 0;  // Undefined, really.
  }
  else{
    sort(scores.begin(), scores.end());
    if (size % 2 == 0){
      return (scores[size / 2 - 1] + scores[size / 2]) / 2;
    }
    else {
      return scores[size / 2];
    }
  }
}

pdsp::PDSPCalPhoton::PDSPCalPhoton(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fRawOpWaveformLabel(p.get<art::InputTag>("RawOpWaveformLabel")),
  fResponseFile(p.get<std::string>("ResponseFile")),
  fFiltFunc(p.get<std::string>("FiltFunc")),
  fPECalib(p.get<std::vector<double>>("PECalib")),
  ntbin(2000)
{
  // Wiener filtered waveform
  produces<std::vector<recob::OpWaveform>>("wiener");
  produces<art::Assns<raw::OpDetWaveform, recob::OpWaveform>>("wiener");
  // Apply Gaussian filter and remove electronics response
  produces<std::vector<recob::OpWaveform>>();
  produces<art::Assns<raw::OpDetWaveform, recob::OpWaveform>>();

  // Get response functions
  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fResponseFile, fullname);
  TFile *resfile;
  if (fullname.empty()) {
    throw cet::exception("PDSPCalPhoton") << "Unable to find response file "  << fResponseFile;
  }
  else{
    resfile = TFile::Open(fullname.c_str());
    wienerfilter = (TH1D*)resfile->Get("wienerfilter");
    elecres = (TH1D*)resfile->Get("elecres");
    if (!wienerfilter) throw cet::exception("PDSPCalPhoton") << "Unable to find wienerfilter";
    if (!elecres) throw cet::exception("PDSPCalPhoton") << "Unable to find elecres";
  }

  // Save wiener filter
  for (int i = 0; i<wienerfilter->GetNbinsX(); ++i){
    vwiener.push_back(wienerfilter->GetBinContent(i+1));
  }

  // Save inverse response function
  TF1 *respfun = new TF1("respfun", elecrespfun, 0, ntbin, 5);
  respfun->SetParameter(0, 10.23);
  respfun->SetParameter(1, 8.101);
  respfun->SetParameter(2, 66.43);
  respfun->SetParameter(3, 72.98);
  respfun->SetParameter(4, 94.63);

  fftw_complex *in, *out;
  fftw_plan pf;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  pf = fftw_plan_dft_1d(ntbin, in, out, FFTW_FORWARD, FFTW_MEASURE);
  double resp_integral = 0;
  for (int i = 0; i<ntbin; ++i){
    in[i][0] = respfun->Eval(i);
    in[i][1] = 0;
    resp_integral += in[i][0];
  }
  for (int i = 0; i<ntbin; ++i){
    in[i][0] *= resp_integral;
  }
  fftw_execute(pf);
  for (int i = 0; i<ntbin; ++i){
    std::complex<double> tmp (1,0);
    std::complex<double> res (out[i][0], out[i][1]);
    vinvres.push_back(tmp/res);
    //std::cout<<i<<" "<<real(vinvres.back())<<" "<<imag(vinvres.back())<<std::endl;
  }
  fftw_destroy_plan(pf);
  fftw_free(in); fftw_free(out);
  delete respfun;

  // Construct filter function
  filtfun = new TF1("filtfun",fFiltFunc.c_str(),0,150);

  //double  PEcalib[300]={1565.51, 1572.83, 1568.89, 1597.42, 1566.34, 1567.61, 1796.34, 1577.88, 1539.07, 1554.94, 1545.95, 1568.67, 1559.3, 1523.57, 1545.92, 1545.83, 1537.22, 1530.41, 1559.09, 1708.33, 1532.2, 1542.95, 1540.17, 1837.27, 1542.26, 1, 1547.04, 1525.85, 1797.64, 1562.35, 1542.65, 1500.22, 1558.41, 1697.22, 1536.21, 1529.25, 1, 1522.8, 1529.15, 1543.78, 1, 1, 1, 1, 1, 1, 1, 1, 1557.63, 1836.62, 1542.79, 1855.54, 1536.01, 1563.32, 1554.35, 1572.27, 1533.53, 1537.15, 1, 1524.63, 1564.09, 1554.8, 1, 1357.43, 1557.81, 1, 1538.47, 1563.21, 1553.47, 1557.65, 1535.63, 1350.47, 1555.16, 1, 1553.14, 1, 1564.94, 1717.32, 1568.56, 1559.69, 1558.74, 1529.36, 1, 1522.59, 1552.92, 1554.34, 1557.97, 1561.54, 1, 1, 1, 1, 1, 1, 1, 1, 1511.98, 1558.58, 1539.22, 1531.83, 1540.52, 1, 1557.68, 1554.41, 1719.46, 1554.15, 1711.39, 1710.53, 1540.22, 1558.55, 1, 1532.33, 1340.25, 1556.37, 1547.44, 1552.58, 1553.36, 1555.81, 1562.98, 1, 1720.11, 1543.24, 1568.4, 1553.57, 1682.39, 1862.54, 1736.57, 1685.79, 1558.13, 1546.75, 1554.38, 1559.44, 792.387, 768.789, 786.324, 811.68, 814.341, 672.574, 736.12, 642.818, 647.111, 636.881, 756.772, 770.625, 1560.33, 1543.5, 1575.76, 1541.66, 1537.91, 1562.15, 1580.87, 1541.14, 1529.42, 1553.19, 1549.07, 1544.29, 1, 1562.64, 1571.92, 1555.76, 1, 1555.4, 1566.6, 1543.83, 1530.93, 1540.15, 1345.81, 1552.64, 1547.64, 1561.12, 1564.48, 1553.25, 1738.38, 1586.85, 1561.54, 1726.41, 1842.11, 1556.9, 1531.88, 1788, 1693.92, 1568.82, 1776.21, 1542.83, 1, 1, 1, 1, 1, 1, 1, 1, 1561.53, 1559.44, 1555.41, 1500.87, 1568.88, 1553.48, 1761.53, 1348.08, 1800.21, 1710.85, 1588.96, 1566.39, 1577.34, 1592.33, 1591.09, 1558.41, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1002.26, 1, 1038.31, 996.86, 1005.85, 983.454, 1018.45, 1039.98, 1020.16, 1011.72, 1033.08, 983.217, 1006.96, 1024.45, 1011.36, 984.045, 1008.55, 1028.3, 1002.73, 992.306, 1029.01, 1013.14, 1070.02, 992.292, 1006.06, 978.098, 963.095, 947.012, 973.652, 1008.66, 982.414, 1026.21, 994.26, 1003.63, 1032.23, 767.735, 788.662, 772.74, 633.259, 813.255, 651.201, 687.14, 754.99, 645.796, 770.245, 737.236, 827.26, 988.743, 980.494, 1001.89, 998.85, 992.995, 1025.81, 989.474, 1032.17, 986.536, 980.31, 1000.97, 978.489, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (int i = 0; i<300; ++i){
    std::cout<<i<<" "<<fPECalib[i]<<std::endl;
  }
}

void pdsp::PDSPCalPhoton::produce(art::Event& e)
{
  // Implementation of required member function here.
  std::vector < art::Ptr < raw::OpDetWaveform > > wfList;
  auto wfHandle = e.getHandle<std::vector<raw::OpDetWaveform>>(fRawOpWaveformLabel);
  if (!wfHandle.isValid()){
    throw cet::exception("PDSPCalPhoton") << "Unable to get waveforms using label "  << fRawOpWaveformLabel;
  }
  else{
    art::fill_ptr_vector(wfList, wfHandle);
  }

  auto out_recowaveforms1 = std::make_unique< std::vector< recob::OpWaveform > >();
  auto out_recowaveforms2 = std::make_unique< std::vector< recob::OpWaveform > >();
  auto recorawassn1 = std::make_unique< art::Assns<raw::OpDetWaveform, recob::OpWaveform> >();
  auto recorawassn2 = std::make_unique< art::Assns<raw::OpDetWaveform, recob::OpWaveform> >();

  fftw_complex *in, *out;
  fftw_plan pf;
  fftw_plan pb;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  pf = fftw_plan_dft_1d(ntbin, in, out, FFTW_FORWARD, FFTW_MEASURE);
  pb = fftw_plan_dft_1d(ntbin, in, out, FFTW_BACKWARD, FFTW_MEASURE);

  // Vector after applying Wiener filter
  std::vector<float> wfwiener(ntbin);
  // Vector after applying Gaussian filter and remove electronics response
  std::vector<float> wfimpulse(ntbin);

  for (auto const& wf : wfList) {
    if (int(wf->Waveform().size()) != ntbin){
      throw cet::exception("PDSPCalPhoton") << "Waveform size does not match "  << wf->Waveform().size()<<" "<<ntbin;
    }
    
    // Get baseline
    //double baseline = CalcMedian(wf->Waveform());
    double baseline = 0;
    TH1D *basehelp= new TH1D("basehelp","basehelp",2000, 1300,1800);
    for(size_t j=0; j<wf->Waveform().size(); j++){
      if(j<1000){
        basehelp->Fill(wf->Waveform()[j]);
      }
    }
    int basebinmax = basehelp->GetMaximumBin();
    baseline = basehelp->GetXaxis()->GetBinCenter(basebinmax);
    basehelp->Delete();   

    for (int i = 0; i<ntbin; ++i){
      in[i][0] = wf->Waveform()[i] - baseline;
      in[i][1] = 0;
    }

    fftw_execute(pf);

    // FFT of waveform
    std::vector<std::complex<double>> wffft(wf->Waveform().size());

    for (int i = 0; i<ntbin; ++i){
      double filt = vwiener[i*vwiener.size()/ntbin];
      wffft[i] = std::complex<double>(out[i][0], out[i][1]);
      // Apply Wiener filter
      in[i][0] = out[i][0]*filt;
      in[i][1] = out[i][1]*filt;
    }

    fftw_execute(pb);

    for (int i = 0; i<ntbin; ++i){
      wfwiener[i] = out[i][0]/ntbin;
    }
    recob::OpWaveform out_recowaveFinal(wf->TimeStamp(), wf->ChannelNumber(), wfwiener);
    out_recowaveforms1->emplace_back(std::move(out_recowaveFinal));
    util::CreateAssn(*this, e, *out_recowaveforms1, wf, *recorawassn1);

    // Apply Gaussian filter and remove electronics response
    for (int i = 0; i<ntbin; ++i){
      double filt = filtfun->Eval(150.*i/ntbin);
      if (i>ntbin/2){
        filt = filtfun->Eval(150.*(ntbin-i)/ntbin);
      }
      std::complex<double> tmp = wffft[i]*vinvres[i]*fPECalib[wf->ChannelNumber()];
      in[i][0] = real(tmp)*filt;
      in[i][1] = imag(tmp)*filt;
      //std::cout<<i<<" "<<wffft[i]<<" "<<vinvres[i*vinvres.size()/ntbin]<<" "<<tmp<<" "<<filt<<" "<<in[i][0]<<" "<<in[i][1]<<std::endl;
    }
    fftw_execute(pb);

    for (int i = 0; i<ntbin; ++i){
      wfimpulse[i] = out[i][0]/ntbin;
      //std::cout<<out[i][0]<<std::endl;
    }
    std::vector<float> wfimpulse_shift(ntbin);
    for (int i = 0; i<ntbin; ++i){
      if (i<95){
        wfimpulse_shift[i] = wfimpulse[ntbin-95+i];
      }
      else{
        wfimpulse_shift[i] = wfimpulse[i-95]; 
      }
    }
    recob::OpWaveform out_recowaveFinal2(wf->TimeStamp(), wf->ChannelNumber(), wfimpulse_shift);
    out_recowaveforms2->emplace_back(std::move(out_recowaveFinal2));
    util::CreateAssn(*this, e, *out_recowaveforms2, wf, *recorawassn2);

  }
  e.put(std::move(out_recowaveforms1),"wiener");
  e.put(std::move(out_recowaveforms2));
  e.put(std::move(recorawassn1),"wiener");
  e.put(std::move(recorawassn2));
}

DEFINE_ART_MODULE(pdsp::PDSPCalPhoton)
