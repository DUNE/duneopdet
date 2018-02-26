//=========================================================
// AlgoSSPLeadingEdge.cc
// This class provides the digital output for the 
// Leading Edge Discrimintator, used in SSP module.
//
// Vitor Luzio, UFABC, 2017
// Based on AlgoSiPM.cxx
//=========================================================



#include "dune/OpticalDetector/AlgoSSPLeadingEdge.h"

namespace pmtana {

  //---------------------------------------------------------------------------
  AlgoSSPLeadingEdge::AlgoSSPLeadingEdge(const fhicl::ParameterSet &pset,const std::string name)
    : PMTPulseRecoBase(name)
  {

    _adc_thres = pset.get< float  >("ADCThreshold"   );
//    _min_width = pset.get< float  >("MinWidth"       );
//    _2nd_thres = pset.get< float  >("SecondThreshold");
    _pedestal  = pset.get< float  >("Pedestal"       );
    _dwindow   = pset.get< int    >("DWindow"        );
    _rdwindow  = pset.get< size_t >("ReadoutWd"      );
    _pretrg    = pset.get< size_t >("PreTrg"         );

 //  std::cout << "VITOR DEBUG" << std::endl;
 //  std::cout << "ADCThreshold = " << _adc_thres << std::endl;
 //  std::cout << "DWindow = " << _dwindow << std::endl;
 //  std::cout << "Readout Window = " << _rdwindow << std::endl;
 //  std::cout << "Pre Trigger = " << _pretrg << std::endl;
 //  std::cout << "VITOR DEBUG" << std::endl;




//    _nsigma = 5;

    Reset();

  }

  //---------------------------------------------------------------------------
  AlgoSSPLeadingEdge::~AlgoSSPLeadingEdge()
  {}

  //---------------------------------------------------------------------------
  void AlgoSSPLeadingEdge::Reset()
  {

    PMTPulseRecoBase::Reset();

  }

  //---------------------------------------------------------------------------
  bool AlgoSSPLeadingEdge::RecoPulse( const pmtana::Waveform_t& wf,
                            const pmtana::PedestalMean_t& ped_mean,
                            const pmtana::PedestalSigma_t& ped_rms )
  {

    bool   fire          = false;
    int    counter       = 0;
    double pedestal      = _pedestal;
    double threshold     = _adc_thres;
    double threshold2    = _adc_thres;
    threshold           += pedestal;
 //   double pre_threshold = _2nd_thres;
 //   pre_threshold       += pedestal;

    int	   d_window	  = _dwindow; 
    int	   timer	  = 0;  //time to delay the wf and compare amplitudes
    int	   readout_window = _rdwindow; //from 0 to 2046 in SSP
    int    readout_pretrigger = _pretrg; //(from 0 to 2047 in SSP)
    double threshold_cmp;
    int    tamanho, cnt2;

    Reset();
    tamanho = wf.size();


    for (short const &value : wf) {

	 cnt2 = value;
         if(cnt2 != value){std::cout << "teste VITOR" << std::endl; }

         if (counter < (tamanho - d_window )) {

                threshold_cmp = wf[counter+d_window]-wf[counter];       
			

		if (threshold_cmp >= threshold2){
			
			if(timer==0){trg_wvf.emplace_back(counter+d_window);}

                        if(!fire && timer==0){
                                fire = true;
				if(readout_pretrigger < (counter + d_window)){
                                	_pulse.t_start = counter + d_window - readout_pretrigger;
				}else{
					_pulse.t_start = 1;
				}
                                timer++;
                        }else if(timer>0 && timer < d_window){
                                timer++;
                        }else if(timer == d_window){
                                timer=0;
                        }

                }else {
                        if(timer>0 && timer < d_window){timer++;}
                        else if(timer==d_window || timer==0 ){timer=0;}

                }
                if (fire==true && counter == (_pulse.t_start + readout_window)){ //we need to configure the _min_widht to be 700 bins
                         fire = false;
                         _pulse.t_end = counter - 1;
                         _pulse_v.push_back(_pulse);
                         _pulse.reset_param();

                }
        }

        counter++;

    }

    if (fire) { //It's not mine, but is a good strategy to handle the pulses at the end of waveform

      // Take care of a pulse that did not finish within the readout window
   
      fire = false;
      _pulse.t_end = tamanho;
      _pulse_v.push_back(_pulse);
      _pulse.reset_param();

    }

    return true;

  }

}
