#include "PDVD_PDMapAlg.hh"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"


//------------------------------------------------------------------------------
//--- opdet::PDVD_PDMapAlg implementation
//------------------------------------------------------------------------------

namespace opdet {

  PDVD_PDMapAlg::PDVD_PDMapAlg(const fhicl::ParameterSet&)
  {
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file("PDVD_PDS_Mapping.json", fname);
    std::ifstream i(fname, std::ifstream::in);
    i >> PDmap;
    i.close();
  }

PDVD_PDMapAlg::~PDVD_PDMapAlg()
  { }

  bool PDVD_PDMapAlg::isPDType(size_t ch, std::string pdname) const
  {
    if(PDmap.at(ch)["pd_type"] == std::string(pdname)) return true;
    return false;
  }

  bool PDVD_PDMapAlg::isSensitiveToAr(size_t ch) const
  {
    return PDmap.at(ch)["sens_Ar"].get<bool>();
  }

  bool PDVD_PDMapAlg::isSensitiveToXe(size_t ch) const
  {
    return PDmap.at(ch)["sens_Xe"].get<bool>();
  }

  std::string PDVD_PDMapAlg::pdType(size_t ch) const
  {
    return PDmap.at(ch)["pd_type"];
  }

  double PDVD_PDMapAlg::Efficiency(size_t ch) const
  {
  return PDmap.at(ch)["eff"];
  }

  std::string PDVD_PDMapAlg::hardwareChannel(size_t ch) const
  {
    return PDmap.at(ch)["HardwareChannel"];
  }

  std::vector<int> PDVD_PDMapAlg::getChannelsOfType(std::string pdname) const
  {
    std::vector<int> out_ch_v;
    for (size_t ch = 0; ch < PDmap.size(); ch++) {
      if (PDmap.at(ch)["pd_type"] == pdname) out_ch_v.push_back(ch);
    }
    return out_ch_v;
  }

  auto PDVD_PDMapAlg::getChannelEntry(size_t ch) const
  {
    return PDmap.at(ch);
  }

  size_t PDVD_PDMapAlg::size() const
  {
    return PDmap.size();
  }

  // template<>
  nlohmann::json PDVD_PDMapAlg::getCollectionWithProperty(std::string property)
  {
    nlohmann::json subSetPDmap;
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  [property](const nlohmann::json e)->bool
                    {return e[property];} );
    return subSetPDmap;
  }





DEFINE_ART_CLASS_TOOL(PDVD_PDMapAlg)
}

