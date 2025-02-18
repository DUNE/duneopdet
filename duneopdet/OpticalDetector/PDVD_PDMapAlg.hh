#ifndef PDVD_PDMAPALG_HH
#define PDVD_PDMAPALG_HH

#include "PDMapAlg.h"
//#include "art/Utilities/ToolMacros.h"
//#include "art/Utilities/make_tool.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <string>

//#include "art_root_io/TFileService.h"

#include "json.hpp"

namespace opdet {

  class PDVD_PDMapAlg : PDMapAlg{

  public:
    //Default constructor
    explicit PDVD_PDMapAlg(const fhicl::ParameterSet& pset);
    PDVD_PDMapAlg() : PDVD_PDMapAlg(fhicl::ParameterSet()) {}
    //Default destructor
    ~PDVD_PDMapAlg();

    nlohmann::json getCollectionWithProperty(std::string property);
    template<typename T> nlohmann::json getCollectionWithProperty(std::string property, T property_value);
    template<typename Func> nlohmann::json getCollectionFromCondition(Func condition);

    // struct Config {};

    //  PDVD_PDMapAlg(Config const&) {}

    // void setup() {}

    bool isPDType(size_t ch, std::string pdname) const override;
    bool isSensitiveToAr(size_t ch) const;
    bool isSensitiveToXe(size_t ch) const;

    std::string pdType(size_t ch) const override;
    double Efficiency(size_t ch) const;
    std::string hardwareChannel(size_t ch) const;
      
    std::vector<int> getChannelsOfType(std::string pdname) const;
    auto getChannelEntry(size_t ch) const;
      
    size_t size() const;

  private:
    nlohmann::json PDmap;

  }; // class PDVD_PDMapAlg

  template<typename T>
  nlohmann::json PDVD_PDMapAlg::getCollectionWithProperty(std::string property, T property_value)
  {
    nlohmann::json subSetPDmap;
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  [property, property_value](auto const& e)->bool
                    {return e[property] == property_value;} );
    return subSetPDmap;
  }

  template<typename Func>
  nlohmann::json PDVD_PDMapAlg::getCollectionFromCondition(Func condition)
  {
    nlohmann::json subSetPDmap;
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  condition);
    return subSetPDmap;
  }

} // namespace

#endif // PDVD_PDMAPALG_HH

