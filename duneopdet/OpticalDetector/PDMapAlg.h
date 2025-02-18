#ifndef PDMapAlg_H
#define PDMapAlg_H

#include "fhiclcpp/ParameterSet.h"
//#include "art/Utilities/ToolMacros.h"

#include <string>

namespace opdet {

  //virtual base class
  class PDMapAlg {

  public:
    //Default destructor
    virtual ~PDMapAlg() noexcept = default;

    //only required element: define string per opdet channel number
    virtual std::string pdType(size_t ch) const = 0;

    //helper for PDType lookup
    virtual bool isPDType(size_t ch, std::string pdname) const
    { return (pdType(ch)==pdname); }

  private:

  }; // class PDMapAlg

  /*
  //a default dummy one
  class PDMapAlgSimple : PDMapAlg
  {
  public:
    explicit PDMapAlgSimple(const fhicl::ParameterSet& pset)
    { fType = pset.get<std::string>("Type","pmt"); }
    
    ~PDMapAlgSimple();

    std::string pdType(size_t ch) const override { return fType; }

  private:
    std::string fType;
    
  };


  DEFINE_ART_CLASS_TOOL(PDMapAlgSimple)
  */
} // namespace


#endif // PDMAPALG_H

