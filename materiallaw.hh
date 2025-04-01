#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {
namespace FluidMatrix {

template <class Scalar>
class RelativPermeability : public BrooksCoreyRegularization<Scalar> 
{
    using ParentType = BrooksCoreyRegularization<Scalar>;
    
public:
    using ParentType::ParentType;

    template<class MaterialLaw>
    void init(const MaterialLaw* m, const std::string& paramGroup)
    {
        ParentType::init(m, paramGroup);
        std::cout << "Material Law Initialized" << std::endl;

        swr_ = getParam<Scalar>(paramGroup + ".Swr");
        
        mode_ = getParam<std::string>("MaterialLaw.Mode");

        if (paramGroup == "Facies1"){
            std::cout << "------ Using Material Law: " << mode_ << " ------" << std::endl;
        }


        if (mode_ != "master1" && mode_ != "master2" && mode_ != "brooks" && mode_ != "genuchten")
        {
            throw std::runtime_error("MaterialLaw.Mode must be either 'master1', 'master2', 'brooks', or 'genuchten'");
        }

        if (mode_ == "genuchten")
        {
            vanGenuchtenN_ = 2.5;
        }
    }

    OptionalScalar<Scalar> krw(const Scalar swe) const
    {

        if (mode_ == "brooks")
            return ParentType::krw(swe);

        if (mode_ == "genuchten")
            return krwGenuchten(swe);
        
        if (mode_ == "master1")
            return krwMaster1(swe);

        if (mode_ == "master2")
            return krwMaster2(swe);        

        return ParentType::krw(swe);
    }

    OptionalScalar<Scalar> krwGenuchten(const Scalar swe) const
    {
        using std::clamp;
        using std::pow;
        using std::sqrt;

        // Use the stored residual saturation
        Scalar swr = swr_;
        // Compute effective saturation: Se = (S - Swr) / (1 - Swr)
        Scalar se = (swe - swr) / (1.0 - swr);        

        Scalar seClamped = clamp(se, 0.0, 1.0);
        
        Scalar m = 1.0 - 1.0/vanGenuchtenN_;
        // van Genuchten water relative permeability:
        // krw = sqrt(Se) * [1 - (1 - Se^(1/m))^m]^2
        Scalar term = 1.0 - pow(seClamped, 1.0/m);
        Scalar krw_val = sqrt(seClamped) * pow(1 - pow(term, m), 2);
        return krw_val;        
    }

    OptionalScalar<Scalar> krwMaster1(const Scalar swe) const
    {
        using std::clamp;
        using std::pow;

        Scalar sweClamped = clamp(swe, 0.0, 1.0);

        Scalar krw_max = 0.98;
        Scalar swr = swr_;
        Scalar cw = 1.8;
        
        Scalar sEff = (sweClamped - swr) / (1.0 - swr);

        Scalar kr = krw_max * pow(sEff, cw);
        return kr;
    }

    OptionalScalar<Scalar> krwMaster2(const Scalar swe) const
    {
        using std::clamp;
        using std::pow;

        Scalar sweClamped = clamp(swe, 0.0, 1.0);

        Scalar krw_max = 0.88;
        Scalar swr = swr_;
        Scalar cw = 3.9;
        
        Scalar sEff = (sweClamped - swr) / (1.0 - swr);

        Scalar kr = krw_max * pow(sEff, cw);
        return kr;      
    }

    OptionalScalar<Scalar> krn(const Scalar swe) const
    {

        if (mode_ == "brooks")
            return ParentType::krn(swe);

        if (mode_ == "genuchten")
            return krnGenuchten(swe);
        
        if (mode_ == "master1")
            return krnMaster1(swe);

        if (mode_ == "master2")
            return krnMaster2(swe);        

        return ParentType::krn(swe);
    }

    OptionalScalar<Scalar> krnGenuchten(const Scalar swe) const
    {
        using std::clamp;
        using std::pow;

        Scalar swr = swr_;
        Scalar se = (swe - swr) / (1.0 - swr);

        Scalar seClamped = clamp(se, 0.0, 1.0);
        Scalar m = 1.0 - 1.0/vanGenuchtenN_;
        // van Genuchten non-wetting relative permeability:
        // Example formulation: krn = (1 - Se)^2 * [1 - Se^(1/m)]^(2*m)
        Scalar krn_val = pow(1 - seClamped, 2) * pow(1 - pow(seClamped, 1.0/m), 2*m);
        return krn_val;
    }

    OptionalScalar<Scalar> krnMaster1(const Scalar swe) const
    {
        using std::clamp;
        using std::pow;

        Scalar sweClamped = clamp(swe, 0.0, 1.0);

        Scalar krCO2_max = 0.56;
        Scalar swr = swr_;
        Scalar cn = 2.2;
        
        Scalar sEff = (sweClamped - swr) / (1.0 - swr);

        Scalar sCO2Eff = 1.0 - sEff;

        Scalar kr = krCO2_max * pow(sCO2Eff, cn);
        return kr;        
    }

    OptionalScalar<Scalar> krnMaster2(const Scalar swe) const
    {
        using std::clamp;
        using std::pow;

        Scalar sweClamped = clamp(swe, 0.0, 1.0);

        Scalar krCO2_max = 0.56;
        Scalar swr = swr_;
        Scalar cn = 2.0;
        
        Scalar sEff = (sweClamped - swr) / (1.0 - swr);

        Scalar sCO2Eff = 1.0 - sEff;

        Scalar kr = krCO2_max * pow(sCO2Eff, cn);
        return kr;       
    }

private:
    std::string mode_;

    Scalar vanGenuchtenN_;

    Scalar swr_ = 0;

};

template<typename Scalar = double>
using RelativPermeabilityDefault = TwoPMaterialLaw<Scalar, BrooksCorey, RelativPermeability<Scalar>, TwoPEffToAbsDefaultPolicy>;

}
}