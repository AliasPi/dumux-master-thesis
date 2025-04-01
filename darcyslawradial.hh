#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm> // for std::max

#include <dumux/common/properties.hh>

namespace Dumux
{

template <class TypeTag>
class DarcysLawRadial{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using CO2 = typename FluidSystem::CO2;

public:
    //--------------------------------------------------------------------------
    // Numerical derivative of viscosity wrt p using central difference:
    //
    //   dmu/dp ~ [ mu(p + dp) - mu(p - dp) ] / (2 * dp).
    //
    // We pick dp adaptively to be the maximum of some fraction of p or 1000 Pa,
    // ensuring dp is not too small or too large.
    //--------------------------------------------------------------------------
    static Scalar dViscositydP_numerical(Scalar T, Scalar p)
    {
        // Safety: if p < 1 Pa, bump it up to avoid negative or near-zero range
        Scalar pSafe = std::max(p, 1.0);

        // Let's choose dp as a fraction of p, but at least 1000 Pa:
        Scalar dp = std::max(1.0e3, 1.0e-5 * pSafe);

        // Evaluate mu at p + dp, p - dp
        Scalar muPlus  = CO2::liquidViscosity(T, pSafe + dp);
        Scalar muMinus = CO2::liquidViscosity(T, (pSafe > dp ? pSafe - dp : 1.0));

        // Central difference
        Scalar derivative = (muPlus - muMinus) / (2.0 * dp);
        return derivative;
    }

    //------------------------------------------------------------------------------
    // Newton's method to solve for pWell in the radial-flow equation:
    //
    //   F(p) = p - pInf - [ q * mu(p ) / (2π k h ) ] * ln(rInf / rWell) = 0.
    //
    // With a numerical derivative for mu(p). We also handle possible issues like
    // near-zero derivatives, negative p, non-convergence, etc.
    //
    // Inputs:
    //   - pInf         : far-field reservoir pressure [Pa]
    //   - rInfinity    : outer radius where pInf applies [m]
    //   - rWell        : well radius [m]
    //   - k            : permeability [m²]
    //   - h            : reservoir thickness [m]
    //   - injectionRate: mass injection rate of CO₂ [kg/s]
    //   - T            : temperature [K]
    //   - tol          : convergence tolerance (for p) [Pa]
    //   - maxIter      : maximum allowed iterations
    //
    // Return: the final pWell (injection pressure).
    //------------------------------------------------------------------------------
    static Scalar computeInjectionPressureNewton(Scalar pInf,
                                        Scalar rInfinity,
                                        Scalar rWell,
                                        Scalar k,
                                        Scalar h,
                                        Scalar injectionRate,
                                        Scalar T,
                                        Scalar tol     = 1.0,
                                        int    maxIter = 50)
    {
        // Pi
        const Scalar PI = 3.14159265358979323846;

        // 1) Initial volumetric flow [m³/s], using density at pInf
        Scalar rho_init = CO2::liquidDensity(T, pInf > 0.0 ? pInf : 1.0);
        Scalar q_volume = injectionRate / rho_init;

        // 2) log term
        if (rInfinity <= rWell) {
            throw std::runtime_error("Error: rInfinity must be > rWell.");
        }
        Scalar logTerm = std::log(rInfinity / rWell);

        // 3) Initial guess for pWell
        //    Start at pInf plus 1 bar => 1e5 Pa, but also avoid pInf < 0
        Scalar pWell = (pInf > 0.0 ? pInf : 1.0e5) + 1.0e5;

        // 4) Newton iteration
        bool converged = false;
        Scalar pWell_new = pWell;

        for (int iter = 0; iter < maxIter; ++iter)
        {
            // i) Evaluate mu at current guess
            Scalar muNow = CO2::liquidViscosity(T, pWell);

            // ii) Evaluate derivative dmu/dp using numeric approximation
            Scalar dMuDpNow = dViscositydP_numerical(T, pWell);

            // iii) Define F(pWell) = pWell - pInf - (q_volume * muNow / (2π k h)) * logTerm
            Scalar factor   = (q_volume * muNow) / (2.0 * PI * k * h);
            Scalar F_value  = pWell - pInf - factor * logTerm;

            // iv) F'(pWell) = 1 - [ (q_volume / (2π k h)) * dMuDpNow * logTerm ]
            Scalar dfactor_dp = (q_volume * dMuDpNow) / (2.0 * PI * k * h);
            Scalar F_prime    = 1.0 - dfactor_dp * logTerm;

            // Safety check: Avoid division by zero or near-zero
            if (std::fabs(F_prime) < 1e-14) {
                std::cerr << "[WARNING] F'(pWell) is near zero. "
                        << "Aborting iteration to prevent numerical blow-up.\n";
                break;
            }

            // v) Newton update
            pWell_new = pWell - F_value / F_prime;

            // Safety: ensure pWell_new doesn't go negative
            if (pWell_new < 1.0) {
                pWell_new = 1.0;
            }

            // vi) Convergence check
            if (std::fabs(pWell_new - pWell) < tol) {
                pWell = pWell_new;
                converged = true;
                break;
            }

            // (Optional) If density changes significantly with p, update q_volume
            // Scalar rho_now = CO2::density(T, pWell_new);
            // q_volume       = injectionRate / rho_now;

            // vii) Update for next iteration
            pWell = pWell_new;
        }

        // 5) Check final status
        if (!converged) {
            std::cerr << "[WARNING] Newton did not converge within " << maxIter 
                    << " iterations. Last pWell: " << pWell << " Pa.\n";
        }

        return pWell;
    }

};
    
}



// //------------------------------------------------------------------------------
// // Demo main()
// //------------------------------------------------------------------------------
// int main()
// {
//     // Example reservoir + well data (replace with real site data for Norway)
//     Scalar pInf          = 10.0e6;     // Far-field pressure [Pa]
//     Scalar rInfinity     = 1000.0;     // Outer boundary radius [m]
//     Scalar rWell         = 0.15;       // Well radius [m] (Paper SPE11)
//     Scalar k             = 1.0e-13;    // Permeability [m²]
//     Scalar h             = 1;          // Reservoir thickness [m] 1m for 2D (Paper SPE11)
//     Scalar injectionRate = 10.0;       // mass injection rate [kg/s] (CO₂)
//     Scalar T             = 318.15;     // Temperature [K]

//     // Tolerance and iterations
//     Scalar tol     = 10.0;   // 10 Pa tolerance
//     int    maxIter = 50;

//     // Compute injection pressure using Newton’s method with numeric derivative
//     Scalar pWell = computeInjectionPressureNewton(
//         pInf, rInfinity, rWell,
//         k, h,
//         injectionRate,
//         T,
//         tol,
//         maxIter
//     );

//     // Print result
//     std::cout << "Final injection well pressure (Newton + numeric dmu/dp): " 
//               << pWell << " Pa\n";

//     return 0;
// }
