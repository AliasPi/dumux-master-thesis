#ifndef INTERPOLATE_PCENTRY_HH
#define INTERPOLATE_PCENTRY_HH

#include <dune/common/fvector.hh>
#include <algorithm>

//! This class provides an interpolation between two capillary entry pressures,
//! one coming from facies1 (pc1) and the other from facies5 (pc5).
template<typename Scalar, int dim>
class InterpolatePcEntry
{
public:
    //! Constructor: supply the capillary entry pressures for facies1 and facies5.
    InterpolatePcEntry(Scalar pcEntryFacies1, Scalar pcEntryFacies5)
        : pc1_(pcEntryFacies1), pc5_(pcEntryFacies5)
    {}

    //! Given a global position, return an interpolated capillary entry pressure.
    /*!
      The interpolation region is defined (in simulation coordinates) as:
      - x ∈ [2450, 3700]
      - y ∈ [2350, 2500]   (note: y = original y + 2000)
      
      A simple average of the normalized distances in x and y is used as weight.
      For positions below the region the weight is 0 (use pc1) and above the region the weight is 1 (use pc5).
    */
    Scalar getPcEntry(const Dune::FieldVector<Scalar, dim>& pos) const
    {
        // Define region boundaries in simulation coordinates.
        const Scalar x_min = 2450.0;
        const Scalar x_max = 3700.0;
        const Scalar y_min = 2350.0; // 350 + 2000
        const Scalar y_max = 2500.0; // 500 + 2000

        // Compute a normalized distance in x.
        Scalar tx = 0.0;
        if (pos[0] <= x_min)
            tx = 0.0;
        else if (pos[0] >= x_max)
            tx = 1.0;
        else
            tx = (pos[0] - x_min) / (x_max - x_min);

        // Compute a normalized distance in y.
        Scalar ty = 0.0;
        if (pos[1] <= y_min)
            ty = 0.0;
        else if (pos[1] >= y_max)
            ty = 1.0;
        else
            ty = (pos[1] - y_min) / (y_max - y_min);

        // Use the average of the two normalized distances as the overall weight.
        Scalar weight = (tx + ty) / 2.0;

        // Return the interpolated capillary entry pressure:
        // weight==0 → use facies1 value, weight==1 → use facies5 value.
        return (1.0 - weight) * pc1_ + weight * pc5_;
    }

private:
    Scalar pc1_, pc5_;
};

//! Helper function that returns true if the given position lies in the interpolation region.
template<typename Scalar, int dim>
bool isWithinInterface(const Dune::FieldVector<Scalar, dim>& pos)
{
    const Scalar x_min = 2450.0;
    const Scalar x_max = 3700.0;
    const Scalar y_min = 2350.0; // 350 + 2000
    const Scalar y_max = 2500.0; // 500 + 2000

    return (pos[0] >= x_min && pos[0] <= x_max &&
            pos[1] >= y_min && pos[1] <= y_max);
}

#endif // INTERPOLATE_PCENTRY_HH
