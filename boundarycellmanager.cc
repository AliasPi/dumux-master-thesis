#include "boundarycellmanager.hh"
#include <stdexcept>

// Meyer's Singleton: function-local static ensures exactly one instance
template <class Scalar>
BoundaryCellManager<Scalar>& BoundaryCellManager<Scalar>::getInstance()
{
    static BoundaryCellManager<Scalar> instance;
    return instance;
}

template <class Scalar>
void BoundaryCellManager<Scalar>::setBoundaryCells(const std::vector<bool>& boundaryCells)
{
    boundaryCells_ = boundaryCells;
}

template <class Scalar>
void BoundaryCellManager<Scalar>::setBoundaryCellVolumes(const std::vector<Scalar>& boundaryCellVolumes)
{
    boundaryCellVolumes_ = boundaryCellVolumes;
}

template <class Scalar>
bool BoundaryCellManager<Scalar>::isBoundaryCell(int index) const
{
    if (index < 0 || index >= static_cast<int>(boundaryCells_.size()))
    {
        // Optionally handle out-of-range index. For example, throw or return false:
        // throw std::out_of_range("Index out of range in isBoundaryCell()");
        return false; 
    }
    return boundaryCells_[index];
}

template <class Scalar>
Scalar BoundaryCellManager<Scalar>::boundaryCellVolume(int index) const
{
    if (index < 0 || index >= static_cast<int>(boundaryCellVolumes_.size()))
    {
        // Optionally handle out-of-range index. For example, throw or return false:
        // throw std::out_of_range("Index out of range in isBoundaryCell()");
        return 0.0; 
    }
    return boundaryCellVolumes_[index];
}

template class BoundaryCellManager<double>;
template class BoundaryCellManager<float>;
