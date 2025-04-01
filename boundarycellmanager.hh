#ifndef BOUNDARYCELLMANAGER_HH
#define BOUNDARYCELLMANAGER_HH

#include <vector>

template <class Scalar>
class BoundaryCellManager
{
public:
    // Singleton access method (Meyer's Singleton)
    static BoundaryCellManager<Scalar>& getInstance();

    // Set the vector of boundary cells
    void setBoundaryCells(const std::vector<bool>& boundaryCells);
    void setBoundaryCellVolumes(const std::vector<Scalar>& boundaryCellVolumes);


    // Check if a given cell index is on the boundary
    bool isBoundaryCell(int index) const;
    Scalar boundaryCellVolume(int index) const;

    // Delete copy constructor and assignment
    BoundaryCellManager(const BoundaryCellManager<Scalar>&) = delete;
    BoundaryCellManager<Scalar>& operator=(const BoundaryCellManager<Scalar>&) = delete;

private:
    // Private constructor and destructor to enforce singleton
    BoundaryCellManager() = default;
    ~BoundaryCellManager() = default;

    // The vector storing which cells are boundary
    std::vector<bool> boundaryCells_;
    std::vector<Scalar> boundaryCellVolumes_;
};

#endif // BOUNDARYCELLMANAGER_HH
