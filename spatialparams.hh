// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CO2Tests
 * \brief Definition of the spatial parameters for the heterogeneous
 *        problem which uses the non-isothermal or isothermal CO2
 *        fully implicit model.
 */

#ifndef DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH
#define DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH

#include <iostream>
#include <dumux/io/grid/griddata.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dune/common/fmatrix.hh>
#include "materiallaw.hh"
#include "boundarycellmanager.hh"

namespace Dumux {

template<class GridGeometry, class Scalar>
class HeterogeneousSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry,
                                       Scalar,
                                       HeterogeneousSpatialParams<GridGeometry, Scalar>>
{
    using Grid = typename GridGeometry::Grid;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, HeterogeneousSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using PcKrSwCurve = FluidMatrix::RelativPermeabilityDefault<Scalar>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using PermeabilityType = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    
    HeterogeneousSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<const GridData<Grid>> gridData)
    : ParentType(gridGeometry)
    , gridData_(gridData)
    {
        std::cout << "Spatial Params Init" << std::endl;
        boundaryVolumeMultiplier_ = getParam<Scalar>("Problem.BoundaryVolumeMultiplier");
        useBoundaryVolume_ = getParam<bool>("Problem.UseBoundaryVolume");
        dispersivity_ = getParam<Scalar>("Problem.Dispersivity");
        readFaciesProperties();         
    }

    // returns the permeabilty to simulation
    // gets called every time step of the simulation
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {        
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return permeability(eIdx);
    }


    PermeabilityType permeability(std::size_t eIdx) const
    {        
        PermeabilityType K;
        // Get horizontal permeability
        const auto kh = permeabilities_[physicalGroupIdx(eIdx)];
        
        // Initialize matrix to 0
        K = 0.0;
        
        // Set horizontal permeability in x and y directions
        K[0][0] = kh;
        if (dimWorld > 1) {
            K[1][1] = kh;
        }
        
        // Set vertical permeability (0.1 * horizontal)
        if (dimWorld > 2) {
            K[2][2] = 0.1 * kh;  // For 3D
        } else if (dimWorld > 1) {
            K[1][1] = 0.1 * kh;  // For 2D
        }
        
        return K;
        }

    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        std::cout << "###############################################DISPER##################################" << dispersivity_ << std::endl;
        return dispersivity_;

    }

    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                          const SubControlVolume& scv,
                                          const ElementSolution& elemSol) const
    {      
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return makeFluidMatrixInteraction(pcKrSwCurves_[physicalGroupIdx(eIdx)]);
    }

    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return inertVolumeFraction(eIdx);         
    }

    Scalar inertVolumeFraction(std::size_t eIdx) const {
        return 1- porosities_[physicalGroupIdx(eIdx)]; 
    }

    template<class ElementSolution, class SolidState>
    Scalar solidThermalConductivity(const Element& element,
                                    const SubControlVolume& scv,
                                    const ElementSolution& elemSol,
                                    const SolidState& solidState) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return solidThermalConductivities_[physicalGroupIdx(eIdx)];
    }

    void readFaciesProperties()
    {
        int num_facies = 7;
        pcKrSwCurves_.reserve(num_facies);
        permeabilities_.resize(num_facies);
        porosities_.resize(num_facies);
        solidThermalConductivities_.resize(num_facies);

        for(int i = 0; i < num_facies; i++)
        {
            int facies = i + 1;
            pcKrSwCurves_.emplace_back("Facies" + std::to_string(facies));
            permeabilities_[i] = getParam<Scalar>("Facies" + std::to_string(facies) + ".Permeability");
            porosities_[i] = getParam<Scalar>("Facies" + std::to_string(facies) + ".Porosity");   
            solidThermalConductivities_[i] = getParam<Scalar>("Facies" + std::to_string(facies) + ".SolidThermalConductivity");
        }     
    }

    void readPhysicalGroups()
    {
        const auto& gg = this->gridGeometry();
        const auto& gridView = gg.gridView();
        physicalGroupIdx_.resize(gridView.size(0));

        for(const auto& element : elements(gridView))
        {          
            const auto eIdx = this->gridGeometry().elementMapper().index(element);    
            physicalGroupIdx_[eIdx] = gridData_->parameters(element)[0];
        }
    }

    void determineBoundaryCells() {
        const auto& gridView = this->gridGeometry().gridView();
        const auto size = gridView.size(0);

        std::vector<bool> boundaryCells(size, false);
        std::vector<Scalar> boundaryCellVolumes(size, 0.0);

        if (!useBoundaryVolume_){
            BoundaryCellManager<Scalar>::getInstance().setBoundaryCells(boundaryCells);
            BoundaryCellManager<Scalar>::getInstance().setBoundaryCellVolumes(boundaryCellVolumes);
            return;
        }

        for (const auto& element : elements(gridView)) {
            auto eIdx = gridView.indexSet().index(element);

            for (const auto& is : intersections(gridView, element)) {
                if (is.boundary()) {                    
                    boundaryCells[eIdx] = true;

                    if (physicalGroupIdx_[eIdx] == 1 || physicalGroupIdx_[eIdx] == 7) // Volume Multiplier does not apply to facies 1 and 7
                        continue;

                    if (!isLateralBoundary(is.geometry().center()[dimWorld-1])) {
                        continue;
                    }

                    boundaryCellVolumes[eIdx] += is.geometry().volume() * boundaryVolumeMultiplier_;
                }
            }            
        }

        BoundaryCellManager<Scalar>::getInstance().setBoundaryCells(boundaryCells);
        BoundaryCellManager<Scalar>::getInstance().setBoundaryCellVolumes(boundaryCellVolumes);
    }

    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::BrineIdx; }

    Scalar temperatureAtPos(const GlobalPosition &globalPos) const {
    return (70.0 + 273.15) - 0.025 * globalPos[dimWorld-1];
}

    int physicalGroupIdx(std::size_t eIdx) const
    { return physicalGroupIdx_[eIdx] - 1; }

    bool isLateralBoundary(Scalar height) const
    { return height > 1e-6 && height + 1e-6 < this->gridGeometry().bBoxMax()[dimWorld-1]; }

private:   
    Scalar boundaryVolumeMultiplier_;
    bool useBoundaryVolume_;
    Scalar dispersivity_;
    std::shared_ptr<const GridData<Grid>> gridData_;

    std::vector<int> physicalGroupIdx_;
    std::vector<Scalar> permeabilities_;
    std::vector<Scalar> porosities_;
    std::vector<Scalar> solidThermalConductivities_;

    std::vector<PcKrSwCurve> pcKrSwCurves_;
};

} // end namespace Dumux

#endif
