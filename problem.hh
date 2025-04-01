// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CO2Tests
 * \brief Definition of a problem, where CO2 is injected into a reservoir.
 */

#ifndef DUMUX_HETEROGENEOUS_PROBLEM_HH
#define DUMUX_HETEROGENEOUS_PROBLEM_HH

#include <cmath>

#include <dune/grid/common/partitionset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/parallel/vectorcommdatahandle.hh>

#include <dumux/discretization/box/scvftoscvboundarytypes.hh>
#include <dumux/common/gridcapabilities.hh>

#include <dumux/porousmediumflow/problem.hh>

#include "darcyslawradial.hh"

namespace Dumux {

template <class TypeTag >
class HeterogeneousProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using PointSource =  GetPropType<TypeTag, Properties::PointSource>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using CO2 = typename FluidSystem::CO2;

    enum
    {
        // Definiert in Volume Variables
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx,
        

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        // phase presence index
        firstPhaseOnly = Indices::firstPhaseOnly,

        // component indices
        BrineIdx = FluidSystem::BrineIdx,
        CO2Idx = FluidSystem::CO2Idx,
        

        // equation indices
        //Massenbilanz Brine
        conti0EqIdx = Indices::conti0EqIdx,
        // Massenbilanz CO2
        contiCO2EqIdx = conti0EqIdx + CO2Idx,
    };

#if !ISOTHERMAL
    enum {
        temperatureIdx = Indices::temperatureIdx,
        // Energiebilanz
        energyEqIdx = Indices::energyEqIdx,
    };
#endif

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

    // the discretization method we are using
    using DiscretizationMethod = DiscretizationMethods::Box;
    static constexpr DiscretizationMethod discMethod{};
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;

    // 2D oder 3D
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    template<class SpatialParams>
    HeterogeneousProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    , bottomBoundaryId_(10)
    , topBoundaryId_(11)
    , isWellInjecting1_(false)
    , isWellInjecting2_(false)
    {
        nTemperature_        = getParam<int>("FluidSystem.NTemperature");
        nPressure_           = getParam<int>("FluidSystem.NPressure");
        pressureLow_         = getParam<Scalar>("FluidSystem.PressureLow");
        pressureHigh_        = getParam<Scalar>("FluidSystem.PressureHigh");
        temperatureLow_      = getParam<Scalar>("FluidSystem.TemperatureLow");
        temperatureHigh_     = getParam<Scalar>("FluidSystem.TemperatureHigh");
        name_                = getParam<std::string>("Problem.Name");

        useWell1_            = getParam<bool>("Well1.UseWell");
        injectionRateWell1_  = getParam<Scalar>("Well1.InjectionRate");
        positionWell1_       = getParam<std::vector<Scalar>>("Well1.Position");

        useWell2_            = getParam<bool>("Well2.UseWell");
        injectionRateWell2_  = getParam<Scalar>("Well2.InjectionRate");
        positionWell2_       = getParam<std::vector<Scalar>>("Well2.Position");

#if !ISOTHERMAL
        injectionPressureWell1_    = getParam<Scalar>("Well1.InjectionPressure");
        injectionTemperatureWell1_ = getParam<Scalar>("Well1.InjectionTemperature");

        injectionPressureWell2_    = getParam<Scalar>("Well2.InjectionPressure");
        injectionTemperatureWell2_ = getParam<Scalar>("Well2.InjectionTemperature");
#endif

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);  

        this->spatialParams().readPhysicalGroups();
        this->spatialParams().determineBoundaryCells();

        // precompute the boundary types for the box method from the cell-centered boundary types
        scvfToScvBoundaryTypes_.computeBoundaryTypes(*this);
    }

    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolume &scv) const
    { 
        return scvfToScvBoundaryTypes_.boundaryTypes(scv); 
    }

    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes bcTypes;

        const auto boundaryId = scvf.boundaryFlag();

        if (boundaryId == bottomBoundaryId_ ) {
            // Massenbilanzen Neumann No-Flow Randbedingungen
            bcTypes.setAllNeumann();
#if !ISOTHERMAL
            // Energiebilanz Dirichlet feste Temperaturen 
            bcTypes.setDirichlet(energyEqIdx);
#endif
        }
        else if (boundaryId == topBoundaryId_) {
            // Massenbilanzen Neumann No-Flow Randbedingungen
            bcTypes.setAllNeumann();
#if !ISOTHERMAL
            // Energiebilanz Dirichlet feste Temperaturen 
            bcTypes.setDirichlet(energyEqIdx);
#endif
        }
        else {
            // Massenbilanzen und Energiebilanz Neumann No-Flow Randbedingungen
            bcTypes.setAllNeumann();            
        }
        

        return bcTypes;
    }


    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector fluxes(0.0);
        return fluxes;
    }

    // Hinzufügen der Wells
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        // TODO: In 3D there should be multiple point sources modeling the horizontal wells
        std::cout << "Adding Point Sources" << std::endl;

        if (useWell1_){
#if !ISOTHERMAL
            auto well1 = PointSource({
                positionWell1_[0], // x
                positionWell1_[1], // y
                positionWell1_[2]  // z
                }, 
                {
                    0, // Brine Massenbilanz
                    injectionRateWell1_, // CO2 Massenbilanz
                    0 // Energiebilanz
                }, 1 // Id
            );
#else
            auto well1 = PointSource({positionWell1_[0], positionWell1_[1], positionWell1_[2]}, {0, injectionRateWell1_}, 1);
#endif
            pointSources.push_back(well1);
        }

        if (useWell2_){
#if !ISOTHERMAL
            auto well2 = PointSource({positionWell2_[0], positionWell2_[1], positionWell2_[2]}, {0, injectionRateWell2_, 0}, 2);
#else
            auto well2 = PointSource({positionWell2_[0], positionWell2_[1], positionWell2_[2]}, {0, injectionRateWell2_}, 2);
#endif
            
            pointSources.push_back(well2);
        }
    }

    // Berechnung während der Simulation
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        const auto& volVars = elemVolVars[scv];  

        if (source.id() == 1) {
            if (isWellInjecting1_) {   
                // TODO: Think about how to handle the pressure
                // The pressure should be atleast be the pressure of the domain
                // https://stacks.stanford.edu/file/druid:xk906gg6034/Shu05.pdf
#if !ISOTHERMAL
                // Druck der Lagerstätte
                const auto pressure = volVars.pressure(0);
                // Energieeintrag
                const auto energySource = injectionRateWell1_ * CO2::liquidEnthalpy(injectionTemperatureWell1_, pressure);
                source = NumEqVector({0, injectionRateWell1_, energySource });  
#endif             
            }
            else {
                // Keine Injektion
                source = typename PointSource::Values(0.0);
            }
        }
        else if (source.id() == 2) {
            if (isWellInjecting2_) {
#if !ISOTHERMAL
                const auto pressure = volVars.pressure(0);
                const auto energySource = injectionRateWell2_ * CO2::liquidEnthalpy(injectionTemperatureWell2_, pressure);
                source = NumEqVector({0, injectionRateWell2_, energySource });
#endif
                
            }
            else {
                source = typename PointSource::Values(0.0);
            }
        }    

        
    }

    template<class VTKWriter>
    void addFieldsToWriter(VTKWriter& vtk)
    {
        const auto numElements = this->gridGeometry().gridView().size(0);
        const auto numDofs = this->gridGeometry().numDofs();

        vtkPermeabilityH_.resize(numElements);   
        vtkPermeabilityV_.resize(numElements);

        vtkPorosity_.resize(numElements);    
        vtkPermeability_.resize(numElements);           

        vtk.addField(vtkPorosity_, "cellwisePorosity");
        
        vtk.addField(vtkPermeabilityH_, "horizontalPermeability");
        vtk.addField(vtkPermeabilityV_, "verticalPermeability");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.permeability()[0][0]; }, "effPermeabilityH");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.permeability()[1][1]; }, "effPermeabilityV");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(phase0Idx); }, "density_liq");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(phase1Idx); }, "density_gas");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.averageMolarMass(phase0Idx); }, "avgMolarMass_liq");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.averageMolarMass(phase1Idx); }, "avgMolarMass_gas");        

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(phase0Idx); }, "saturation_liq");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(phase1Idx); }, "saturation_gas");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.moleFraction(phase0Idx, BrineIdx); }, "moleFractionW_liq");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.moleFraction(phase0Idx, CO2Idx); }, "moleFractionN_liq");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.moleFraction(phase1Idx, BrineIdx); }, "moleFractionW_gas");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.moleFraction(phase1Idx, CO2Idx); }, "moleFractionN_gas");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.viscosity(phase0Idx); }, "viscosity_liq");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.viscosity(phase1Idx); }, "viscosity_gas");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.molarDensity(phase0Idx); }, "molarDensit_liq");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.molarDensity(phase1Idx); }, "molarDensity_gas");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.relativePermeability(phase0Idx); }, "relativePermeability_liq");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.relativePermeability(phase1Idx); }, "relativePermeability_gas");




#if !ISOTHERMAL
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.enthalpy(BrineIdx); }, "enthalpyW");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.enthalpy(CO2Idx); }, "enthalpyN");
#else
        vtkTemperature_.resize(numDofs, 0.0);
        vtk.addField(vtkTemperature_, "T");
#endif

        const auto& gridView = this->gridGeometry().gridView();
        auto fvGeometry = localView(this->gridGeometry());
        Scalar maxTemp = -std::numeric_limits<Scalar>::max();
        Scalar minTemp = std::numeric_limits<Scalar>::max();
        Scalar maxPress = -std::numeric_limits<Scalar>::max();
        Scalar minPress = std::numeric_limits<Scalar>::max();
        GlobalPosition maxTempPos, minTempPos, maxPressPos, minPressPos;
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();
                const auto& pos = scv.dofPosition();
                const auto pressure = approximateInitialReservoirPressure_(pos, 20);
                const auto temp = this->spatialParams().temperatureAtPos(pos);

                if(temp > maxTemp) {
                    maxTemp = temp;
                    maxTempPos = pos;
                }
                if(temp < minTemp) {
                    minTemp = temp;
                    minTempPos = pos;
                }
                if(pressure > maxPress) {
                    maxPress = pressure;
                    maxPressPos = pos;
                }
                if(pressure < minPress) {
                    minPress = pressure;
                    minPressPos = pos;
                }

#if ISOTHERMAL
                vtkTemperature_[dofIdxGlobal] = this->spatialParams().temperatureAtPos(scv.dofPosition());
#endif
            }
           
            vtkPorosity_[eIdx] = 1- this->spatialParams().inertVolumeFraction(eIdx);
            auto perm = this->spatialParams().permeability(eIdx);
            vtkPermeabilityH_[eIdx] = perm[0][0];
            vtkPermeabilityV_[eIdx] = perm[dimWorld-1][dimWorld-1];
        }    
    }

    void setWellInjecting1(bool isWellInjecting)
    { isWellInjecting1_ = isWellInjecting; }

    void setWellInjecting2(bool isWellInjecting)
    { isWellInjecting2_ = isWellInjecting; }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }
      
    const std::string& name() const
    { return name_; }    
    
private:
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(firstPhaseOnly);

        Scalar pressure = approximateInitialReservoirPressure_(globalPos, 20);
        
        const Scalar moleFracLiquidCO2 = 0.0;
        const Scalar moleFracLiquidBrine = 1.0 - moleFracLiquidCO2;

        const Scalar meanM = FluidSystem::molarMass(BrineIdx)*moleFracLiquidBrine
                             + FluidSystem::molarMass(CO2Idx)*moleFracLiquidCO2;

        if(useMoles) // mole-fraction formulation
            values[switchIdx] = moleFracLiquidCO2;
        else // mass-fraction formulation
            values[switchIdx] = moleFracLiquidCO2*FluidSystem::molarMass(CO2Idx)/meanM;

        values[pressureIdx] = pressure;

#if !ISOTHERMAL
        values[temperatureIdx] = this->spatialParams().temperatureAtPos(globalPos);
#endif

        return values;
    }

    float approximateInitialReservoirPressure_(const GlobalPosition &globalPos, int iterations) const
    {
        const Scalar temp = this->spatialParams().temperatureAtPos(globalPos);
        Scalar densityW = 1000; // initial guess
        const Scalar overburdenDepth = 2000.0; // [m] extra depth
        Scalar pressure = 1.0e5; // atmospheric pressure at the “surface”
        
        // Iterative refinement of the pressure.
        // Add the pressure due to overburden to the computed hydrostatic column.
        for (int i = 0; i < iterations; i++){
            densityW = FluidSystem::Brine::liquidDensity(temp, pressure);
            pressure = 1.0e5 + (overburdenDepth + (this->gridGeometry().bBoxMax()[dimWorld-1] - globalPos[dimWorld-1])) * densityW * 9.81;
        }
        return pressure;
    }


    int bottomBoundaryId_, topBoundaryId_;

    std::string name_ ;

    int nTemperature_, nPressure_;
    Scalar pressureLow_, pressureHigh_, temperatureLow_, temperatureHigh_;

    std::vector<Scalar> vtkPorosity_, vtkTemperature_, vtkPermeability_, vtkPermeabilityH_, vtkPermeabilityV_;

    bool isWellInjecting1_, isWellInjecting2_;
    bool useWell1_, useWell2_;
    Scalar injectionRateWell1_, injectionRateWell2_;        
    std::vector<Scalar> positionWell1_, positionWell2_;
    ScvfToScvBoundaryTypes<BoundaryTypes, DiscretizationMethod> scvfToScvBoundaryTypes_;

#if !ISOTHERMAL
    Scalar injectionPressureWell1_, injectionTemperatureWell1_, injectionPressureWell2_, injectionTemperatureWell2_;
#endif

};

}

#endif