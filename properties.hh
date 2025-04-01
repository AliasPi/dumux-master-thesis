#ifndef DUMUX_HETEROGENEOUS_PROPERTIES_HH
#define DUMUX_HETEROGENEOUS_PROPERTIES_HH

#include <dune/alugrid/grid.hh>
#include "discretization/box/box.hh"
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/defaultco2table.hh>
#include <dumux/material/components/h2o.hh>

// Fluid System
#include <dumux/material/fluidsystems/brineco2.hh>
// Model
#include "model/model.hh"

#ifndef ISOTHERMAL
#define ISOTHERMAL 0
#endif

#ifndef THREEDIMENSIONAL
#define THREEDIMENSIONAL 0
#endif

#include <dumux/common/pointsource.hh>

#include "spatialparams.hh"
#include "problem.hh"

// Aufbau der Simulation (was soll benutzt werden)

namespace Dumux::Properties {

// Set the model and the discretization method
namespace TTag {
struct Heterogeneous { using InheritsFrom = std::tuple<TwoPTwoCCO2, BoxModel>; };
}

#if !THREEDIMENSIONAL
template<class TypeTag>
struct Grid<TypeTag, TTag::Heterogeneous> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::Heterogeneous> { using type = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>; };
#endif

// template<class TypeTag>
// struct Grid<TypeTag, TTag::Heterogeneous> { using type = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Heterogeneous> { using type = HeterogeneousProblem<TypeTag>; };

// // Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Heterogeneous>
{
    using type = HeterogeneousSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                            GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Heterogeneous>
{
    using type = FluidSystems::BrineCO2<GetPropType<TypeTag, Properties::Scalar>,
                                        Components::CO2<GetPropType<TypeTag, Properties::Scalar>, GeneratedCO2Tables::CO2Tables>,
                                        Components::TabulatedComponent<Components::H2O<GetPropType<TypeTag, Properties::Scalar>>>,
                                        FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/false>>;
};

// Use Moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Heterogeneous> { static constexpr bool value = true; };

template<class TypeTag>
struct PointSource<TypeTag, TTag::Heterogeneous> { 
    using type = IdPointSource<Dune::FieldVector<typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::ctype, GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld>,
                               Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>, 
                               int>; 
};

#if !ISOTHERMAL
namespace TTag {
struct HeterogeneousNI { using InheritsFrom = std::tuple<TwoPTwoCCO2NI, BoxModel>; };
}

#if !THREEDIMENSIONAL
template<class TypeTag>
struct Grid<TypeTag, TTag::HeterogeneousNI> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::HeterogeneousNI> { using type = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>; };
#endif

// template<class TypeTag>
// struct Grid<TypeTag, TTag::HeterogeneousNI> { using type = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::HeterogeneousNI> { using type = HeterogeneousProblem<TypeTag>; };

// // Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HeterogeneousNI>
{
    using type = HeterogeneousSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                            GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::HeterogeneousNI>
{
    using type = FluidSystems::BrineCO2<GetPropType<TypeTag, Properties::Scalar>,
                                        Components::CO2<GetPropType<TypeTag, Properties::Scalar>, GeneratedCO2Tables::CO2Tables>,
                                        Components::TabulatedComponent<Components::H2O<GetPropType<TypeTag, Properties::Scalar>>>,
                                        FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/true>>;
};

// // We set the solidSystem  used for our simulation
// template<class TypeTag>
// struct SolidSystem<TypeTag, TTag::HeterogeneousNI>
// {
//     using Scalar = GetPropType<TypeTag, Properties::Scalar>;
//     using type = SolidSystems::FaciesSolidSystem<Scalar, Components::Facies1<Scalar>, Components::Facies2<Scalar>>;
// };

// Use Moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::HeterogeneousNI> { static constexpr bool value = false; };

template<class TypeTag>
struct PointSource<TypeTag, TTag::HeterogeneousNI> { 
    using type = IdPointSource<Dune::FieldVector<typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::ctype, GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld>,
                                                Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>, 
                                                int>; 
};
#endif

}

#endif