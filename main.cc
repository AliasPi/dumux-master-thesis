#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <filesystem>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include "properties.hh"
#include "masteroftime.hh"
#include "shiftedgridgeometry.hh"
int main(int argc, char** argv){
    using namespace Dumux;

    Dumux::initialize(argc, argv);

    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    // use the type tag defined in properties.hh (model and discretization)
    using TypeTag = Properties::TTag::TYPETAG;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();
    using BaseGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    // Define GridGeometry as the shifted version.
    using GridGeometry = ShiftedGridGeometry<BaseGridGeometry>;

    // Now create the shifted grid geometry:
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, 2000.0);

    std::cout << "Initialize Spatial Params" << std::endl;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto gridData = gridManager.getGridData();
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry, gridData);

    std::cout << "Initialize Problem" << std::endl;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);
    problem->computePointSourceMap();

    std::cout << "Solution Vector" << std::endl;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    std::cout << "Grid Variables" << std::endl;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    std::cout << "Time Loop Setup" << std::endl;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    std::cout << "VTK Writer" << std::endl;
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    std::filesystem::path fullPath(problem->name());
    std::filesystem::path dirPath;
    std::string baseName = fullPath.filename().string();
    if(fullPath.is_absolute()){
        dirPath = fullPath.parent_path();
    } else {
        dirPath = std::filesystem::current_path() / fullPath.parent_path();
    }
    std::filesystem::create_directories(dirPath);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, (dirPath / baseName).string());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    problem->addFieldsToWriter(vtkWriter); // Add custom fieldds to writer
    vtkWriter.write(0.0);

    std::cout << "Initialize Time Loop" << std::endl;
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    auto masterOfTime = std::make_shared<MasterOfTime<TypeTag, Scalar>>(*timeLoop, *problem);

    std::cout << "Assembler" << std::endl;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    std::cout << "Linear Solver" << std::endl;
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    std::cout << "Newton Solver" << std::endl;
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    

    std::cout << "Starting Time Loop" << std::endl;
    timeLoop->start(); do
    {       
        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        masterOfTime->advanceTimeStep();
        Scalar time = masterOfTime->time();

        vtkWriter.write(time);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(false);
    }

    return 0;
}