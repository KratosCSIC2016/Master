/*
==============================================================================
KratosMultiphaseApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

//Father of config
#include "custom_strategies/solver_strategy_configuration.h"

//Flow & transport
#include "custom_strategies/fractional_iterative_strategy.h"
//New Multiphase
#include "custom_strategies/sequential_multiphase_fem_strategy.h"

//configuration files
#include "custom_strategies/fractional_iterative_configuration.h"
//New Multiphase
#include "custom_strategies/sequential_multiphase_fem_configuration.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    //typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************
    //
    
    class_< SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
            boost::noncopyable >
            ("SolverStrategyConfiguration", init< ModelPart&, unsigned int>())
            ;

    class_< FractionalIterativeConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType > >,
            boost::noncopyable >
            ("FractionalIterativeConfiguration", init< ModelPart&, LinearSolverType::Pointer, LinearSolverType::Pointer,
             unsigned int >());
			 
	class_< SequentialMultiphaseFEMConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType > >,
            boost::noncopyable >
            ("SequentialMultiphaseFEMConfiguration", init< ModelPart&, LinearSolverType::Pointer, LinearSolverType::Pointer, LinearSolverType::Pointer, LinearSolverType::Pointer,
             unsigned int >()); 

    class_< FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("FractionalIterativeStrategy",
             init < ModelPart&,
             SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >&,
             bool,
             double, double,
             int, int,
             unsigned int, unsigned int,
             bool, bool, unsigned int
             >())
            .def("iterationPressure", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::iterationPressure)
            .def("iterationConcentration", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::iterationConcentration)
            .def("StorageOldSPressIteration", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::StorageOldSPressIteration)
            .def("StorageOldSConcIteration", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::StorageOldSConcIteration)
            .def("isConverged", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::isConverged)
            .def("CalculateFlowBalances", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::CalculateFlowBalances)
            .def("calculateDensityNodes", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateDensityNodes)
            .def("SetEchoLevel", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetEchoLevel)
            .def("Clear", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            .def("AssignInitialStepValues", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssignInitialStepValues)
            .def("PredictSV", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::PredictSV)
            .def("InitializeFractionalStep", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
            .def("CalculateDarcyFlowPressure", &FractionalIterativeStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::CalculateDarcyFlowPressure)
        ;
    
	class_< SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("SequentialMultiphaseFEMStrategy",
             init < ModelPart&,
             SolverStrategyConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >&,
             bool,
             unsigned int, unsigned int,
             bool, bool, unsigned int
             >())
            .def("systemIteration", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::systemIteration)
            //.def("StorageOldIteration", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::StorageOldIteration)   
            //.def("isLocalConverged", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::isLocalConverged)
            .def("checkNumberIterations", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::checkNumberIterations)
            .def("calculateLocalEquations", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateLocalEquations)
            .def("calculateMultiphaseProperties", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateMultiphaseProperties)
            .def("calculateThermoNonWetProperties", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateThermoNonWetProperties)
            .def("calculateThermoWetProperties", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateThermoWetProperties) 
            .def("calculate_Pc", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculate_Pc)
            .def("calculate_kr", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculate_kr)
            .def("calculatePhasesDarcyFlow", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculatePhasesDarcyFlow)
            .def("SetEchoLevel", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetEchoLevel)
            .def("Clear", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            .def("assignInitialStepValues", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::assignInitialStepValues)
            .def("predictSV", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::predictSV)
            //.def("inicializateSV", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::inicializateSV)
            .def("calculateInitialConditions", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateInitialConditions)
            .def("calculateEnthalpyNode", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::calculateEnthalpyNode)    
            .def("InitializeFractionalStep", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
            .def("writeResults", &SequentialMultiphaseFEMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::writeResults)
        ;
}

}  // namespace Python.

} // Namespace Kratos

