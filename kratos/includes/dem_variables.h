//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//











#if !defined(KRATOS_DEM_VARIABLES_H_INCLUDED )
#define  KRATOS_DEM_VARIABLES_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "includes/kratos_components.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "containers/weak_pointer_vector.h"
#include "containers/periodic_variables_container.h"

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

//TODO: move to the Kratos DEM_application or eventually to the FluidDynamicsAsNeeded
namespace Kratos
{
     //for DEM Application:
    KRATOS_DEFINE_VARIABLE( double, PARTICLE_MASS )
    KRATOS_DEFINE_VARIABLE( double, RADIUS )
    KRATOS_DEFINE_VARIABLE( double, SEARCH_TOLERANCE )
    KRATOS_DEFINE_VARIABLE( double, AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION )
    KRATOS_DEFINE_VARIABLE( double, DEM_DELTA_TIME )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( TOTAL_FORCES )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DAMP_FORCES )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_VEL )

    KRATOS_DEFINE_VARIABLE( Vector, NEIGHBOURS_IDS_DOUBLE )
    KRATOS_DEFINE_VARIABLE( vector<double>, BASSET_HISTORIC_INTEGRANDS )
    KRATOS_DEFINE_VARIABLE( vector<double>, HINSBERG_TAIL_CONTRIBUTIONS )
    KRATOS_DEFINE_VARIABLE( Vector, PARTICLE_ROTATE_SPRING_FAILURE_TYPE )

    KRATOS_DEFINE_VARIABLE( vector<int>, OLD_NEIGHBOURS_IDS )
    KRATOS_DEFINE_VARIABLE( vector<int>, INI_NEIGHBOURS_IDS )
    KRATOS_DEFINE_VARIABLE( vector<int>, CONTINUUM_INI_NEIGHBOURS_IDS )
    KRATOS_DEFINE_VARIABLE( vector<int>, NEIGHBOURS_IDS )
    KRATOS_DEFINE_VARIABLE( vector<int>, PARTICLE_INITIAL_FAILURE_ID )
    KRATOS_DEFINE_VARIABLE( vector<int>, CONTINUUM_PARTICLE_INITIAL_FAILURE_ID )

    KRATOS_DEFINE_VARIABLE( int, FIXED_MESH_OPTION)
    KRATOS_DEFINE_VARIABLE( int, PARTICLE_MATERIAL )

    KRATOS_DEFINE_VARIABLE( std::string, ELEMENT_TYPE )

    typedef vector<array_1d<double,3> > VectorArray3Double;
    KRATOS_DEFINE_VARIABLE( VectorArray3Double, PARTICLE_ROTATE_SPRING_MOMENT )

    // Swimming DEM Application BEGINNING
    KRATOS_DEFINE_VARIABLE( int, NUMBER_OF_INIT_BASSET_STEPS )
    KRATOS_DEFINE_VARIABLE( int, COUPLING_TYPE)
    KRATOS_DEFINE_VARIABLE( int, NON_NEWTONIAN_OPTION )
    KRATOS_DEFINE_VARIABLE( int, MANUALLY_IMPOSED_DRAG_LAW_OPTION )
    KRATOS_DEFINE_VARIABLE( int, DRAG_MODIFIER_TYPE )
    KRATOS_DEFINE_VARIABLE( int, BUOYANCY_FORCE_TYPE )
    KRATOS_DEFINE_VARIABLE( int, DRAG_FORCE_TYPE )
    KRATOS_DEFINE_VARIABLE( int, VIRTUAL_MASS_FORCE_TYPE )
    KRATOS_DEFINE_VARIABLE( int, BASSET_FORCE_TYPE )
    KRATOS_DEFINE_VARIABLE( int, LIFT_FORCE_TYPE )
    KRATOS_DEFINE_VARIABLE( int, MAGNUS_FORCE_TYPE )
    KRATOS_DEFINE_VARIABLE( int, HYDRO_TORQUE_TYPE )
    KRATOS_DEFINE_VARIABLE( int, FLUID_MODEL_TYPE )
    KRATOS_DEFINE_VARIABLE( int, DRAG_POROSITY_CORRECTION_TYPE )
    KRATOS_DEFINE_VARIABLE( int, TIME_STEPS_PER_QUADRATURE_STEP )
    KRATOS_DEFINE_VARIABLE( int, QUADRATURE_ORDER )
    KRATOS_DEFINE_VARIABLE( double, POWER_LAW_TOLERANCE )
    KRATOS_DEFINE_VARIABLE( double, PARTICLE_SPHERICITY )
    KRATOS_DEFINE_VARIABLE( double, INIT_DRAG_FORCE )
    KRATOS_DEFINE_VARIABLE( double, DRAG_LAW_SLOPE )
    KRATOS_DEFINE_VARIABLE( double, SOLID_FRACTION )
    KRATOS_DEFINE_VARIABLE( double, SOLID_FRACTION_RATE )
    KRATOS_DEFINE_VARIABLE( double, FLUID_FRACTION )
    KRATOS_DEFINE_VARIABLE( double, FLUID_FRACTION_RATE )
    KRATOS_DEFINE_VARIABLE( double, PHASE_FRACTION )
    KRATOS_DEFINE_VARIABLE( double, PHASE_FRACTION_RATE )
    KRATOS_DEFINE_VARIABLE( double, SOLID_FRACTION_PROJECTED )
    KRATOS_DEFINE_VARIABLE( double, FLUID_FRACTION_PROJECTED )
    KRATOS_DEFINE_VARIABLE( double, FLUID_DENSITY_PROJECTED )
    KRATOS_DEFINE_VARIABLE( double, FLUID_VISCOSITY_PROJECTED )
    KRATOS_DEFINE_VARIABLE( double, REYNOLDS_NUMBER )
    KRATOS_DEFINE_VARIABLE( double, DRAG_COEFFICIENT )
    KRATOS_DEFINE_VARIABLE( double, SHEAR_RATE_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_OLD )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_OLD_OLD )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SLIP_VELOCITY )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ACCELERATION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ADDITIONAL_FORCE )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ADDITIONAL_FORCE_OLD )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ADDITIONAL_FORCE_OLD_OLD )
    KRATOS_DEFINE_VARIABLE( double, DELTA_TIME_QUADRATURE )
    KRATOS_DEFINE_VARIABLE( double, LAST_TIME_APPENDING )


            
    //SINTERING VARIABLES
    KRATOS_DEFINE_VARIABLE( double, ATOMIC_VOLUME )
    KRATOS_DEFINE_VARIABLE( double, SURFACE_ENERGY )
    KRATOS_DEFINE_VARIABLE( double, DIHEDRAL_ANGLE )
    KRATOS_DEFINE_VARIABLE( double, SINTERING_START_TEMPERATURE )
    KRATOS_DEFINE_VARIABLE( double, RELAXATION_TIME )
    KRATOS_DEFINE_VARIABLE( double, LARGE_VISCOSITY_COEFFICIENT )
    KRATOS_DEFINE_VARIABLE( double, THERMAL_EXPANSION_COEFFICIENT )
    KRATOS_DEFINE_VARIABLE( double, PRE_EXP_DIFFUSION_COEFFICIENT )
    KRATOS_DEFINE_VARIABLE( double, GB_WIDTH )
    KRATOS_DEFINE_VARIABLE( double, ENTHAPLY_ACTIVATION )
    ///
            
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( HYDRODYNAMIC_FORCE )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( HYDRODYNAMIC_MOMENT )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_VEL_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_VEL_PROJECTED_RATE )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_VEL_LAPL_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_VEL_LAPL_RATE_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_ACCEL_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_FLUID_ACCEL_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_VORTICITY_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( HYDRODYNAMIC_REACTION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MEAN_HYDRODYNAMIC_REACTION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DRAG_REACTION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LIFT_FORCE )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VIRTUAL_MASS_FORCE )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( BASSET_FORCE )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( BUOYANCY )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PRESSURE_GRADIENT )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PRESSURE_GRAD_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SOLID_FRACTION_GRADIENT )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SOLID_FRACTION_GRADIENT_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_FRACTION_GRADIENT )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_FRACTION_GRADIENT_PROJECTED )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PHASE_FRACTION_GRADIENT )
    // Swimming DEM Application END

}  // namespace Kratos.

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_DEM_VARIABLES_H_INCLUDED  defined