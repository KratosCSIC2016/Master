/* 
 * File:   DEM_application_variables.h
 * Author: Salva Latorre
 *
 * Created on October 9, 2014, 10:54 AM
 */

#ifndef KRATOS_DEM_APPLICATION_VARIABLES_H
#define	KRATOS_DEM_APPLICATION_VARIABLES_H

#include "includes/define.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"
#include "includes/condition.h"
#include "utilities/quaternion.h"
#include "custom_utilities/cluster_information.h"

namespace Kratos
{

#define DEM_COPY_SECOND_TO_FIRST_3(a, b)            a[0]  = b[0]; a[1]  = b[1]; a[2]  = b[2];
#define DEM_COPY_SECOND_TO_FIRST_4(a, b)            a[0]  = b[0]; a[1]  = b[1]; a[2]  = b[2]; a[3]  = b[3];
#define DEM_ADD_SECOND_TO_FIRST(a, b)               a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
#define DEM_SET_COMPONENTS_TO_ZERO_3(a)             a[0]  = 0.0;  a[1]  = 0.0;  a[2]  = 0.0;
#define DEM_SET_COMPONENTS_TO_ZERO_3x3(a)           a[0][0] = 0.0; a[0][1] = 0.0; a[0][2] = 0.0; a[1][0] = 0.0; a[1][1] = 0.0; a[1][2] = 0.0; a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 0.0;
#define DEM_MULTIPLY_BY_SCALAR_3(a, b)              a[0] = b * a[0]; a[1] = b * a[1]; a[2] = b * a[2];
#define DEM_MODULUS_3(a)                            sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
#define DEM_MODULUS_2(a)                            sqrt(a[0] * a[0] + a[1] * a[1])
#define DEM_INNER_PRODUCT_3(a, b)                       (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
#define DEM_SET_TO_CROSS_OF_FIRST_TWO_3(a, b, c)    c[0] = a[1] * b[2] - a[2] * b[1]; c[1] = a[2] * b[0] - a[0] * b[2]; c[2] = a[0] * b[1] - a[1] * b[0];
#define DEM_COPY_SECOND_TO_FIRST_3x3(a, b)          a[0][0] = b[0][0]; a[0][1] = b[0][1]; a[0][2] = b[0][2]; \
                                                    a[1][0] = b[1][0]; a[1][1] = b[1][1]; a[1][2] = b[1][2]; \
                                                    a[2][0] = b[2][0]; a[2][1] = b[2][1]; a[2][2] = b[2][2];
    
  KRATOS_DEFINE_VARIABLE(WeakPointerVector< Element >, CONTINUUM_INI_NEIGHBOUR_ELEMENTS)
  KRATOS_DEFINE_VARIABLE(WeakPointerVector< Element >, NODE_TO_NEIGH_ELEMENT_POINTER)

  //constitutive law
  KRATOS_DEFINE_VARIABLE(std::string, DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME)
  KRATOS_DEFINE_VARIABLE(std::string, DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME)
  //KRATOS_DEFINE_VARIABLE( DEMConstitutiveLaw::Pointer, DEM_CONSTITUTIVE_LAW_POINTER )
   
  KRATOS_DEFINE_VARIABLE(std::string, PROBABILITY_DISTRIBUTION)
  KRATOS_DEFINE_VARIABLE(std::string, EXCENTRICITY_PROBABILITY_DISTRIBUTION)
  
  // OPTIONS AND FLAGS
  KRATOS_DEFINE_VARIABLE(int, TOP)
  KRATOS_DEFINE_VARIABLE(int, BOTTOM)
  KRATOS_DEFINE_VARIABLE(int, FORCE_INTEGRATION_GROUP)
  KRATOS_DEFINE_VARIABLE(int, BOUNDING_BOX_OPTION)
  KRATOS_DEFINE_VARIABLE(int, ROTATION_OPTION)
  KRATOS_DEFINE_VARIABLE(int, CRITICAL_TIME_OPTION)
  KRATOS_DEFINE_VARIABLE(int, VIRTUAL_MASS_OPTION)
  KRATOS_DEFINE_VARIABLE(int, SEARCH_CONTROL) 
  KRATOS_DEFINE_VARIABLE(double, COORDINATION_NUMBER)
  KRATOS_DEFINE_VARIABLE(double, MAX_AMPLIFICATION_RATIO_OF_THE_SEARCH_RADIUS)
  KRATOS_DEFINE_VARIABLE(vector<int>, SEARCH_CONTROL_VECTOR)
  KRATOS_DEFINE_VARIABLE(int, CLEAN_INDENT_OPTION)
  KRATOS_DEFINE_VARIABLE(int, TRIHEDRON_OPTION)
  KRATOS_DEFINE_VARIABLE(int, ROLLING_FRICTION_OPTION)
  KRATOS_DEFINE_VARIABLE(int, POISSON_EFFECT_OPTION)
  KRATOS_DEFINE_VARIABLE(int, SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION)
  KRATOS_DEFINE_VARIABLE(int, NEIGH_INITIALIZED)
  KRATOS_DEFINE_VARIABLE(int, TRIAXIAL_TEST_OPTION)
  KRATOS_DEFINE_VARIABLE(int, FIX_VELOCITIES_FLAG)
  KRATOS_DEFINE_VARIABLE(int, COMPUTE_STRESS_TENSOR_OPTION)
  KRATOS_DEFINE_VARIABLE(int, PARTICLE_ID)
  KRATOS_DEFINE_VARIABLE(bool, CONTAINS_CLUSTERS)
  KRATOS_DEFINE_VARIABLE(bool, RANDOM_ORIENTATION)
  KRATOS_DEFINE_VARIABLE(int, LOCAL_RESOLUTION_METHOD)
  KRATOS_DEFINE_VARIABLE(int, COMPUTE_FEM_RESULTS_OPTION)
  KRATOS_DEFINE_VARIABLE(int, BREAKABLE_CLUSTER)
  KRATOS_DEFINE_VARIABLE(ClusterInformation, CLUSTER_INFORMATION)
  KRATOS_DEFINE_VARIABLE(std::string, CLUSTER_FILE_NAME)
  
  KRATOS_DEFINE_VARIABLE(double, INITIAL_VELOCITY_X_VALUE)
  KRATOS_DEFINE_VARIABLE(double, INITIAL_VELOCITY_Y_VALUE)
  KRATOS_DEFINE_VARIABLE(double, INITIAL_VELOCITY_Z_VALUE)
  KRATOS_DEFINE_VARIABLE(double, INITIAL_ANGULAR_VELOCITY_X_VALUE)
  KRATOS_DEFINE_VARIABLE(double, INITIAL_ANGULAR_VELOCITY_Y_VALUE)
  KRATOS_DEFINE_VARIABLE(double, INITIAL_ANGULAR_VELOCITY_Z_VALUE)

  // *************** Continuum only BEGIN *************
  KRATOS_DEFINE_VARIABLE(bool, DELTA_OPTION)
  KRATOS_DEFINE_VARIABLE(int, CASE_OPTION)
  KRATOS_DEFINE_VARIABLE(double, SKIN_SPHERE)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_COHESION)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_TENSION)
  
  KRATOS_DEFINE_VARIABLE(int, PROPERTIES_ID)
  KRATOS_DEFINE_VARIABLE(int, CONTACT_MESH_OPTION)
  //KRATOS_DEFINE_VARIABLE(int, FAILURE_CRITERION_OPTION)
  KRATOS_DEFINE_VARIABLE(int, CONCRETE_TEST_OPTION)
  KRATOS_DEFINE_VARIABLE(int, COHESIVE_GROUP)
  KRATOS_DEFINE_VARIABLE(int, IF_BOUNDARY_ELEMENT)
  KRATOS_DEFINE_VARIABLE(Vector, IF_BOUNDARY_FACE)
  KRATOS_DEFINE_VARIABLE(vector<int>, PARTICLE_CONTACT_FAILURE_ID)
  // *************** Continuum only END ***************

  // MATERIAL PARAMETERS

  KRATOS_DEFINE_VARIABLE(double, NODAL_MASS_COEFF)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_MOMENT_OF_INERTIA)
  KRATOS_DEFINE_VARIABLE(double, ROLLING_FRICTION)
  KRATOS_DEFINE_VARIABLE(double, HISTORICAL_MIN_K)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_INERTIA)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_DENSITY)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_FRICTION)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_STATIC_FRICTION_COEF)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_DYNAMIC_FRICTION_COEF)
  KRATOS_DEFINE_VARIABLE(double, COEFFICIENT_OF_RESTITUTION)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_ROTATION_DAMP_RATIO)
  KRATOS_DEFINE_VARIABLE(double, DAMPING_GAMMA)
  KRATOS_DEFINE_VARIABLE(double, K_NORMAL)
  KRATOS_DEFINE_VARIABLE(double, K_TANGENTIAL)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_RADIUS)
  KRATOS_DEFINE_VARIABLE(double, MAX_STRESS)
  KRATOS_DEFINE_VARIABLE(double, ALPHA)
  KRATOS_DEFINE_VARIABLE(double, ALPHA_FUNCTION)
  KRATOS_DEFINE_VARIABLE(double, GAMMA)
  KRATOS_DEFINE_VARIABLE(double, EXCENTRICITY)
  KRATOS_DEFINE_VARIABLE(double, EXCENTRICITY_STANDARD_DEVIATION)
  KRATOS_DEFINE_VARIABLE(double, FABRIC_COEFFICIENT)
  KRATOS_DEFINE_VARIABLE(double, POISSON_VALUE)
  KRATOS_DEFINE_VARIABLE(double, KT_FACTOR)
  KRATOS_DEFINE_VARIABLE(double, INTERNAL_COHESION)
  
  // *************** Nano-particle only BEGIN *************
  KRATOS_DEFINE_VARIABLE(double, CATION_CONCENTRATION)
  // *************** Nano-particle only END *************

  // *************** Continuum only BEGIN *************
  KRATOS_DEFINE_VARIABLE(double, SLOPE_FRACTION_N1)
  KRATOS_DEFINE_VARIABLE(double, SLOPE_FRACTION_N2)
  KRATOS_DEFINE_VARIABLE(double, SLOPE_FRACTION_N3)
  KRATOS_DEFINE_VARIABLE(double, SLOPE_LIMIT_COEFF_C1)
  KRATOS_DEFINE_VARIABLE(double, SLOPE_LIMIT_COEFF_C2)
  KRATOS_DEFINE_VARIABLE(double, SLOPE_LIMIT_COEFF_C3)
  KRATOS_DEFINE_VARIABLE(double, YOUNG_MODULUS_PLASTIC)
  KRATOS_DEFINE_VARIABLE(double, PLASTIC_YIELD_STRESS)
  KRATOS_DEFINE_VARIABLE(double, DAMAGE_FACTOR)
  KRATOS_DEFINE_VARIABLE(double, SHEAR_ENERGY_COEF)
  KRATOS_DEFINE_VARIABLE(double, DONZE_G1)
  KRATOS_DEFINE_VARIABLE(double, DONZE_G2)
  KRATOS_DEFINE_VARIABLE(double, DONZE_G3)
  KRATOS_DEFINE_VARIABLE(double, DONZE_MAX_DEF)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_FAILURE)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_ORIENTATION)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_SIGMA)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_TAU)
  KRATOS_DEFINE_VARIABLE(double, FAILURE_CRITERION_STATE)
  KRATOS_DEFINE_VARIABLE(double, UNIDIMENSIONAL_DAMAGE)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_SIGMA_MIN)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_TAU_ZERO)
  KRATOS_DEFINE_VARIABLE(double, CONTACT_INTERNAL_FRICC)

  // *************** Continuum only END *************

  // GEOMETRIC PARAMETERS

  // *************** Continuum only BEGIN *************
  KRATOS_DEFINE_VARIABLE(double, LOCAL_CONTACT_AREA_HIGH)
  KRATOS_DEFINE_VARIABLE(double, LOCAL_CONTACT_AREA_LOW)
  KRATOS_DEFINE_VARIABLE(double, MEAN_CONTACT_AREA)
  KRATOS_DEFINE_VARIABLE(double, REPRESENTATIVE_VOLUME)
  KRATOS_DEFINE_VARIABLE(boost::numeric::ublas::vector<int>,  NEIGHBOUR_IDS)
  KRATOS_DEFINE_VARIABLE(Vector,  NEIGHBOURS_CONTACT_AREAS)
  // *************** Continuum only END ***************
  
  // INLET PARAMETERS
    
  KRATOS_DEFINE_VARIABLE(double,INLET_START_TIME)    
  KRATOS_DEFINE_VARIABLE(double,INLET_STOP_TIME)
  KRATOS_DEFINE_VARIABLE(double,INLET_NUMBER_OF_PARTICLES)
  KRATOS_DEFINE_VARIABLE(double,STANDARD_DEVIATION)
  KRATOS_DEFINE_VARIABLE(double,MAX_RAND_DEVIATION_ANGLE)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LINEAR_VELOCITY)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(INLET_INITIAL_VELOCITY)

  // KINEMATICS

  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_ROTATION_ANGLE)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(EULER_ANGLES)
  KRATOS_DEFINE_VARIABLE(double,ORIENTATION_REAL) // JIG: SHOULD BE REMOVED IN THE FUTURE
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ORIENTATION_IMAG) // JIG: SHOULD BE REMOVED IN THE FUTURE
  KRATOS_DEFINE_VARIABLE(Quaternion<double>, ORIENTATION)
  KRATOS_DEFINE_VARIABLE(Quaternion<double>, AUX_ORIENTATION)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_AUX_ANGULAR_VELOCITY)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DELTA_DISPLACEMENT)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DELTA_ROTA_DISPLACEMENT)
  KRATOS_DEFINE_VARIABLE(double, VELOCITY_START_TIME) 
  KRATOS_DEFINE_VARIABLE(double, VELOCITY_STOP_TIME)
  KRATOS_DEFINE_VARIABLE(double, ANGULAR_VELOCITY_START_TIME) 
  KRATOS_DEFINE_VARIABLE(double, ANGULAR_VELOCITY_STOP_TIME)
  KRATOS_DEFINE_VARIABLE(int, RIGID_BODY_MOTION)


  // FORCE AND MOMENTUM

  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_MOMENT)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ROLLING_RESISTANCE_MOMENT)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAX_ROTA_MOMENT)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ELASTIC_FORCES)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONTACT_FORCES)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RIGID_ELEMENT_FORCE)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(TANGENTIAL_ELASTIC_FORCES)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(FORCE_REACTION)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MOMENT_REACTION)
  KRATOS_DEFINE_VARIABLE(double, DEM_PRESSURE) 
  KRATOS_DEFINE_VARIABLE(double, DEM_NODAL_AREA)
    
  // ENERGY

  KRATOS_DEFINE_VARIABLE(double, PARTICLE_ELASTIC_ENERGY)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_ROTATIONAL_KINEMATIC_ENERGY)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_GRAVITATIONAL_ENERGY)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_INELASTIC_VISCODAMPING_ENERGY)
  KRATOS_DEFINE_VARIABLE(double, PARTICLE_INELASTIC_FRICTIONAL_ENERGY)
  KRATOS_DEFINE_VARIABLE(int, COMPUTE_ENERGY_OPTION)
  
  // *************** Continuum only BEGIN *************
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(INITIAL_ROTA_MOMENT)
  KRATOS_DEFINE_VARIABLE(Vector, PARTICLE_BLOCK_CONTACT_FORCE)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(EXTERNAL_APPLIED_FORCE)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(EXTERNAL_APPLIED_MOMENT)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_CONTACT_FORCE)
  KRATOS_DEFINE_VARIABLE(VectorArray3Double, PARTICLE_CONTACT_FORCES )
  KRATOS_DEFINE_VARIABLE(double, NEIGHBOUR_SIZE)


  // CONCRETE TEST

  KRATOS_DEFINE_VARIABLE(double, FIXED_VEL_TOP)
  KRATOS_DEFINE_VARIABLE(double, FIXED_VEL_BOT)
  KRATOS_DEFINE_VARIABLE(double, AREA_VERTICAL_TAPA)
  KRATOS_DEFINE_VARIABLE(double, AREA_VERTICAL_CENTRE)

  // TENSION

  KRATOS_DEFINE_APPLICATION_VARIABLE( DEM_APPLICATION, Matrix, DEM_STRESS_TENSOR )

  // APPLIED LOADS

  KRATOS_DEFINE_VARIABLE(double, BLAST_RADIUS)
  KRATOS_DEFINE_VARIABLE(int   , BLAST_CURVE)
  KRATOS_DEFINE_VARIABLE(double, BLAST_PRESSURE_MAX)
  KRATOS_DEFINE_VARIABLE(double, BLAST_TIME_PRESSURE_MAX)
  KRATOS_DEFINE_VARIABLE(double, BLAST_SHAPE_FACTOR)
  KRATOS_DEFINE_VARIABLE(double, BLAST_TIME_DELAY)
  KRATOS_DEFINE_VARIABLE(int   , BLAST_BOREHOLE)
  KRATOS_DEFINE_VARIABLE(int   , BLAST_NPOINTS)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_1)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_2)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_3)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_4)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_5)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_6)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_7)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BLAST_COORDINATES_8)
  // *************** Continuum only END *************

  // Possible future blocks (no FEM) interaction

  KRATOS_DEFINE_VARIABLE(Vector, PARTICLE_BLOCK_CONTACT_FAILURE_ID)
  KRATOS_DEFINE_VARIABLE(Vector, PARTICLE_BLOCK_IF_INITIAL_CONTACT)
  KRATOS_DEFINE_VARIABLE(WeakPointerVector<Element >,    NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)
  KRATOS_DEFINE_VARIABLE(WeakPointerVector<Condition >,  NEIGHBOUR_RIGID_FACES)
  KRATOS_DEFINE_VARIABLE(WeakPointerVector<Element >,    NEIGHBOUR_PARTICLE_OF_RIGID_FACE)
  KRATOS_DEFINE_VARIABLE(Vector,  NEIGHBOUR_RIGID_FACES_PRAM)
  KRATOS_DEFINE_VARIABLE(Vector,  NEIGHBOUR_RIGID_FACES_ELASTIC_CONTACT_FORCE)
  KRATOS_DEFINE_VARIABLE(Vector,  NEIGHBOUR_RIGID_FACES_TOTAL_CONTACT_FORCE)

  // DUMMIES INT AND DOUBLE VARIABLES
  KRATOS_DEFINE_VARIABLE(int, DUMMY_SWITCH)
  
  // EXPORTS
  
  KRATOS_DEFINE_VARIABLE(double, EXPORT_ID)
  KRATOS_DEFINE_VARIABLE(double, EXPORT_PARTICLE_FAILURE_ID)
  KRATOS_DEFINE_VARIABLE(int, PRINT_EXPORT_ID)
  KRATOS_DEFINE_VARIABLE(int, PRINT_STRESS_TENSOR_OPTION)
  
  // For DEM_FEM Element
  
  KRATOS_DEFINE_VARIABLE(double, LOCAL_DAMP_RATIO)
  
  // For the DEM_Clusters Element
  
  KRATOS_DEFINE_VARIABLE(double, CLUSTER_VOLUME)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRINCIPAL_MOMENTS_OF_INERTIA)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_ANGULAR_VELOCITY)
  KRATOS_DEFINE_VARIABLE(double, CHARACTERISTIC_LENGTH)
  KRATOS_DEFINE_VARIABLE(double, SPRAYED_MATERIAL)

  // DUMMY VARIABLES FOR CALCULATE
  KRATOS_DEFINE_VARIABLE(double, CALCULATE_COMPUTE_NEW_NEIGHBOURS_HISTORICAL_DATA)
  KRATOS_DEFINE_VARIABLE(double, CALCULATE_COMPUTE_NEW_RIGID_FACE_NEIGHBOURS_HISTORICAL_DATA)
  KRATOS_DEFINE_VARIABLE(double, CALCULATE_SET_INITIAL_DEM_CONTACTS)
  KRATOS_DEFINE_VARIABLE(double, CALCULATE_SET_INITIAL_FEM_CONTACTS)

  //Cfeng,131013,RigidFace
  
  KRATOS_DEFINE_VARIABLE(double, RIGID_FACE_ROTA_SPEED)
  KRATOS_DEFINE_VARIABLE(double, RIGID_FACE_AXIAL_SPEED)
  KRATOS_DEFINE_VARIABLE(int,    RIGID_FACE_PROP_ID)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RIGID_FACE_ROTA_ORIGIN_COORD)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RIGID_FACE_ROTA_AXIAL_DIR)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RIGID_FACE_ROTA_GLOBAL_VELOCITY)
  KRATOS_DEFINE_VARIABLE(double, RIGID_FACE_BEGIN_TIME)
  KRATOS_DEFINE_VARIABLE(double, RIGID_FACE_END_TIME)
  KRATOS_DEFINE_VARIABLE(int, RIGID_FACE_FLAG)
  KRATOS_DEFINE_VARIABLE(Vector, RIGID_FACE_COMPUTE_MOVEMENT)
  
  //SLS DEM-FEM
  KRATOS_DEFINE_VARIABLE(double, WALL_FRICTION)
  KRATOS_DEFINE_VARIABLE(double, SHEAR_STRESS)
  KRATOS_DEFINE_VARIABLE(double, NON_DIMENSIONAL_VOLUME_WEAR)
  KRATOS_DEFINE_VARIABLE(double, IMPACT_WEAR)
  KRATOS_DEFINE_VARIABLE(double, SEVERITY_OF_WEAR)
  KRATOS_DEFINE_VARIABLE(double, BRINELL_HARDNESS)
  KRATOS_DEFINE_VARIABLE(bool  , COMPUTE_WEAR)
  KRATOS_DEFINE_VARIABLE(double, IMPACT_WEAR_SEVERITY)
  KRATOS_DEFINE_VARIABLE(double, WALL_COHESION)
  //BOUNDING BOX
  KRATOS_DEFINE_VARIABLE(double, BOUNDING_BOX_START_TIME)
  KRATOS_DEFINE_VARIABLE(double, BOUNDING_BOX_STOP_TIME)
  
  //OPTIMIZATION 
  KRATOS_DEFINE_VARIABLE(double, TOTAL_CONTACT_DISTANCES)

  // *************** Thermal only BEGIN *************
  KRATOS_DEFINE_VARIABLE(double, HEATFLUX)
  KRATOS_DEFINE_VARIABLE(double, THERMAL_CONDUCTIVITY)
  // *************** Thermal only END ***************  
  
class DEMFlags
  {
  public:
    KRATOS_DEFINE_LOCAL_FLAG(HAS_ROTATION);
    KRATOS_DEFINE_LOCAL_FLAG(IS_SINTERING);
    KRATOS_DEFINE_LOCAL_FLAG(HAS_ROLLING_FRICTION);
    KRATOS_DEFINE_LOCAL_FLAG(HAS_CRITICAL_TIME);  
    KRATOS_DEFINE_LOCAL_FLAG(FIXED_VEL_X); 
    KRATOS_DEFINE_LOCAL_FLAG(FIXED_VEL_Y); 
    KRATOS_DEFINE_LOCAL_FLAG(FIXED_VEL_Z); 
    KRATOS_DEFINE_LOCAL_FLAG(FIXED_ANG_VEL_X); 
    KRATOS_DEFINE_LOCAL_FLAG(FIXED_ANG_VEL_Y); 
    KRATOS_DEFINE_LOCAL_FLAG(FIXED_ANG_VEL_Z);
    KRATOS_DEFINE_LOCAL_FLAG(BELONGS_TO_A_CLUSTER);
    KRATOS_DEFINE_LOCAL_FLAG(HAS_STRESS_TENSOR);	
    KRATOS_DEFINE_LOCAL_FLAG(COPIED_STRESS_TENSOR);
    KRATOS_DEFINE_LOCAL_FLAG(COPIED_STRESS_TENSOR2);
    KRATOS_DEFINE_LOCAL_FLAG(PRINT_STRESS_TENSOR);
  };  
}

#endif	/* KRATOS_DEM_APPLICATION_VARIABLES_H */
