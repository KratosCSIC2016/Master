//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  VictorBez  $
//   Date:                $Date:  2015         $
//   Revision:            $Revision: 1.0       $
//
//

 
#if !defined(KRATOS_MULTIPHASE_APPLICATION_H_INCLUDED )
#define  KRATOS_MULTIPHASE_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes:none in this case 

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "porous_media_application_variables.h"

#include "custom_elements/flowpressuretrans_2d.h" 
#include "custom_elements/multiphaseve.h"
#include "custom_elements/multiphasefem.h"
#include "custom_conditions/sinksourcenonve.h"
#include "custom_conditions/sinksourcepressure.h"
#include "custom_conditions/sinksourcenonfem.h"
#include "custom_conditions/massflow.h"
#include "custom_conditions/massflowve.h"
#include "custom_conditions/interfaceve_2.h"

#include "includes/variables.h"
#include "includes/condition.h"         //we'll also need conditions for the point heat loads

#include "includes/ublas_interface.h"


using namespace std;

namespace Kratos
{

 
	///@name Kratos Globals

	///@{ 


        //General transport & flow
        KRATOS_DEFINE_VARIABLE(int, IS_FLOW_STATIONARY)
        KRATOS_DEFINE_VARIABLE(int, IS_TRANSPORT_STATIONARY)
        KRATOS_DEFINE_VARIABLE(double, HEAD_LEVEL)
        KRATOS_DEFINE_VARIABLE(double, SPECIFIC_STORAGE) //Eliminar de flowpressTransElem (cambiar por _elem)
        KRATOS_DEFINE_VARIABLE(double, DENSITY_ELEM)
        
        //General transport & flow, multiphase
        //KRATOS_DEFINE_VARIABLE(int, IS_BUOYANCY)
 
        //KRATOS_DEFINE_VARIABLE(double, SATURATION_WET)
        //KRATOS_DEFINE_VARIABLE(double, SATURATION_NON)
        //KRATOS_DEFINE_VARIABLE(double, PRESSURE_WET)
        //KRATOS_DEFINE_VARIABLE(double, PRESSURE_NON)
        //KRATOS_DEFINE_VARIABLE(double, POROSITY_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, PERMEABILITY_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, PERMEABILITY_WET_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, PERMEABILITY_NON_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, VISCOSITY_WET_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, VISCOSITY_NON_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, DENSITY_WET_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, DENSITY_NON_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, DIFFUSION_COEFFICIENT_ELEM)
        KRATOS_DEFINE_VARIABLE(double, PERMEABILITY) //Eliminar de flowpressTransElem (cambiar por _elem)
        KRATOS_DEFINE_VARIABLE(double, DIFFUSION_COEFFICIENT) //Eliminar de flowpressTransElem (cambiar por _elem)
        
        //General multiphase
        //KRATOS_DEFINE_VARIABLE(int, IS_CAPILLARITY_NEGLECTED)
        //KRATOS_DEFINE_VARIABLE(string, CAPILLARITY_PRESSURE_MODEL)
        //KRATOS_DEFINE_VARIABLE(string, RELATIVE_PERMEABILITY_MODEL)
        //KRATOS_DEFINE_VARIABLE(double, RESIDUAL_SW)
        //KRATOS_DEFINE_VARIABLE(double, RESIDUAL_SN)
        //KRATOS_DEFINE_VARIABLE(double, EXPONENT_KR_LAW)
        //KRATOS_DEFINE_VARIABLE(double, LAMBDA_LEVERETT)
        //KRATOS_DEFINE_VARIABLE(double, ENTRY_PRESSURE)       
                
        //KRATOS_DEFINE_VARIABLE(double, SATURATION_WET_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, CAPILLARITY_PRESSURE_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, DERIV_PC_SW_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, DERIV_2_PC_SW_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, DERIV_KRW_SW_ELEM)
        //KRATOS_DEFINE_VARIABLE(double, DERIV_KRN_SW_ELEM) 
        
        //Conditions flow & transport
        KRATOS_DEFINE_VARIABLE(double, LEAKAGE_COEFFICIENT)
        KRATOS_DEFINE_VARIABLE(double, LEVEL)
        KRATOS_DEFINE_VARIABLE(double, PRESCRIBED_VALUE)
        KRATOS_DEFINE_VARIABLE(double, SINK_SOURCE_PRESS)
        KRATOS_DEFINE_VARIABLE(double, SINK_SOURCE)
        
        //Conditions Multiphase
        KRATOS_DEFINE_VARIABLE(double, SINK_SOURCE_PHASE_PRESS)
        KRATOS_DEFINE_VARIABLE(double, SINK_SOURCE_PHASE_SAT)
        KRATOS_DEFINE_VARIABLE(double, SINK_SOURCE_PRESS_Z)
        KRATOS_DEFINE_VARIABLE(double, SINK_SOURCE_HEIGHT)
        
        //Strategy flow & transport
        KRATOS_DEFINE_VARIABLE(double, DARCY_FLOW_X)
        KRATOS_DEFINE_VARIABLE(double, DARCY_FLOW_Y)
        //KRATOS_DEFINE_VARIABLE(double, STORAGE_BALANCE)
        //KRATOS_DEFINE_VARIABLE(double, DARCY_FLOW_BALANCE) 
        //KRATOS_DEFINE_VARIABLE(double, SINKSOURCE_BALANCE)
		KRATOS_DEFINE_VARIABLE(double, PRESSURE_OLD_IT)
        KRATOS_DEFINE_VARIABLE(double, CONCENTRATION_OLD_IT)
		KRATOS_DEFINE_VARIABLE(double, CONCENTRATION)
        
        
        //Strategy multiphase
        //KRATOS_DEFINE_VARIABLE(double, ENTHALPY_WET_NODE)
        //KRATOS_DEFINE_VARIABLE(double, ENTHALPY_NON_NODE)
        
        KRATOS_DEFINE_VARIABLE(double, FIRST_SV_OLD_IT)
        KRATOS_DEFINE_VARIABLE(double, SECOND_SV_OLD_IT)
        
        KRATOS_DEFINE_VARIABLE(double, MAXDIFF_SV_1)
        KRATOS_DEFINE_VARIABLE(double, MAXDIFF_SV_2)
        KRATOS_DEFINE_VARIABLE(double, MAXITER_SV_1)
        KRATOS_DEFINE_VARIABLE(double, MAXITER_SV_2)
        KRATOS_DEFINE_VARIABLE(double, MAXNORMREL_NR)
        KRATOS_DEFINE_VARIABLE(double, MAXITER_NR)
        
        KRATOS_DEFINE_VARIABLE(int, TIMES_TO_PRINT)
        
        // sinksource.cpp condition. Generalization of input (testing)
        KRATOS_DEFINE_VARIABLE(string, SV)
        KRATOS_DEFINE_VARIABLE(string, SV_X)
        KRATOS_DEFINE_VARIABLE(string, SV_Y)
        KRATOS_DEFINE_VARIABLE(string, SV_Z)
        
	class KratosMultiphaseApplication : public KratosApplication
	{
	public:

		/// Pointer definition of KratosMultiphaseApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosMultiphaseApplication);


		/// Default constructor.
		KratosMultiphaseApplication();


		/// Destructor.
		virtual ~KratosMultiphaseApplication(){} 


		virtual void Register();


		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosMultiphaseApplication";
		}


		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}


		///// Print object's data.
      		virtual void PrintData(std::ostream& rOStream) const
      		{
      			KRATOS_WATCH("in the custom Multiphase application");
     		 	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
			rOStream << "Variables:" << std::endl;
			KratosComponents<VariableData>().PrintData(rOStream);
			rOStream << std::endl;
			rOStream << "Elements:" << std::endl;
			KratosComponents<Element>().PrintData(rOStream);
			rOStream << std::endl;
			rOStream << "Conditions:" << std::endl;
			KratosComponents<Condition>().PrintData(rOStream);
      		}	



	protected:


	private:

                const FlowPressureTrans2D  mFlowPressureTrans2D; //and here is our element.
                const MultiphaseVE<1> mMultiphaseVE1D;
                const MultiphaseVE<2> mMultiphaseVE2D;
				const MultiphaseFEM<1> mMultiphaseFEM1D;
				const MultiphaseFEM<2> mMultiphaseFEM2D;
                
                const SinkSourcePressure mPointSinkSourcePressure;
                const SinkSourcePressure mLineSinkSourcePressure;
                //const SinkSourcePressure mTriangleSinkSourcePressure;
               
                const SinkSourceNonVE  mPointSinkSourceNonVE;
                const SinkSourceNonVE  mLineSinkSourceNonVE;
                
                const MassFlow mMassFlow;
                
                const MassFlowVE  mMassFlowVE;
                
                const InterfaceVE_2  mInterfaceVE_2;
				
				const SinkSourceNonFEM  mPointSinkSourceNonFEM;
				const SinkSourceNonFEM  mLineSinkSourceNonFEM;
 

		/// Assignment operator.
		KratosMultiphaseApplication& operator=(KratosMultiphaseApplication  const& rOther);


		/// Copy constructor.
		KratosMultiphaseApplication(KratosMultiphaseApplication const& rOther);

	}; // Class KratosMultiphaseApplication 


}  // namespace Kratos.

#endif // KRATOS_MULTIPHASE_APPLICATION_H_INCLUDED  defined 


