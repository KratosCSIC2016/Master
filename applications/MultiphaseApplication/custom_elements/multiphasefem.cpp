#include "multiphasefem.h"
#include "utilities/geometry_utilities.h"

namespace Kratos {

/*
 * public Test<TDim,TNumNodes> functions
 */

template< unsigned int TDim,
          unsigned int TNumNodes >
void MultiphaseFEM<TDim,TNumNodes>::Initialize()
{
//    mSpanWagnerEOS = new SpanWagnerEOS();
//    this->mSpanWagnerEOS->FillTables();
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
			//After calculateLocalEquations and before solving elemental matrix
            this->InterpolateElemFields(rCurrentProcessInfo);
            //this->CalculatePhasesDarcyFlow(rCurrentProcessInfo);
            break;
        }
        case 2:
        {
			//Hace falta??
            //this->InterpolateElemFields(rCurrentProcessInfo);
            //this->CalculatePhasesDarcyFlow(rCurrentProcessInfo);
            break;
        }
        case 3:
        {
			//After calculateLocalEquations in the first step in order to solve 
			//initial CO2 phase injected in the reservoir. NO COMPLETE (only first step).
            //this->InterpolateElemFields(rCurrentProcessInfo);
            //this->CalculatePhasesDarcyFlow(rCurrentProcessInfo);
            break;
        }
        case 4:
        {
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::InterpolateElemFields(ProcessInfo &rCurrentProcessInfo)
{
    //Parameters
    const unsigned int isNonPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_NON_PHASE];
    const unsigned int isWetPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_WET_PHASE];
    const unsigned int isCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
    
    double saturationNonElem = 0.0;
    double permRelWetElem= 0.0;
    double permRelNonElem= 0.0;
    
    double PcElem= 0.0;
    double derivPc_SnElem= 0.0;

    double densityNon_Elem = 0.0;
    double viscosityNon_Elem = 0.0;
    double compressibilityNon_Elem = 0.0;
    double specificHeatNon_Elem = 0.0;
    
    double densityWet_Elem = 0.0;
    double viscosityWet_Elem = 0.0;
    double compressibilityWet_Elem = 0.0;
    double specificHeatWet_Elem = 0.0;
    
    unsigned int numberOfNodes = this->GetGeometry().PointsNumber();
    double weight = 1.0 / ((double) numberOfNodes);

    for (unsigned int node = 0; node < numberOfNodes; node++)
    {  
       const double saturationNonNode = this->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON);
       const double permRelWetNode = this->GetGeometry()[node].FastGetSolutionStepValue(PERMEABILITY_WET_NODE);
       const double permRelNonNode = this->GetGeometry()[node].FastGetSolutionStepValue(PERMEABILITY_NON_NODE);

       saturationNonElem += (weight * saturationNonNode);
       permRelWetElem += (weight * permRelWetNode);
       permRelNonElem += (weight * permRelNonNode);

       if(isCapillarityNeglected !=1)
       {
            const double PcNode = this->GetGeometry()[node].FastGetSolutionStepValue(CAPILLARITY_PRESSURE_NODE);
            const double derivPc_SnNode = this->GetGeometry()[node].FastGetSolutionStepValue(DERIV_PC_SN_NODE);

            PcElem += (weight * PcNode);
            derivPc_SnElem += (weight * derivPc_SnNode);   
       }
       
       //Calculate effective saturation by node
       const double SwrNode  = rCurrentProcessInfo[RESIDUAL_SW];
       const double SnrNode  = rCurrentProcessInfo[RESIDUAL_SN];                       
       const double SeNode = (1.0 - saturationNonNode-SwrNode)/(1.0-SwrNode-SnrNode);
           
       if(isNonPhaseNonConstant !=1)
       {     
           double densityNonNode = 0.0;
           double viscosityNonNode = 0.0;
           double compressibilityNonNode = 0.0;
           double specificHeatNonNode = 0.0;
           
           if( SeNode >= 1.0 ) //Wet Phase
            {
               densityNonNode = 0.0; 
               viscosityNonNode = 0.0; 
               compressibilityNonNode = 0.0;
               specificHeatNonNode = 0.0; 
            }
             else //if( Se_Elem <= 0.00000000001 ) NonWet Phase
            {
               densityNonNode = this->GetGeometry()[node].FastGetSolutionStepValue(DENSITY_NON_NODE); 
               viscosityNonNode = this->GetGeometry()[node].FastGetSolutionStepValue(VISCOSITY_NON_NODE); 
               compressibilityNonNode = this->GetGeometry()[node].FastGetSolutionStepValue(COMPRESSIBLITY_NON_NODE);
               specificHeatNonNode = this->GetGeometry()[node].FastGetSolutionStepValue(SPECIFIC_HEAT_NON_NODE); 
            }
                      
            densityNon_Elem += (weight * densityNonNode);
            viscosityNon_Elem += (weight * viscosityNonNode); 
            compressibilityNon_Elem += (weight * compressibilityNonNode); 
            specificHeatNon_Elem += (weight * specificHeatNonNode);  
       }
       
       if(isWetPhaseNonConstant !=1)
       {     
           double densityWetNode = 0.0;
           double viscosityWetNode = 0.0;
           double compressibilityWetNode = 0.0;
           double specificHeatWetNode = 0.0;
           
           if( SeNode <= 0.0 ) //Wet Phase
            {
               densityWetNode = 0.0; 
               viscosityWetNode = 0.0; 
               compressibilityWetNode = 0.0;
               specificHeatWetNode = 0.0; 
            }
             else //if( Se_Elem <= 0.00000000001 ) NonWet Phase
            {
               densityWetNode = this->GetGeometry()[node].FastGetSolutionStepValue(DENSITY_WET_NODE); 
               viscosityWetNode = this->GetGeometry()[node].FastGetSolutionStepValue(VISCOSITY_WET_NODE); 
               compressibilityWetNode = this->GetGeometry()[node].FastGetSolutionStepValue(COMPRESSIBLITY_WET_NODE);
               specificHeatWetNode = this->GetGeometry()[node].FastGetSolutionStepValue(SPECIFIC_HEAT_WET_NODE); 
            }
                      
            densityWet_Elem += (weight * densityWetNode);
            viscosityWet_Elem += (weight * viscosityWetNode); 
            compressibilityWet_Elem += (weight * compressibilityWetNode); 
            specificHeatWet_Elem += (weight * specificHeatWetNode);  
       }
       
    }

    this->SetValue(SATURATION_NON_ELEM,saturationNonElem);
	this->SetValue(SATURATION_WET_ELEM, (1.0 - saturationNonElem));
    this->SetValue(PERMEABILITY_WET_ELEM,permRelWetElem);
    this->SetValue(PERMEABILITY_NON_ELEM,permRelNonElem);

    if(isCapillarityNeglected !=1)
    {

        this->SetValue(CAPILLARITY_PRESSURE_ELEM,PcElem);
        this->SetValue(DERIV_PC_SN_ELEM,derivPc_SnElem);
    }
    
    if(isNonPhaseNonConstant !=1)
    {
        this->SetValue(DENSITY_NON_ELEM,densityNon_Elem);
        this->SetValue(VISCOSITY_NON_ELEM,viscosityNon_Elem);
        this->SetValue(COMPRESSIBLITY_NON_ELEM,compressibilityNon_Elem);
        this->SetValue(SPECIFIC_HEAT_NON_ELEM,specificHeatNon_Elem);
    }
    
    if(isWetPhaseNonConstant !=1)
    {
        this->SetValue(DENSITY_WET_ELEM,densityWet_Elem);
        this->SetValue(VISCOSITY_WET_ELEM,viscosityWet_Elem);
        this->SetValue(COMPRESSIBLITY_WET_ELEM,compressibilityWet_Elem);
        this->SetValue(SPECIFIC_HEAT_WET_ELEM,specificHeatWet_Elem);
    }
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                            ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
                    case 1:
                    {
                            this->PwEquationIdVector(rResult,rCurrentProcessInfo);
                            break;
                    }
                    case 2:
                    {
                            this->SnEquationIdVector(rResult,rCurrentProcessInfo);
                            break;
                    }
                    case 3:
                    {
                            this->HnEquationIdVector(rResult,rCurrentProcessInfo);
                            break;
                    }
                    case 4:
                    {
                            this->TEquationIdVector(rResult,rCurrentProcessInfo);
                            break;
                    }
                    default:
                    {
                            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index AT THE EquationIdVector: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                    }
		}

    KRATOS_CATCH("");
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::PwEquationIdVector(EquationIdVectorType& rResult,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes,false);	

    for (unsigned int node=0;node< number_of_nodes ;node++)
    {	
        rResult[node] = GetGeometry()[node].GetDof(PRESSURE_WET).EquationId();
    } 
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::SnEquationIdVector(EquationIdVectorType& rResult,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes,false);	

    for (unsigned int node=0;node<number_of_nodes;node++)
    {		
        rResult[node] = GetGeometry()[node].GetDof(SATURATION_NON).EquationId();
    }  
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::HnEquationIdVector(EquationIdVectorType& rResult,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes,false);	

    for (unsigned int node=0;node<number_of_nodes;node++)
    {		
        rResult[node] = GetGeometry()[node].GetDof(ENTHALPY_NON_NODE).EquationId();
    }  
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::TEquationIdVector(EquationIdVectorType& rResult,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes,false);	

    for (unsigned int node=0;node<number_of_nodes;node++)
    {		
        rResult[node] = GetGeometry()[node].GetDof(TEMPERATURE_NODE).EquationId();
    }  
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                      ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
                this->GetPwEquationDofList(rElementalDofList,rCurrentProcessInfo);
                break;
        }
        case 2:
        {
                this->GetSnEquationDofList(rElementalDofList,rCurrentProcessInfo);
                break;
        }
        case 3:
        {
                this->GetHnEquationDofList(rElementalDofList,rCurrentProcessInfo);
                break;
        }
        case 4:
        {
                this->GetTEquationDofList(rElementalDofList,rCurrentProcessInfo);
                break;
        }
        default:
        {
                KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index at the GetDofList: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetPwEquationDofList(DofsVectorType& rElementalDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(rElementalDofList.size() != number_of_nodes)
            rElementalDofList.resize(number_of_nodes);	

    for (unsigned int node=0;node<number_of_nodes;node++)
            rElementalDofList[node] = GetGeometry()[node].pGetDof(PRESSURE_WET);
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetSnEquationDofList(DofsVectorType& rElementalDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();	

    if(rElementalDofList.size() != number_of_nodes)
            rElementalDofList.resize(number_of_nodes);	

    for (unsigned int node=0;node<number_of_nodes;node++)
            rElementalDofList[node] = GetGeometry()[node].pGetDof(SATURATION_NON);
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetHnEquationDofList(DofsVectorType& rElementalDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();	

    if(rElementalDofList.size() != number_of_nodes)
            rElementalDofList.resize(number_of_nodes);	

    for (unsigned int node=0;node<number_of_nodes;node++)
            rElementalDofList[node] = GetGeometry()[node].pGetDof(ENTHALPY_NON_NODE);
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetTEquationDofList(DofsVectorType& rElementalDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();	

    if(rElementalDofList.size() != number_of_nodes)
            rElementalDofList.resize(number_of_nodes);	

    for (unsigned int node=0;node<number_of_nodes;node++)
            rElementalDofList[node] = GetGeometry()[node].pGetDof(TEMPERATURE_NODE);
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    
    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
                this->CalculateSumPhasesEquationPw(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
                break;
        }
        case 2:
        {
                this->CalculateNonWettingEquationSn(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
				break;
        }
        case 3:
        {
                this->CalculateNonWettingEnergyEquationHn(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
                break;
        }
        case 4:
        {
                this->CalculateWettingEnergyEquationT(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
                break;
        }
        default:
        {
                KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index at the CalculateLocalSystem: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    
    MatrixType temp(0,0);

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
                this->CalculateSumPhasesEquationPw(temp,rRightHandSideVector,rCurrentProcessInfo);
                break;
        }
        case 2:
        {
                this->CalculateNonWettingEquationSn(temp,rRightHandSideVector,rCurrentProcessInfo);
				break;
        }
        case 3:
        {
                this->CalculateNonWettingEnergyEquationHn(temp,rRightHandSideVector,rCurrentProcessInfo);
                break;
        }
        case 4:
        {
                this->CalculateWettingEnergyEquationT(temp,rRightHandSideVector,rCurrentProcessInfo);
                break;
        }
        default:
        {
                KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index at the CalculateRightHandSide: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}	

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateNiNj(MatrixType &rNiNj)
{
    rNiNj = ZeroMatrix(TNumNodes, TNumNodes);

    const double lumpedFactor = 1.0 / double(TNumNodes);
    rNiNj = lumpedFactor * IdentityMatrix(TNumNodes, TNumNodes);

}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateNiNj(MatrixType &rNiNj, const VectorType& aN)
{
    rNiNj = ZeroMatrix(TNumNodes, TNumNodes);

    for(SizeType node_i = 0; node_i < TNumNodes; node_i++)
        for (SizeType node_j = 0; node_j < TNumNodes; node_j++)
            rNiNj(node_i,node_j) += aN[node_i] * aN[node_j];
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateBiBj(MatrixType& rBiBj, const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > aDN_DX)
{
    rBiBj = ZeroMatrix(TNumNodes, TNumNodes);

    noalias( rBiBj ) = prod(aDN_DX, trans(aDN_DX));
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateSumPhasesEquationPw(MatrixType& rLeftHandSideMatrix,
                                                                  VectorType& rRightHandSideVector,
                                                                  ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NDof = 1; //Pw
    const SizeType localSize = TNumNodes * NDof;
    const double theta = 1.0;
        
    //resizing as needed the LHS	
    if(rLeftHandSideMatrix.size1() != localSize)	
        rLeftHandSideMatrix.resize(localSize,localSize,false);	
    noalias(rLeftHandSideMatrix) = ZeroMatrix(localSize,localSize); //resetting LHS	
    
    //resizing as needed the RHS
    if(rRightHandSideVector.size() != localSize)	
        rRightHandSideVector.resize(localSize,false);	
    rRightHandSideVector = ZeroVector(localSize); //resetting RHS
    
    //Parameters
    const unsigned int isTransient = rCurrentProcessInfo[IS_TRANSIENT];
    const unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
    const unsigned int isCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
    const unsigned int isNonPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_NON_PHASE];
    const unsigned int isWetPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_WET_PHASE];
    double invDt = 0.0;
    if(isTransient !=1)
    {
        const double deltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
        invDt = 1.0/ deltaTime;
    }
     
    VectorType gravity = ZeroVector(3);
    if(isBuoyancy !=1) gravity = rCurrentProcessInfo[GRAVITY]; 
    
    //Constant
    const double porosity_Elem = GetProperties()[POROSITY_ELEM];
    const double  intrinsecPermeability_Elem = GetProperties()[PERMEABILITY_ELEM];
    const double compressibilityMedium_Elem = GetProperties()[COMPRESSIBILITY_MEDIUM_ELEM];
	const double thermalExpansionWet_Elem = GetProperties()[THERMAL_COEFF_WET_ELEM];
	const double thermalExpansionNon_Elem = GetProperties()[THERMAL_COEFF_NON_ELEM];

    //Non-Constant
    const double saturationNon_Elem = GetValue(SATURATION_NON_ELEM);
    const double saturationWet_Elem = 1.0 - saturationNon_Elem;
    const double permeabilityNon_Elem = GetValue(PERMEABILITY_NON_ELEM);
    const double permeabilityWet_Elem = GetValue(PERMEABILITY_WET_ELEM);
    
    //Calculate capillaryPressureModel (Pc, derivPc_Sw, dderivPc_Sw)
    double derivPc_Sn_Elem = 0.0;
    if(isCapillarityNeglected !=1)
        derivPc_Sn_Elem = GetValue(DERIV_PC_SN_ELEM);

    //Compressible NonWetting phase case
    double densityNon_Elem = 0.0;
    double viscosityNon_Elem = 0.0;
    double compressibilityNon_Elem = 0.0;
    if(isNonPhaseNonConstant !=1)
    {
        densityNon_Elem = GetValue(DENSITY_NON_ELEM);
        viscosityNon_Elem = GetValue(VISCOSITY_NON_ELEM);
        if(viscosityNon_Elem == 0.0) viscosityNon_Elem = 1.0; //Case of wetting phase
        compressibilityNon_Elem = GetValue(COMPRESSIBLITY_NON_ELEM);  
    }
    else
    {
        densityNon_Elem = GetProperties()[DENSITY_NON_ELEM];
        viscosityNon_Elem = GetProperties()[VISCOSITY_NON_ELEM];
		compressibilityNon_Elem = GetProperties()[COMPRESSIBLITY_NON_ELEM];
    }
    
    //Compressible Wetting phase case
    double densityWet_Elem = 0.0; 
    double viscosityWet_Elem = 0.0; 
    double compressibilityWet_Elem = 0.0;
    if(isWetPhaseNonConstant !=1)
    {
        densityWet_Elem = GetValue(DENSITY_WET_ELEM);
        viscosityWet_Elem = GetValue(VISCOSITY_WET_ELEM);
        if(viscosityWet_Elem == 0.0) viscosityWet_Elem = 1.0; //Case of nonwetting phase
        compressibilityWet_Elem = GetValue(COMPRESSIBLITY_WET_ELEM);  
    }
    else
    {
        densityWet_Elem = GetProperties()[DENSITY_WET_ELEM];
        viscosityWet_Elem = GetProperties()[VISCOSITY_WET_ELEM];  
		compressibilityWet_Elem = GetProperties()[COMPRESSIBLITY_WET_ELEM];
    }
    
    //Element variable non-constant (no properties)
    const double mobility_wet_Elem = permeabilityWet_Elem/viscosityWet_Elem;
    const double mobility_non_Elem = permeabilityNon_Elem/viscosityNon_Elem;
    const double totalMobility_Elem = mobility_wet_Elem + mobility_non_Elem;
    
    //State variables: Pwk, Pwk+1 & IncrPw, Snk, Snk+1 & IncrSw, Hnk, Hnk+1 & IncrHn, Tk, Tk+1 & IncrT,
    //array_1d<double,TNumNodes>  
    VectorType pressureWetStepK_1 = ZeroVector(TNumNodes);
    VectorType pressureWetStepK = ZeroVector(TNumNodes);
    VectorType saturationNonStepK_1 = ZeroVector(TNumNodes);
    VectorType saturationNonStepK = ZeroVector(TNumNodes);
    VectorType saturationNonIncrement = ZeroVector(TNumNodes);
    VectorType enthalpyNonIncrement = ZeroVector(TNumNodes);
    VectorType temperatureIncrement = ZeroVector(TNumNodes);
    for(unsigned int node = 0; node< TNumNodes; node++)
    {
        const double Pw_k_1 = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET);
        const double Pw_k = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET,1);
        const double Sn_k_1 = GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON);
        const double Sn_k = GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON, 1);
        const double IncrementSn = Sn_k_1 - Sn_k;
        const double Hn_k_1 = GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE);
        const double Hn_k = GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE, 1);
        const double IncrementHn = Hn_k_1 - Hn_k; 
        const double T_k_1 = GetGeometry()[node].FastGetSolutionStepValue(TEMPERATURE_NODE);
        const double T_k = GetGeometry()[node].FastGetSolutionStepValue(TEMPERATURE_NODE, 1);
        const double IncrementT = T_k_1 - T_k; 

        pressureWetStepK_1[node] = Pw_k_1;
        pressureWetStepK[node] = Pw_k;
        saturationNonStepK_1[node] = Sn_k_1;
        saturationNonStepK[node] = Sn_k;
        saturationNonIncrement[node] = IncrementSn;
        enthalpyNonIncrement[node] = IncrementHn;
        temperatureIncrement[node] = IncrementT;
    }

    ////Element Geometry (N,gradN)
    //area
    double domain_L_A_V = 0.0;
    //GradN
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX = ZeroMatrix(TNumNodes,TDim);
    //N
    array_1d<double, TNumNodes > N = ZeroVector(TNumNodes); 
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, domain_L_A_V); 
    
	//Case of 1D problems
	if (TNumNodes == 2)
		domain_L_A_V = 1.0;

    //Elemental matrix
    //boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes > 
    MatrixType NiNj = ZeroMatrix(TNumNodes,TNumNodes);
    this->CalculateNiNj(NiNj);
    MatrixType BiBj= ZeroMatrix(TNumNodes,TNumNodes);
    this->CalculateBiBj(BiBj, DN_DX);
    
    ////LHS
	//Solid & fluids compressibility storage LHS. Is medium Compressibility negative ??
	const double wet_comp_coeff = saturationWet_Elem*compressibilityWet_Elem;
	const double non_comp_coeff = saturationNon_Elem*compressibilityNon_Elem;
	const double total_comp_coeff = (porosity_Elem*(wet_comp_coeff + wet_comp_coeff) - compressibilityMedium_Elem);
    noalias(rLeftHandSideMatrix) +=     domain_L_A_V * invDt * total_comp_coeff * NiNj;

    //*Total mobility darcyFlow LHS.
    noalias(rLeftHandSideMatrix) +=     theta * domain_L_A_V * intrinsecPermeability_Elem * totalMobility_Elem * BiBj;

    //RHS
	//Solid & fluids compressibility storage RHS. Is medium Compressibility negative ??
	noalias(rRightHandSideVector) +=	domain_L_A_V * invDt * total_comp_coeff  * prod(NiNj, pressureWetStepK);

	//Capillarity storage ONLY RHS.
	const double cap_storage_coeff = ((porosity_Elem * compressibilityNon_Elem) - compressibilityMedium_Elem) * derivPc_Sn_Elem * saturationNon_Elem;
	noalias(rRightHandSideVector) -=	domain_L_A_V * invDt * cap_storage_coeff * prod(NiNj, saturationNonIncrement);

	//Wetting thermal expansion RHS.
    noalias(rRightHandSideVector) -=    domain_L_A_V * invDt * porosity_Elem 
                                        * saturationWet_Elem * thermalExpansionWet_Elem 
                                        * prod(NiNj, temperatureIncrement);
    
	//Non-wetting thermal expansion RHS.
    noalias(rRightHandSideVector) -=    domain_L_A_V * invDt * porosity_Elem 
                                        * saturationNon_Elem * thermalExpansionNon_Elem 
                                        * prod(NiNj, enthalpyNonIncrement);
    
	//Total mobility darcyFlow RHS.
    noalias(rRightHandSideVector) -=    (1.0 - theta) * domain_L_A_V * intrinsecPermeability_Elem * totalMobility_Elem * prod(BiBj,pressureWetStepK);
    
	//*Capillarity darcyFlow RHS.
    noalias(rRightHandSideVector)-=     domain_L_A_V * intrinsecPermeability_Elem * mobility_non_Elem * derivPc_Sn_Elem * prod(BiBj,saturationNonStepK_1) ;

	//Buoyancy darcyFlow RHS.
	const double buoyancy_coeff = intrinsecPermeability_Elem * (mobility_wet_Elem * densityWet_Elem + mobility_non_Elem * densityNon_Elem);
    noalias(rRightHandSideVector) +=    domain_L_A_V *buoyancy_coeff
										* prod(DN_DX, gravity);

    //*Residual RHS -= LHS*DUMMY_UNKNOWNs
    noalias(rRightHandSideVector) -=    prod(rLeftHandSideMatrix,pressureWetStepK_1);
    
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateNonWettingEquationh(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int NDof = 1; //h
    const SizeType localSize = TNumNodes * NDof;
    const double theta = rCurrentProcessInfo.GetValue(THETA); //IMP = 1.0 & EXP = 0.0. Should be always explicit
        
    //resizing as needed the LHS	
    if(rLeftHandSideMatrix.size1() != localSize)	
        rLeftHandSideMatrix.resize(localSize,localSize,false);	
    noalias(rLeftHandSideMatrix) = ZeroMatrix(localSize,localSize); 	
    
    //resizing as needed the RHS
    if(rRightHandSideVector.size() != localSize)	
        rRightHandSideVector.resize(localSize,false);	
    rRightHandSideVector = ZeroVector(localSize);
    
    //Parameters
    const unsigned int isTransient = rCurrentProcessInfo[IS_TRANSIENT];
    const unsigned int isNonPhaseCompressible = rCurrentProcessInfo[IS_COMPRESSIBLE_NON_PHASE];
    const unsigned int isWetPhaseCompressible = rCurrentProcessInfo[IS_COMPRESSIBLE_WET_PHASE];
    const unsigned int isDissolutionOn = rCurrentProcessInfo[DISSOLUTION_ON];
    const double thicknessReservoir = rCurrentProcessInfo[THICKNESS_RESERVOIR];
    const double Swr = rCurrentProcessInfo[RESIDUAL_SW];
    
    double invDt = 0.0;
    if(isTransient !=1)
    {
        const double deltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
        invDt = 1.0/ deltaTime;
    }

    VectorType gravity = ZeroVector(3);
    gravity = rCurrentProcessInfo[GRAVITY]; 
    
    //Constant
    const double porosity_Elem = this->GetProperties()[POROSITY_ELEM];
    const double  intrinsecPermeability_Elem = this->GetProperties()[PERMEABILITY_ELEM];
    
    //Non-Constant
    const double height_Elem =   this->GetValue(HEIGHT_ELEM); 
    const double permeabilityNon_Elem = this->GetValue(PERMEABILITY_NON_ELEM);
    const double permeabilityWet_Elem = this->GetValue(PERMEABILITY_WET_ELEM);
    
    //Compressible NonWetting phase case
    double densityNon_Elem = 0.0;
    double viscosityNon_Elem = 0.0;
    double compressibilityNon_Elem = 0.0;
    //double thermalExpansionNon_Elem = 0.0;
    if(isNonPhaseCompressible !=1)
    {
        densityNon_Elem = GetValue(DENSITY_NON_ELEM);
        viscosityNon_Elem = GetValue(VISCOSITY_NON_ELEM);
        if(viscosityNon_Elem == 0.0) viscosityNon_Elem = 1.0; //Case of wetting phase
        compressibilityNon_Elem = GetValue(COMPRESSIBLITY_NON_ELEM);  
        //thermalExpansionNon_Elem = GetValue(THERMAL_COEFF_NON_ELEM); 
    }
    else
    {
        densityNon_Elem = this->GetProperties()[DENSITY_NON_ELEM];
        viscosityNon_Elem = this->GetProperties()[VISCOSITY_NON_ELEM];
    }
        
    //Compressible Wetting phase case
    double densityWet_Elem = 0.0; 
    double viscosityWet_Elem = 0.0; 
    if(isWetPhaseCompressible !=1)
    {
        densityWet_Elem = this->GetValue(DENSITY_WET_ELEM);
        viscosityWet_Elem = this->GetValue(VISCOSITY_WET_ELEM);
        if(viscosityWet_Elem == 0.0) viscosityWet_Elem = 1.0; //Case of nonwetting phase 
    }
    else
    {
        if(isDissolutionOn !=1)
            densityWet_Elem = GetValue(DENSITY_WET_ELEM);
        else
            densityWet_Elem = GetProperties()[DENSITY_WET_ELEM];
            
        viscosityWet_Elem = this->GetProperties()[VISCOSITY_WET_ELEM]; 
    }
    
    //Element variable non-constant (no properties)
    const double mobility_wet_Elem = (permeabilityWet_Elem/viscosityWet_Elem);
    const double mobility_non_Elem = (permeabilityNon_Elem/viscosityNon_Elem);
    const double mobility_total_Elem = mobility_non_Elem + mobility_wet_Elem;
    const double fractional_non_Elem = mobility_non_Elem/mobility_total_Elem;
    const double fractional_wet_Elem = mobility_wet_Elem/mobility_total_Elem;
    
    VectorType darcyFlowNon = ZeroVector(TDim);
    darcyFlowNon = this->GetValue(DARCYFLOW_NON);
    VectorType darcyFlowWet = ZeroVector(TDim);
    darcyFlowWet = this->GetValue(DARCYFLOW_WET);
    VectorType darcyFlowTotal = ZeroVector(TDim);
    darcyFlowTotal = this->GetValue(DARCYFLOW_NON) + this->GetValue(DARCYFLOW_WET);
    
    //State variables: Pk+1, Swk, Swk+1 & IncrSw
    VectorType pressureZStepK_1 = ZeroVector(TNumNodes);
    VectorType pressureZStepK = ZeroVector(TNumNodes);
    VectorType pressureZIncrement = ZeroVector(TNumNodes);
    MatrixType pressureZIncrement_Matrix= ZeroMatrix(TNumNodes,TNumNodes);
    VectorType heightStepK_1 = ZeroVector(TNumNodes);
    VectorType heightStepK = ZeroVector(TNumNodes);
    VectorType sinkSourceBalance = ZeroVector(TNumNodes);
    VectorType zTop = ZeroVector(TNumNodes);
    for(unsigned int node = 0; node< TNumNodes; node++)
    {
        const double Pz_k_1 = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_Z);
        const double Pz_k = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_Z,1);
        const double IncrementPz = Pz_k_1 - Pz_k; 
        const double h_k_1 = GetGeometry()[node].FastGetSolutionStepValue(HEIGHT);
        const double h_k = GetGeometry()[node].FastGetSolutionStepValue(HEIGHT, 1);
        const double zTop_node = GetGeometry()[node].Z();
        const double skBalance = -GetGeometry()[node].FastGetSolutionStepValue(DARCY_FLOW_SUM_BALANCE);
        const double nElementsByNode = GetGeometry()[node].FastGetSolutionStepValue(N_NODES);
        
        pressureZStepK_1[node] = Pz_k_1;
        pressureZStepK[node] = Pz_k;
        pressureZIncrement[node] = IncrementPz; 
        pressureZIncrement_Matrix(node,node)= IncrementPz;
        heightStepK_1[node] = h_k_1;
        heightStepK[node] = h_k;
        zTop[node] = zTop_node;
        sinkSourceBalance[node] = skBalance/nElementsByNode;
    }

    ////Element Geometry (N,gradN)
    //area
    double domain_L_A_V = 0.0;
    //GradN
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX = ZeroMatrix(TNumNodes,TDim);
    //N
    array_1d<double, TNumNodes > N = ZeroVector(TNumNodes); 
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, domain_L_A_V); 
    
    //Elemental matrix
    //boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes > 
    MatrixType NiNj = ZeroMatrix(TNumNodes,TNumNodes);
    this->CalculateNiNj(NiNj);
    MatrixType BiBj= ZeroMatrix(TNumNodes,TNumNodes);
    this->CalculateBiBj(BiBj, DN_DX);

    ////LHS
    //term STORAGE ETA
    noalias(rLeftHandSideMatrix) +=     domain_L_A_V * porosity_Elem * ( 1.0 - Swr)* invDt * NiNj;
    
    //term STORAGE compressibility ETA (en k+1, Pz_k+1) 
    noalias(rLeftHandSideMatrix) +=     domain_L_A_V * porosity_Elem * ( 1.0 - Swr) *invDt * compressibilityNon_Elem 
                                        *prod(NiNj,pressureZIncrement_Matrix);
    
    //term ADVECTION (qn+qw)*deriv{fn}{h}*h_k+1, NO BFECC
    VectorType q_DNDX =                 ZeroVector(TNumNodes);
    noalias(q_DNDX) =                   prod(DN_DX, darcyFlowTotal);
    const double numerator =            (1.0-Swr)*viscosityWet_Elem*viscosityNon_Elem*thicknessReservoir;
    const double denominator =          pow(thicknessReservoir*viscosityNon_Elem + height_Elem*((1.0-Swr)*viscosityWet_Elem-viscosityNon_Elem),2.0);
    const double deriv_fn_h =           (numerator / denominator);
    noalias(rLeftHandSideMatrix) +=     domain_L_A_V * deriv_fn_h * outer_prod(N, q_DNDX);

    //term MOVEMENT BETWEEN PHASES (Should be zero: EXPLICIT)
    noalias(rLeftHandSideMatrix) +=     theta * domain_L_A_V * thicknessReservoir * intrinsecPermeability_Elem 
                                        * mobility_non_Elem * fractional_wet_Elem * (densityWet_Elem - densityNon_Elem) 
                                        * (-1)* gravity(1)* BiBj;
    
    ////RHS
    //term STORAGE ETA
    noalias(rRightHandSideVector) +=    domain_L_A_V * porosity_Elem * ( 1.0 - Swr)* invDt * prod(NiNj, heightStepK);
    
    //term FRACTIONAL FLOW
    noalias(rRightHandSideVector) -=    fractional_non_Elem * sinkSourceBalance;
    
    //term MOVEMENT BETWEEN PHASES (Should not be zero: EXPLICIT)
    noalias(rRightHandSideVector) -=    (1.0 - theta) * domain_L_A_V * thicknessReservoir * intrinsecPermeability_Elem 
                                        * mobility_non_Elem * fractional_wet_Elem * (densityWet_Elem - densityNon_Elem)
                                        * (-1)* gravity(1)* prod(BiBj,heightStepK);

    //term SLOPE 
    noalias(rRightHandSideVector) +=    domain_L_A_V * thicknessReservoir * intrinsecPermeability_Elem * mobility_non_Elem 
                                        * (densityNon_Elem*(fractional_non_Elem - 1.0)+densityWet_Elem* fractional_wet_Elem)
                                        * (-1)*gravity(1)* prod(BiBj,zTop);

    //term STORAGE compressibility ETA (en k, Pz_k),put it at the Pz_+1, 
//    noalias(rRightHandSideVector) -=    domain_L_A_V * invDt * porosity_Elem 
//                                        *( 1.0 - Swr) * compressibilityNon_Elem * height_Elem
//                                        * prod(NiNj, pressureZIncrement);
    
    //Residual RHS -= LHS*DUMMY_UNKNOWNs
    noalias(rRightHandSideVector) -=    prod(rLeftHandSideMatrix,heightStepK_1); 
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateNonWettingEquationSn(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{
	const unsigned int NDof = 1; //Sn
	const SizeType localSize = TNumNodes * NDof;
	const double theta = rCurrentProcessInfo.GetValue(THETA); //IMP = 1.0 & EXP = 0.0. Should be always explicit

	//resizing as needed the LHS	
	if (rLeftHandSideMatrix.size1() != localSize)
		rLeftHandSideMatrix.resize(localSize, localSize, false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(localSize, localSize); //resetting LHS	

																	 //resizing as needed the RHS
	if (rRightHandSideVector.size() != localSize)
		rRightHandSideVector.resize(localSize, false);
	rRightHandSideVector = ZeroVector(localSize); //resetting RHS

	//Parameters
	const unsigned int isTransient = rCurrentProcessInfo[IS_TRANSIENT];
	const unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
	const unsigned int isCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
	const unsigned int isWetPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_WET_PHASE];
	const unsigned int isNonPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_NON_PHASE];
	const unsigned int lambda = rCurrentProcessInfo[LAMBDA];
	const unsigned int Snr = rCurrentProcessInfo[RESIDUAL_SN];
	const unsigned int Swr = rCurrentProcessInfo[RESIDUAL_SW];
	
	double invDt = 0.0;
	if (isTransient != 1)
	{
		const double deltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
		invDt = 1.0 / deltaTime;
	}

	VectorType gravity = ZeroVector(3);
	if (isBuoyancy != 1) gravity = rCurrentProcessInfo[GRAVITY];

	//Constant
	const double porosity_Elem = GetProperties()[POROSITY_ELEM];
	const double  intrinsecPermeability_Elem = GetProperties()[PERMEABILITY_ELEM];
	const double compressibilityMedium_Elem = GetProperties()[COMPRESSIBILITY_MEDIUM_ELEM];
	const double thermalExpansionNon_Elem = GetProperties()[THERMAL_COEFF_NON_ELEM];

	//Non-Constant
	const double saturationNon_Elem = GetValue(SATURATION_NON_ELEM);
	const double saturationWet_Elem = 1.0 - saturationNon_Elem;
	const double permeabilityNon_Elem = GetValue(PERMEABILITY_NON_ELEM);
	const double permeabilityWet_Elem = GetValue(PERMEABILITY_WET_ELEM);

	//Calculate capillaryPressureModel (Pc, derivPc_Sw, dderivPc_Sw)
	double derivPc_Sn_Elem = 0.0;
	if (isCapillarityNeglected != 1)
		derivPc_Sn_Elem = GetValue(DERIV_PC_SN_ELEM);

	//Compressible NonWetting phase case
	double densityNon_Elem = 0.0;
	double viscosityNon_Elem = 0.0;
	double compressibilityNon_Elem = 0.0;
	if (isNonPhaseNonConstant != 1)
	{
		densityNon_Elem = GetValue(DENSITY_NON_ELEM);
		viscosityNon_Elem = GetValue(VISCOSITY_NON_ELEM);
		if (viscosityNon_Elem == 0.0) viscosityNon_Elem = 1.0; //Case of wetting phase
		compressibilityNon_Elem = GetValue(COMPRESSIBLITY_NON_ELEM);
	}
	else
	{
		densityNon_Elem = GetProperties()[DENSITY_NON_ELEM];
		viscosityNon_Elem = GetProperties()[VISCOSITY_NON_ELEM];
		compressibilityNon_Elem = GetProperties()[COMPRESSIBLITY_NON_ELEM];
	}

	//Compressible Wetting phase case
	double densityWet_Elem = 0.0;
	double viscosityWet_Elem = 0.0;
	if (isWetPhaseNonConstant != 1)
	{
		densityWet_Elem = GetValue(DENSITY_WET_ELEM);
		viscosityWet_Elem = GetValue(VISCOSITY_WET_ELEM);
		if (viscosityWet_Elem == 0.0) viscosityWet_Elem = 1.0; //Case of nonwetting phase
	}
	else
	{
		densityWet_Elem = GetProperties()[DENSITY_WET_ELEM];
		viscosityWet_Elem = GetProperties()[VISCOSITY_WET_ELEM];
	}

	//Element variable non-constant (no properties)
	const double mobility_wet_Elem = permeabilityWet_Elem / viscosityWet_Elem;
	const double mobility_non_Elem = permeabilityNon_Elem / viscosityNon_Elem;
	const double totalMobility_Elem = mobility_wet_Elem + mobility_non_Elem;
	const double fractionalFlow_wet_Elem = mobility_wet_Elem / totalMobility_Elem;
	const double fractionalFlow_non_Elem = mobility_non_Elem / totalMobility_Elem;

	VectorType darcyFlowNon = ZeroVector(TDim);
	darcyFlowNon = this->GetValue(DARCYFLOW_NON);
	VectorType darcyFlowWet = ZeroVector(TDim);
	darcyFlowWet = this->GetValue(DARCYFLOW_WET);
	VectorType darcyFlowTotal = ZeroVector(TDim);
	darcyFlowTotal = this->GetValue(DARCYFLOW_NON) + this->GetValue(DARCYFLOW_WET);

	//State variables: Pwk, Pwk+1, IncrPw & IncrPw_Mat, Snk, Snk+1 & IncrSn, IncrSn & IncrSn_Mat, skSourceBalancek+1
	//array_1d<double,TNumNodes>  
	VectorType pressureWetStepK_1 = ZeroVector(TNumNodes);
	VectorType pressureWetStepK = ZeroVector(TNumNodes);
	VectorType pressureWetIncrement = ZeroVector(TNumNodes);  
	MatrixType pressureWetIncrement_Matrix = ZeroMatrix(TNumNodes, TNumNodes);
	VectorType saturationNonStepK_1 = ZeroVector(TNumNodes);
	VectorType saturationNonStepK = ZeroVector(TNumNodes);
	VectorType saturationNonIncrement = ZeroVector(TNumNodes);
	MatrixType saturationNonIncrement_Matrix = ZeroMatrix(TNumNodes, TNumNodes);
	VectorType enthalpyNonIncrement = ZeroVector(TNumNodes);
	MatrixType enthalpyNonIncrement_Matrix = ZeroMatrix(TNumNodes, TNumNodes);
	VectorType sinkSourceBalance = ZeroVector(TNumNodes); 
	for (unsigned int node = 0; node< TNumNodes; node++)
	{
		const double Pw_k_1 = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET);
		const double Pw_k = GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET, 1);
		const double IncrementPw = Pw_k_1 - Pw_k;
		const double Sn_k_1 = GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON);
		const double Sn_k = GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON, 1);
		const double IncrementSn = Sn_k_1 - Sn_k;
		const double Hn_k_1 = GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE);
		const double Hn_k = GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE, 1);
		const double IncrementHn = Hn_k_1 - Hn_k;
		const double skBalance = -GetGeometry()[node].FastGetSolutionStepValue(DARCY_FLOW_SUM_BALANCE); 
		const double nElementsByNode = GetGeometry()[node].FastGetSolutionStepValue(N_NODES); //1D need comented and put to 2

		pressureWetStepK_1[node] = Pw_k_1;
		pressureWetStepK[node] = Pw_k;
		pressureWetIncrement[node] = IncrementPw;
		pressureWetIncrement_Matrix(node, node) = IncrementPw;
		saturationNonStepK_1[node] = Sn_k_1;
		saturationNonStepK[node] = Sn_k;
		saturationNonIncrement[node] = IncrementSn;
		saturationNonIncrement_Matrix(node, node) = IncrementSn;
		enthalpyNonIncrement[node] = IncrementHn;
		enthalpyNonIncrement_Matrix(node, node) = IncrementHn;
		sinkSourceBalance[node] = skBalance / nElementsByNode;
	}

	////Element Geometry (N,gradN)
	//area
	double domain_L_A_V = 0.0;
	//GradN
	boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX = ZeroMatrix(TNumNodes, TDim);
	//N
	array_1d<double, TNumNodes > N = ZeroVector(TNumNodes);
	GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, domain_L_A_V);

	//Case of 1D problems
	//if (TNumNodes == 2)
	//	domain_L_A_V = 1.0;

	//Elemental matrix
	//boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes > 
	MatrixType NiNj = ZeroMatrix(TNumNodes, TNumNodes);
	this->CalculateNiNj(NiNj);
	MatrixType BiBj = ZeroMatrix(TNumNodes, TNumNodes);
	this->CalculateBiBj(BiBj, DN_DX);

	////LHS
	//Compressibility storage Sn ONLY LHS. (en k+1, Sn_k+1) 
	const double nonWetting_comp_coeff = ( (porosity_Elem*compressibilityNon_Elem) - compressibilityMedium_Elem );
	noalias(rLeftHandSideMatrix) +=		domain_L_A_V * invDt * nonWetting_comp_coeff * prod(NiNj, pressureWetIncrement_Matrix);

	//Capillarity storage LHS.
	const double cap_storage_coeff =	nonWetting_comp_coeff * derivPc_Sn_Elem * saturationNon_Elem;
	noalias(rLeftHandSideMatrix) +=		domain_L_A_V * invDt * cap_storage_coeff * NiNj;
	//Way 2. Rwmove: Capillarity storage RHS.
	//noalias(rLeftHandSideMatrix) +=     domain_L_A_V * invDt * cap_storage_coeff * prod(NiNj, saturationNonIncrement_Matrix);

	//Storage Sn LHS.
	noalias(rLeftHandSideMatrix) +=		domain_L_A_V * invDt * porosity_Elem * NiNj;
	
	//Non-wetting thermal expansion ONLY LHS. (en k+1, Sn_k+1) 
	noalias(rLeftHandSideMatrix) +=		domain_L_A_V * invDt 
										* porosity_Elem * thermalExpansionNon_Elem
										* prod(NiNj, enthalpyNonIncrement_Matrix);

	//*term ADVECTION (qn+qw)*deriv{fn}{Sn}*Sn_k+1, LHS. (Should be zero: EXPLICIT) DOUBT!!! theta added!!
	VectorType q_DNDX = ZeroVector(TNumNodes);
	noalias(q_DNDX) = prod(DN_DX, darcyFlowTotal);
	const double numerator1 = lambda*viscosityWet_Elem*pow((saturationNon_Elem - Snr), (double)(lambda - 1.0));
	const double numerator2 = viscosityWet_Elem*pow((saturationNon_Elem - Snr), (double)(lambda))*(lambda*viscosityWet_Elem*pow((saturationNon_Elem - Snr), (double)(lambda - 1.0))-lambda*viscosityNon_Elem*pow((1.0 - saturationNon_Elem - Swr), (double)(lambda - 1.0)));
	const double denominator1 = (pow((saturationNon_Elem - Snr), (double)(lambda))*viscosityWet_Elem)+(pow((1.0 - saturationNon_Elem - Swr), (double)(lambda))*viscosityNon_Elem);
	const double denominator2 = pow(((pow((saturationNon_Elem - Snr), (double)(lambda))*viscosityWet_Elem) + (pow((1.0 - saturationNon_Elem - Swr), (double)(lambda))*viscosityNon_Elem)),2.0);
	const double deriv_fn_Sn = (numerator1 / denominator1) - (numerator2 / denominator2);
	noalias(rLeftHandSideMatrix) += theta * domain_L_A_V * deriv_fn_Sn * outer_prod(N, q_DNDX);

	//*Capillarity darcyFlow LHS. (Should be zero: EXPLICIT)
	noalias(rLeftHandSideMatrix) +=		theta * domain_L_A_V * intrinsecPermeability_Elem
										* mobility_non_Elem * (1.0 - fractionalFlow_non_Elem) * derivPc_Sn_Elem
										* BiBj;

	//RHS
	//Capillarity storage RHS.
	noalias(rRightHandSideVector) +=	domain_L_A_V * invDt * cap_storage_coeff  * prod(NiNj, saturationNonStepK);

	//*Storage Sn RHS.
	noalias(rRightHandSideVector) +=	domain_L_A_V * invDt * porosity_Elem  * prod(NiNj, saturationNonStepK);


	//*term ADVECTION (qn+qw)*deriv{fn}{Sn}*Sn_k+1, RHS. (Should be not zero: EXPLICIT) DOUBT!!!
	noalias(rRightHandSideVector) -= (1.0 - theta) * domain_L_A_V * deriv_fn_Sn * prod(outer_prod(N, q_DNDX), saturationNonStepK);

	//Capillarity darcyFlow RHS. (Should be not zero: EXPLICIT)
	noalias(rRightHandSideVector) -=	(1.0 - theta) * domain_L_A_V * intrinsecPermeability_Elem
										* mobility_non_Elem * (1.0 - fractionalFlow_non_Elem) * derivPc_Sn_Elem
										* prod(BiBj, saturationNonStepK);

	//*term FRACTIONAL FLOW
	noalias(rRightHandSideVector) -= fractionalFlow_non_Elem * sinkSourceBalance;

	//Buoyancy darcyFlow RHS.
	const double buoyancy_coeff = intrinsecPermeability_Elem * mobility_non_Elem *(((1.0 - fractionalFlow_non_Elem)*densityNon_Elem)  - (fractionalFlow_wet_Elem*densityWet_Elem));
	noalias(rRightHandSideVector) +=	domain_L_A_V *buoyancy_coeff
										* prod(DN_DX, gravity);

	//Compressibility storage Sn RHS. (en k, Pw_k) ,put it at the Pw_+1, 
	//    noalias(rRightHandSideVector) -=    domain_L_A_V * invDt
	//									      *saturationNon_Elem * nonWetting_comp_coeff * prod(NiNj, pressureWetIncrement);

	//Residual RHS -= LHS*DUMMY_UNKNOWNs
	noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, saturationNonStepK_1);
    
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateNonWettingEnergyEquationHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NDof = 1; //Hn
    const SizeType localSize = TNumNodes * NDof;
        
    //resizing as needed the LHS	
    if(rLeftHandSideMatrix.size1() != localSize)	
        rLeftHandSideMatrix.resize(localSize,localSize,false);	
    noalias(rLeftHandSideMatrix) = ZeroMatrix(localSize,localSize); //resetting LHS	
    
    //resizing as needed the RHS
    if(rRightHandSideVector.size() != localSize)	
        rRightHandSideVector.resize(localSize,false);	
    rRightHandSideVector = ZeroVector(localSize); //resetting RHS
    
    //Parameters
    const unsigned int isTransient = rCurrentProcessInfo[IS_TRANSIENT];
    const unsigned int isNonPhaseCompressible = rCurrentProcessInfo[IS_COMPRESSIBLE_NON_PHASE];
    double invDt = 0.0;
    if(isTransient !=1)
    {
        const double deltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
        invDt = 1.0/ deltaTime;
    }
    
    //Constant
    const double porosity_Elem = GetProperties()[POROSITY_ELEM];
    const double thermalConductivityNon_Elem = GetProperties()[THERMALCONDUCTIVITY_NON_ELEM];
    
    //Non-Constant
    const double saturationWet_Elem = GetValue(SATURATION_WET_ELEM);
    const double saturationNon_Elem = 1.0 - saturationWet_Elem;

    //Compressible NonWetting phase case
    double densityNon_Elem = 0.0;
    double specificHeatNon_Elem = 0.0;
    if(isNonPhaseCompressible !=1)
    {        
        densityNon_Elem = GetValue(DENSITY_NON_ELEM);
        specificHeatNon_Elem = GetValue(SPECIFIC_HEAT_NON_ELEM);
        if(densityNon_Elem == 0.0) densityNon_Elem = 1.0; //Case of wetting phase
        if(specificHeatNon_Elem == 0.0) specificHeatNon_Elem = 1.0; //Case of wetting phase
    }
    else
    {
        densityNon_Elem = GetProperties()[DENSITY_NON_ELEM];
        specificHeatNon_Elem = GetProperties()[SPECIFIC_HEAT_NON_ELEM];
    }
    
    VectorType darcyFlowNon = ZeroVector(TDim);
    darcyFlowNon = this->GetValue(DARCYFLOW_NON);
    
    //State variables: Hnk, Hnk+1 & IncrHn
    VectorType enthalpyNonStepK_1 = ZeroVector(TNumNodes);
    VectorType enthalpyNonStepK = ZeroVector(TNumNodes);
    VectorType temperatureNonStepK_1 = ZeroVector(TNumNodes);
    for(unsigned int node = 0; node< TNumNodes; node++)
    {
        const double hn_k_1 = GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE);
        const double hn_k = GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE,1);
        const double T_k_1 = GetGeometry()[node].FastGetSolutionStepValue(TEMPERATURE_NODE);
        enthalpyNonStepK_1[node] = hn_k_1;
        enthalpyNonStepK[node] = hn_k;
        temperatureNonStepK_1[node] = T_k_1;
    }

    ////Element Geometry (N,gradN)
    //area
    double domain_L_A_V = 0.0;
    //GradN
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX = ZeroMatrix(TNumNodes,TDim);
    //N
    array_1d<double, TNumNodes > N = ZeroVector(TNumNodes); 
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, domain_L_A_V); 
    
	//Case of 1D problems
	if (TNumNodes == 2)
		domain_L_A_V = 1.0;

    //Elemental matrix
    MatrixType NiNj = ZeroMatrix(TNumNodes,TNumNodes);
    this->CalculateNiNj(NiNj);
    MatrixType BiBj= ZeroMatrix(TNumNodes,TNumNodes);
    this->CalculateBiBj(BiBj, DN_DX);

    ////LHS
    //Storage
    noalias(rLeftHandSideMatrix) +=     domain_L_A_V * invDt * porosity_Elem 
                                        *saturationNon_Elem*NiNj;
    
    //Advective contribution to the stiffness matrix 
    array_1d<double, TDim > q_DN = ZeroVector(TDim);
    noalias(q_DN) = prod(DN_DX, darcyFlowNon);
    noalias(rLeftHandSideMatrix) +=     domain_L_A_V*outer_prod(N, q_DN);
    
    //RHS
    //Diffusion (thermal conduction)
    //noalias(rRightHandSideVector) +=    domain_L_A_V*(thermalConductivityNon_Elem/densityNon_Elem)
    //                                    * prod( BiBj, temperatureNonStepK_1);
    
    //Diffusion (thermal conduction)
    noalias(rRightHandSideVector) +=    domain_L_A_V * saturationNon_Elem
                                        * (thermalConductivityNon_Elem/(densityNon_Elem*specificHeatNon_Elem))
                                        * prod( BiBj, enthalpyNonStepK_1);
    
    //Storage
    noalias(rRightHandSideVector) +=    domain_L_A_V * invDt * porosity_Elem 
                                        *saturationNon_Elem* prod(NiNj, enthalpyNonStepK);
                
    
    //Residual RHS -= LHS*DUMMY_UNKNOWNs
    noalias(rRightHandSideVector) -=    prod(rLeftHandSideMatrix,enthalpyNonStepK_1);
    
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculateWettingEnergyEquationT(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::CalculatePhasesDarcyFlow(ProcessInfo& rCurrentProcessInfo)
{
	//Parameters
	const unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
	const unsigned int isCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
	const unsigned int isWetPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_WET_PHASE];
	const unsigned int isNonPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_NON_PHASE];
	const unsigned int Snr = rCurrentProcessInfo[RESIDUAL_SN];
	const unsigned int Swr = rCurrentProcessInfo[RESIDUAL_SW];

	VectorType gravity = ZeroVector(3);
	if (isBuoyancy != 1) gravity = rCurrentProcessInfo[GRAVITY];

	//Constant
	const double porosity_Elem = GetProperties()[POROSITY_ELEM];
	const double  intrinsecPermeability_Elem = GetProperties()[PERMEABILITY_ELEM];
	const double compressibilityMedium_Elem = GetProperties()[COMPRESSIBILITY_MEDIUM_ELEM];
	const double thermalExpansionNon_Elem = GetProperties()[THERMAL_COEFF_NON_ELEM];

	//Non-Constant
	const double saturationNon_Elem = GetValue(SATURATION_NON_ELEM);
	const double saturationWet_Elem = 1.0 - saturationNon_Elem;
	const double Se_Elem = (saturationWet_Elem - Swr) / (1.0 - Snr -Swr);
	const double permeabilityNon_Elem = GetValue(PERMEABILITY_NON_ELEM);
	const double permeabilityWet_Elem = GetValue(PERMEABILITY_WET_ELEM);

	//Calculate capillaryPressureModel (Pc, derivPc_Sw, dderivPc_Sw)
	double derivPc_Sn_Elem = 0.0;
	if (isCapillarityNeglected != 1)
		derivPc_Sn_Elem = GetValue(DERIV_PC_SN_ELEM);

	//Compressible NonWetting phase case
	double densityNon_Elem = 0.0;
	double viscosityNon_Elem = 0.0;
	double compressibilityNon_Elem = 0.0;
	if (isNonPhaseNonConstant != 1)
	{
		densityNon_Elem = GetValue(DENSITY_NON_ELEM);
		viscosityNon_Elem = GetValue(VISCOSITY_NON_ELEM);
		if (viscosityNon_Elem == 0.0) viscosityNon_Elem = 1.0; //Case of wetting phase
		compressibilityNon_Elem = GetValue(COMPRESSIBLITY_NON_ELEM);
	}
	else
	{
		densityNon_Elem = GetProperties()[DENSITY_NON_ELEM];
		viscosityNon_Elem = GetProperties()[VISCOSITY_NON_ELEM];
		compressibilityNon_Elem = GetProperties()[COMPRESSIBLITY_NON_ELEM];
	}

	//Compressible Wetting phase case
	double densityWet_Elem = 0.0;
	double viscosityWet_Elem = 0.0;
	if (isWetPhaseNonConstant != 1)
	{
		densityWet_Elem = GetValue(DENSITY_WET_ELEM);
		viscosityWet_Elem = GetValue(VISCOSITY_WET_ELEM);
		if (viscosityWet_Elem == 0.0) viscosityWet_Elem = 1.0; //Case of nonwetting phase
	}
	else
	{
		densityWet_Elem = GetProperties()[DENSITY_WET_ELEM];
		viscosityWet_Elem = GetProperties()[VISCOSITY_WET_ELEM];
	}

	//Element variable non-constant (no properties)
	const double mobility_wet_Elem = permeabilityWet_Elem / viscosityWet_Elem;
	const double mobility_non_Elem = permeabilityNon_Elem / viscosityNon_Elem;
	const double totalMobility_Elem = mobility_wet_Elem + mobility_non_Elem;
	const double fractionalFlow_wet_Elem = mobility_wet_Elem / totalMobility_Elem;
	const double fractionalFlow_non_Elem = mobility_non_Elem / totalMobility_Elem;

    ////Element Geometry (N,gradN)
    //area
    double domain_L_A_V = 0.0;
    //GradN
    //Matrix 
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > DN_DX = ZeroMatrix(TNumNodes,TDim);
    //N
    //Vector
    array_1d<double, TNumNodes > N = ZeroVector(TNumNodes); 
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, domain_L_A_V); 
    
	//Case of 1D problems
	if (TNumNodes == 2)
		domain_L_A_V = 1.0;

    //Compute gradient of head
	VectorType gradPw = ZeroVector(TDim);
	VectorType gradPw_NonBuoy = ZeroVector(TDim);
	VectorType gradPw_WetBuoy = ZeroVector(TDim);
	VectorType saturationNonGrad = ZeroVector(TDim);

    for (unsigned int node = 0; node < TNumNodes ; node++)
    {
       const double currentPressureWet_Node =  this->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET);
       const double currentSaturationNon_Node =  this->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON);
       for (unsigned int dim = 0; dim < TDim; dim++)
        {
           gradPw[dim]+=currentPressureWet_Node*DN_DX.at_element(node,dim);
           saturationNonGrad[dim]+=currentSaturationNon_Node*DN_DX.at_element(node,dim);
        }

    }

    gradPw_NonBuoy = gradPw;
    gradPw_WetBuoy = gradPw;

    for (unsigned int dim = 0; dim < TDim; dim++)
    {
       gradPw_NonBuoy[dim] -= gravity[dim] * (densityNon_Elem) ;
       gradPw_WetBuoy[dim] -= gravity[dim] * (densityWet_Elem) ;
    }

    //K*gradh, isotropic
    Matrix Kintrinsec_Matrix = ZeroMatrix(TDim,TDim);
    for (unsigned int dim = 0; dim < TDim; dim++)
    {
        Kintrinsec_Matrix(dim,dim)= intrinsecPermeability_Elem;
    }

	VectorType nonDarcyFlow_Elem = ZeroVector(TDim);
	VectorType wetDarcyFlow_Elem = ZeroVector(TDim);

    if( Se_Elem >= 1.0 /*0.9999*/ )
    {
		wetDarcyFlow_Elem = -mobility_wet_Elem*prod(Kintrinsec_Matrix, trans(gradPw_WetBuoy));

    }
    else if( Se_Elem <= 0.0 )
    {
        nonDarcyFlow_Elem =     -mobility_non_Elem*prod(Kintrinsec_Matrix,trans(gradPw_NonBuoy));
								-mobility_non_Elem*derivPc_Sn_Elem*prod(Kintrinsec_Matrix, trans(saturationNonGrad));
    }
    else
    {
		nonDarcyFlow_Elem = -mobility_non_Elem*prod(Kintrinsec_Matrix, trans(gradPw_NonBuoy));
							-mobility_non_Elem*derivPc_Sn_Elem*prod(Kintrinsec_Matrix, trans(saturationNonGrad));
		wetDarcyFlow_Elem = -mobility_wet_Elem*prod(Kintrinsec_Matrix, trans(gradPw_WetBuoy));
    }

    this->SetValue(DARCYFLOW_NON,nonDarcyFlow_Elem);
    this->SetValue(DARCYFLOW_WET,wetDarcyFlow_Elem);
            
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetPwValues(Vector& rValues, const int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i)
        rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_WET,Step);
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetSnValues(Vector& rValues,
                                          const int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rValues[i] = rGeom[i].FastGetSolutionStepValue(SATURATION_NON,Step);
    }
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetHnValues(Vector& rValues,
                                          const int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rValues[i] = rGeom[i].FastGetSolutionStepValue(ENTHALPY_NON_NODE,Step);
    }
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
void MultiphaseFEM<TDim,TNumNodes>::GetTValues(Vector& rValues,
                                          const int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rValues[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE_NODE,Step);
    }
}

template< unsigned int TDim,
          unsigned int TNumNodes  >
double MultiphaseFEM<TDim,TNumNodes>::ElementSize(/*ShapeFunctionDerivativesType &rDN_DX*/)
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    // calculate minimum element length (used in stabilization Tau)
    array_1d<double,3> Edge(3,0.0);
    Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    double ElemSize = Edge[0]*Edge[0];
    for (SizeType d = 1; d < TDim; d++)
        ElemSize += Edge[d]*Edge[d];

    for (SizeType i = 2; i < NumNodes; i++)
        for(SizeType j = 0; j < i; j++)
        {
            Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
            double Length = Edge[0]*Edge[0];
            for (SizeType d = 1; d < TDim; d++)
                Length += Edge[d]*Edge[d];
            if (Length < ElemSize) ElemSize = Length;
        }
    return sqrt(ElemSize);
}

/*
 * Template class definition (this should allow us to compile the desired template instantiations)
 */

template class MultiphaseFEM<1>;
template class MultiphaseFEM<2>;

}
