/* *********************************************************
 *
 *   Last Modified by:    $Author: VictorBez $
 *   Date:                $Date: 2015        $
 *   Revision:            $Revision: 1.0     $
 *
 * ***********************************************************/


#if !defined(KRATOS_SEQUENTIAL_MULTIPHASE_FEM_STRATEGY)
#define  KRATOS_SEQUENTIAL_MULTIPHASE_FEM_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"

//Inheritance
#include "solving_strategies/strategies/solving_strategy.h"

//Configuration files of Strategy (to load objects of configuration attributes to deal with Matrix)
#include "custom_strategies/solver_strategy_configuration.h"
#include "custom_strategies/sequential_multiphase_fem_configuration.h"


//Application
#include "multiphase_application.h"
#include "porous_media_application_variables.h"
#include "custom_utilities/spanwagnereos.h"

namespace Kratos
{
    
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
	class SequentialMultiphaseFEMStrategy
	: public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

	/** Counted pointer of ClassName */
	KRATOS_CLASS_POINTER_DEFINITION(SequentialMultiphaseFEMStrategy);

	typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

	typedef typename BaseType::TDataType TDataType;

	typedef typename BaseType::DofsArrayType DofsArrayType;

	typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

	typedef typename BaseType::TSystemVectorType TSystemVectorType;

	typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

	typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;


	/*
	 * From incompresible_fluid_application/ strategies/ custom_strategies/ fractional_step_streategy.h
	 */

	 //Doc-> constructor
	SequentialMultiphaseFEMStrategy(
		ModelPart& model_part,
		SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& aSolverStrategyConfiguration,
		bool aReformDofAtEachIteration = true,
		unsigned int aTimeOrder = 1,
		unsigned int aDomainSize = 1,
		bool aPredictorCorrector = false,
		bool aMoveMeshFlag = false,
		unsigned int aEchoLevel = 3)
		: SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, aMoveMeshFlag), mSolverStrategyConfiguration(aSolverStrategyConfiguration)
	{
		KRATOS_TRY

		this->mPredictorOrder = aTimeOrder;
		this->mTimeOrder = aTimeOrder;
		this->mDomainSize = aDomainSize;
		this->mPredictorCorrector = aPredictorCorrector;
		this->mReformDofAtEachIteration = aReformDofAtEachIteration;
		this->mEchoLevel = aEchoLevel;

		//performs checks to verify the quality of the input
		this->Check();

		ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

		this->mNoIterative = rCurrentProcessInfo[NO_ITERATIVE];

		//Tolerance & max iter of last two system moved herein
		//due to troubles with python & boost in argument constructor in Kratos
		this->mTolerance_Pw = rCurrentProcessInfo[MAXDIFF_PW];
		this->mMaxIter_Pw = rCurrentProcessInfo[MAXITER_PW];
		this->mTolerance_Sn = rCurrentProcessInfo[MAXDIFF_SN];
		this->mMaxIter_Sn = rCurrentProcessInfo[MAXITER_SN];
		this->mTolerance_Hn = rCurrentProcessInfo[MAXDIFF_HN];
		this->mMaxIter_Hn = rCurrentProcessInfo[MAXITER_HN];
		this->mTolerance_T = rCurrentProcessInfo[MAXDIFF_T];
		this->mMaxIter_T = rCurrentProcessInfo[MAXITER_T];

		this->mRelativePermeabilityModel = rCurrentProcessInfo[RELATIVE_PERMEABILITY_MODEL];

		this->mIsCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
		this->mCapillarityPressureModel = rCurrentProcessInfo[CAPILLARITY_PRESSURE_MODEL];

		this->mIsNonPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_NON_PHASE];
		this->mIsWetPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_WET_PHASE];

		this->mIsNonIsoThermal = rCurrentProcessInfo[IS_NONISOTHERMAL_NON];
		this->mIsWetIsoThermal = rCurrentProcessInfo[IS_NONISOTHERMAL_WET];

		mSpanWagnerEOS = new SpanWagnerEOS();
		this->mEOS = rCurrentProcessInfo[CO2_EOS];
		if (this->mEOS == "SpanWagner")
			this->mSpanWagnerEOS->FillTables();

		rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
		this->mSumPhasesEqStrategy_Pw = mSolverStrategyConfiguration.pGetStrategy(std::string("SumPhasesEqStrategy_Pw"));
		rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
		this->mNonPhaseEqStrategy_Sn = mSolverStrategyConfiguration.pGetStrategy(std::string("NonPhaseEqStrategy_Sn"));
		if (this->mIsNonIsoThermal != 1)
		{
			rCurrentProcessInfo[FRACTIONAL_STEP] = 3;
			this->mNonWetPhaseEnergyEqStrategy_Hn = mSolverStrategyConfiguration.pGetStrategy(std::string("NonWetPhaseEnergyEqStrategy_Hn"));
			if (this->mIsWetIsoThermal != 1)
			{
				rCurrentProcessInfo[FRACTIONAL_STEP] = 4;
				this->mWetPhaseEnergyEqStrategy_T = mSolverStrategyConfiguration.pGetStrategy(std::string("WetPhaseEnergyEqStrategy_T"));
			}
		}

		//Set to initial values state variable buffers 
		this->inicializateSV(PRESSURE_WET, PRESSURE_WET_OLD_ITER);
		this->inicializateSV(SATURATION_NON, SATURATION_NON_OLD_ITER);
		this->inicializateSV(ENTHALPY_NON_NODE, ENTHALPY_NON_OLD_ITER);
		this->inicializateSV(TEMPERATURE_NODE, TEMPERATURE_OLD_ITER);

		//OJOOOOOOOOO... 1D
		this->CalculateNodalArea();

		//Inicializate balance totals to storage & print
		const int numberOfVariables = 8;
		this->mBalanceTotalsToPrint = ZeroVector(numberOfVariables);

		this->mStep = 1;

		KRATOS_CATCH("")
	}

	virtual ~SequentialMultiphaseFEMStrategy()
	{
		delete this->mSpanWagnerEOS;
	}

	void CalculateNodalArea()
	{
		const unsigned int NDim = mSolverStrategyConfiguration.GetDomainSize();

		for (int node = 0; node < static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
		{
			ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;
			itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
			itNode->FastGetSolutionStepValue(N_NODES) = 0.0;
		}

		for (int elem = 0; elem < static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
		{
			ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin() + elem;

			unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
			const double weightElem = 1.0 / numberOfNodes;

			////Element Geometry (N,gradN)
			//area
			double domain_L_A_V = 0.0;
			//GradN
			Matrix DN_DX = ZeroMatrix(numberOfNodes, NDim);
			//N
			Vector N = ZeroVector(numberOfNodes);
			GeometryUtils::CalculateGeometryData(itElem->GetGeometry(), DN_DX, N, domain_L_A_V);

			//Case of 1D problems
			if (numberOfNodes == 2)
				domain_L_A_V = 1.0;

			for (int inode = 0; inode < numberOfNodes; inode++)
			{
				itElem->GetGeometry()[inode].FastGetSolutionStepValue(NODAL_AREA) += domain_L_A_V*weightElem;
				itElem->GetGeometry()[inode].FastGetSolutionStepValue(N_NODES) += 1.0;
			}

			//Calculate proper_side of element neccesary to stability condition (Only triangles!!!))
			// X
			double currentSideMeasure_X = 0.0;
			const double x0 = itElem->GetGeometry()[0].X();
			const double x1 = itElem->GetGeometry()[1].X();
			const double x2 = itElem->GetGeometry()[2].X();
			const double y0 = itElem->GetGeometry()[0].Y();
			const double y1 = itElem->GetGeometry()[1].Y();
			const double y2 = itElem->GetGeometry()[2].Y();
			const double z2 = itElem->GetGeometry()[2].Z();

			double sideMeasure_X_1 = abs(x0 - x1);
			double sideMeasure_X_2 = abs(x1 - x2);
			double sideMeasure_X_3 = abs(x2 - x0);
			if (sideMeasure_X_1 == 0.0) sideMeasure_X_1 = 1e+9;
			if (sideMeasure_X_2 == 0.0) sideMeasure_X_2 = 1e+9;
			if (sideMeasure_X_3 == 0.0) sideMeasure_X_3 = 1e+9;
			currentSideMeasure_X = min(sideMeasure_X_1, sideMeasure_X_2);
			if (currentSideMeasure_X > min(sideMeasure_X_2, sideMeasure_X_3))
				currentSideMeasure_X = min(sideMeasure_X_2, sideMeasure_X_3);

			const double minSideElem_X = currentSideMeasure_X;
			itElem->SetValue(PROPER_SIDE_X, minSideElem_X);

			// Y
			double currentSideMeasure_Y = 0.0;
			double sideMeasure_Y_1 = abs(y0 - y1);
			double sideMeasure_Y_2 = abs(y1 - y2);
			double sideMeasure_Y_3 = abs(y2 - y0);
			if (sideMeasure_Y_1 == 0.0) sideMeasure_Y_1 = 1e+9;
			if (sideMeasure_Y_2 == 0.0) sideMeasure_Y_2 = 1e+9;
			if (sideMeasure_Y_3 == 0.0) sideMeasure_Y_3 = 1e+9;
			currentSideMeasure_Y = min(sideMeasure_Y_1, sideMeasure_Y_2);
			if (currentSideMeasure_Y > min(sideMeasure_Y_2, sideMeasure_Y_3))
				currentSideMeasure_Y = min(sideMeasure_Y_2, sideMeasure_Y_3);

			const double minSideElem_Y = currentSideMeasure_Y;
			itElem->SetValue(PROPER_SIDE_Y, minSideElem_Y);
		}
	}

	void Initialize()
	{
		//Calculate enthalpy from Pressure & Temperature... 1D comented
		this->calculateInitialConditions();
	}

	void InitializeTimeStep()
	{
		//Only get in if its nonIsothermal (enthalpy transport not working)... 1D comented
		if (this->mIsNonIsoThermal == 1)
			this->calculateEnthalpyNode();
	}

	double Solve()
	{
		KRATOS_TRY

		Timer::Start("solve");

		//Initial Guess, if not is by default method
		//Predict();
		//PredictSV(this->mStep, this->mPredictorOrder);

		//Calculate Non-wetting enthalpy: Set Non-Wetting enthalpy: hn = (P,T)
		if (this->mStep == 1)
			this->Initialize();
		else
			this->InitializeTimeStep();

		//Initial Guess: Set to initial iteration value the previous time SV value
		this->SetPreviousTimeToOldIterationSV(SATURATION_NON_OLD_ITER, SATURATION_NON);
		this->SetPreviousTimeToOldIterationSV(PRESSURE_WET_OLD_ITER, PRESSURE_WET);
		this->SetPreviousTimeToOldIterationSV(ENTHALPY_NON_OLD_ITER, ENTHALPY_NON_NODE);
		this->SetPreviousTimeToOldIterationSV(TEMPERATURE_OLD_ITER, TEMPERATURE_NODE);

		//Iteration loop
		unsigned int globalIteration = 0;
		bool isGlobalConverged = false;
		while (!isGlobalConverged)
		{
			std::ofstream myfile;
			myfile.open("ConvergedManager.txt", std::ios::app);
			myfile << "Global iteration = " << globalIteration << std::endl;

			this->calculateLocalEquations();

			//this->writeResults("one", this->mStep);

			std::cout << "Solving... Sum Phases eq. (Pw system): " << std::endl;
			if (globalIteration != 0) this->StorageOldIteration(PRESSURE_WET, PRESSURE_WET_OLD_ITER);

			const double norm_Pw = systemIteration(1);

			double maxDiff_Pw = 0.0;
			const bool isConverged_Pw = isLocalConverged("Pressure Wet", PRESSURE_WET, PRESSURE_WET_OLD_ITER, this->mTolerance_Pw, maxDiff_Pw);
			myfile << "Pw maxDiff =  " << maxDiff_Pw << std::endl;
			std::cout << "Pw maxDiff = " << maxDiff_Pw << std::endl;

			//this->writeResults("two", this->mStep);

			//ESTABA EN CLASICO
			//this->calculateLocalEquations();

			//ESTABA EN VE, EN CLÁSICO EN EL ELEMENTO
			this->calculatePhasesDarcyFlow();

			//Calculate balances (massFlow)
			this->CalculateSumBalances();

			std::cout << "Solving... Non Phase eq. (Sn system): " << std::endl;
			if (globalIteration != 0) this->StorageOldIteration(SATURATION_NON, SATURATION_NON_OLD_ITER);
			const double norm_Sn = systemIteration(2);
			double maxDiff_Sn = 0.0;
			const bool isConverged_Sn = isLocalConverged("Saturation Non", SATURATION_NON, SATURATION_NON_OLD_ITER, this->mTolerance_Sn, maxDiff_Sn);
			myfile << "Sn maxDiff = " << maxDiff_Sn << std::endl;
			std::cout << "Sn maxDiff = " << maxDiff_Sn << std::endl;

			//this->writeResults("three", this->mStep);

			//Calculate balances (massFlow)
			//this->CalculateNonWettingBalances();
			//this->CalculateWettingBalances();

			bool isConverged_Hn = true;
			bool isConverged_T = true;
			if (this->mIsNonIsoThermal != 1)
			{
				std::cout << "Solving... Non Wetting Phase Energy eq. (Hn system): " << std::endl;
				if (globalIteration != 0) this->StorageOldIteration(ENTHALPY_NON_NODE, ENTHALPY_NON_OLD_ITER);
				const double norm_Hn = systemIteration(3);
				double maxDiff_Hn = 0.0;
				isConverged_Hn = isLocalConverged("Enthalpy", ENTHALPY_NON_NODE, ENTHALPY_NON_OLD_ITER, this->mTolerance_Hn, maxDiff_Hn);
				myfile << "Hn maxDiff = " << maxDiff_Hn << std::endl;
				std::cout << "Hn maxDiff = " << maxDiff_Hn << std::endl;

				//this->CalculateNonWettingEnthalpyBalances();

				if (this->mIsWetIsoThermal != 1)
				{
					std::cout << "Solving... Wetting Phase energy eq. (T system): " << std::endl;
					if (globalIteration != 0) this->StorageOldIteration(TEMPERATURE_NODE, TEMPERATURE_OLD_ITER);
					const double norm_T = systemIteration(4);
					double maxDiff_T = 0.0;
					isConverged_T = isLocalConverged("Temperature", TEMPERATURE_NODE, TEMPERATURE_OLD_ITER, this->mTolerance_T, maxDiff_T);
					myfile << "T maxDiff = " << maxDiff_T << std::endl;
					std::cout << "T maxDiff = " << maxDiff_T << std::endl;

					//this->CalculateWettingTemperatureBalances();
				}
			}

			globalIteration++;

			//Only iterative first time & first iteration
			bool firstIteration_t_0 = true;
			if (this->mStep == 1)
				if (globalIteration == 1)
					firstIteration_t_0 = false;
				else
					firstIteration_t_0 = true;

			//If not iterative (IMPES... Always), Incompressible or converged:
			if (this->mNoIterative != 1 && firstIteration_t_0 ||
				this->mIsNonPhaseNonConstant == 1 && this->mIsWetPhaseNonConstant == 1 ||
				isConverged_Pw && isConverged_Sn && isConverged_Hn && isConverged_T)
			{
				isGlobalConverged = true;

				myfile << "GLOBAL convergence ACHIEVED. Time: " << this->mStep << " .Iteration: " << globalIteration << std::endl;
				std::cout << "GLOBAL convergence ACHIEVED. Time: " << this->mStep << " .Iteration: " << globalIteration << std::endl;
				//this->writeResults("Converged", this->mStep);

				//Calculate All balances
				this->CalculateSumBalances();
				//this->CalculateNonWettingBalances();
				//this->CalculateWettingBalances();
				//if (this->mIsNonIsoThermal != 1)
				//{
					//this->CalculateNonWettingEnthalpyBalances();
					//if (this->mIsNonWetIsoThermal != 1)
						//this->CalculateWettingTemperatureBalances();
				//}

				//Parameters & properties
				this->calculateLocalEquations();

				//write results in a custom output
				//this->writeBalances("Converged", this->mStep, globalIteration);
				this->writeResults("Converged", this->mStep);
			}
			else
			{
				myfile << "GLOBAL convergence NOT achieved. " << "Iteration: " << globalIteration << std::endl;
				std::cout << "GLOBAL convergence NOT achieved. " << "Iteration: " << globalIteration << std::endl;

				bool isMaxNumberOfIterationsReached = false;
				isMaxNumberOfIterationsReached = this->checkNumberIterations(globalIteration, this->mMaxIter_Pw);
				isMaxNumberOfIterationsReached = this->checkNumberIterations(globalIteration, this->mMaxIter_Sn);
				if (this->mIsNonIsoThermal != 1)
				{
					isMaxNumberOfIterationsReached = this->checkNumberIterations(globalIteration, this->mMaxIter_Hn);
					if (this->mIsWetIsoThermal != 1)
						isMaxNumberOfIterationsReached = this->checkNumberIterations(globalIteration, this->mMaxIter_T);
				}

				if (isMaxNumberOfIterationsReached) exit(1);
			}

			myfile.close();

		}

		if (this->mReformDofAtEachIteration == true)
			this->Clear();

		this->mStep += 1;
		Timer::Stop("solve");

		return 0.0;

		KRATOS_CATCH("")
	}

	void calculateInitialConditions()
	{
		this->calculateEnthalpyNode();
	}

	void calculateEnthalpyNode()
	{

		ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

		//Calculate effective saturation
		const double Swr = rCurrentProcessInfo[RESIDUAL_SW];
		const double Snr = rCurrentProcessInfo[RESIDUAL_SN];

		for (int node = 0; node < static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
		{
			ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;

			const double saturationNon_Node = itNode->FastGetSolutionStepValue(SATURATION_NON);
			const double saturationWet_Node = 1.0 - saturationNon_Node;
			const double SwEfective_Node = (saturationWet_Node - Swr) / (1.0 - Swr - Snr);

			if (SwEfective_Node < 1.0)
			{
				const double pressureWetNode = itNode->FastGetSolutionStepValue(PRESSURE_WET);
				const double pressureCapillarNode = itNode->FastGetSolutionStepValue(CAPILLARITY_PRESSURE_NODE);
				double pressureNonNode = 0.0;
				if (this->mIsCapillarityNeglected != 1)
					pressureNonNode = pressureWetNode + pressureCapillarNode;
				else
					pressureNonNode = pressureWetNode;

				const double temperatureNode = itNode->FastGetSolutionStepValue(TEMPERATURE_NODE);

				double enthalpyNonNode = 0.0;
				this->mSpanWagnerEOS->interpolate_calculate_variable_PT("enthalpy", enthalpyNonNode, pressureNonNode, temperatureNode);
				itNode->FastGetSolutionStepValue(ENTHALPY_NON_NODE) = enthalpyNonNode;

			}
		}
	}

	void calculateLocalEquations()
	{
		//Kr & derivatives
		ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

		for (int node = 0; node < static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
		{
			ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;

			this->calculateMultiphaseProperties(itNode, rCurrentProcessInfo);

			if (this->mIsNonPhaseNonConstant != 1)
				this->calculateThermoNonWetProperties(itNode, rCurrentProcessInfo);

			if (this->mIsWetPhaseNonConstant != 1)
				this->calculateThermoWetProperties(itNode, rCurrentProcessInfo);
		}
	}

	void calculateMultiphaseProperties(ModelPart::NodesContainerType::iterator& aNode, ProcessInfo& rCurrentProcessInfo)
	{
		//Const parameters
		const double Swr = rCurrentProcessInfo[RESIDUAL_SW];
		const double Snr = rCurrentProcessInfo[RESIDUAL_SN];

		//Calculate Saturations
		const double saturationNon_Node = aNode->FastGetSolutionStepValue(SATURATION_NON);
		const double saturationWet_Node = 1.0 - saturationNon_Node;
		aNode->FastGetSolutionStepValue(SATURATION_WET) = saturationWet_Node;
		const double Se_Node = (saturationWet_Node - Swr) / (1.0 - Swr - Snr);

		//Set capillarity to zero (only it will modify in case of capillar problem)
		aNode->FastGetSolutionStepValue(CAPILLARITY_PRESSURE_NODE) = 0.0;
		const double pressureWetNode = aNode->FastGetSolutionStepValue(PRESSURE_WET);
		aNode->FastGetSolutionStepValue(PRESSURE_NON) = aNode->FastGetSolutionStepValue(PRESSURE_WET);

		//Calculate Krn, Krw
		double permRelNon_Node = 0.0;
		double permRelWet_Node = 0.0;
		if (this->mRelativePermeabilityModel == "Simple")
			this->calculate_kr_simple(permRelWet_Node, permRelNon_Node, Se_Node, rCurrentProcessInfo);
		else
			this->calculate_kr(permRelWet_Node, permRelNon_Node, Se_Node, Swr, Snr, rCurrentProcessInfo);

		//Set variables
		aNode->FastGetSolutionStepValue(PERMEABILITY_NON_NODE) = permRelNon_Node;
		aNode->FastGetSolutionStepValue(PERMEABILITY_WET_NODE) = permRelWet_Node;

		//Calculate capillaryPressureModel (Pc, derivPc_Sn, dderivPc_Sn)
		if (this->mIsCapillarityNeglected != 1)
		{
			double Pc_Node = 0.0;
			double derivPc_Sn_Node = 0.0;
			if (this->mCapillarityPressureModel == "BrooksCorey")
				this->calculate_Pc_BrooksCorey(Pc_Node, derivPc_Sn_Node, Se_Node, Swr, Snr, rCurrentProcessInfo);
			else
				this->calculate_Pc(Pc_Node, derivPc_Sn_Node, Se_Node, Swr, Snr, rCurrentProcessInfo);

			aNode->FastGetSolutionStepValue(CAPILLARITY_PRESSURE_NODE) = Pc_Node;
			aNode->FastGetSolutionStepValue(DERIV_PC_SN_NODE) = derivPc_Sn_Node;
			if (this->mIsCapillarityNeglected != 1)
				aNode->FastGetSolutionStepValue(PRESSURE_NON) = pressureWetNode + Pc_Node;

		}

	}

	void calculate_kr_simple(double& rKrw, double& rKrn, const double aSe, ProcessInfo& rCurrentProcessInfo)
	{
		const unsigned int lambda = rCurrentProcessInfo[LAMBDA];
		rKrw = pow(aSe, (double)lambda);
		rKrn = pow(1.0 - aSe, (double)lambda);
	}

	void calculate_kr(double& rKrw, double& rKrn, const double aSe, const double aSwr, const double aSnr, ProcessInfo& rCurrentProcessInfo)
	{

		if (this->mRelativePermeabilityModel == "SimpleConstitutiveLaw")
		{
			const int exp = rCurrentProcessInfo[EXPONENT_KR_LAW];

			if (aSe <  0.0 || aSe > 1.0)
			{
				rKrw = 0.0;
				rKrn = 0.0;
			}
			else
			{
				rKrw = pow(aSe, double(exp));
				rKrn = pow((1.0 - aSe), double(exp));
			}
			//          rDerivKrw_Sw = ( exp/(1.0-aSwr-aSnr) )*pow(aSe, double(exp-1));
			//          rDerivKrn_Sw = (-1.0) * ( exp/(1.0-aSwr-aSnr) )*pow(1.0-aSe, double(exp-1));
		}
		else if (this->mRelativePermeabilityModel == "LeverettJ-Function")
		{
			const int lambda = rCurrentProcessInfo[LAMBDA_LEVERETT];
			//const double Pd = rCurrentProcessInfo[ENTRY_PRESSURE]; 

			if (aSe <  0.0 || aSe > 1.0)
			{
				rKrw = 0.0;
				rKrn = 0.0;
			}
			else
			{
				rKrw = pow(aSe, ((2.0 + 3.0*lambda) / lambda));
				rKrn = pow((1.0 - aSe), 2.0) * (1.0 - pow(aSe, (2.0 + lambda) / lambda));
			}
			//          rDerivKrw_Sw = ((2.0+3.0*lambda)/lambda) * (1.0/(1.0-aSwr-aSnr)) * pow(aSe, (2.0+2.0*lambda) );
			//          rDerivKrn_Sw = (-1.0) * pow((1.0-aSe), 2.0) * ( (2.0+lambda)/lambda ) * pow(aSe, (2.0/lambda) ) -
			//                           2*(1.0-aSe) * (1.0- pow(aSe,(2.0+lambda)/lambda));
		}
		else
		{
			KRATOS_THROW_ERROR(std::logic_error, "Dont exist that relative permeability model", "");
		}

	}
	
	void calculate_Pc_BrooksCorey(double& rPc, double& rdPc_dSn, const double aSe, const double aSwr, const double aSnr, ProcessInfo& rCurrentProcessInfo)
	{

		const int lambda = rCurrentProcessInfo[LAMBDA];
		const double Pd = rCurrentProcessInfo[ENTRY_PRESSURE];

		rPc = Pd* pow(aSe, (-1.0 / lambda));

		if (rPc < 0.0)
			rdPc_dSn = 0.0;
		else
			rdPc_dSn = ((Pd) / (lambda*(1.0 - aSwr - aSnr))) * (pow(aSe, (-(1.0) / (lambda))));

		//            rddPc_ddSw = ( (Pd * (1+lambda))/( pow(lambda, 2.0) *(1.0-aSwr-aSnr)) ) * ( pow(aSe, (-(1.0+2.0*lambda)/lambda ) ) );

	}

    void calculate_Pc(double& rPc, double& rdPc_dSw, const double aSe, const double aSwr, const double aSnr, ProcessInfo& rCurrentProcessInfo)
    {

		if (this->mCapillarityPressureModel == "VanGenuchten") 
        {

        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Dont exist that capillarity pressure model", "");
        }
    }
    
    void calculateThermoNonWetProperties(ModelPart::NodesContainerType::iterator& aNode, ProcessInfo& rCurrentProcessInfo)
    {
		const double Swr = rCurrentProcessInfo[RESIDUAL_SW];
		const double Snr = rCurrentProcessInfo[RESIDUAL_SN];

        const double saturationNon_Node = aNode->FastGetSolutionStepValue(SATURATION_NON);
		const double saturationWet_Node = 1.0 - saturationNon_Node;
		const double pressureWet_Node = aNode->FastGetSolutionStepValue(PRESSURE_WET);
		const double pressureCapillar_Node = aNode->FastGetSolutionStepValue(CAPILLARITY_PRESSURE_NODE);
        const double pressureNon_Node = pressureWet_Node + pressureCapillar_Node;
        const double enthalpyNon_Node = aNode->FastGetSolutionStepValue(ENTHALPY_NON_NODE);
 
		const double SwEfective_Node = (saturationWet_Node - Swr) / (1.0 - Swr - Snr);

		if (SwEfective_Node < 1.0) //NonWet Phase
		{
			std::vector <double> propertiesVec;
			this->mSpanWagnerEOS->computeProperties(propertiesVec, pressureNon_Node, enthalpyNon_Node);

			const double temperature_Node = propertiesVec[0];
			const double densityNon_Node = propertiesVec[1];
			const double specificHeatNon_Node = propertiesVec[2];
			const double viscosityNon_Node = propertiesVec[3];
			const double compressibilityNon_Node = propertiesVec[4];

			if (this->mIsNonIsoThermal != 1 && saturationNon_Node >= 1e-03)
				aNode->FastGetSolutionStepValue(TEMPERATURE_NON_NODE) = temperature_Node;

			aNode->FastGetSolutionStepValue(DENSITY_NON_NODE) = densityNon_Node;
			aNode->FastGetSolutionStepValue(SPECIFIC_HEAT_NON_NODE) = specificHeatNon_Node;
			aNode->FastGetSolutionStepValue(VISCOSITY_NON_NODE) = viscosityNon_Node;
			aNode->FastGetSolutionStepValue(COMPRESSIBLITY_NON_NODE) = compressibilityNon_Node;
		}
		else //Wet Phase
		{
			aNode->FastGetSolutionStepValue(DENSITY_NON_NODE) = 0.0;
			aNode->FastGetSolutionStepValue(SPECIFIC_HEAT_NON_NODE) = 0.0;
			aNode->FastGetSolutionStepValue(VISCOSITY_NON_NODE) = 0.0;
			aNode->FastGetSolutionStepValue(COMPRESSIBLITY_NON_NODE) = 0.0;
		}

    }
    
    void calculateThermoWetProperties(ModelPart::NodesContainerType::iterator& aNode, ProcessInfo& rCurrentProcessInfo)
    {
        
    }
    
   // void calculatePhasesDarcyFlow()
   // {
   //     ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();    
   //     
   //     //Parameters
   //     const unsigned int isCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
   //     const unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
   //     const unsigned int NDim = mSolverStrategyConfiguration.GetDomainSize();
   //     
   //     for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
   //     {

   //         ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;

   //         unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();
   //     
   //         Vector gravity = ZeroVector(3);
   //         if(isBuoyancy !=1) gravity = rCurrentProcessInfo[GRAVITY]; 
   //         
   //         //Calculate effective saturation            
   //         const double saturationWet_Elem = itElem->GetValue(SATURATION_WET_ELEM);
   //         const double Swr_Elem  = rCurrentProcessInfo[RESIDUAL_SW];
   //         const double Snr_Elem  = rCurrentProcessInfo[RESIDUAL_SN];    
   //         const double Se_Elem = (saturationWet_Elem-Swr_Elem)/(1.0-Swr_Elem-Snr_Elem);
   //          double Id = itElem->Id();
   //          double AASw_Elem = itElem->GetValue(SATURATION_WET_ELEM);
   //          double AAPn_Elem = itElem->GetValue(PRESSURE_NON_ELEM);
   //          double AAHn_Elem = itElem->GetValue(ENTHALPY_NON_ELEM);
   //         cout << "Strategy: "<< " Id " << Id << " Sw " << AASw_Elem << " Pn " << AAPn_Elem << " Hn " << AAHn_Elem << endl;
   //         
   //         //Constant
   //         const double  intrinsecPermeability_Elem = itElem->GetProperties()[PERMEABILITY_ELEM];
   //         
   //         //Non-Constant                   
   //         const double permeabilityNon_Elem = itElem->GetValue(PERMEABILITY_NON_ELEM);
   //         const double permeabilityWet_Elem = itElem->GetValue(PERMEABILITY_WET_ELEM);
   //         const double densityWet_Elem = itElem->GetProperties()[DENSITY_WET_ELEM]; 
   //         const double viscosityWet_Elem = itElem->GetProperties()[VISCOSITY_WET_ELEM]; 
   //         double densityNon_Elem = 0.0;
   //         double viscosityNon_Elem = 0.0;

   //         //Calculate capillaryPressureModel (derivPc_Sw)
   //         double derivPc_Sw_Elem = 0.0;
   //         if(isCapillarityNeglected !=1)
   //             derivPc_Sw_Elem = itElem->GetValue(DERIV_PC_SW_ELEM);
   //         
   //         if(this->mIsNonPhaseNonConstant !=1)
   //         {
   //             viscosityNon_Elem = itElem->GetValue(VISCOSITY_NON_ELEM);
   //             densityNon_Elem = itElem->GetValue(DENSITY_NON_ELEM);
   //         }
   //         else
   //         {
   //             viscosityNon_Elem = itElem->GetProperties()[VISCOSITY_NON_ELEM];
   //             densityNon_Elem = itElem->GetProperties()[DENSITY_NON_ELEM];
   //         }

   //         //Element variable non-constant (no properties)
   //         const double mobility_wet_Elem = permeabilityWet_Elem/viscosityWet_Elem;
   //         const double mobility_non_Elem = permeabilityNon_Elem/viscosityNon_Elem;

   //         ////Element Geometry (N,gradN)
   //         //area
   //         double domain_L_A_V = 0.0;
   //         //GradN
   //         Matrix DN_DX = ZeroMatrix(numberOfNodes,NDim);
   //         //N
   //         Vector N = ZeroVector(numberOfNodes); 
   //         GeometryUtils::CalculateGeometryData(itElem->GetGeometry(), DN_DX, N, domain_L_A_V); 
   //         
			////Case of 1D problems
			//if (numberOfNodes == 2)
			//	domain_L_A_V = 1.0;

   //         //Compute gradient of head
   //         Vector gradPn = ZeroVector(NDim);
   //         Vector gradPn_NonBuoy = ZeroVector(NDim);
   //         Vector gradPn_WetBuoy = ZeroVector(NDim);
   //         Vector saturationWetGrad = ZeroVector(NDim);

   //         for (unsigned int node = 0; node < numberOfNodes ; node++)
   //         {
   //            const double currentPressureNon_Node =  itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_NON);
   //            const double currentSaturationWet_Node =  itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_WET);
   //            for (unsigned int dim = 0; dim < NDim; dim++)
   //             {
   //                gradPn[dim]+=currentPressureNon_Node*DN_DX.at_element(node,dim);
   //                saturationWetGrad[dim]+=currentSaturationWet_Node*DN_DX.at_element(node,dim);
   //             }

   //         }

   //         gradPn_NonBuoy = gradPn;
   //         gradPn_WetBuoy = gradPn;
   //         
   //         for (unsigned int dim = 0; dim < NDim; dim++)
   //         {
   //            gradPn_NonBuoy[dim] -= gravity[dim] * (densityNon_Elem) ;
   //            gradPn_WetBuoy[dim] -= gravity[dim] * (densityWet_Elem) ;
   //         }

   //         //K*gradh, isotropic
   //         Matrix Kintrinsec_Matrix = ZeroMatrix(NDim,NDim);
   //         for (unsigned int dim = 0; dim < NDim; dim++)
   //         {
   //             Kintrinsec_Matrix(dim,dim)= intrinsecPermeability_Elem;
   //         }

   //         Vector nonDarcyFlow_Elem = ZeroVector(NDim);
   //         Vector wetDarcyFlow_Elem = ZeroVector(NDim);
   //         
   //         if( Se_Elem >= 0.99999999999 )
   //         {
   //             wetDarcyFlow_Elem =     -mobility_wet_Elem*prod(Kintrinsec_Matrix,trans(gradPn_WetBuoy))
   //                                     +mobility_wet_Elem*derivPc_Sw_Elem*prod(Kintrinsec_Matrix,trans(saturationWetGrad));
   //                             
   //         }
   //         else if( Se_Elem <= 0.00000000001 )
   //         {
   //             nonDarcyFlow_Elem =     -mobility_non_Elem*prod(Kintrinsec_Matrix,trans(gradPn_NonBuoy));
   //         }
   //         else
   //         {
   //             nonDarcyFlow_Elem =     -mobility_non_Elem*prod(Kintrinsec_Matrix,trans(gradPn_NonBuoy));
   //             wetDarcyFlow_Elem =     -mobility_wet_Elem*prod(Kintrinsec_Matrix,trans(gradPn_WetBuoy))
   //                                     +mobility_wet_Elem*derivPc_Sw_Elem*prod(Kintrinsec_Matrix,trans(saturationWetGrad));
   //         }
   //         
   //         itElem->SetValue(DARCYFLOW_NON,nonDarcyFlow_Elem);
   //         itElem->SetValue(DARCYFLOW_WET,wetDarcyFlow_Elem);
   //         
   //     }
   // }
    

	void calculatePhasesDarcyFlow()
	{
		ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
		const unsigned int NDim = mSolverStrategyConfiguration.GetDomainSize();

		//Parameters
		const unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
		const unsigned int isCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
		const unsigned int isWetPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_WET_PHASE];
		const unsigned int isNonPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_NON_PHASE];
		const unsigned int Snr = rCurrentProcessInfo[RESIDUAL_SN];
		const unsigned int Swr = rCurrentProcessInfo[RESIDUAL_SW];

		Vector gravity = ZeroVector(3);
		if (isBuoyancy != 1) gravity = rCurrentProcessInfo[GRAVITY];

		for (auto itElem = BaseType::GetModelPart().ElementsBegin();
			itElem != BaseType::GetModelPart().ElementsEnd(); ++itElem)
		{
			//ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;

			unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();

			//Constant
			const double porosity_Elem = itElem->GetProperties()[POROSITY_ELEM];
			const double  intrinsecPermeability_Elem = itElem->GetProperties()[PERMEABILITY_ELEM];
			const double compressibilityMedium_Elem = itElem->GetProperties()[COMPRESSIBILITY_MEDIUM_ELEM];
			const double thermalExpansionNon_Elem = itElem->GetProperties()[THERMAL_COEFF_NON_ELEM];

			//Non-Constant
			const double saturationNon_Elem = itElem->GetValue(SATURATION_NON_ELEM);
			const double saturationWet_Elem = 1.0 - saturationNon_Elem;
			const double Se_Elem = (saturationWet_Elem - Swr) / (1.0 - Snr - Swr);
			const double permeabilityNon_Elem = itElem->GetValue(PERMEABILITY_NON_ELEM);
			const double permeabilityWet_Elem = itElem->GetValue(PERMEABILITY_WET_ELEM);

			//Calculate capillaryPressureModel (Pc, derivPc_Sw, dderivPc_Sw)
			double derivPc_Sn_Elem = 0.0;
			if (isCapillarityNeglected != 1)
				derivPc_Sn_Elem = itElem->GetValue(DERIV_PC_SN_ELEM);

			//Compressible NonWetting phase case
			double densityNon_Elem = 0.0;
			double viscosityNon_Elem = 0.0;
			double compressibilityNon_Elem = 0.0;
			if (isNonPhaseNonConstant != 1)
			{
				densityNon_Elem = itElem->GetValue(DENSITY_NON_ELEM);
				viscosityNon_Elem = itElem->GetValue(VISCOSITY_NON_ELEM);
				if (viscosityNon_Elem == 0.0) viscosityNon_Elem = 1.0; //Case of wetting phase
				compressibilityNon_Elem = itElem->GetValue(COMPRESSIBLITY_NON_ELEM);
			}
			else
			{
				densityNon_Elem = itElem->GetProperties()[DENSITY_NON_ELEM];
				viscosityNon_Elem = itElem->GetProperties()[VISCOSITY_NON_ELEM];
				compressibilityNon_Elem = itElem->GetProperties()[COMPRESSIBLITY_NON_ELEM];
			}

			//Compressible Wetting phase case
			double densityWet_Elem = 0.0;
			double viscosityWet_Elem = 0.0;
			if (isWetPhaseNonConstant != 1)
			{
				densityWet_Elem = itElem->GetValue(DENSITY_WET_ELEM);
				viscosityWet_Elem = itElem->GetValue(VISCOSITY_WET_ELEM);
				if (viscosityWet_Elem == 0.0) viscosityWet_Elem = 1.0; //Case of nonwetting phase
			}
			else
			{
				densityWet_Elem = itElem->GetProperties()[DENSITY_WET_ELEM];
				viscosityWet_Elem = itElem->GetProperties()[VISCOSITY_WET_ELEM];
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
			Matrix DN_DX = ZeroMatrix(numberOfNodes, NDim);
			//N
			Vector N = ZeroVector(numberOfNodes);
			GeometryUtils::CalculateGeometryData(itElem->GetGeometry(), DN_DX, N, domain_L_A_V);

			//Compute gradient of head
			Vector gradPw = ZeroVector(NDim);
			Vector gradPw_NonBuoy = ZeroVector(NDim);
			Vector gradPw_WetBuoy = ZeroVector(NDim);
			Vector saturationNonGrad = ZeroVector(NDim);

			for (unsigned int node = 0; node < numberOfNodes; node++)
			{
				const double currentPressureWet_Node = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET);
				const double currentSaturationNon_Node = itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON);
				for (unsigned int dim = 0; dim < NDim; dim++)
				{
					gradPw[dim] += currentPressureWet_Node*DN_DX.at_element(node, dim);
					saturationNonGrad[dim] += currentSaturationNon_Node*DN_DX.at_element(node, dim);
				}

			}

			gradPw_NonBuoy = gradPw;
			gradPw_WetBuoy = gradPw;

			for (unsigned int dim = 0; dim < NDim; dim++)
			{
				gradPw_NonBuoy[dim] -= gravity[dim] * (densityNon_Elem);
				gradPw_WetBuoy[dim] -= gravity[dim] * (densityWet_Elem);
			}

			//K*gradh, isotropic
			Matrix Kintrinsec_Matrix = ZeroMatrix(NDim, NDim);
			for (unsigned int dim = 0; dim < NDim; dim++)
			{
				Kintrinsec_Matrix(dim, dim) = intrinsecPermeability_Elem;
			}

			Vector nonDarcyFlow_Elem = ZeroVector(NDim);
			Vector wetDarcyFlow_Elem = ZeroVector(NDim);

			if (Se_Elem >= 1.0 /*0.9999*/)
			{
				wetDarcyFlow_Elem = -mobility_wet_Elem*prod(Kintrinsec_Matrix, trans(gradPw_WetBuoy));

			}
			else if (Se_Elem <= 0.0)
			{
				nonDarcyFlow_Elem = -mobility_non_Elem*prod(Kintrinsec_Matrix, trans(gradPw_NonBuoy));
									-mobility_non_Elem*derivPc_Sn_Elem*prod(Kintrinsec_Matrix, trans(saturationNonGrad));
			}
			else
			{
				nonDarcyFlow_Elem = -mobility_non_Elem*prod(Kintrinsec_Matrix, trans(gradPw_NonBuoy));
									-mobility_non_Elem*derivPc_Sn_Elem*prod(Kintrinsec_Matrix, trans(saturationNonGrad));
				wetDarcyFlow_Elem = -mobility_wet_Elem*prod(Kintrinsec_Matrix, trans(gradPw_WetBuoy));
			}

			itElem->SetValue(DARCYFLOW_NON, nonDarcyFlow_Elem);
			itElem->SetValue(DARCYFLOW_WET, wetDarcyFlow_Elem);
		}
	}

	void CalculateSumBalances()
	{
		//Set to zero balance fields
		for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
			node != BaseType::GetModelPart().NodesEnd(); ++node)
		{
			node->FastGetSolutionStepValue(STORAGE_SUM_BALANCE) = 0.0;
			node->FastGetSolutionStepValue(DARCY_FLOW_SUM_BALANCE) = 0.0;
			node->FastGetSolutionStepValue(DIRICHLET_SUM_BALANCE) = 0.0;
			node->FastGetSolutionStepValue(NEUMANN_SUM_BALANCE) = 0.0;
		}

		//Process info parameters
		ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

		//Parameters
		const unsigned int NDim = mSolverStrategyConfiguration.GetDomainSize();
		//Parameters
		const unsigned int isTransient = rCurrentProcessInfo[IS_TRANSIENT];
		const unsigned int isBuoyancy = rCurrentProcessInfo[IS_BUOYANCY];
		const unsigned int isCapillarityNeglected = rCurrentProcessInfo[IS_CAPILLARITY_NEGLECTED];
		const unsigned int isNonPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_NON_PHASE];
		const unsigned int isWetPhaseNonConstant = rCurrentProcessInfo[IS_NONCONSTANT_WET_PHASE];
		double invDt = 0.0;
		if (isTransient != 1)
		{
			const double deltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
			invDt = 1.0 / deltaTime;
		}

		Vector gravity = ZeroVector(3);
		gravity = rCurrentProcessInfo[GRAVITY];

		//Loop over elements
		for (int elem = 0; elem < static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
		{
			ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin() + elem;

			unsigned int numberOfNodes = itElem->GetGeometry().PointsNumber();

			//Constant
			const double porosity_Elem = itElem->GetProperties()[POROSITY_ELEM];
			const double  intrinsecPermeability_Elem = itElem->GetProperties()[PERMEABILITY_ELEM];
			const double compressibilityMedium_Elem = itElem->GetProperties()[COMPRESSIBILITY_MEDIUM_ELEM];
			const double thermalExpansionWet_Elem = itElem->GetProperties()[THERMAL_COEFF_WET_ELEM];
			const double thermalExpansionNon_Elem = itElem->GetProperties()[THERMAL_COEFF_NON_ELEM];

			//Non-Constant
			const double saturationNon_Elem = itElem->GetValue(SATURATION_NON_ELEM);
			const double saturationWet_Elem = 1.0 - saturationNon_Elem;
			const double permeabilityNon_Elem = itElem->GetValue(PERMEABILITY_NON_ELEM);
			const double permeabilityWet_Elem = itElem->GetValue(PERMEABILITY_WET_ELEM);

			//Calculate capillaryPressureModel (Pc, derivPc_Sw, dderivPc_Sw)
			double derivPc_Sn_Elem = 0.0;
			if (isCapillarityNeglected != 1)
				derivPc_Sn_Elem = itElem->GetValue(DERIV_PC_SN_ELEM);

			//Compressible NonWetting phase case
			double densityNon_Elem = 0.0;
			double viscosityNon_Elem = 0.0;
			double compressibilityNon_Elem = 0.0;
			if (isNonPhaseNonConstant != 1)
			{
				densityNon_Elem = itElem->GetValue(DENSITY_NON_ELEM);
				viscosityNon_Elem = itElem->GetValue(VISCOSITY_NON_ELEM);
				if (viscosityNon_Elem == 0.0) viscosityNon_Elem = 1.0; //Case of wetting phase
				compressibilityNon_Elem = itElem->GetValue(COMPRESSIBLITY_NON_ELEM);
			}
			else
			{
				densityNon_Elem = itElem->GetProperties()[DENSITY_NON_ELEM];
				viscosityNon_Elem = itElem->GetProperties()[VISCOSITY_NON_ELEM];
				compressibilityNon_Elem = itElem->GetProperties()[COMPRESSIBLITY_NON_ELEM];
			}

			//Compressible Wetting phase case
			double densityWet_Elem = 0.0;
			double viscosityWet_Elem = 0.0;
			double compressibilityWet_Elem = 0.0;
			if (isWetPhaseNonConstant != 1)
			{
				densityWet_Elem = itElem->GetValue(DENSITY_WET_ELEM);
				viscosityWet_Elem = itElem->GetValue(VISCOSITY_WET_ELEM);
				if (viscosityWet_Elem == 0.0) viscosityWet_Elem = 1.0; //Case of nonwetting phase
				compressibilityWet_Elem = itElem->GetValue(COMPRESSIBLITY_WET_ELEM);
			}
			else
			{
				densityWet_Elem = itElem->GetProperties()[DENSITY_WET_ELEM];
				viscosityWet_Elem = itElem->GetProperties()[VISCOSITY_WET_ELEM];
				compressibilityWet_Elem = itElem->GetProperties()[COMPRESSIBLITY_WET_ELEM];
			}

			//Element variable non-constant (no properties)
			const double mobility_wet_Elem = permeabilityWet_Elem / viscosityWet_Elem;
			const double mobility_non_Elem = permeabilityNon_Elem / viscosityNon_Elem;
			const double totalMobility_Elem = mobility_wet_Elem + mobility_non_Elem;

			////Element Geometry (N,gradN, gradgradN)    
			//area
			double domain_L_A_V = 0.0;

			//GradGradN
			boost::numeric::ublas::bounded_matrix<double, 3, 3 > DN_DXX = ZeroMatrix(3, 3); // 1D : 1,1
			//GradN
			boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX = ZeroMatrix(3, 2); // 1D :  2,1
			//N
			array_1d<double, 3 > N = ZeroVector(3);// 1D : 2

			//Elemental matrix
			Matrix NiNj = ZeroMatrix(numberOfNodes, numberOfNodes);
			this->CalculateNiNj(NiNj, numberOfNodes, NDim);

			//Geometry staff
			GeometryUtils::CalculateGeometryData(itElem->GetGeometry(), DN_DXX, DN_DX, N, domain_L_A_V);

			//Case of 1D problems
			//if (numberOfNodes == 2)
			//	domain_L_A_V = 1.0;

			//State variables: Pwk, Pwk+1 & IncrPw, Snk, Snk+1 & IncrSw, Hnk, Hnk+1 & IncrHn, Tk, Tk+1 & IncrT,
			//array_1d<double,TNumNodes>  
			Vector pressureWetStepK_1 = ZeroVector(numberOfNodes);
			Vector pressureWetStepK = ZeroVector(numberOfNodes);
			Vector diffPressureWet = ZeroVector(numberOfNodes);
			Vector saturationNonStepK_1 = ZeroVector(numberOfNodes);
			Vector saturationNonStepK = ZeroVector(numberOfNodes);
			Vector diffSaturationNon = ZeroVector(numberOfNodes);
			Vector enthalpyNonIncrement = ZeroVector(numberOfNodes);
			Vector temperatureIncrement = ZeroVector(numberOfNodes);
			for (unsigned int node = 0; node< numberOfNodes; node++)
			{
				const double Pw_k_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET);
				const double Pw_k = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET, 1);
				const double IncrementPw = Pw_k_1 - Pw_k;
				const double Sn_k_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON);
				const double Sn_k = itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON, 1);
				const double IncrementSn = Sn_k_1 - Sn_k;
				const double Hn_k_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE);
				const double Hn_k = itElem->GetGeometry()[node].FastGetSolutionStepValue(ENTHALPY_NON_NODE, 1);
				const double IncrementHn = Hn_k_1 - Hn_k;
				const double T_k_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(TEMPERATURE_NODE);
				const double T_k = itElem->GetGeometry()[node].FastGetSolutionStepValue(TEMPERATURE_NODE, 1);
				const double IncrementT = T_k_1 - T_k;

				pressureWetStepK_1[node] = Pw_k_1;
				pressureWetStepK[node] = Pw_k;
				saturationNonStepK_1[node] = Sn_k_1;
				saturationNonStepK[node] = Sn_k;
				diffPressureWet[node] = IncrementPw;
				diffSaturationNon[node] = IncrementSn;
				enthalpyNonIncrement[node] = IncrementHn;
				temperatureIncrement[node] = IncrementT;
			}

			////DARCY_FLOW_SUM_BALANCE
			Matrix MatElemPressure = ZeroMatrix(numberOfNodes, numberOfNodes);
			Matrix MatElemCapillarity = ZeroMatrix(numberOfNodes, numberOfNodes);

			const int ntensor =  numberOfNodes;//1D... 1
			//K_darcyTerm... 1D 
			Vector K_darcyTerm = ZeroVector(numberOfNodes);
			K_darcyTerm[0] = intrinsecPermeability_Elem * totalMobility_Elem;
			K_darcyTerm[1] = 0.0;
			K_darcyTerm[2] = intrinsecPermeability_Elem * totalMobility_Elem;

			//K_diffCapillarity... 1D
			Vector K_Capillarity = ZeroVector(numberOfNodes);
			K_Capillarity[0] = intrinsecPermeability_Elem * derivPc_Sn_Elem * mobility_non_Elem;
			K_Capillarity[1] = 0.0;
			K_Capillarity[2] = intrinsecPermeability_Elem * derivPc_Sn_Elem * mobility_non_Elem;

			//K_diffBuoyancy... 1D
			Vector K_Buoyancy = ZeroVector(numberOfNodes);
			K_Buoyancy[0] = intrinsecPermeability_Elem*(mobility_wet_Elem*densityWet_Elem + mobility_non_Elem * densityNon_Elem);
			K_Buoyancy[1] = 0.0;
			K_Buoyancy[2] = intrinsecPermeability_Elem*(mobility_wet_Elem*densityWet_Elem + mobility_non_Elem * densityNon_Elem);

			//Loop over nodes            
			int k = 0;    // Row of GradGrad matrix
			for (int i = 0; i< numberOfNodes - 1; i++)    // Loop of columns of local matrix
			{
				for (int j = i + 1; j < numberOfNodes; j++)	// Loop of columns of local matrix
				{
					for (int itensor = 0; itensor < ntensor; itensor++)	// Loop columns of GradGrad matrix (or of tensor values)
					{
						double elem_ij_Pressure = -DN_DXX.at_element(k, itensor)* K_darcyTerm[itensor];
						double elem_ij_Capillarity = -DN_DXX.at_element(k, itensor)* K_Capillarity[itensor];

						MatElemPressure(i, j) += elem_ij_Pressure;
						MatElemPressure(j, i) += elem_ij_Pressure;
						MatElemPressure(i, i) -= elem_ij_Pressure;
						MatElemPressure(j, j) -= elem_ij_Pressure;

						MatElemCapillarity(i, j) += elem_ij_Capillarity;
						MatElemCapillarity(j, i) += elem_ij_Capillarity;
						MatElemCapillarity(i, i) -= elem_ij_Capillarity;
						MatElemCapillarity(j, j) -= elem_ij_Capillarity;

					}
					k++;
				}
			}
			// Set DARCY_FLOW_SUM_BALANCE
			Vector balanceDarcyFlowSum = ZeroVector(numberOfNodes);

			Vector balancePressure = ZeroVector(numberOfNodes);
			noalias(balancePressure) = prod(MatElemPressure, pressureWetStepK_1);
			Vector balanceCapillarity = ZeroVector(numberOfNodes);
			noalias(balanceCapillarity) = prod(MatElemCapillarity, saturationNonStepK_1);
			Vector balanceBuoyancy = ZeroVector(numberOfNodes);
			//loop over nodes of element
			for (int inode = 0; inode < numberOfNodes; inode++)
			{
				double result = 0;
				//double densityNode = itElem->GetGeometry()[inode].FastGetSolutionStepValue(DENSITY);
				//loop over spatial dimensions
				for (int idim = 0; idim < NDim; idim++)
				{
					const int k_index = idim*NDim; //Only 2D           
					result += K_Buoyancy[k_index] * gravity[idim] * DN_DX.at_element(inode, idim)*domain_L_A_V;
				}

				//store result 
				balanceBuoyancy[inode] = result;

			}

			balanceDarcyFlowSum = balancePressure + balanceCapillarity + balanceBuoyancy;

			////STORAGE_SUM_BALANCE
			Vector storageCompressibilityBalance = ZeroVector(numberOfNodes);
			Vector storageCapillarityBalance = ZeroVector(numberOfNodes);
			if (isTransient != 1)
			{
				//Solid & fluids compressibility storage LHS. Is medium Compressibility negative ??
				const double wet_comp_coeff = saturationWet_Elem*compressibilityWet_Elem;
				const double non_comp_coeff = saturationNon_Elem*compressibilityNon_Elem;
				const double total_comp_coeff = (porosity_Elem*(wet_comp_coeff + wet_comp_coeff) - compressibilityMedium_Elem);
				noalias(storageCompressibilityBalance) = domain_L_A_V * invDt * total_comp_coeff *prod(NiNj, diffPressureWet);

				//Capillarity storage ONLY RHS. Is medium Compressibility negative ??
				const double cap_storage_coeff = ((porosity_Elem * compressibilityNon_Elem) - compressibilityMedium_Elem) * derivPc_Sn_Elem * saturationNon_Elem;
				noalias(storageCapillarityBalance) = domain_L_A_V * invDt * cap_storage_coeff * prod(NiNj, diffSaturationNon);
			}

			////DIRICHLET_SUM_BALANCE & set BALANCES
			//Loop over nodes 
			for (int inode = 0; inode < numberOfNodes; inode++)
			{
				double storageBalanceNode = 0.0;
				if (isTransient != 1)
				{
					storageBalanceNode = storageCompressibilityBalance[inode] + storageCapillarityBalance[inode];
					itElem->GetGeometry()[inode].FastGetSolutionStepValue(STORAGE_SUM_BALANCE) += storageBalanceNode;
				}

				itElem->GetGeometry()[inode].FastGetSolutionStepValue(DARCY_FLOW_SUM_BALANCE) += balanceDarcyFlowSum[inode];

				if (itElem->GetGeometry()[inode].IsFixed(PRESSURE_WET))
				{
					const double dirichletBalance = storageBalanceNode - balanceDarcyFlowSum[inode];
					itElem->GetGeometry()[inode].FastGetSolutionStepValue(DIRICHLET_SUM_BALANCE) += dirichletBalance;
				}
			}
		}

		////NEUMANN_SUM_BALANCE
		//Loop over conditions 
		for (int cond = 0; cond<static_cast<int>(BaseType::GetModelPart().Conditions().size()); cond++)
		{
			ModelPart::ConditionsContainerType::iterator itCond = BaseType::GetModelPart().ConditionsBegin() + cond;

			unsigned int numberOfNodes = itCond->GetGeometry().PointsNumber();

			const char ConditionName1[] = "PointSinkSourceNonFEM";
			Condition const& ref_cond1 = KratosComponents<Condition>::Get(ConditionName1);

			if (typeid (ref_cond1) == typeid (*itCond))
			{
				//Loop over nodes 
				for (int inode = 0; inode < numberOfNodes; inode++)
				{
					const double sinkSourceNon = itCond->GetProperties()[SINK_SOURCE_NON];
					itCond->GetGeometry()[inode].FastGetSolutionStepValue(NEUMANN_SUM_BALANCE) = sinkSourceNon;
				}
			}
		}
	}


	void CalculateNiNj(Matrix &rNiNj, unsigned int nNodes, unsigned int nDim)
	{
		rNiNj = ZeroMatrix(nNodes, nNodes);

		const double lumpedFactor = 1.0 / double(nNodes);
		rNiNj = lumpedFactor * IdentityMatrix(nNodes, nNodes);

	}

	void CalculateNiNj(Matrix &rNiNj, const Vector& aN, unsigned int nNodes, unsigned int nDim)
	{
		rNiNj = ZeroMatrix(nNodes, nNodes);

		for (unsigned int node_i = 0; node_i < nNodes; node_i++)
			for (unsigned int node_j = 0; node_j < nNodes; node_j++)
				rNiNj(node_i, node_j) += aN[node_i] * aN[node_j];
	}

	void CalculateBiBj(Matrix& rBiBj, Matrix aDN_DX, unsigned int nNodes, unsigned int nDim)
	{
		rBiBj = ZeroMatrix(nNodes, nNodes);

		noalias(rBiBj) = prod(aDN_DX, trans(aDN_DX));
	}

    template< class TVariableType >
    void StorageOldIteration(const Kratos::Variable<TVariableType> aVariable_Current,
                             const Kratos::Variable<TVariableType> aVariable_Old)
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
             node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            //setting the old value of the pressure & concentration to the current one
            const double variable1 = (node)->FastGetSolutionStepValue(aVariable_Current);
            (node)->FastGetSolutionStepValue(aVariable_Old) = variable1;

            const double variable2 = (node)->FastGetSolutionStepValue(aVariable_Current, 1);
            (node)->FastGetSolutionStepValue(aVariable_Old,1) = variable2;

            const double variable3 = (node)->GetValue(aVariable_Current);
            (node)->SetValue(aVariable_Old,variable3);
        }
       
        KRATOS_CATCH("")
    }

	template< class TVariableType >
	void SetPreviousTimeToOldIterationSV(const Kratos::Variable<TVariableType> aCorrectedVariable,
		const Kratos::Variable<TVariableType> aStateVariable)
	{
		KRATOS_TRY

			for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
				node != BaseType::GetModelPart().NodesEnd(); ++node)
		{
			//setting the old value of the pressure & concentration to the current one
			const double correctedVariable_k = (node)->FastGetSolutionStepValue(aStateVariable, 1);
			(node)->FastGetSolutionStepValue(aCorrectedVariable) = correctedVariable_k;
		}

		KRATOS_CATCH("")
	}

   double systemIteration(const int aFractionalStep)
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        rCurrentProcessInfo[FRACTIONAL_STEP] = aFractionalStep;
        
        double normDx = 0.0;
        switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
        {
            case 1:
            {
                normDx = mSumPhasesEqStrategy_Pw->Solve();
                break;
            }
            case 2:
            {
                normDx = mNonPhaseEqStrategy_Sn->Solve();
                break;
            }
            case 3:
            {
                normDx = mNonWetPhaseEnergyEqStrategy_Hn->Solve();
                break;
            }
            case 4:
            {
                normDx = mWetPhaseEnergyEqStrategy_T->Solve();
                break;
            }
            default:
            {
                KRATOS_THROW_ERROR(std::logic_error,"SequentialStrategy, systemIteration: Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
            }
        }
            
        return normDx;

        KRATOS_CATCH("");
    }
   
    template< class TVariableType >    
    bool isLocalConverged(const std::string aName_SV,
                             const Kratos::Variable<TVariableType> aVariable_Current,
                             const Kratos::Variable<TVariableType> aVariable_Old, 
                             const double aTolerance,
                             double& aMaxDiff)
    {
        KRATOS_TRY;
        
        double maxUpdate_SV = 0.0;
        bool isLocalConverged = false;
        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
               node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            const double oldIter_SV = (node)->FastGetSolutionStepValue(aVariable_Old);
            const double newIter_SV = (node)->FastGetSolutionStepValue(aVariable_Current);
            const double currentDiffSV_Node = std::abs (double (oldIter_SV - newIter_SV));

            if( currentDiffSV_Node > maxUpdate_SV )
                maxUpdate_SV = currentDiffSV_Node;
        }
        
        std::cout << aName_SV << " maxDiff: " << maxUpdate_SV <<std::endl;
        aMaxDiff = maxUpdate_SV;
        
        if(maxUpdate_SV < aTolerance)
        {
            std::cout << aName_SV << " convergence achieved " << std::endl;
            isLocalConverged= true;
        }
       
        return isLocalConverged;
        
         KRATOS_CATCH("")
        
    }
    
    bool checkNumberIterations(const int aGlobalIteration, const int aMaxIteration_SV)
    {
        bool isReachedMaxIter =false;
        
        if(aGlobalIteration > aMaxIteration_SV)
            isReachedMaxIter= true;
        
        return isReachedMaxIter;
    }
    
	void writeResults(std::string aTitol, const int aTime, const unsigned int aIteration = 0)
	{
		ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
		const int printInterval = rCurrentProcessInfo[TIMES_TO_PRINT];

		if (/*this->mStep%printInterval == 0*/ true)
		{
			std::ofstream myfile;
			myfile.open("propertiesModel.txt", std::ios::app);
			myfile.precision(10);

			myfile << "Concepte: " << aTitol << std::endl << std::endl;

			myfile << "Time: " << aTime << " Iteration: " << aIteration << std::endl << std::endl;

			for (int elem = 0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
			{
				ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin() + elem;

				const double permNonElem1 = itElem->GetValue(PERMEABILITY_NON_ELEM);
				const double permWetElem1 = itElem->GetValue(PERMEABILITY_WET_ELEM);
				const double viscosityNonElem = itElem->GetValue(VISCOSITY_NON_ELEM);
				const double densityNonElem = itElem->GetValue(DENSITY_NON_ELEM);
				const double compresssiblityNon = itElem->GetValue(COMPRESSIBLITY_NON_ELEM);
				Vector darcyNon = itElem->GetValue(DARCYFLOW_NON);
				Vector darcyWet = itElem->GetValue(DARCYFLOW_WET);
				Vector darcyTot = itElem->GetValue(DARCYFLOW_TOTAL);
				const double  SwElem = itElem->GetValue(SATURATION_WET_ELEM);
				const double  SnElem = itElem->GetValue(SATURATION_NON_ELEM);

				myfile << "Elem:                     " << itElem->Id() << std::endl;
				myfile << "SATURATION_NON_ELEM       " << SnElem << std::endl;
				myfile << "SATURATION_WET_ELEM               " << SwElem << std::endl;
				myfile << std::endl;
				myfile << "PERMEABILITY_NON_ELEM     " << permNonElem1 << std::endl;
				myfile << "PERMEABILITY_WET_ELEM     " << permWetElem1 << std::endl;
				myfile << std::endl;
				myfile << "VISCOSITY_NON_ELEM        " << viscosityNonElem << std::endl;
				myfile << std::endl;
				myfile << "DENSITY_NON_ELEM          " << densityNonElem << std::endl;
				myfile << std::endl;
				myfile << "COMPRESSIBLITY_NON_ELEM   " << compresssiblityNon << std::endl;
				myfile << std::endl;
				myfile << "DARCYFLOW_NON             " << darcyNon << std::endl;
				myfile << "DARCYFLOW_WET             " << darcyWet << std::endl;
				myfile << "DARCYFLOW_TOTAL           " << darcyTot << std::endl;
				myfile << std::endl << std::endl << std::endl;

				for (unsigned int node = 0; node < itElem->GetGeometry().PointsNumber(); node++)
				{
					const double PwNode = itElem->GetGeometry()[node].FastGetSolutionStepValue(PRESSURE_WET);
					const double SnNode = itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON);
					const double SnNode_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON, 1);
					const double SnNode_old = itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON_OLD_ITER);
					const double SnNode_old_1 = itElem->GetGeometry()[node].FastGetSolutionStepValue(SATURATION_NON_OLD_ITER, 1);
					const double derivFn = itElem->GetGeometry()[node].FastGetSolutionStepValue(DERIV_PC_SN_NODE);


					myfile << "                           Node                 " << itElem->GetGeometry()[node].Id() << std::endl;
					myfile << "                           PRESSURE_Z           " << PwNode << std::endl;
					myfile << "                           SATURATION_NON       " << SnNode << std::endl;
					myfile << "                           SATURATION_NON_K     " << SnNode_1 << std::endl;
					myfile << std::endl;
					myfile << "                           SATURATION_NON_OLD   " << SnNode_old << std::endl;
					myfile << "                           SATURATION_NON_OLD_K " << SnNode_old_1 << std::endl;
					myfile << std::endl;
				}

				myfile << std::endl << std::endl;
			}

			myfile << std::endl << std::endl;

			myfile.close();
		}

		if (/*this->mStep%printInterval == 0*/ true)
		{
			std::ofstream myfile2;
			myfile2.open("Solution.txt", std::ios::app);
			myfile2.precision(10);

			myfile2 << "Concepte: " << aTitol << std::endl << std::endl;

			myfile2 << "Time: " << aTime << " Iteration: " << aIteration << std::endl << std::endl;

			for (int node = 0; node<static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
			{
				ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;

				const double hSaturationPerNode = itNode->FastGetSolutionStepValue(SATURATION_NON);
				myfile2 << node << " " << hSaturationPerNode << std::endl;
			}

			myfile2 << std::endl << std::endl;

			for (int node = 0; node<static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
			{
				ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;

				const double zPressurePerNode = itNode->FastGetSolutionStepValue(PRESSURE_WET);
				myfile2 << node << " " << zPressurePerNode << std::endl;
			}

			myfile2 << std::endl << std::endl;

			for (int node = 0; node<static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
			{
				ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;

				const double nonEnthalpyPerNode = itNode->FastGetSolutionStepValue(ENTHALPY_NON_NODE);
				myfile2 << node << " " << nonEnthalpyPerNode << std::endl;
			}

			myfile2 << std::endl << std::endl;

			for (int node = 0; node<static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
			{
				ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;

				const double nonTemperatureNode = itNode->FastGetSolutionStepValue(TEMPERATURE_NODE);
				myfile2 << node << " " << nonTemperatureNode << std::endl;
			}

			myfile2 << std::endl << std::endl;


			myfile2 << std::endl << std::endl;
			myfile2.close();
		}

		if (false)
		{
			if (this->mStep%printInterval == 0)
			{
				std::ofstream myfileMovie_h;
				std::ofstream myfileMovie_Pz;
				std::ofstream myfileMovie_Tw;
				myfileMovie_h.open("h_times.txt", std::ios::app);
				myfileMovie_Pz.open("Pz_times.txt", std::ios::app);
				myfileMovie_h.precision(10);
				myfileMovie_Pz.precision(10);

				myfileMovie_h << "Start" << std::endl;
				myfileMovie_Pz << "Start" << std::endl;

				if (this->mIsNonIsoThermal != 1)
				{
					myfileMovie_Tw.open("Tw_times.txt", std::ios::app);
					myfileMovie_Tw.precision(10);

					myfileMovie_Tw << "Start" << std::endl;
				}

				for (int node = 0; node < static_cast<int>(BaseType::GetModelPart().Nodes().size()); node++)
				{
					ModelPart::NodesContainerType::iterator itNode = BaseType::GetModelPart().NodesBegin() + node;

					const double x = itNode->X();
					const double y = itNode->Y();
					const double heightCO2 = itNode->FastGetSolutionStepValue(SATURATION_NON);
					const double Pz = (itNode->FastGetSolutionStepValue(PRESSURE_WET)) / 1e5;
					myfileMovie_h << x << "\t" << y << "\t" << heightCO2 << std::endl;
					myfileMovie_Pz << x << "\t" << y << "\t" << Pz << std::endl;

					if (this->mIsNonIsoThermal != 1)
					{
						const double Tw = itNode->FastGetSolutionStepValue(ENTHALPY_NON_NODE);
						myfileMovie_Tw << x << "\t" << y << "\t" << Tw << std::endl;
					}

				}

				myfileMovie_h << "End" << std::endl << std::endl;
				myfileMovie_Pz << "End" << std::endl << std::endl;
				if (this->mIsNonIsoThermal != 1)
					myfileMovie_Tw << "End" << std::endl << std::endl;

			}
		}

	}

    template< class TVariableType >
    void inicializateSV(const Kratos::Variable<TVariableType> aVariable_Current,
                             const Kratos::Variable<TVariableType> aVariable_Old)
    {
        KRATOS_TRY
 
        for (ModelPart::NodeIterator node = BaseType::GetModelPart().NodesBegin();
        node != BaseType::GetModelPart().NodesEnd(); ++node)
        {
            const double pressInitialGuess = (node)->FastGetSolutionStepValue(aVariable_Current);
            (node)->FastGetSolutionStepValue(aVariable_Old,1) = pressInitialGuess;//0.0;
            (node)->FastGetSolutionStepValue(aVariable_Old) = pressInitialGuess;//0.0;
            (node)->SetValue(aVariable_Old, pressInitialGuess /*0.0*/);
        }
        
        KRATOS_CATCH("");
    }
    
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
        this->mSumPhasesEqStrategy_Pw->SetEchoLevel(Level);
        this->mNonPhaseEqStrategy_Sn->SetEchoLevel(Level);
        if(this->mIsNonIsoThermal !=1)
        {
            this->mNonWetPhaseEnergyEqStrategy_Hn->SetEchoLevel(Level);
            this->mWetPhaseEnergyEqStrategy_T->SetEchoLevel(Level);
        }
  
    }
    
    virtual void Clear()
    {
        KRATOS_WATCH("FractionalIterativeStrategy Clear Function called");  
        this->mSumPhasesEqStrategy_Pw->Clear();
        this->mNonPhaseEqStrategy_Sn->Clear();
        if(this->mIsNonIsoThermal !=1)
        {
            this->mNonWetPhaseEnergyEqStrategy_Hn->Clear();
            this->mWetPhaseEnergyEqStrategy_T->Clear();
        }
    }
   
    virtual int Check()
    {
        KRATOS_TRY

        
        //veryfying that the model part has all the variables needed
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE_WET) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----PRESSURE_WET or PRESSURE_NON---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(SATURATION_NON) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----SATURATION_WET or SATURATION_NON---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE_WET_OLD_ITER) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----FIRST_SV_OLD_IT---- variable!!!!!! ERROR", "");
        if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(SATURATION_NON_OLD_ITER) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----SECOND_SV_OLD_IT---- variable!!!!!! ERROR", "");
        
        
        //check that the domain size is correctly prescribed
        if (this->mDomainSize != mSolverStrategyConfiguration.GetDomainSize())
            KRATOS_THROW_ERROR(std::logic_error, "domain size not coinciding", "")

            //verify buffer size
            if (BaseType::GetModelPart().GetBufferSize() < mTimeOrder + 1)
                KRATOS_THROW_ERROR(std::logic_error, "insufficient buffer size. Buffer size should be >= time_order+1", "");

        //check that, in the 2D case, the xy plane is used.
        if (this->mDomainSize == 2)
        {
            double zmin = BaseType::GetModelPart().NodesBegin()->Z();
            double zmax = zmin;
            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                if (i->Z() < zmin) zmin = i->Z();
                else if (i->Z() > zmax) zmax = i->Z();
            }
            if (fabs(zmax - zmin) > 1e-20)
                KRATOS_THROW_ERROR(std::logic_error, "2D model is not in the XY plane!", "")
        }
            
            const char ElementName1[] = "MultiphaseFEM1D";
            Element const& ref_el1 = KratosComponents<Element>::Get(ElementName1);
            const char ElementName2[] = "MultiphaseFEM2D";
            Element const& ref_el2 = KratosComponents<Element>::Get(ElementName2);
            const char ElementName3[] = "MultiphaseIsoparametric2D";
            Element const& ref_el3 = KratosComponents<Element>::Get(ElementName3);
            const char ElementName4[] = "MultiphaseIsoparametric3D";
            Element const& ref_el4 = KratosComponents<Element>::Get(ElementName4);
            const char ElementName5[] = "MultiphaseIsoparametric2D4N";
            Element const& ref_el5 = KratosComponents<Element>::Get(ElementName5);
            const char ElementName6[] = "MultiphaseIsoparametric3D8N";
            Element const& ref_el6 = KratosComponents<Element>::Get(ElementName6);
            for (ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin();
                    it != BaseType::GetModelPart().ElementsEnd(); ++it)
            {
                if (it->Id() < 1)
                    KRATOS_THROW_ERROR(std::logic_error, "Element Id can not be lesser than 1 (0 is not allowed as Id)", "");
                if (typeid (ref_el1) != typeid (*it) && typeid (ref_el2) != typeid (*it) && 
                    typeid (ref_el3) != typeid (*it) && typeid (ref_el4) != typeid (*it) &&
                    typeid (ref_el5) != typeid (*it) && typeid (ref_el6) != typeid (*it))
                {
                    std::cout << "wrong element found --> " << it->Id() << std::endl;
                    KRATOS_THROW_ERROR(std::logic_error, "Fractional step strategy requires Multiphase2D element for the 2D case", "");
                }
                //Chequeo area del elemento, en 1D es 0
                it->Check(BaseType::GetModelPart().GetProcessInfo());
            }

        return 0;
        
        
        KRATOS_CATCH("")
    }
    
    void InitializeFractionalStep(const int aStep, const int aTimeOrder)
    {
        KRATOS_TRY;
        
        KRATOS_CATCH("");
    }
    
    void assignInitialStepValues() 
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }
    
    void predictSV(int step, int prediction_order)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    //////////////////////////////////////////// 
     
    //virtual double GetStageResidualNorm(unsigned int step)
    //{
    //  
    //    if (step <= 3)
    //        return this->mPressureStrategy->GetResidualNorm();
    //    if (step == 4)
    //        return this->mConcentrationStrategy->GetResidualNorm();
    //    else
    //     
    //    
    //    return 0.0;
    //}

    //////////////////////////////////////////// 
        /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // This part needed a specialization from element in order to execute a method from my inherit element
      
    /*
    void calculateDensityDomain()
    {
        //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            itElem->calculateDensity();
        }
    }
    
    void CalculateFlowBalancesDomain()
    {
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            itElem->CalculateFlowBalances(CurrentProcessInfo);
        }
    }
        
    void CalculateDarcyFlowPressureDomain()
    {
        ProcessInfo& CurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        
        //#pragma omp parallel for firstprivate(rhs,lhs)
        for(int elem=0; elem<static_cast<int>(BaseType::GetModelPart().Elements().size()); elem++)
        {
            ModelPart::ElementsContainerType::iterator itElem = BaseType::GetModelPart().ElementsBegin()+elem;
            itElem->CalculateDarcyFlowPressure(CurrentProcessInfo);
        }
    }
 */   
    
protected:
    
    typename BaseType::Pointer mSumPhasesEqStrategy_Pw; 
    typename BaseType::Pointer mNonPhaseEqStrategy_Sn; 
    typename BaseType::Pointer mNonWetPhaseEnergyEqStrategy_Hn; 
    typename BaseType::Pointer mWetPhaseEnergyEqStrategy_T; 
    
    double mTolerance_Pw; 
    double mTolerance_Sn;
    double mTolerance_Hn;
    double mTolerance_T;
    int mMaxIter_Pw; 
    int mMaxIter_Sn; 
    int mMaxIter_Hn; 
    int mMaxIter_T; 
    unsigned int mTimeOrder;
    unsigned int mPredictorOrder; 
    bool mPredictorCorrector; 
    bool mReformDofAtEachIteration;
    int mEchoLevel;
    
private:    

    unsigned int mStep;
    unsigned int mDomainSize; 
    
    bool mIntegrationByGaussPoints;
    
    std::string mRelativePermeabilityModel;
    
    unsigned int mIsCapillarityNeglected;
    std::string mCapillarityPressureModel;
    
    unsigned int mIsNonPhaseNonConstant;
    unsigned int mIsWetPhaseNonConstant;
    unsigned int mIsNonIsoThermal;
	unsigned int mIsWetIsoThermal;
    std::string mEOS;
    SpanWagnerEOS* mSpanWagnerEOS;
	unsigned int mNoIterative;

	Vector mBalanceTotalsToPrint;

    SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& mSolverStrategyConfiguration;

    SequentialMultiphaseFEMStrategy(const SequentialMultiphaseFEMStrategy& Other);

}; /* Class SequentialMultiphaseFEMStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_SEQUENTIAL_MULTIPHASE_FEM_STRATEGY  defined */

