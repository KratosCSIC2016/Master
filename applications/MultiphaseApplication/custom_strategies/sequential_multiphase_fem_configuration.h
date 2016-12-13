#if !defined(KRATOS_SEQUENTIAL_MULTIPHASE_FEM_CONFIGURATION )
#define  KRATOS_SEQUENTIAL_MULTIPHASE_FEM_CONFIGURATION

//External includes
#include "boost/smart_ptr.hpp"

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
//Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
//Builder & solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_slip.h"
//Scheme
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class SequentialMultiphaseFEMConfiguration : public SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::DofsArrayType DofsArrayType; 
    typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
    typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
    
    /*
     * From incompresible_fluid_application/ strategies/ custom_strategies/ fractional_iterative_configuration.h 
     */
    
    SequentialMultiphaseFEMConfiguration(ModelPart& model_part,
                                typename TLinearSolver::Pointer pSystemLinearSolver_Pw,
                                typename TLinearSolver::Pointer pSystemLinearSolver_Sn,
                                typename TLinearSolver::Pointer pSystemLinearSolver_Hn,
                                typename TLinearSolver::Pointer pSystemLinearSolver_T,
                                unsigned int mDomainSize
                               )
        : SolverStrategyConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, mDomainSize)
    {

        const bool CalculateReactions = false;
        const bool CalculateNormDxFlag = true;
        const bool ReformDofAtEachIteration = true;
        
        //SumPhases_Pn      
        this->mSumPhasesEqScheme_Pw = typename SchemeType::Pointer
            (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
        this->mSumPhasesEqSystemBuild_Pw = BuilderSolverTypePointer
            (new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver > (pSystemLinearSolver_Pw) ); 
           //new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pSystemLinearSolver_Pn, PRESSURE_NON));
        this->mSumPhasesEqStrategy_Pw = typename BaseType::Pointer
            (new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > (model_part,mSumPhasesEqScheme_Pw,pSystemLinearSolver_Pw,mSumPhasesEqSystemBuild_Pw,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        this->mSumPhasesEqStrategy_Pw->SetEchoLevel(2);
    
        //WetPhase_Sw
        this->mNonPhaseEqScheme_Sn = typename SchemeType::Pointer
            (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
        this->mNonPhaseEqSystemBuild_Sn = BuilderSolverTypePointer
            (new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver > (pSystemLinearSolver_Sn) ); 
           //new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pSystemLinearSolver_Sn, SATURATION_NON));
        this->mNonPhaseEqStrategy_Sn = typename BaseType::Pointer
            (new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > (model_part,mNonPhaseEqScheme_Sn,pSystemLinearSolver_Sn,mNonPhaseEqSystemBuild_Sn,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        this->mNonPhaseEqStrategy_Sn->SetEchoLevel(2);    
        
        //NonWetPhaseEnergy_Hn
        this->mNonWetPhaseEnergyEqScheme_Hn = typename SchemeType::Pointer
            (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
        this->mNonWetPhaseEnergyEqSystemBuild_Hn = BuilderSolverTypePointer
            (new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver > (pSystemLinearSolver_Hn) ); 
           //new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pSystemLinearSolver_Hn, ENTHALPY_NON_NODE));
        this->mNonWetPhaseEnergyEqStrategy_Hn = typename BaseType::Pointer
            (new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > (model_part,mNonWetPhaseEnergyEqScheme_Hn,pSystemLinearSolver_Hn,mNonWetPhaseEnergyEqSystemBuild_Hn,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        this->mNonWetPhaseEnergyEqStrategy_Hn->SetEchoLevel(2);    
        
        //WetPhaseEnergy_T
        this->mWetPhaseEnergyEqScheme_T = typename SchemeType::Pointer
            (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
        this->mWetPhaseEnergyEqSystemBuild_T = BuilderSolverTypePointer
            (new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver > (pSystemLinearSolver_T) ); 
           //new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pSystemLinearSolver_T, TEMPERATURE_NODE));
        this->mWetPhaseEnergyEqStrategy_T = typename BaseType::Pointer
            (new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > (model_part,mWetPhaseEnergyEqScheme_T,pSystemLinearSolver_T,mWetPhaseEnergyEqSystemBuild_T,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
        this->mWetPhaseEnergyEqStrategy_T->SetEchoLevel(2);    

    }


    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer pGetStrategy(const std::string& strategy_name)
    {
        KRATOS_TRY

        if (strategy_name == std::string("SumPhasesEqStrategy_Pw"))
            return mSumPhasesEqStrategy_Pw;
        else if (strategy_name == std::string("NonPhaseEqStrategy_Sn"))
            return mNonPhaseEqStrategy_Sn;
        else if (strategy_name == std::string("NonWetPhaseEnergyEqStrategy_Hn"))
            return mNonWetPhaseEnergyEqStrategy_Hn;
        else if (strategy_name == std::string("WetPhaseEnergyEqStrategy_T"))
            return mWetPhaseEnergyEqStrategy_T;
        else
            KRATOS_THROW_ERROR(std::invalid_argument, "SequentialMultiphaseFEMConfiguration: Trying to get an inexisting strategy", "");
    
        KRATOS_CATCH("")
    }

    void pSetUpDof(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
	{
	   case 1:
            {
                this->mSumPhasesEqSystemBuild_Pw->SetUpDofSet(mSumPhasesEqScheme_Pw, mSumPhasesEqStrategy_Pw->GetModelPart());
                //DofsArrayType& rDofSetWetting =this->mSumPhasesEqSystemBuild_Pn->GetDofSet();
                break;
            }
	    case 2:
            {
                this->mNonPhaseEqSystemBuild_Sn->SetUpDofSet(mNonPhaseEqScheme_Sn, mNonPhaseEqStrategy_Sn->GetModelPart());
                //DofsArrayType& rDofSetNonWetting =this->mWetPhaseEqSystemBuild_Sw->GetDofSet();
                break;
            }
            case 3:
            {
                this->mNonWetPhaseEnergyEqSystemBuild_Hn->SetUpDofSet(mNonWetPhaseEnergyEqScheme_Hn, mNonWetPhaseEnergyEqStrategy_Hn->GetModelPart());
                //DofsArrayType& rDofSetNonWetting =this->mNonWetPhaseEnergyEqSystemBuild_Hn->GetDofSet();
                break;
            }
            case 4:
            {
                this->mWetPhaseEnergyEqSystemBuild_T->SetUpDofSet(mWetPhaseEnergyEqScheme_T, mWetPhaseEnergyEqStrategy_T->GetModelPart());
                //DofsArrayType& rDofSetNonWetting =this->mWetPhaseEnergyEqSystemBuild_T->GetDofSet();
                break;
            }
	    default:
            {
                    KRATOS_THROW_ERROR(std::logic_error,"SequentialMultiphaseFEMConfiguration: Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
            }
        }
        
        KRATOS_CATCH("")
                
    }

protected:


private:

    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mSumPhasesEqStrategy_Pw;
    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mNonPhaseEqStrategy_Sn;
    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mNonWetPhaseEnergyEqStrategy_Hn;
    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mWetPhaseEnergyEqStrategy_T;

    BuilderSolverTypePointer mSumPhasesEqSystemBuild_Pw; 
    BuilderSolverTypePointer mNonPhaseEqSystemBuild_Sn; 
    BuilderSolverTypePointer mNonWetPhaseEnergyEqSystemBuild_Hn; 
    BuilderSolverTypePointer mWetPhaseEnergyEqSystemBuild_T; 
   
    typename SchemeType::Pointer mSumPhasesEqScheme_Pw;
    typename SchemeType::Pointer mNonPhaseEqScheme_Sn;
    typename SchemeType::Pointer mNonWetPhaseEnergyEqScheme_Hn;
    typename SchemeType::Pointer mWetPhaseEnergyEqScheme_T;

}; /* Class SequentialMultiphaseFEMConfiguration */

} /* namespace Kratos.*/

#endif /* KRATOS_SEQUENTIAL_MULTIPHASE_FEM_CONFIGURATION  defined */
