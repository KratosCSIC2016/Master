#if !defined(KRATOS_SINKSOURCENON_FEM_CONDITION_H_INCLUDED )
#define  KRATOS_SINKSOURCENON_FEM_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "custom_utilities/spanwagnereos.h"

namespace Kratos
{
 class SinkSourceNonFEM : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of PointForce2D
       KRATOS_CLASS_POINTER_DEFINITION(SinkSourceNonFEM);
      
      /// Default constructor. 
	  SinkSourceNonFEM(IndexType NewId, GeometryType::Pointer pGeometry);

	  SinkSourceNonFEM(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~SinkSourceNonFEM();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo);
      
   protected:
 
   private:
       void GetDofListPressureWetPw(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo);
       void GetDofListSaturationNonSn(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo);
       void GetDofListEnthalpyNonHn(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo);

       void EquationIdVectorPressureWetPw(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
       void EquationIdVectorSaturationNonSn(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
       void EquationIdVectorEnthalpyNonHn(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
       
       void CalculateRightHandSideSinkSourcePressureWetPw(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       void CalculateRightHandSideSinkSourceSaturationNonSn(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       void CalculateRightHandSideSinkSourceEnthalpyNonHn(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

       void CalculateLocalSystemSinkSourcePressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       void CalculateLocalSystemSinkSourceSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       void CalculateLocalSystemSinkSourceEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
       
       void CalculateAllPressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag);
       void CalculateAllSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag);
       void CalculateAllEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag);

       void CalculateLocalSystemTrianglePressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       void CalculateLocalSystemLinePressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       void CalculateLocalSystemPointPressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       
       void CalculateLocalSystemTriangleSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       void CalculateLocalSystemLineSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);
       void CalculateLocalSystemPointSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource);

       void CalculateLocalSystemTriangleEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource, const double aNonEnthalpy,const double aNonDensity);
       void CalculateLocalSystemLineEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource, const double aNonEnthalpy,const double aNonDensity);
       void CalculateLocalSystemPointEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource, const double aNonEnthalpy,const double aNonDensity);

       SpanWagnerEOS* mSpanWagnerEOS;
       
       friend class Serializer;

	// A private default constructor necessary for serialization  
        SinkSourceNonFEM() : Condition()
       {
       }
        

  }; // Class SinkSourceNonFEM 

} //namespace kratos 
#endif
