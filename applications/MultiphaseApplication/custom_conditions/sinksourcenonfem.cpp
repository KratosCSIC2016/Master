// Project includes 
#include "includes/define.h"
#include "custom_conditions/sinksourcenonfem.h"
#include "multiphase_application.h"
#include "utilities/math_utils.h"


namespace Kratos
{

	typedef GeometryData::KratosGeometryType KratosGeometryType;

	//************************************************************************************
	//************************************************************************************
	SinkSourceNonFEM::SinkSourceNonFEM(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
		mSpanWagnerEOS = new SpanWagnerEOS();
		//if( EOS == "SpanWagner")
		this->mSpanWagnerEOS->FillTables();
	}

	//************************************************************************************
	//************************************************************************************
	SinkSourceNonFEM::SinkSourceNonFEM(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{

		mSpanWagnerEOS = new SpanWagnerEOS();
		//if( EOS == "SpanWagner")
		this->mSpanWagnerEOS->FillTables();
	}
	Condition::Pointer SinkSourceNonFEM::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new SinkSourceNonFEM(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	SinkSourceNonFEM::~SinkSourceNonFEM()
	{
		delete this->mSpanWagnerEOS;
	}

	//************************************************************************************
	//************************************************************************************
	void SinkSourceNonFEM::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
			switch (rCurrentProcessInfo[FRACTIONAL_STEP])
			{
			case 1:
			{
				this->CalculateRightHandSideSinkSourcePressureWetPw(rRightHandSideVector, rCurrentProcessInfo);
				break;
			}
			case 2:
			{
				this->CalculateRightHandSideSinkSourceSaturationNonSn(rRightHandSideVector, rCurrentProcessInfo);
				break;
			}
			case 3:
			{
				this->CalculateRightHandSideSinkSourceEnthalpyNonHn(rRightHandSideVector, rCurrentProcessInfo);
				break;
			}
			case 4:
			{
				if (rRightHandSideVector.size() != 0)
					rRightHandSideVector.resize(0, false);
				break;
			}
			case 5:
			{
				if (rRightHandSideVector.size() != 0)
					rRightHandSideVector.resize(0, false);
				break;
			}
			default:
			{
				KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
			}
			}
		KRATOS_CATCH("")
	}

	void SinkSourceNonFEM::CalculateRightHandSideSinkSourcePressureWetPw(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			bool CalculateStiffnessMatrixFlag = false;
		MatrixType temp = Matrix();

		this->CalculateAllPressureWetPw(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);

		KRATOS_CATCH("")

	}

	void SinkSourceNonFEM::CalculateRightHandSideSinkSourceSaturationNonSn(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			bool CalculateStiffnessMatrixFlag = false;
		MatrixType temp = Matrix();

		this->CalculateAllSaturationNonSn(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);

		KRATOS_CATCH("")

	}

	void SinkSourceNonFEM::CalculateRightHandSideSinkSourceEnthalpyNonHn(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			bool CalculateStiffnessMatrixFlag = false;
		MatrixType temp = Matrix();

		this->CalculateAllEnthalpyNonHn(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);

		KRATOS_CATCH("")

	}

	void SinkSourceNonFEM::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			switch (rCurrentProcessInfo[FRACTIONAL_STEP])
			{
			case 1:
			{
				this->CalculateLocalSystemSinkSourcePressureWetPw(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
				//rRightHandSideVector[0]=0.0;
				break;
			}
			case 2:
			{
				this->CalculateLocalSystemSinkSourceSaturationNonSn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
				//rRightHandSideVector[0]=0.0;
				break;
			}
			case 3:
			{
				this->CalculateLocalSystemSinkSourceEnthalpyNonHn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
				break;
			}
			case 4:
			{
				if (rLeftHandSideMatrix.size1() != 0)
				{
					rLeftHandSideMatrix.resize(0, 0, false);
					rRightHandSideVector.resize(0, false);
				}
				break;
			}
			case 5:
			{
				if (rLeftHandSideMatrix.size1() != 0)
				{
					rLeftHandSideMatrix.resize(0, 0, false);
					rRightHandSideVector.resize(0, false);
				}
				break;
			}
			default:
			{
				KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
			}
			}
		KRATOS_CATCH("")

	}

	void SinkSourceNonFEM::CalculateLocalSystemSinkSourcePressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		KRATOS_TRY

			bool CalculateStiffnessMatrixFlag = true;

		this->CalculateAllPressureWetPw(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);

		KRATOS_CATCH("")

	}

	void SinkSourceNonFEM::CalculateLocalSystemSinkSourceSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		KRATOS_TRY

			bool CalculateStiffnessMatrixFlag = true;

		this->CalculateAllSaturationNonSn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);

		KRATOS_CATCH("")

	}

	void SinkSourceNonFEM::CalculateLocalSystemSinkSourceEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

		KRATOS_TRY

			bool CalculateStiffnessMatrixFlag = true;

		this->CalculateAllEnthalpyNonHn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);

		KRATOS_CATCH("")

	}
	void SinkSourceNonFEM::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		switch (rCurrentProcessInfo[FRACTIONAL_STEP])
		{
		case 1:
		{
			this->EquationIdVectorPressureWetPw(rResult, rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->EquationIdVectorSaturationNonSn(rResult, rCurrentProcessInfo);
			break;
		}
		case 3:
		{
			this->EquationIdVectorEnthalpyNonHn(rResult, rCurrentProcessInfo);
			break;
		}
		case 4:
		{
			if (rResult.size() != 0)
				rResult.resize(0, false);
			break;
		}
		case 5:
		{
			if (rResult.size() != 0)
				rResult.resize(0, false);
			break;
		}
		default:
		{
			KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

	}

	void SinkSourceNonFEM::EquationIdVectorPressureWetPw(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i = 0; i<number_of_nodes; i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(PRESSURE_WET).EquationId());
		}
	}

	void SinkSourceNonFEM::EquationIdVectorSaturationNonSn(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i = 0; i<number_of_nodes; i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(SATURATION_NON).EquationId());
		}
	}

	void SinkSourceNonFEM::EquationIdVectorEnthalpyNonHn(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i = 0; i<number_of_nodes; i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(ENTHALPY_NON_NODE).EquationId());
		}
	}

	//************************************************************************************
	//************************************************************************************
	void SinkSourceNonFEM::GetDofList(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		switch (rCurrentProcessInfo[FRACTIONAL_STEP])
		{
		case 1:
		{
			this->GetDofListPressureWetPw(rConditionalDofList, rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->GetDofListSaturationNonSn(rConditionalDofList, rCurrentProcessInfo);
			break;
		}
		case 3:
		{
			this->GetDofListEnthalpyNonHn(rConditionalDofList, rCurrentProcessInfo);
			break;
		}
		case 4:
		{
			if (rConditionalDofList.size() != 0)
				rConditionalDofList.resize(0);
			break;
		}
		case 5:
		{
			if (rConditionalDofList.size() != 0)
				rConditionalDofList.resize(0);
			break;
		}
		default:
		{
			KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

	}

	void SinkSourceNonFEM::GetDofListPressureWetPw(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int dim = 1;
		rConditionalDofList.resize(GetGeometry().size()*dim);
		unsigned int index;
		for (unsigned int i = 0; i<GetGeometry().size(); i++)
		{

			index = i*dim;
			rConditionalDofList[index] = (GetGeometry()[i].pGetDof(PRESSURE_WET));
		}
	}

	void SinkSourceNonFEM::GetDofListSaturationNonSn(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int dim = 1;
		rConditionalDofList.resize(GetGeometry().size()*dim);
		unsigned int index;
		for (unsigned int i = 0; i<GetGeometry().size(); i++)
		{

			index = i*dim;
			rConditionalDofList[index] = (GetGeometry()[i].pGetDof(SATURATION_NON));
		}
	}

	void SinkSourceNonFEM::GetDofListEnthalpyNonHn(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int dim = 1;
		rConditionalDofList.resize(GetGeometry().size()*dim);
		unsigned int index;
		for (unsigned int i = 0; i<GetGeometry().size(); i++)
		{

			index = i*dim;
			rConditionalDofList[index] = (GetGeometry()[i].pGetDof(ENTHALPY_NON_NODE));
		}
	}

	//************************************************************************************
	//************************************************************************************

	void SinkSourceNonFEM::CalculateAllPressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag)
	{

		KRATOS_TRY

			//Set Matrix & vector size
			unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize = number_of_nodes;

		double sinkSource = GetProperties()[SINK_SOURCE_NON];

		if (aCalculateStiffnessMatrixFlag == true)
		{
			if (rLeftHandSideMatrix.size1() != MatSize)
				rLeftHandSideMatrix.resize(MatSize, MatSize, false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);
		}

		if (rRightHandSideVector.size() != MatSize)
			rRightHandSideVector.resize(MatSize, false);


		KratosGeometryType typeOfElement = GetGeometry().GetGeometryType();

		switch (typeOfElement) {
		case  GeometryData::Kratos_Point2D:
			this->CalculateLocalSystemPointPressureWetPw(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
			break;
		case GeometryData::Kratos_Line2D2:
			this->CalculateLocalSystemLinePressureWetPw(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
			break;
		case GeometryData::Kratos_Triangle2D3:
			this->CalculateLocalSystemTrianglePressureWetPw(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
			break;
		default:
			KRATOS_THROW_ERROR(std::logic_error, "This element is not yet implement in order to solve source term!!!!! ERROR", "");
			break;
		}

		KRATOS_CATCH("");


	}

	void SinkSourceNonFEM::CalculateAllSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag)
	{

		KRATOS_TRY

		//Set Matrix & vector size
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize = number_of_nodes;

		double sinkSource = GetProperties()[SINK_SOURCE_NON];

		if (aCalculateStiffnessMatrixFlag == true)
		{
			if (rLeftHandSideMatrix.size1() != MatSize)
				rLeftHandSideMatrix.resize(MatSize, MatSize, false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);
		}

		if (rRightHandSideVector.size() != MatSize)
			rRightHandSideVector.resize(MatSize, false);


		KratosGeometryType typeOfElement = GetGeometry().GetGeometryType();

		switch (typeOfElement) {
		case  GeometryData::Kratos_Point2D:
			this->CalculateLocalSystemPointSaturationNonSn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
			break;
		case GeometryData::Kratos_Line2D2:
			this->CalculateLocalSystemLineSaturationNonSn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
			break;
		case GeometryData::Kratos_Triangle2D3:
			this->CalculateLocalSystemTriangleSaturationNonSn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
			break;
		default:
			KRATOS_THROW_ERROR(std::logic_error, "This element is not yet implement in order to solve source term!!!!! ERROR", "");
			break;
		}

		KRATOS_CATCH("");


	}

	void SinkSourceNonFEM::CalculateAllEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag)
	{

		KRATOS_TRY

			//Set Matrix & vector size
			unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize = number_of_nodes;

		double sinkSource = GetProperties()[SINK_SOURCE_NON];
		const double currentHeadPressureWet = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_WET);
		const double currentCapillarPressure = GetGeometry()[0].FastGetSolutionStepValue(CAPILLARITY_PRESSURE_NODE);
		const double currentHeadPressureNon = currentHeadPressureWet + currentCapillarPressure;
		const double injTemperature = GetProperties()[PRESCRIBED_VALUE];
		double injEnthalpy = 0.0;//-92956.4313
		this->mSpanWagnerEOS->interpolate_calculate_variable_PT("enthalpy", injEnthalpy, currentHeadPressureNon, injTemperature);
		double injDensity = 0.0;
		this->mSpanWagnerEOS->interpolate_calculate_variable_PT("density", injDensity, currentHeadPressureNon, injTemperature);

		if (aCalculateStiffnessMatrixFlag == true)
		{
			if (rLeftHandSideMatrix.size1() != MatSize)
				rLeftHandSideMatrix.resize(MatSize, MatSize, false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);
		}

		if (rRightHandSideVector.size() != MatSize)
			rRightHandSideVector.resize(MatSize, false);


		KratosGeometryType typeOfElement = GetGeometry().GetGeometryType();

		switch (typeOfElement) {
		case  GeometryData::Kratos_Point2D:
			this->CalculateLocalSystemTriangleEnthalpyNonHn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource, injEnthalpy, injDensity);
			break;
		case GeometryData::Kratos_Line2D2:
			this->CalculateLocalSystemLineEnthalpyNonHn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource, injEnthalpy, injDensity);
			break;
		case GeometryData::Kratos_Triangle2D3:
			this->CalculateLocalSystemTriangleEnthalpyNonHn(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource, injEnthalpy, injDensity);
			break;
		default:
			KRATOS_THROW_ERROR(std::logic_error, "This element is not yet implement in order to solve source term!!!!! ERROR", "");
			break;
		}

		KRATOS_CATCH("");


	}

	void SinkSourceNonFEM::CalculateLocalSystemPointPressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
	{
		rRightHandSideVector[0] = aSinkSource;
	}

	void SinkSourceNonFEM::CalculateLocalSystemLinePressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
	{
		//calculate lenght
		double xlenght = GetGeometry()[1].X() - GetGeometry()[0].X();
		double ylenght = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double lenght = xlenght*xlenght + ylenght*ylenght;
		lenght = sqrt(lenght);

		rRightHandSideVector[0] = (aSinkSource*lenght) / 2;
		rRightHandSideVector[1] = (aSinkSource*lenght) / 2;

	}

	void SinkSourceNonFEM::CalculateLocalSystemTrianglePressureWetPw(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
	{
		//calculate area
		double Area = GetGeometry().Area();

		rRightHandSideVector[0] = (aSinkSource*Area) / 3;
		rRightHandSideVector[1] = (aSinkSource*Area) / 3;
		rRightHandSideVector[2] = (aSinkSource*Area) / 3;


	}

	void SinkSourceNonFEM::CalculateLocalSystemPointSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
	{
		rRightHandSideVector[0] = aSinkSource;
	}

	void SinkSourceNonFEM::CalculateLocalSystemLineSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
	{
		//calculate lenght
		double xlenght = GetGeometry()[1].X() - GetGeometry()[0].X();
		double ylenght = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double lenght = xlenght*xlenght + ylenght*ylenght;
		lenght = sqrt(lenght);

		rRightHandSideVector[0] = (aSinkSource*lenght) / 2;
		rRightHandSideVector[1] = (aSinkSource*lenght) / 2;

	}

	void SinkSourceNonFEM::CalculateLocalSystemTriangleSaturationNonSn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
	{
		//calculate area
		double Area = GetGeometry().Area();

		rRightHandSideVector[0] = (aSinkSource*Area) / 3;
		rRightHandSideVector[1] = (aSinkSource*Area) / 3;
		rRightHandSideVector[2] = (aSinkSource*Area) / 3;


	}


	void SinkSourceNonFEM::CalculateLocalSystemPointEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, const double aSinkSource, const double aNonEnthalpy, const double aNonDensity)
	{
		const double currentHeadEnthalpy = GetGeometry()[0].FastGetSolutionStepValue(ENTHALPY_NON_NODE);
		const double currentHeadDensity = GetGeometry()[0].FastGetSolutionStepValue(DENSITY_NON_NODE);

		rLeftHandSideMatrix(0, 0) = aSinkSource;

		rRightHandSideVector[0] = (aSinkSource / currentHeadDensity)*((aNonEnthalpy*aNonDensity) - (currentHeadEnthalpy*currentHeadDensity));
	}

	void SinkSourceNonFEM::CalculateLocalSystemLineEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource, const double aNonEnthalpy, const double aNonDensity)
	{
		//calculate lenght
		double xlenght = GetGeometry()[1].X() - GetGeometry()[0].X();
		double ylenght = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double lenght = xlenght*xlenght + ylenght*ylenght;
		lenght = sqrt(lenght);

		rRightHandSideVector[0] = (aSinkSource*aNonEnthalpy*lenght) / 2;
		rRightHandSideVector[1] = (aSinkSource*aNonEnthalpy*lenght) / 2;

	}

	void SinkSourceNonFEM::CalculateLocalSystemTriangleEnthalpyNonHn(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource, const double aNonEnthalpy, const double aNonDensity)
	{
		//calculate area
		double Area = GetGeometry().Area();

		rRightHandSideVector[0] = (aSinkSource*aNonEnthalpy*Area) / 3;
		rRightHandSideVector[1] = (aSinkSource*aNonEnthalpy*Area) / 3;
		rRightHandSideVector[2] = (aSinkSource*aNonEnthalpy*Area) / 3;


	}







	//************************************************************************************
	//************************************************************************************
} // Namespace Kratos
