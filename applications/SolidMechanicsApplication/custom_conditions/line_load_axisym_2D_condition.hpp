//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LINE_LOAD_AXISYM_2D_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_LOAD_AXISYM_2D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/line_load_2D_condition.hpp"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

/// Force Load Condition for 3D and 2D geometries. (base class)

/**
 * Implements a Force Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) LineLoadAxisym2DCondition
    : public LineLoad2DCondition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of LineLoadAxisym2DCondition
    KRATOS_CLASS_POINTER_DEFINITION( LineLoadAxisym2DCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LineLoadAxisym2DCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    LineLoadAxisym2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    LineLoadAxisym2DCondition( LineLoadAxisym2DCondition const& rOther);

    /// Destructor
    virtual ~LineLoadAxisym2DCondition();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
			      NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties ) const;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const;



    //************* COMPUTING  METHODS


    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    LineLoadAxisym2DCondition() {};
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, 
				     const double& rPointNumber);

    /**
     * Calculation and addition of the matrices of the LHS
     */
    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */
    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight);

    /**
     * Calculation of the contidion radius (axisymmetry)
     */
    void CalculateRadius(double & rCurrentRadius,
			 double & rReferenceRadius,
			 const Vector& rN);

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


}; // class LineLoadAxisym2DCondition.

} // namespace Kratos.

#endif // KRATOS_LINE_LOAD_AXISYM_2D_CONDITION_H_INCLUDED defined 
