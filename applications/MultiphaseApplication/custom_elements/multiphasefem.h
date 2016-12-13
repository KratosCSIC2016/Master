/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-10-09 10:34:00 $
//   Revision:            $Revision: 0.1 $
//
//


#if !defined(KRATOS_MULTIPHASEFEM_H_INCLUDED )
#define  KRATOS_MULTIPHASEFEM_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"

// Application includes
#include "porous_media_application_variables.h"
#include "custom_utilities/spanwagnereos.h"

namespace Kratos
{

    ///@addtogroup MultiphaseApplication
    ///@{

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

    /// A stabilized element for the incompressible multiphase equations Pn-Sw.
    /**
     */
    template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
    class MultiphaseFEM : public Element
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Multiphase
        KRATOS_CLASS_POINTER_DEFINITION(MultiphaseFEM);

        /// Node type (default is: Node<3>)
        typedef Node <3> NodeType;

        /// Geometry type (using with given NodeType)
        typedef Geometry<NodeType> GeometryType;

        /// Definition of nodes container type, redefined from GeometryType
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

        /// Vector type for local contributions to the linear system
        typedef Vector VectorType;

        /// Matrix type for local contributions to the linear system
        typedef Matrix MatrixType;

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef std::vector<std::size_t> EquationIdVectorType;

        typedef std::vector< Dof<double>::Pointer > DofsVectorType;

        typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

        //Constructors.

        /// Default constuctor.
        /**
         * @param NewId Index number of the new element (optional)
         */
        MultiphaseFEM(IndexType NewId = 0) :
            Element(NewId)
        {}

        /// Constructor using an array of nodes.
        /**
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        MultiphaseFEM(IndexType NewId, const NodesArrayType& ThisNodes) :
            Element(NewId, ThisNodes)
        {}

        /// Constructor using a geometry object.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        MultiphaseFEM(IndexType NewId, GeometryType::Pointer pGeometry) :
            Element(NewId, pGeometry)
        {}

        /// Constuctor using geometry and properties.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        MultiphaseFEM(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
            Element(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~MultiphaseFEM()
        {}

        /// Create a new element of this type
        /**
         * Returns a pointer to a new Multiphase element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            return Element::Pointer(new MultiphaseFEM(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

        virtual void Initialize();

        /// Initializes the element and all geometric information required for the problem.
        virtual void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo);

        /// Initialize viscosity, adding Smagorinsky eddy viscosity if it is active.
        virtual void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo);

        /// Calculate the element's local contribution to the system for the current step.
        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY;
            KRATOS_THROW_ERROR(std::logic_error,"MultiphaseFEM::CalculateLeftHandSide not implemented","");
            KRATOS_CATCH("");
        }

        // The following methods have different implementations depending on TDim
        /// Provides the global indices for each one of this element's local rows
        /**
         * this determines the elemental equation ID vector for all elemental
         * DOFs
         * @param rResult A vector containing the global Id of each row
         * @param rCurrentProcessInfo the current process info object (unused)
         */
        virtual void EquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo);

        /// Returns a list of the element's Dofs
        /**
         * @param ElementalDofList the list of DOFs
         * @param rCurrentProcessInfo the current process info instance
         */
        virtual void GetDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo);

        /// Turn back information as a string.
        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "MultiphaseFEM #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "MultiphaseFEM" << TDim << "D";
        }

    private:
        
        //Thermodynamics EOS     
        SpanWagnerEOS* mSpanWagnerEOS;
        
        void InterpolateElemFields(ProcessInfo &rCurrentProcessInfo);
        
        void CalculatePhasesDarcyFlow(ProcessInfo &rCurrentProcessInfo);
        
        //Lumped
        void CalculateNiNj(MatrixType &rNiNj);

        //Consistent
        void CalculateNiNj(MatrixType &rNiNj, const VectorType& aN);

        void CalculateBiBj(MatrixType& rBiBj, 
                            const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > aDN_DX);
        
        void CalculateSumPhasesEquationPw(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo);

        void CalculateNonWettingEquationh(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);		
		
        void CalculateNonWettingEquationSn(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);
        
        void CalculateNonWettingEnergyEquationHn(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo);

        void CalculateWettingEnergyEquationT(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);
        
        void PwEquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo);

        void SnEquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo);

        void HnEquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo);
                
        void TEquationIdVector(EquationIdVectorType& rResult,
                      ProcessInfo& rCurrentProcessInfo);
        
        void GetPwEquationDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo);

        void GetSnEquationDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo);

        void GetHnEquationDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo);

        void GetTEquationDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo);
        
        void GetPwValues(Vector& rValues,
                               const int Step = 0);

        void GetSnValues(Vector& rValues,
                               const int Step = 0);

        void GetHnValues(Vector& rValues,
                       const int Step = 0);

        void GetTValues(Vector& rValues,
                               const int Step = 0);
        
        double ElementSize(/*ShapeFunctionDerivativesType& rDN_DX*/);
        
                virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable, 
        std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == CAPILLARITY_PRESSURE_ELEM)
            {
                    // Set output vector (for a single integration point)
                    rValues.resize(1, false);
                    rValues[0]=this->GetValue(CAPILLARITY_PRESSURE_ELEM);
            }
            else if (rVariable == DENSITY_NON_ELEM)
            {
                    // Set output vector (for a single integration point)
                    rValues.resize(1, false);
                    rValues[0]=this->GetValue(DENSITY_NON_ELEM);
            }
            else if (rVariable == VISCOSITY_NON_ELEM)
            {
                    // Set output vector (for a single integration point)
                    rValues.resize(1, false);
                    rValues[0]=this->GetValue(VISCOSITY_NON_ELEM);
            }
            else if (rVariable == COMPRESSIBLITY_NON_ELEM)
            {
                    // Set output vector (for a single integration point)
                    rValues.resize(1, false);
                    rValues[0]=this->GetValue(COMPRESSIBLITY_NON_ELEM);
            }
            else if (rVariable == DERIV_MUN_PN_ELEM)
            {
                    // Set output vector (for a single integration point)
                    rValues.resize(1, false);
                    rValues[0]=this->GetValue(DERIV_MUN_PN_ELEM);
            }
            else if (rVariable == TEMPERATURE_ELEM)
            {
                    // Set output vector (for a single integration point)
                    rValues.resize(1, false);
                    rValues[0]=this->GetValue(TEMPERATURE_ELEM);
            }
            else // Default behaviour (returns elemental data)
            {
                    rValues.resize(1, false);
            }
        }
        
        virtual void GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
        std::vector<array_1d<double,3> >& rValues,const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == DARCYFLOW_NON)
            {
                // Set output vector (for a single integration point)
                rValues.resize(1);
                array_1d<double, 3 > & rDarcyFlowNon_Elem = rValues[0];
                rDarcyFlowNon_Elem[0] = this->GetValue(DARCYFLOW_NON_X);
                rDarcyFlowNon_Elem[1] = this->GetValue(DARCYFLOW_NON_Y);
                rDarcyFlowNon_Elem[2] = 0.0;//this->GetValue(DARCYFLOW_NON_Z);
            }
            else if (rVariable == DARCYFLOW_WET)
            {
                rValues.resize(1);
                array_1d<double, 3 > & rDarcyFlowWet_Elem = rValues[0];
                rDarcyFlowWet_Elem[0] = this->GetValue(DARCYFLOW_WET_X);
                rDarcyFlowWet_Elem[1] = this->GetValue(DARCYFLOW_WET_Y);
                rDarcyFlowWet_Elem[2] = 0.0;//this->GetValue(DARCYFLOW_NON_Z);
            }
            else // Default behaviour (returns elemental data)
            {
                rValues.resize(1);
            }
        }
        
        ///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        }

        /// Assignment operator.
        MultiphaseFEM & operator=(MultiphaseFEM const& rOther);

        /// Copy constructor.
        MultiphaseFEM(MultiphaseFEM const& rOther);

        ///@}

    }; // Class MultiphaseFEM

    /// input stream function
    template< unsigned int TDim,
              unsigned int TNumNodes  >
    inline std::istream& operator >>(std::istream& rIStream,
                                     MultiphaseFEM<TDim,TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim,
          unsigned int TNumNodes >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const MultiphaseFEM<TDim,TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

    ///@} // PorousMedia Application group

} // namespace Kratos.

#endif // KRATOS_MULTIPHASEFEM_H_INCLUDED  defined
