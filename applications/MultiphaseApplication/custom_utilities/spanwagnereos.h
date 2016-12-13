#if !defined(KRATOS_SPANWAGNEREOS_UTILITY_INCLUDED )
#define  KRATOS_SPANWAGNEREOS_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>
#include <utility>
#include <QtXml>

// Project includes 
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"

//Application
#include "porous_media_application_variables.h"

// Agus et al. Span & Wagner et al.
#define Pcriticalpoint 7377300  //Pa
#define Tcriticalpoint 304.1282 //Kelvin
#define Ptriplepoint 517950  //Pa
#define Ttriplepoint 216.592 //Kelvin
#define hcriticalPoint -174534.2750  //J/Kg
#define RhocriticalPoint 467.60 //Kg/m3
#define CpcriticalPoint 22722108.72 //J/Kg
#define MucriticalPoint 0.00003573	//Pa.s

using namespace std;

namespace Kratos
{
    //this class is to be modified by the user to customize the interpolation process
    //template< unsigned int TDim>
    class SpanWagnerEOS 
    {
    public:
	
        KRATOS_CLASS_POINTER_DEFINITION(SpanWagnerEOS);

        typedef std::vector <double> Vec;
        typedef std::vector <double> &rVec;
        typedef std::pair <double,double> Pair;
        typedef std::vector < Pair > VecPair;  
        typedef std::pair < std::pair < Vec, Vec >, std::pair < Vec, Vec > > InterpolationPoints;
        
        
        SpanWagnerEOS(/*ModelPart& model_part*/)
                /*: mr_model_part(model_part)*/              //mr_model_part is saved as private variable (declared at the end of the file)  
        {
                KRATOS_TRY	

                QDir homeDirectory = QDir::homePath();
				//This two lines has to be commented in case of use linux
				homeDirectory.cdUp();
				homeDirectory.cdUp();
				//
                homeDirectory.cd("kratos");
                homeDirectory.cd("applications");
                homeDirectory.cd("MultiphaseApplication");
                homeDirectory.cd("custom_utilities");
                QString fileNameRaw = "spanWagner";
                QString fileName = QString ("%1.xml").arg(fileNameRaw);

                this->mXmlAbsolutePath = homeDirectory.absoluteFilePath(fileName);

                KRATOS_CATCH("")	
        }
		
        void FillTables()
        {
//            this->CreateSaturationTables(this->mSaturationLinesMap);
            
            this->CreatePropertiesTables(this->mPropertiesMap);
            
            this->CreatePropertiesTables_2();
        }

        void CreateSaturationTables(std::map <QString, std::vector <std::pair <double, double> > > rSaturationLinesMap)
        {                    
            QFile file(this->mXmlAbsolutePath);
            QDomDocument *domDocument=new QDomDocument();
            domDocument->setContent(&file,true);

            QDomElement root=domDocument->firstChild().toElement();

            delete domDocument;

            Vec pressure, enthalpy, temperature, density, specificHeat;

            QList <QString> namesTables;
            namesTables <<"saturationLine" << "enthalpyGasSaturation"  << "enthalpyLiquidSaturation"  << "densityGasSaturation"  << "densityLiquidSaturation"  << "specificHeatGasSaturation"  << "specificHeatLiquidSaturation" << "viscosityGasSaturation"  << "viscosityLiquidSaturation";

            VecPair container;

            for (unsigned i=0;i< namesTables.size();i++)
            {

                    QDomNodeList nodeList = root.elementsByTagName(namesTables[i]);

                    for(unsigned nNode=0; nNode < nodeList.count(); nNode++)
                    {

                            QDomElement node =nodeList.at(nNode).toElement();
                            QString textNode=node.text().simplified();
                            QStringList listTable = textNode.split(" ", QString::SkipEmptyParts);

                            for(unsigned i=0;i < listTable.size(); i++)
                                    {
                                                    i++;

                                                    double P_hg=listTable.at(i).toDouble();
                                                    double hg=listTable.at(i-1).toDouble();
                                                    Pair Phg = std::make_pair(P_hg, hg);

                                                    container.push_back(Phg);
                                    }
                     }

                    rSaturationLinesMap.insert( std::pair< QString, std::vector < Pair > > (namesTables[i],container) );
                    container.clear();

            }

        }
                
        void CreatePropertiesTables(std::map < QString, Vec > &rPropertiesMap)
        {


            QFile file(mXmlAbsolutePath);
            QDomDocument *domDocument=new QDomDocument();
            domDocument->setContent(&file,true);

            QDomElement root=domDocument->firstChild().toElement();

            delete domDocument;

            Vec pressure, enthalpy, temperature, density, specificHeat, pressureMu, enthalpyMu, temperatureMu, viscosity;

            QDomNodeList nodeListEnthalpy = root.elementsByTagName("enthalpy");
            QDomNodeList nodeListDensity = root.elementsByTagName("density");
            QDomNodeList nodeListSpecificHeat = root.elementsByTagName("specificHeat");
            QDomNodeList nodeListViscosity = root.elementsByTagName("viscosity");

            std::vector < std::pair < QString, std::vector < std::pair <double, Pair  > > > > dataContainer;

            // Fill P, T, h, rho, Cp
            for(unsigned nNode=0; nNode < nodeListEnthalpy.count(); nNode++)
            {

                QDomElement nodeEnthalpyTable =nodeListEnthalpy.at(nNode).toElement();
                QDomElement nodeDensityTable =nodeListDensity.at(nNode).toElement();
                QDomElement nodeSpecificHeatTable =nodeListSpecificHeat.at(nNode).toElement();
                QDomElement nodeViscosityTable =nodeListViscosity.at(nNode).toElement();

                QString textNodeEnthalpyTable=nodeEnthalpyTable.text().simplified();
                QString textNodeDensityTable=nodeDensityTable.text().simplified();
                QString textNodeSpecificHeatTable=nodeSpecificHeatTable.text().simplified();
                QString textNodeViscosityTable=nodeViscosityTable.text().simplified();

                double P=nodeEnthalpyTable.attribute("data").toDouble();

                QStringList listEnthalpyTable = textNodeEnthalpyTable.split(" ", QString::SkipEmptyParts);
                QStringList listDensityTable = textNodeDensityTable.split(" ", QString::SkipEmptyParts);
                QStringList listSpecificHeatTable = textNodeSpecificHeatTable.split(" ", QString::SkipEmptyParts);
                QStringList listViscosityTable = textNodeViscosityTable.split(" ", QString::SkipEmptyParts);

                for(unsigned i=0;i < listEnthalpyTable.size(); i++)
                        {
                                        i++;

                                        double T=listEnthalpyTable.at(i-1).toDouble();
                                        double h=listEnthalpyTable.at(i).toDouble();
                                        double rho=listDensityTable.at(i).toDouble();
                                        double Cp=listSpecificHeatTable.at(i).toDouble();

                                        pressure.push_back(P);
                                        temperature.push_back(T);
                                        enthalpy.push_back(h);
                                        density.push_back(rho);
                                        specificHeat.push_back(Cp);
                        }

            }

            // Fill viscosity tables
            for(unsigned nNode=0; nNode < nodeListViscosity.count(); nNode++)
            {

                    QDomElement nodeViscosityTable =nodeListViscosity.at(nNode).toElement();

                    QString textNodeViscosityTable=nodeViscosityTable.text().simplified();

                    double PMu=nodeViscosityTable.attribute("data").toDouble();

                    QStringList listViscosityTable = textNodeViscosityTable.split(" ", QString::SkipEmptyParts);

                    for(unsigned i=0;i < listViscosityTable.size(); i++)
                            {
                                            i+=2;

                                            double TMu=listViscosityTable.at(i-2).toDouble();
                                            double hMu=listViscosityTable.at(i-1).toDouble();
                                            double Mu=listViscosityTable.at(i).toDouble();

                                            pressureMu.push_back(PMu);
                                            enthalpyMu.push_back(hMu);
                                            temperatureMu.push_back(TMu);
                                            viscosity.push_back(Mu);
                            }

            }

             rPropertiesMap.insert( std::pair< QString, Vec > ("pressure",pressure));
             rPropertiesMap.insert( std::pair< QString, Vec > ("enthalpy",enthalpy));
             rPropertiesMap.insert( std::pair< QString, Vec > ("temperature",temperature));
             rPropertiesMap.insert( std::pair< QString, Vec > ("density",density));
             rPropertiesMap.insert( std::pair< QString, Vec > ("specificHeat",specificHeat));
             rPropertiesMap.insert( std::pair< QString, Vec > ("pressureMu",pressureMu));
             rPropertiesMap.insert( std::pair< QString, Vec > ("enthalpyMu",enthalpyMu));
             rPropertiesMap.insert( std::pair< QString, Vec > ("temperatureMu",temperatureMu));
             rPropertiesMap.insert( std::pair< QString, Vec > ("viscosity",viscosity));
        }
        
        void CreatePropertiesTables_2()
        {

            QFile file(this->mXmlAbsolutePath);
            QDomDocument *domDocument=new QDomDocument();
            domDocument->setContent(&file,true);

            QDomElement root=domDocument->firstChild().toElement();

            delete domDocument;

            QDomNodeList VariableNodeList = root.elementsByTagName("variablesTable");
            QStringList variables = VariableNodeList.at(0).toElement().text().simplified().split(" ", QString::SkipEmptyParts);

            for(unsigned i=0; i < variables.size(); i++)
            {

            std::vector < std::pair <double, Pair > >   tableAux;

            QDomNodeList nodeList = root.elementsByTagName(variables.at(i));

            if(variables.at(i)=="viscosity")
            {

                            for(unsigned nNode=0; nNode < nodeList.count(); nNode++)
                            {

                            QDomElement node =nodeList.at(nNode).toElement();

                            QString textNode=node.text().simplified();
                            double data=node.attribute("data").toDouble();

                            QStringList list = textNode.split(" ", QString::SkipEmptyParts);

                            for(unsigned i=0;i < list.size(); i++)
                                    {
                                                    i+=2;

                                                    double varY=list.at(i-2).toDouble();
                                                    double varX=list.at(i).toDouble();

                                                    Pair varX_varY = std::make_pair(varY,varX);
                                                    std::pair <double, Pair > currrentline = std::make_pair(data,varX_varY);

                                                    tableAux.push_back(currrentline);
                                    }

                            }

            }
            else
            {
                            for(unsigned nNode=0; nNode < nodeList.count(); nNode++)
                            {

                            QDomElement node =nodeList.at(nNode).toElement();

                            QString textNode=node.text().simplified();
                            double data=node.attribute("data").toDouble();

                            QStringList list = textNode.split(" ", QString::SkipEmptyParts);

                            for(unsigned i=0;i < list.size(); i++)
                                    {
                                                    i++;

                                                    double varY=list.at(i-1).toDouble();
                                                    double varX=list.at(i).toDouble();

                                                    Pair varX_varY = std::make_pair(varY,varX);
                                                    pair <double, Pair > currrentline = std::make_pair(data,varX_varY);

                                                    tableAux.push_back(currrentline);
                                    }

                            }
            }

            mPropertiesMap_2.insert ( std::pair< QString, std::vector < std::pair <double, Pair > > >(variables.at(i),tableAux) );
            }

        }
        
        void computeProperties(rVec rProperties, const double aPressure, const double aEnthalpy)
        {

           const double P = aPressure;
           const double h = aEnthalpy;

           //Properties
           double T,rho,Cp,Mu,deriv_rho_P,beta,deriv_2_rho_P,deriv_beta_P,deriv_mu_P;

           //look for 4 points to interpolate in tables 
           InterpolationPoints points;
           this->LookforInterpolationPoints_PhTrhoCp(points, P,h);
           InterpolationPoints pointsMu;
           this->lookforInterpolationPoints_Mu(pointsMu, P,h);

           //interpolate temperature
           this->interpolateVariable(T, "temperature",P,h,points); /////Valores de los bordes mal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           //interpolate density
           this->interpolateVariable(rho, "density",P,h,points);
           //interpolate specificHeat
           this->interpolateVariable(Cp, "specificHeat",P,h,points);
           //interpolate viscosity
           this->interpolateVariable(Mu, "viscosity",P,h,pointsMu);
           
           //interpolate deriv_Rho_Press
           this->interpolateVariableDerivPressure(deriv_rho_P,"density",P,h,points);
           //calculate beta 
           beta = (1.0/rho)*deriv_rho_P;
           
           //interpolate deriv_2_Rho_Press 
           this->interpolateVariableSecondDerivPressure(deriv_2_rho_P,"density",P,h,points);  
           //calculate derivBeta_P
           deriv_beta_P = (1.0/rho)*(deriv_2_rho_P - beta);
           
           //interpolate deriv_mu_P 
           this->interpolateVariableSecondDerivPressure(deriv_mu_P,"viscosity",P,h,points);  

            rProperties.push_back(T);
            rProperties.push_back(rho);
            rProperties.push_back(Cp);
            rProperties.push_back(Mu);
            rProperties.push_back(beta);
            rProperties.push_back(deriv_rho_P);
            rProperties.push_back(deriv_2_rho_P);
            rProperties.push_back(deriv_beta_P);
            rProperties.push_back(deriv_mu_P);
//                std::cout << "rho = " << rho  << std::endl; 
//        std::cout << "Mu = " << Mu << std::endl; 
//        std::cout << "beta = " << beta  << std::endl;
	}
   
        void lookforValueOnSaturationLine (double& rInterpolatedValue, const QString aSaturationLineName, const double aPressure)
        {

            std::vector < Pair > saturationLineTable = this->mSaturationLinesMap.find(aSaturationLineName)->second;

            Pair  storagePressure; //P1,P2
            Pair  storageEnthalpy; //h1,h2

            for(unsigned i=1;i < saturationLineTable.size(); i++)
                if(aPressure <= saturationLineTable[i].first && aPressure >= saturationLineTable[i-1].first)
                {
                    storagePressure=std::make_pair(saturationLineTable[i-1].first,saturationLineTable[i].first);//P
                    storageEnthalpy=std::make_pair(saturationLineTable[i-1].second,saturationLineTable[i].second);//h
                     break;
                }

            rInterpolatedValue =storageEnthalpy.first
                                +(( storageEnthalpy.second - storageEnthalpy.first )/( storagePressure.first - storagePressure.second ))
                                *(storagePressure.first-aPressure);
        }
        
        ~SpanWagnerEOS()
        {
            
        }

        void Calculate()  
        {
            KRATOS_TRY
//            double area;                 
//            for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); //looping the elements
//                ielem!=mr_model_part.ElementsEnd(); ielem++)
//            {
//                    Geometry<Node<3> >& geom = ielem->GetGeometry(); 
//                    area=CalculateArea(geom);                      
//            }

            KRATOS_CATCH("")
        } 
	

        void LookforInterpolationPoints_PhTrhoCp (InterpolationPoints& rInterpolationPoints, const double aPressure, const double aEnthalpy)
        {
	// Interpolation tables out of transitionPhase
	Vec pressureTable = this->mPropertiesMap.find("pressure")->second;
	Vec enthalpyTable = this->mPropertiesMap.find("enthalpy")->second;
	Vec temperatureTable = this->mPropertiesMap.find("temperature")->second;
	Vec densityTable = this->mPropertiesMap.find("density")->second;
	Vec specificHeatTable = this->mPropertiesMap.find("specificHeat")->second;

	Vec highPressure,highPressureEnthalpy,highPressureTemperature,highPressureDensity,highPressureSpecificHeat;
	Vec lowPressure,lowPressureEnthalpy,lowPressureTemperature,lowPressureDensity,lowPressureSpecificHeat;
	double highPressureNr=0;
	double lowPressureNr=0;
	for(unsigned i=0;i < pressureTable.size()-1; i++)
		if(aPressure <= pressureTable[i+1] && aPressure >= pressureTable[i])
		{
			highPressureNr = pressureTable[i+1];
			lowPressureNr = pressureTable[i];

			break;

		}

	if(highPressureNr==0 || lowPressureNr==0)
            KRATOS_THROW_ERROR(std::logic_error, "Pressure is out of tables range during interpolation, value of pressure is and Enthalpy is ", "");
//            ("Pressure is out of tables range during interpolation, value of pressure is %0 and Enthalpy is %1").arg(aPressure).arg(aEnthalpy)
        
	for(unsigned i=0;i < pressureTable.size()-1; i++)
	{
		if( lowPressureNr == pressureTable[i] )
		{
					lowPressure.push_back(pressureTable[i]);
					lowPressureEnthalpy.push_back(enthalpyTable[i]);
					lowPressureTemperature.push_back(temperatureTable[i]);
					lowPressureDensity.push_back(densityTable[i]);
					lowPressureSpecificHeat.push_back(specificHeatTable[i]);
		}
		else if (highPressureNr == pressureTable[i])
		{
					highPressure.push_back(pressureTable[i+1]);
					highPressureEnthalpy.push_back(enthalpyTable[i+1]);
					highPressureTemperature.push_back(temperatureTable[i+1]);
					highPressureDensity.push_back(densityTable[i+1]);
					highPressureSpecificHeat.push_back(specificHeatTable[i+1]);
		}
	}

	this->increasingOrder5 (lowPressureEnthalpy, lowPressure, lowPressureTemperature,lowPressureDensity,lowPressureSpecificHeat);
	this->increasingOrder5 (highPressureEnthalpy, highPressure, highPressureTemperature,highPressureDensity,highPressureSpecificHeat);

	pair < Vec, Vec > pointsLowerPressure, pointsHigherPressure;

	for(unsigned i=0;i < lowPressureEnthalpy.size()-1; i++)
		if(aEnthalpy <= lowPressureEnthalpy[i+1] && aEnthalpy >= lowPressureEnthalpy[i])
		{
			Vec TempRhoCpi;
			TempRhoCpi.push_back(lowPressure[i]);
			TempRhoCpi.push_back(lowPressureEnthalpy[i]);
			TempRhoCpi.push_back(lowPressureTemperature[i]);
			TempRhoCpi.push_back(lowPressureDensity[i]);
			TempRhoCpi.push_back(lowPressureSpecificHeat[i]);

			Vec TempRhoCpiPlus1;
			TempRhoCpiPlus1.push_back(lowPressure[i+1]);
			TempRhoCpiPlus1.push_back(lowPressureEnthalpy[i+1]);
			TempRhoCpiPlus1.push_back(lowPressureTemperature[i+1]);
			TempRhoCpiPlus1.push_back(lowPressureDensity[i+1]);
			TempRhoCpiPlus1.push_back(lowPressureSpecificHeat[i+1]);

			pointsLowerPressure = std::make_pair(TempRhoCpi,TempRhoCpiPlus1);

			break;
		}

	 for(unsigned i=0;i < highPressureEnthalpy.size()-1; i++)
		if(aEnthalpy <= highPressureEnthalpy[i+1] && aEnthalpy >= highPressureEnthalpy[i])
		{
			Vec TempRhoCpi;
			TempRhoCpi.push_back(highPressure[i]);
			TempRhoCpi.push_back(highPressureEnthalpy[i]);
			TempRhoCpi.push_back(highPressureTemperature[i]);
			TempRhoCpi.push_back(highPressureDensity[i]);
			TempRhoCpi.push_back(highPressureSpecificHeat[i]);

			Vec TempRhoCpiPlus1;
			TempRhoCpiPlus1.push_back(highPressure[i+1]);
			TempRhoCpiPlus1.push_back(highPressureEnthalpy[i+1]);
			TempRhoCpiPlus1.push_back(highPressureTemperature[i+1]);
			TempRhoCpiPlus1.push_back(highPressureDensity[i+1]);
			TempRhoCpiPlus1.push_back(highPressureSpecificHeat[i+1]);

			pointsHigherPressure = std::make_pair(TempRhoCpi,TempRhoCpiPlus1);

			break;
		}

		if (pointsLowerPressure.first.empty() || pointsLowerPressure.second.empty() || pointsHigherPressure.first.empty() || pointsHigherPressure.second.empty())
		{

//			QString msg= QString("Enthalpy is out of tables range during interpolation, value of pressure is %0 and Enthalpy is %1").arg(aPressure).arg(aEnthalpy);
//			GeneralException e(msg);
//			throw(e);
		}

                rInterpolationPoints = std::make_pair(pointsLowerPressure,pointsHigherPressure);
        }	

        void lookforInterpolationPoints_Mu (InterpolationPoints& rInterpolationPoints,const double aPressure, const double aEnthalpy)
        {
            // Interpolation tables out of transitionPhase
            Vec pressureTable = this->mPropertiesMap.find("pressureMu")->second;
            Vec  enthalpyTable = this->mPropertiesMap.find("enthalpyMu")->second;
            Vec  temperatureTable = this->mPropertiesMap.find("temperatureMu")->second;
            Vec  viscosityTable = this->mPropertiesMap.find("viscosity")->second;

            Vec highPressure,highPressureEnthalpy,highPressureTemperature,highPressureViscosity;
            Vec lowPressure,lowPressureEnthalpy,lowPressureTemperature,lowPressureViscosity;
            double highPressureNr=0;
            double lowPressureNr=0;
            for(unsigned i=0;i < pressureTable.size()-1; i++)
                    if(aPressure <= pressureTable[i+1] && aPressure >= pressureTable[i])
                    {
                            highPressureNr = pressureTable[i+1];
                            lowPressureNr = pressureTable[i];

                            break;

                    }

            if(highPressureNr==0 || lowPressureNr==0)
            {
    //		QString msg= QString("Pressure is out of tables range during interpolation, value of pressure is %0 and Enthalpy is %1").arg(aPressure).arg(aEnthalpy);
    //		GeneralException e(msg);
    //		throw(e);
            }

            for(unsigned i=0;i < pressureTable.size()-1; i++)
            {
                    if( lowPressureNr == pressureTable[i] )
                    {
                                            lowPressure.push_back(pressureTable[i]);
                                            lowPressureEnthalpy.push_back(enthalpyTable[i]);
                                            lowPressureTemperature.push_back(temperatureTable[i]);
                                            lowPressureViscosity.push_back(viscosityTable[i]);
                    }
                    else if (highPressureNr == pressureTable[i])
                    {
                                            highPressure.push_back(pressureTable[i+1]);
                                            highPressureEnthalpy.push_back(enthalpyTable[i+1]);
                                            highPressureTemperature.push_back(temperatureTable[i+1]);
                                            highPressureViscosity.push_back(viscosityTable[i]);
                    }
            }

            this->increasingOrder4 (lowPressureEnthalpy, lowPressure, lowPressureTemperature,lowPressureViscosity);
            this->increasingOrder4 (highPressureEnthalpy, highPressure, highPressureTemperature,highPressureViscosity);

            pair < Vec , Vec > pointsLowerPressure, pointsHigherPressure;

            for(unsigned i=0;i < lowPressureEnthalpy.size()-1; i++)
                    if(aEnthalpy <= lowPressureEnthalpy[i+1] && aEnthalpy >= lowPressureEnthalpy[i])
                    {
                            Vec TempRhoCpi;
                            TempRhoCpi.push_back(lowPressure[i]);
                            TempRhoCpi.push_back(lowPressureEnthalpy[i]);
                            TempRhoCpi.push_back(lowPressureTemperature[i]);
                            TempRhoCpi.push_back(lowPressureViscosity[i]);

                            Vec TempRhoCpiPlus1;
                            TempRhoCpiPlus1.push_back(lowPressure[i+1]);
                            TempRhoCpiPlus1.push_back(lowPressureEnthalpy[i+1]);
                            TempRhoCpiPlus1.push_back(lowPressureTemperature[i+1]);
                            TempRhoCpiPlus1.push_back(lowPressureViscosity[i+1]);

                            pointsLowerPressure = std::make_pair(TempRhoCpi,TempRhoCpiPlus1);

                            break;
                    }

             for(unsigned i=0;i < highPressureEnthalpy.size()-1; i++)
                    if(aEnthalpy <= highPressureEnthalpy[i+1] && aEnthalpy >= highPressureEnthalpy[i])
                    {
                            Vec TempRhoCpi;
                            TempRhoCpi.push_back(highPressure[i]);
                            TempRhoCpi.push_back(highPressureEnthalpy[i]);
                            TempRhoCpi.push_back(highPressureTemperature[i]);
                            TempRhoCpi.push_back(highPressureViscosity[i]);

                            Vec TempRhoCpiPlus1;
                            TempRhoCpiPlus1.push_back(highPressure[i+1]);
                            TempRhoCpiPlus1.push_back(highPressureEnthalpy[i+1]);
                            TempRhoCpiPlus1.push_back(highPressureTemperature[i+1]);
                            TempRhoCpiPlus1.push_back(highPressureViscosity[i+1]);

                            pointsHigherPressure = std::make_pair(TempRhoCpi,TempRhoCpiPlus1);

                            break;
                    }

                    if (pointsLowerPressure.first.empty() || pointsLowerPressure.second.empty() || pointsHigherPressure.first.empty() || pointsHigherPressure.second.empty())
                    {

    //			QString msg= QString("Enthalpy is out of tables range during interpolation, value of pressure is %0 and Enthalpy is %1").arg(aPressure).arg(aEnthalpy);
    //			GeneralException e(msg);
    //			throw(e);
                    }

            rInterpolationPoints = std::make_pair(pointsLowerPressure,pointsHigherPressure);
		
        }	

        void interpolateVariable (double& rVariableValue, const QString aVariable, const double aPressure, const double aEnthalpy, const  InterpolationPoints aPoints)
        {
            double variableID;
            if(aVariable=="temperature")
                    variableID=2;
            else if(aVariable=="density")
                    variableID=3;
            else if(aVariable=="specificHeat")
                    variableID=4;
            else if(aVariable=="viscosity")
                    variableID=3;
            else
            {
    //		GeneralException e("This variable to interpolate desnt exist");
    //		throw(e);
            }

            pair < Vec, Vec > pointsLowerPressure, pointsHigherPressure;

            pointsLowerPressure=aPoints.first;
            pointsHigherPressure=aPoints.second;

            const double variableHigherPressure = pointsHigherPressure.first[variableID] + ( (pointsHigherPressure.second[variableID] -  pointsHigherPressure.first[variableID]) / (pointsHigherPressure.second[1] -  pointsHigherPressure.first[1]) )*(aEnthalpy - pointsHigherPressure.first[1]);
            const double variableLowerPressure  = pointsLowerPressure.first[variableID] + ( (pointsLowerPressure.second[variableID] -  pointsLowerPressure.first[variableID]) / (pointsLowerPressure.second[1] -  pointsLowerPressure.first[1]) )*(aEnthalpy - pointsLowerPressure.first[1]);

            rVariableValue = variableLowerPressure + ( (variableHigherPressure -  variableLowerPressure) / (pointsHigherPressure.second[0] -  pointsLowerPressure.second[0]) )*(aPressure - pointsLowerPressure.second[0]);
//              std::cout <<  rVariableValue  << std::endl;
        }
          
        void interpolateVariableDerivPressure (double& rDerivateValue, const QString aVariable, const double aPressure, const double aEnthalpy, const  InterpolationPoints aPoints )
        {
                double variableID;
                if(aVariable=="temperature")
                        variableID=2;
                else if(aVariable=="density")
                        variableID=3;
                else if(aVariable=="specificHeat")
                        variableID=4;
                else if(aVariable=="viscosity")
                        variableID=3;
                else
                {
//                        GeneralException e("This variable to interpolate desnt exist");
//                        throw(e);
                }

                pair < Vec, Vec > pointsLowerPressure, pointsHigherPressure;

                pointsLowerPressure=aPoints.first;
                pointsHigherPressure=aPoints.second;

                const double variableHigherPressure = pointsHigherPressure.first[variableID] + ( (pointsHigherPressure.second[variableID] -  pointsHigherPressure.first[variableID]) / (pointsHigherPressure.second[1] -  pointsHigherPressure.first[1]) )*(aEnthalpy - pointsHigherPressure.first[1]);
                const double variableLowerPressure  = pointsLowerPressure.first[variableID] + ( (pointsLowerPressure.second[variableID] -  pointsLowerPressure.first[variableID]) / (pointsLowerPressure.second[1] -  pointsLowerPressure.first[1]) )*(aEnthalpy - pointsLowerPressure.first[1]);

                rDerivateValue = (variableHigherPressure -  variableLowerPressure) / (pointsHigherPressure.second[0] -  pointsLowerPressure.second[0]); 
       
        }
        
        void interpolateVariableSecondDerivPressure(double& rDerivateValue, const QString aVariable, const double aPressure, const double aEnthalpy, const  InterpolationPoints aPoints )
        {
            double variableID;
            if(aVariable=="temperature")
                    variableID=2;
            else if(aVariable=="density")
                    variableID=3;
            else if(aVariable=="specificHeat")
                    variableID=4;
            else if(aVariable=="viscosity")
                    variableID=3;
            else
            {
//                        GeneralException e("This variable to interpolate desnt exist");
//                        throw(e);
            }
            
            pair < Vec, Vec > pointsLowerPressure, pointsHigherPressure;

            pointsLowerPressure=aPoints.first;
            pointsHigherPressure=aPoints.second;

            const double variableHigherPressure = pointsHigherPressure.first[variableID] + ( (pointsHigherPressure.second[variableID] -  pointsHigherPressure.first[variableID]) / (pointsHigherPressure.second[1] -  pointsHigherPressure.first[1]) )*(aEnthalpy - pointsHigherPressure.first[1]);
            const double variableLowerPressure  = pointsLowerPressure.first[variableID] + ( (pointsLowerPressure.second[variableID] -  pointsLowerPressure.first[variableID]) / (pointsLowerPressure.second[1] -  pointsLowerPressure.first[1]) )*(aEnthalpy - pointsLowerPressure.first[1]);
            const double variablePointValue = variableLowerPressure + ( (variableHigherPressure -  variableLowerPressure) / (pointsHigherPressure.second[0] -  pointsLowerPressure.second[0]) )*(aPressure - pointsLowerPressure.second[0]);
            
            rDerivateValue = (variableHigherPressure - 2*variablePointValue + variableLowerPressure)/( pow( (pointsHigherPressure.second[0] -  pointsLowerPressure.second[0]), 2.0 ) );
            
        }
        
        void increasingOrder4 (Vec &mainVectToOrdering, Vec  &secondVecToOrdering, Vec  &thirdVecToOrdering,Vec  &fourthVecToOrdering)
        {
            double aux;
            for(unsigned c1=0;c1<=mainVectToOrdering.size()-1;c1++)
                for(unsigned c2=0;c2<mainVectToOrdering.size()-1;c2++)
                    if(mainVectToOrdering[c2]>mainVectToOrdering[c2+1])
                    {

                    aux=mainVectToOrdering[c2];
                    mainVectToOrdering[c2]=mainVectToOrdering[c2+1];
                    mainVectToOrdering[c2+1]=aux;

                    aux=secondVecToOrdering[c2];
                    secondVecToOrdering[c2]=secondVecToOrdering[c2+1];
                    secondVecToOrdering[c2+1]=aux;

                    aux=thirdVecToOrdering[c2];
                    thirdVecToOrdering[c2]=thirdVecToOrdering[c2+1];
                    thirdVecToOrdering[c2+1]=aux;

                    aux=fourthVecToOrdering[c2];
                    fourthVecToOrdering[c2]=fourthVecToOrdering[c2+1];
                    fourthVecToOrdering[c2+1]=aux;

                    }

        }

        void increasingOrder5 (Vec &mainVectToOrdering, Vec  &secondVecToOrdering, Vec  &thirdVecToOrdering,Vec  &fourthVecToOrdering,Vec  &fifthVecToOrdering)
        {
            double aux;
            for(unsigned c1=0;c1<=mainVectToOrdering.size()-1;c1++)
                for(unsigned c2=0;c2<mainVectToOrdering.size()-1;c2++)
                        if(mainVectToOrdering[c2]>mainVectToOrdering[c2+1])
                        {

                        aux=mainVectToOrdering[c2];
                        mainVectToOrdering[c2]=mainVectToOrdering[c2+1];
                        mainVectToOrdering[c2+1]=aux;

                        aux=secondVecToOrdering[c2];
                        secondVecToOrdering[c2]=secondVecToOrdering[c2+1];
                        secondVecToOrdering[c2+1]=aux;

                        aux=thirdVecToOrdering[c2];
                        thirdVecToOrdering[c2]=thirdVecToOrdering[c2+1];
                        thirdVecToOrdering[c2+1]=aux;

                        aux=fourthVecToOrdering[c2];
                        fourthVecToOrdering[c2]=fourthVecToOrdering[c2+1];
                        fourthVecToOrdering[c2+1]=aux;

                        aux=fifthVecToOrdering[c2];
                        fifthVecToOrdering[c2]=fifthVecToOrdering[c2+1];
                        fifthVecToOrdering[c2+1]=aux;

                        }
        }
    
        void interpolate_calculate_variable_PT(const string aVariableName, double& rVariableValue, const double aResult, const double aVariableA)
        {
            /*
            "enthalpy";
            "density";
            "specificHeat";
            "viscosity";
            */

            std::vector < std::pair <double, Pair > > tableAux=mPropertiesMap_2.find(QString::fromStdString(aVariableName))->second;
            std::vector < std::pair <double, Pair > > pointsLowerResult, pointsHigherResult;
            double resulti=0;
            double resultiplus1=0;
            double variableBLowerResult,variableBHigherResult;

            for(unsigned i=0;i < tableAux.size()-1; i++)
                    if(aResult <= tableAux[i+1].first && aResult >= tableAux[i].first  && tableAux[i+1].first != tableAux[i].first)
                    {
                            resulti=tableAux[i].first;
                            resultiplus1=tableAux[i+1].first;
                    }
                    /*else if(result == tableAux[i].first)
                    {

                            double variableBHigherResult = pointsHigherResult[0].second.second + ( (pointsHigherResult[1].second.second -  pointsHigherResult[0].second.second) / (pointsHigherResult[1].second.first -  pointsHigherResult[0].second.first) )*(variableA - pointsHigherResult[0].second.first);

                            aResult[0]=variableBHigherResult;
                            return;
                    }*/

                    if (!resulti || !resultiplus1)
                    {
//                            GeneralException e("result number not in table on  interpolateVariableB  "+variable+"   "+QString::number(result));
//                            throw(e);
                    }

            for(unsigned i=0;i < tableAux.size()-1; i++)
                    if(tableAux[i].first == resulti && aVariableA <= tableAux[i+1].second.first &&  aVariableA >= tableAux[i].second.first)
                            {
                                    pointsLowerResult.push_back(std::make_pair(tableAux[i].first,std::make_pair(tableAux[i].second.first,tableAux[i].second.second)));
                                    pointsLowerResult.push_back(std::make_pair(tableAux[i+1].first,std::make_pair(tableAux[i+1].second.first,tableAux[i+1].second.second)));
                            }
                    else if(tableAux[i].first == resultiplus1 && aVariableA <= tableAux[i+1].second.first &&  aVariableA >= tableAux[i].second.first)
                            {
                                    pointsHigherResult.push_back(std::make_pair(tableAux[i].first,std::make_pair(tableAux[i].second.first,tableAux[i].second.second)));
                                    pointsHigherResult.push_back(std::make_pair(tableAux[i+1].first,std::make_pair(tableAux[i+1].second.first,tableAux[i+1].second.second)));
                            }
                    else if(!pointsHigherResult.empty() && !pointsLowerResult.empty())
                            break;

                    if (pointsHigherResult.empty() || pointsLowerResult.empty())
                    {
//                            GeneralException e("variableA number not in table on  interpolateVariableB  "+variable+"   "+QString::number(variableA));
//                            throw(e);
                    }

                    variableBHigherResult = pointsHigherResult[0].second.second + ( (pointsHigherResult[1].second.second -  pointsHigherResult[0].second.second) / (pointsHigherResult[1].second.first -  pointsHigherResult[0].second.first) )*(aVariableA - pointsHigherResult[0].second.first);
                    variableBLowerResult = pointsLowerResult[0].second.second + ( (pointsLowerResult[1].second.second -  pointsLowerResult[0].second.second) / (pointsLowerResult[1].second.first -  pointsLowerResult[0].second.first) )*(aVariableA - pointsLowerResult[0].second.first);

                    double variableBFinal = variableBLowerResult + ( (variableBHigherResult -  variableBLowerResult) / (pointsHigherResult[0].first -  pointsLowerResult[0].first) )*(aResult - pointsLowerResult[0].first);

                    rVariableValue=variableBFinal;
        }
        
        void calculateTemperature(double& rTemperatureValue, const double aPressureNode, const double aEnthalpyNode)
        {
           //look for 4 points to interpolate in tables 
           InterpolationPoints points;
           this->LookforInterpolationPoints_PhTrhoCp(points, aPressureNode,aEnthalpyNode);
           //interpolate temperature
           this->interpolateVariable(rTemperatureValue, "temperature", aPressureNode, aEnthalpyNode, points);
        }
        
        void calculateDensity(double& rDensityValue, const double aPressureNode, const double aEnthalpyNode)
        {
           //look for 4 points to interpolate in tables 
           InterpolationPoints points;
           this->LookforInterpolationPoints_PhTrhoCp(points, aPressureNode,aEnthalpyNode);
           //interpolate temperature
           this->interpolateVariable(rDensityValue, "density", aPressureNode, aEnthalpyNode, points);
        }
        
	private:
            
        /// \brief Table to storage data from zone of out of transition phase
        std::map < QString, std::vector <double> > mPropertiesMap;

        std::map < QString, std::vector < std::pair < double, Pair > > >  mPropertiesMap_2;
        
        /// \brief Table to storage salturation values from all variables function of Pressure -> P,hl  P,hg	P,rhol	P,rhol	P,rhog  P,Cpl	P,Cpg
        std::map <QString, std::vector <std::pair <double, double> > > mSaturationLinesMap;

        double CalculateArea(Element::GeometryType& geom)
        {
                return 0.5 * ((geom[1].X() - geom[0].X())*(geom[2].Y() - geom[0].Y())- (geom[1].Y() - geom[0].Y())*(geom[2].X() - geom[0].X()));
        }

//        ModelPart& mr_model_part;

        QString mXmlAbsolutePath;
        
	};

}  // namespace Kratos.

#endif // KRATOS_SPANWAGNEREOS_UTILITY_INCLUDED  defined
