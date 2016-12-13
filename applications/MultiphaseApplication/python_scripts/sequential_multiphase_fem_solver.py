from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.MultiphaseApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def AddVariables(model_part):
	model_part.AddNodalSolutionStepVariable(PRESSURE_WET)
	model_part.AddNodalSolutionStepVariable(SATURATION_NON)
	model_part.AddNodalSolutionStepVariable(ENTHALPY_NON_NODE)
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_NODE)
	model_part.AddNodalSolutionStepVariable(PRESSURE_WET_OLD_ITER)
	model_part.AddNodalSolutionStepVariable(SATURATION_NON_OLD_ITER)
	model_part.AddNodalSolutionStepVariable(ENTHALPY_NON_OLD_ITER)
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_OLD_ITER)
	model_part.AddNodalSolutionStepVariable(PRESSURE_NON)
	model_part.AddNodalSolutionStepVariable(SATURATION_WET)
	model_part.AddNodalSolutionStepVariable(PERMEABILITY_WET_NODE)
	model_part.AddNodalSolutionStepVariable(PERMEABILITY_NON_NODE)
	model_part.AddNodalSolutionStepVariable(CAPILLARITY_PRESSURE_NODE)
	model_part.AddNodalSolutionStepVariable(DERIV_PC_SN_NODE)
	model_part.AddNodalSolutionStepVariable(DENSITY_NON_NODE)
	model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT_NON_NODE)
	model_part.AddNodalSolutionStepVariable(VISCOSITY_NON_NODE)
	model_part.AddNodalSolutionStepVariable(COMPRESSIBLITY_NON_NODE)
	model_part.AddNodalSolutionStepVariable(NODAL_AREA)
	model_part.AddNodalSolutionStepVariable(N_NODES)
	model_part.AddNodalSolutionStepVariable(STORAGE_SUM_BALANCE)
	model_part.AddNodalSolutionStepVariable(DARCY_FLOW_SUM_BALANCE)
	model_part.AddNodalSolutionStepVariable(DIRICHLET_SUM_BALANCE)
	model_part.AddNodalSolutionStepVariable(NEUMANN_SUM_BALANCE)
	print("variables for the fractional iterative solver added correctly")

def AddDofs(model_part):

	for node in model_part.Nodes:
		node.AddDof(PRESSURE_WET)        
		node.AddDof(SATURATION_NON)
		node.AddDof(ENTHALPY_NON_NODE)
		node.AddDof(TEMPERATURE_NODE)

	print("dofs for the fractional iterative solver added correctly")

class SequentialMultiphaseFEMSolver:

	def __init__(self, model_part, domain_size):

		self.model_part = model_part
		self.domain_size = domain_size

		#General parameters
		#Picard configuration parameters
		self.timeOrder = 1
		self.predictorCorrector = False
		self.CalculateReactions = False
		self.ReformDofAtEachIteration = True
		self.CalculateNormDxFlag = True
		self.MoveMeshFlag = False
		self.echoLevel = 2
        
		# definition of the solvers
		self.first_linear_solver = SkylineLUFactorizationSolver()
		self.second_linear_solver = SkylineLUFactorizationSolver()
		self.third_linear_solver = SkylineLUFactorizationSolver()
		self.fourth_linear_solver = SkylineLUFactorizationSolver()

		#others

		#pDiagPrecond = DiagonalPreconditioner() X
		#pILUPrecond = ILU0Preconditioner()
		#self.pressure_linear_solver = BICGSTABSolver(1e-3, 5000, pDiagPrecond) X
		#self.concentration_linear_solver = BICGSTABSolver(1e-3, 5000, pDiagPrecond) X
		#self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
		#self.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)

		#tol = 1e-5 #stopping tolerance...  1e-6 
		#max_it = 30 #max number of  iterations.... 500
		#precond = DiagonalPreconditioner()
		#self.NR_solver = CGSolver(tol,max_it, precond)

		#spalart_allmaras_linear_solver 
		#pPrecond = DiagonalPreconditioner()
		#self.NR_solver = BICGSTABSolver(1e-9, 5000, pPrecond)
    
	def Initialize(self):

		self.domain_size = int(self.domain_size)
		self.solverConfiguration = SequentialMultiphaseFEMConfiguration(
			self.model_part,
			self.first_linear_solver,
			self.second_linear_solver,
			self.third_linear_solver,
			self.fourth_linear_solver, 
			self.domain_size)

		self.timeOrder = int(self.timeOrder)
		self.predictorCorrector = bool(self.predictorCorrector)
		self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
		self.MoveMeshFlag = bool(self.MoveMeshFlag)
		self.echoLevel = int(self.echoLevel)

		self.solver = SequentialMultiphaseFEMStrategy(
			self.model_part,
			self.solverConfiguration,
			self.ReformDofAtEachIteration,
			self.timeOrder,
			self.domain_size,
			self.predictorCorrector,
			self.MoveMeshFlag,
			self.echoLevel)

		self.solver.Check()

		(self.solver).SetEchoLevel(self.echoLevel)
		print("finished initialization of the Sequential Multiphase Strategy")

	def Solve(self):

		print("just before solve")
		print(self.model_part)
		(self.solver).Solve()

	def Clear(self):
		(self.solver).Clear()
