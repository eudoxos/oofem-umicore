/*
 * DarcyFlow.C
 *
 *  Created on: Feb 25, 2010
 *      Author: carl
 */

#include "darcyflow.h"
#include "element.h"
#include "inputrecord.h"
#include "timestep.h"
#include "classfactory.h"
#include "sparselinsystemnm.h"
#include "mathfem.h"
#include "tr1darcy.h"
#include "sparsemtrx.h"
#include "nrsolver.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "processcomm.h"
#endif

#include <iostream>
#include <fstream>
#include <cstdio>

namespace oofem {
REGISTER_EngngModel(DarcyFlow);

DarcyFlow :: DarcyFlow(int i, EngngModel *_master) : EngngModel(i, _master)
{
    this->nMethod = NULL;
    this->ndomains = 1;
    this->hasAdvanced = false;
    this->stiffnessMatrix = NULL;
}

DarcyFlow :: ~DarcyFlow()
{
    delete PressureField;
    delete nMethod;
}

IRResultType DarcyFlow :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    // Create solution space for EID_ConservationEquation
    PressureField = new PrimaryField(this, 1, FT_Pressure, EID_ConservationEquation, 1);
#if 0
 #ifdef __PARALLEL_MODE


    printf("Parallel mode!\n");
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                               this->giveNumberOfProcesses(),
                                               this->commMode);
    }

 #endif
#endif
    return IRRT_OK;
}

void DarcyFlow :: solveYourselfAt(TimeStep *tStep)
{
    /*
     * Assemble and solve system of equations as given timestep tStep.
     *
     */

    OOFEM_LOG_INFO("Parallelflag: %u\n", parallelFlag);

    FloatArray *solutionVector = NULL;
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    // Move solution space to current timestep
    if ( !hasAdvanced ) {
        PressureField->advanceSolution(tStep);
        hasAdvanced = true;
    }

    // Point pointer SolutionVector to current solution in VelocityPressureField
    solutionVector = PressureField->giveSolutionVector(tStep);
    solutionVector->resize(neq);
    solutionVector->zero();

    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        this->stiffnessMatrix->buildInternalStructure( this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering() );
    }


    // Build initial/external load (LoadVector)
    this->externalForces.resize(neq);
    this->externalForces.zero();
    this->assembleVectorFromElements( this->externalForces, tStep, EID_ConservationEquation, ExternalForcesVector, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );

    this->incrementOfSolution.resize(neq);
    this->internalForces.resize(neq);

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    double loadLevel;
    int currentIterations;
    this->updateComponent( tStep, InternalRhs, this->giveDomain(1) );
    this->updateComponent( tStep, NonLinearLhs, this->giveDomain(1) );

    FloatArray incrementalLoadVector(0); // Should be allowed to be null
    NM_Status status = this->nMethod->solve(this->stiffnessMatrix,
                                            & ( this->externalForces ),
                                            NULL,
                                            solutionVector,
                                            & ( this->incrementOfSolution ),
                                            & ( this->internalForces ),
                                            this->ebeNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total, // Why this naming scheme? Should be RLM_Total, and ReferenceLoadInputModeType
                                            currentIterations,
                                            tStep);

    if ( status & NM_NoSuccess ) {
        OOFEM_ERROR2( "DarcyFlow :: couldn't solve for time step %d\n", tStep->giveNumber() );
    }

#define DUMPMATRICES 0
#if DUMPMATRICES
    FloatMatrix LHS_backup;
    lhs->toFloatMatrix(LHS_backup);
    DumpMatricesToFile(& LHS_backup, & rhs, NULL);
#endif

    this->updateYourself(tStep);
}

void DarcyFlow :: DumpMatricesToFile(FloatMatrix *LHS, FloatArray *RHS, FloatArray *SolutionVector)
{
    FloatMatrix K;

    FILE *rhsFile = fopen("RHS.txt", "w");
    // rhs.printYourself();

    for ( int i = 1; i <= RHS->giveSize(); i++ ) {
        fprintf( rhsFile, "%0.15e\n", RHS->at(i) );
    }
    fclose(rhsFile);

    FILE *lhsFile = fopen("LHS.txt", "w");

    for ( int i = 1; i <= this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ); i++ ) {
        for ( int j = 1; j <= this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ); j++ ) {
            fprintf( lhsFile, "%0.15e\t", LHS->at(i, j) );
        }
        fprintf(lhsFile, "\n");
    }

    fclose(lhsFile);

    if ( SolutionVector == NULL ) {
        return;
    }

    FILE *SolutionFile = fopen("SolutionVector.txt", "w");
    for ( int i = 1; i <= SolutionVector->giveSize(); i++ ) {
        fprintf( SolutionFile, "%0.15e\n", SolutionVector->at(i) );
    }
    fclose(SolutionFile);
}
void DarcyFlow :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    DofIDItem type = iDof->giveDofID();
    if ( type == P_f ) {
        iDof->printSingleOutputAt(stream, tStep, 'p', VM_Total, 1);
    } else {
        _error("printDofOutputAt: unsupported dof type");
    }
}

void DarcyFlow :: updateYourself(TimeStep *tStep)
{
    EngngModel :: updateYourself(tStep);
}

double DarcyFlow :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return PressureField->giveUnknownValue(dof, mode, tStep);
}

void DarcyFlow :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    switch ( cmpn ) {
    case InternalRhs:
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, EID_ConservationEquation,  InternalForcesVector, VM_Total,
                             EModelDefaultEquationNumbering(), d, & this->ebeNorm);
        break;

    case NonLinearLhs:

        this->stiffnessMatrix->zero();
        this->assemble( this->stiffnessMatrix, tStep, EID_ConservationEquation, StiffnessMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        break;

    default:
        _error2("updateComponent: Unknown component id (%d)", ( int ) cmpn);
    }
}

int DarcyFlow :: forceEquationNumbering(int id)  // Is this really needed???!?
{
    int neq = EngngModel :: forceEquationNumbering(id);

    this->equationNumberingCompleted = false;
    if ( this->stiffnessMatrix ) {
        delete this->stiffnessMatrix;
        this->stiffnessMatrix = NULL;
    }

    return neq;
}

NumericalMethod *DarcyFlow :: giveNumericalMethod(MetaStep *mStep)
{
    /*
     * Returns pointer to NumericalMethod object. The solver used for StokesFlow is SparseLinearSystemNM.
     * If no solver has bee initialized, create one, otherwise, return the existing solver.
     */

    if ( this->nMethod ) {
        return this->nMethod;
    }

    this->nMethod = new NRSolver(this->giveDomain(1), this);
    if ( !nMethod ) {
        OOFEM_ERROR("giveNumericalMethod: numerical method creation failed");
    }
    return this->nMethod;
}

TimeStep *DarcyFlow :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();

    StateCounterType counter = 1;
    delete previousStep;

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, ( double ) istep, 0., counter);
    // time and dt variables are set eq to 0 for statics - has no meaning
    return currentStep;
}

#ifdef __PETSC_MODULE
void DarcyFlow :: initPetscContexts()
{
    PetscContext *petscContext;
    petscContextList->growTo(ndomains);
    for ( int i = 1; i <= this->ndomains; i++ ) {
        petscContext =  new PetscContext(this);
        petscContextList->put(i, petscContext);
    }
}
#endif
}
