/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "linearstatic.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "sparsemtrx.h"
#include "verbose.h"
#include "structuralelement.h"
#include "structuralelementevaluator.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

#ifdef __PARALLEL_MODE
 #include "fetisolver.h"
 #include "sparsemtrx.h"
 #include "problemcomm.h"
 #include "communicator.h"
#endif

namespace oofem {
REGISTER_EngngModel(LinearStatic);

LinearStatic :: LinearStatic(int i, EngngModel *_master) : StructuralEngngModel(i, _master), loadVector(), displacementVector()
{
    stiffnessMatrix = NULL;
    ndomains = 1;
    nMethod = NULL;
    initFlag = 1;
    solverType = ST_Direct;

#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__NODE_CUT;
    nonlocalExt = 0;
    communicator = nonlocCommunicator = NULL;
    commBuff = NULL;
#endif
}


LinearStatic :: ~LinearStatic()
{
    if ( stiffnessMatrix ) {
        delete stiffnessMatrix;
    }
    if ( nMethod ) {
        delete nMethod;
    }
}


NumericalMethod *LinearStatic :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }

    if ( isParallel() ) {
        if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
            nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
        }
    } else {
        nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
    }

    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed");
    }

    return nMethod;
}

IRResultType
LinearStatic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                               this->giveNumberOfProcesses(),
                                               this->commMode);
    }

#endif


    return IRRT_OK;
}


double LinearStatic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }
#endif

    if ( tStep != this->giveCurrentStep() ) {
        _error("giveUnknownComponent: unknown time step encountered");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:
    case VM_Incremental:
        if ( displacementVector.isNotEmpty() ) {
            return displacementVector.at(eq);
        } else {
            return 0.;
        }

    default:
        _error("giveUnknownComponent: Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *LinearStatic :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    //int mstep = 1;
    StateCounterType counter = 1;

    if ( previousStep != NULL ) {
        delete previousStep;
    }

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, ( double ) istep, 0., counter);
    // time and dt variables are set eq to 0 for staics - has no meaning
    return currentStep;
}


void LinearStatic :: solveYourself()
{
#ifdef __PARALLEL_MODE
    if ( this->isParallel() ) {
 #ifdef __VERBOSE_PARALLEL
        // force equation numbering before setting up comm maps
        int neq = this->giveNumberOfDomainEquations(EID_MomentumBalance);
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif

        // set up communication patterns
        // needed only for correct shared rection computation
        communicator->setUpCommunicationMaps(this, true);
        if ( nonlocalExt ) {
            nonlocCommunicator->setUpCommunicationMaps(this, true);
        }
    }
#endif

    StructuralEngngModel :: solveYourself();
}



void LinearStatic :: solveYourselfAt(TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif

        //
        // first step  assemble stiffness Matrix
        //
        stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( stiffnessMatrix == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );

        this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, StiffnessMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );

        initFlag = 0;
    }

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling load\n");
#endif

    //
    // allocate space for displacementVector
    //
    displacementVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    displacementVector.zero();

    //
    // assembling the load vector
    //
    loadVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    loadVector.zero();
    this->assembleVector( loadVector, tStep, EID_MomentumBalance, ExternalForcesVector, VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // internal forces (from Dirichlet b.c's, or thermal expansion, etc.)
    //
    FloatArray internalForces( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    internalForces.zero();
    this->assembleVector( internalForces, tStep, EID_MomentumBalance, InternalForcesVector, VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1) );

    loadVector.subtract(internalForces);

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetatStepumber() ) );

    //
    // call numerical model to solve arose problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("\n\nSolving ...\n\n");
#endif
    NM_Status s = nMethod->solve(stiffnessMatrix, & loadVector, & displacementVector);
    if ( !( s & NM_Success ) ) {
        OOFEM_ERROR("LinearStatic :: solverYourselfAt - No success in solving system.");
    }

    tStep->incrementStateCounter();            // update solution state counter
}


contextIOResultType LinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    FILE *file = NULL;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = StructuralEngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = displacementVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}


contextIOResultType LinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    int istep, iversion;
    FILE *file = NULL;

    this->resolveCorrespondingtStepumber(istep, iversion, obj);

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = StructuralEngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = displacementVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}


void
LinearStatic :: terminate(TimeStep *tStep)
{
    StructuralEngngModel :: terminate(tStep);
    this->printReactionForces(tStep, 1);
    fflush( this->giveOutputStream() );
}


void
LinearStatic :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}



void
LinearStatic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total);
}


#ifdef __PARALLEL_MODE
int
LinearStatic :: estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType)
{
    int mapSize = commMap.giveSize();
    int ndofs, count = 0, pcount = 0;
    IntArray locationArray;
    Domain *domain = this->giveDomain(1);
    DofManager *dman;
    Dof *jdof;

    if ( packUnpackType == ProblemCommMode__ELEMENT_CUT ) {
        for ( int i = 1; i <= mapSize; i++ ) {
            count += domain->giveDofManager( commMap.at(i) )->giveNumberOfDofs();
        }

        return ( buff.givePackSize(MPI_DOUBLE, 1) * count );
    } else if ( packUnpackType == ProblemCommMode__NODE_CUT ) {
        for ( int i = 1; i <= mapSize; i++ ) {
            ndofs = ( dman = domain->giveDofManager( commMap.at(i) ) )->giveNumberOfDofs();
            for ( int j = 1; j <= ndofs; j++ ) {
                jdof = dman->giveDof(j);
                if ( jdof->isPrimaryDof() && ( jdof->__giveEquationNumber() ) ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        // --------------------------------------------------------------------------------
        // only pcount is relevant here, since only prescribed components are exchanged !!!!
        // --------------------------------------------------------------------------------

        return ( buff.givePackSize(MPI_DOUBLE, 1) * pcount );
    } else if ( packUnpackType == ProblemCommMode__REMOTE_ELEMENT_MODE ) {
        for ( int i = 1; i <= mapSize; i++ ) {
            count += domain->giveElement( commMap.at(i) )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}
#endif
} // end namespace oofem
