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

#include "petscsolver.h"

#include "petscsparsemtrx.h"
#include "engngm.h"
#include "floatarray.h"
#include "verbose.h"
#include "timer.h"
#include "error.h"
#include "classfactory.h"

#include <petscksp.h>

namespace oofem {
REGISTER_SparseLinSolver(PetscSolver, ST_Petsc);

PetscSolver :: PetscSolver(Domain *d, EngngModel *m) : SparseLinearSystemNM(d, m) { }

PetscSolver :: ~PetscSolver() { }

NM_Status PetscSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x)
{
    int neqs;

    // first check whether Lhs is defined
    if ( !A ) {
        OOFEM_ERROR("PETScSolver :: solve: unknown Lhs");
    }

    // and whether Rhs
    if ( !b ) {
        OOFEM_ERROR("PETScSolver :: solve: unknown Rhs");
    }

    // and whether previous Solution exist
    if ( !x ) {
        OOFEM_ERROR("PETScSolver :: solve: unknown solution array");
    }

    if ( x->giveSize() != ( neqs = b->giveSize() ) ) {
        OOFEM_ERROR("PETScSolver :: solve: size mismatch");
    }

    if ( A->giveType() != SMT_PetscMtrx ) {
        OOFEM_ERROR("PETScSolver :: solve: PetscSparseMtrx Expected");
    }

    PetscSparseMtrx *Lhs = ( PetscSparseMtrx * ) A;

    Vec globRhsVec;
    Vec globSolVec;

    // "Parallel" context automatically uses sequential alternative if the engineering problem is sequential.
    PetscContext *context = engngModel->givePetscContext( Lhs->giveDomainIndex() );

    /*
     * scatter and gather rhs to global representation
     */
    context->createVecGlobal(& globRhsVec);
    context->scatter2G(b, globRhsVec, ADD_VALUES);

    VecDuplicate(globRhsVec, & globSolVec);

    //VecView(globRhsVec,PETSC_VIEWER_STDOUT_WORLD);
    NM_Status s = this->petsc_solve(Lhs, globRhsVec, globSolVec);
    //VecView(globSolVec,PETSC_VIEWER_STDOUT_WORLD);

    context->scatterG2N(globSolVec, x, INSERT_VALUES);

    VecDestroy(& globSolVec);
    VecDestroy(& globRhsVec);

    return s;
}

NM_Status
PetscSolver :: petsc_solve(PetscSparseMtrx *Lhs, Vec b, Vec x)
{
    int nite;
    PetscErrorCode err;
    KSPConvergedReason reason;
    if ( Lhs->giveType() != SMT_PetscMtrx ) {
        OOFEM_ERROR("PETScSolver :: petsc_solve: PetscSparseMtrx Expected");
    }

    Timer timer;
    timer.startTimer();

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *  Create the linear solver and set various options
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if ( !Lhs->kspInit ) {
#ifdef __PARALLEL_MODE
        MPI_Comm comm = engngModel->giveParallelComm();
#else
        MPI_Comm comm = PETSC_COMM_SELF;
#endif
        KSPCreate(comm, & Lhs->ksp);
        Lhs->kspInit = true;
    }

    /*
     * Set operators. Here the matrix that defines the linear system
     * also serves as the preconditioning matrix.
     */
    if ( Lhs->newValues ) { // Optimization for successive solves
        ///@todo I'm not 100% on the choice MatStructure. SAME_NONZERO_PATTERN should be safe.
        if ( this->engngModel->requiresUnknownsDictionaryUpdate() ) {
            KSPSetOperators(Lhs->ksp, * Lhs->giveMtrx(), * Lhs->giveMtrx(), DIFFERENT_NONZERO_PATTERN);
        } else {
            //KSPSetOperators(Lhs->ksp, * Lhs->giveMtrx(), * Lhs->giveMtrx(), SAME_PRECONDITIONER);
            KSPSetOperators(Lhs->ksp, * Lhs->giveMtrx(), * Lhs->giveMtrx(), SAME_NONZERO_PATTERN);
        }
        Lhs->newValues = false;

        /*
         * Set linear solver defaults for this problem (optional).
         * - The following two statements are optional; all of these
         * parameters could alternatively be specified at runtime via
         * KSPSetFromOptions().  All of these defaults can be
         * overridden at runtime, as indicated below.
         */
        KSPSetTolerances(Lhs->ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

        /*
         * Set runtime options, e.g.,
         * -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
         * These options will override those specified above as long as
         * KSPSetFromOptions() is called _after_ any other customization
         * routines.
         */
        KSPSetFromOptions(Lhs->ksp);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *  Solve the linear system
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //MatView(*Lhs->giveMtrx(),PETSC_VIEWER_STDOUT_SELF);
    err = KSPSolve(Lhs->ksp, b, x);
    if ( err != 0 ) {
        OOFEM_ERROR2("PetscSolver:  Error when solving: %d\n", err);
    }
    KSPGetConvergedReason(Lhs->ksp, & reason);
    KSPGetIterationNumber(Lhs->ksp, & nite);

    if ( reason >= 0 ) {
        //OOFEM_LOG_INFO("PetscSolver:  Converged. KSPConvergedReason: %d, number of iterations: %d\n", reason, nite);
    } else {
        OOFEM_WARNING3("PetscSolver:  Diverged! KSPConvergedReason: %d, number of iterations: %d\n", reason, nite);
    }

    timer.stopTimer();
    OOFEM_LOG_INFO( "PetscSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );

    if ( reason < 0 ) {
        return NM_NoSuccess;
    } else {
        return NM_Success;
    }
}


//#ifndef __PARALLEL_MODE
#if 0
///@todo Parallel mode of this.
NM_Status PetscSolver :: solve(SparseMtrx *A, FloatMatrix &B, FloatMatrix &X)
{
    if ( !A ) {
        OOFEM_ERROR("PETScSolver :: solve: Unknown Lhs");
    }
    if ( A->giveType() != SMT_PetscMtrx ) {
        OOFEM_ERROR("PETScSolver :: solve: PetscSparseMtrx Expected");
    }

    PetscSparseMtrx *Lhs = ( PetscSparseMtrx * ) A;

    Vec globRhsVec;
    Vec globSolVec;

    bool newLhs = true;
    int rows = B.giveNumberOfRows();
    int cols = B.giveNumberOfColumns();
    NM_Status s;
    X.resize(rows, cols);
    double *Xptr = X.givePointer();

    for ( int i = 0; i < cols; ++i ) {
        VecCreateSeqWithArray(PETSC_COMM_SELF, rows, B.givePointer() + rows * i, & globRhsVec);
        VecDuplicate(globRhsVec, & globSolVec);
        s = this->petsc_solve(Lhs, globRhsVec, globSolVec, newLhs);
        if ( !( s & NM_Success ) ) {
            OOFEM_WARNING2("PetscSolver :: solve - No success at solving column %d", i + 1);
            return s;
        }
        newLhs = false;
        double *ptr;
        VecGetArray(globSolVec, & ptr);
        for ( int j = 0; j < rows; ++j ) {
            Xptr [ j + rows * i ] = ptr [ j ];
        }
        VecRestoreArray(globSolVec, & ptr);
    }
    VecDestroy(globSolVec);
    VecDestroy(globRhsVec);
    return s;
}
#endif
} // end namespace oofem
