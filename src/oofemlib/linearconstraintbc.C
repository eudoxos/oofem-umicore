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

#include "linearconstraintbc.h"
#include "classfactory.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"
#include "loadtimefunction.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
#include "node.h"

namespace oofem {
REGISTER_BoundaryCondition(LinearConstraintBC);


LinearConstraintBC :: LinearConstraintBC(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
    this->md = new Node(0, domain);
    // this is internal lagrange multiplier used to enforce the receiver constrain
    // this allocates a new equation related to this constraint
    this->md->appendDof( new MasterDof( 0, this->md, ( DofIDItem ) ( d->giveNextFreeDofID() ) ) );
    this->lhsType.resize(0);
    this->rhsType.resize(0);
}


IRResultType LinearConstraintBC :: initializeFrom(InputRecord *ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    const char *__proc = "initializeFrom";
    IRResultType result;
    rhsLtf = 0;

    IR_GIVE_FIELD(ir, weights, _IFT_LinearConstraintBC_weights);
    IR_GIVE_FIELD(ir, rhs, _IFT_LinearConstraintBC_rhs);
    IR_GIVE_FIELD(ir, dofmans, _IFT_LinearConstraintBC_dofmans);
    IR_GIVE_FIELD(ir, dofs, _IFT_LinearConstraintBC_dofs);
    if ( weights.giveSize() != dofmans.giveSize() ) {
        OOFEM_ERROR3( "Size mismatch, weights %d and dofmans %d", weights.giveSize(), dofmans.giveSize() );
    }
    IR_GIVE_OPTIONAL_FIELD(ir, weightsLtf, _IFT_LinearConstraintBC_weightsltf);
    IR_GIVE_OPTIONAL_FIELD(ir, rhsLtf, _IFT_LinearConstraintBC_rhsltf);

    IR_GIVE_FIELD(ir, lhsType, _IFT_LinearConstraintBC_lhstype);
    IR_GIVE_FIELD(ir, rhsType, _IFT_LinearConstraintBC_rhstype);

    return IRRT_OK;
}


void LinearConstraintBC :: giveLocArray(const UnknownNumberingScheme &r_s,  IntArray &locr, int &lambda_eq)
{
    int size = this->weights.giveSize();
    Dof *idof;

    locr.resize(size);
    // assemble location array
    for ( int _i = 1; _i <= size; _i++ ) {
        idof = this->domain->giveDofManager( this->dofmans.at(_i) )->giveDof( this->dofs.at(_i) );
        locr.at(_i) = r_s.giveDofEquationNumber(idof);
    }

    lambda_eq = r_s.giveDofEquationNumber( md->giveDof(1) );
}


void LinearConstraintBC :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                                    CharType type, const UnknownNumberingScheme &r_s,
                                    const UnknownNumberingScheme &c_s)
{
    int size = this->weights.giveSize();
    IntArray lambdaeq(1);
    FloatMatrix contrib(size, 1), contribt;
    IntArray locr(size), locc(size);

    if ( !this->lhsType.contains( ( int ) type ) ) {
        return;
    }
    this->giveLocArray( r_s, locr, lambdaeq.at(1) );

    if ( this->isImposed(tStep) ) {
        for ( int _i = 1; _i <= size; _i++ ) { // loop over dofs
            double factor = 1.;
            if ( weightsLtf.giveSize() ) {
                factor = domain->giveLoadTimeFunction( weightsLtf.at(_i) )->__at( tStep->giveIntrinsicTime() );
            }
            contrib.at(_i, 1) = this->weights.at(_i) * factor;
        }
        contribt.beTranspositionOf(contrib);

        answer->assemble(lambdaeq, locr, contribt);
        answer->assemble(locr, lambdaeq, contrib);
    } else {
        // the bc is not imposed at specific time step, however in order to make the equation system regular
        // we initialize the allocated equation to the following form 1*labmda = 0, forcing lagrange multiplier
        // of inactive condition to be zero.
        FloatMatrix help(1, 1);
        help.at(1, 1) = 1.0;
        answer->assemble(lambdaeq, lambdaeq, help);
    }
}

void LinearConstraintBC :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                          CharType type, ValueModeType mode,
                                          const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    IntArray loc, lambdaeq(1);
    FloatArray vec(1);
    double factor = 1.;

    if ( !this->rhsType.contains( ( int ) type ) ) {
        return;
    }
    if ( !this->isImposed(tStep) ) {
        return;
    }

    if ( type == InternalForcesVector ) {
        // compute true residual
        int size = this->weights.giveSize();
        Dof *idof;

        // assemble location array
        for ( int _i = 1; _i <= size; _i++ ) {
            factor = 1.;
            if ( weightsLtf.giveSize() ) {
                factor = domain->giveLoadTimeFunction( weightsLtf.at(_i) )->__at( tStep->giveIntrinsicTime() );
            }
            idof = this->domain->giveDofManager( this->dofmans.at(_i) )->giveDof( this->dofs.at(_i) );
            if ( s.giveDofEquationNumber(idof) ) {
                answer.at( s.giveDofEquationNumber(idof) ) += md->giveDof(1)->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
            if ( s.giveDofEquationNumber( md->giveDof(1) ) ) {
                answer.at( s.giveDofEquationNumber( md->giveDof(1) ) ) += idof->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
        }
    } else {
        // use rhs value

        if ( rhsLtf ) {
            factor = domain->giveLoadTimeFunction(rhsLtf)->__at( tStep->giveIntrinsicTime() );
        }
        this->giveLocArray( s, loc, lambdaeq.at(1) );
        vec.at(1) = rhs * factor;
        answer.assemble(vec, lambdaeq);
    }
}

void LinearConstraintBC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) {
    rows.resize(3);
    cols.resize(3);

    IntArray loc, lambdaeq(1);
    this->giveLocArray( r_s, loc, lambdaeq.at(1) );
    // column block
    rows [ 0 ] = loc;
    cols [ 0 ] = lambdaeq;
    // row block
    cols [ 1 ] = loc;
    rows [ 1 ] = lambdaeq;
    // diagonal enry (some sparse mtrx implementation requaire this)
    rows [ 2 ] = lambdaeq;
    cols [ 2 ] = lambdaeq;
}


contextIOResultType
LinearConstraintBC :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsLtf.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofmans.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream->write(& rhs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream->write(& rhsLtf, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
LinearConstraintBC :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsLtf.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofmans.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream->read(& rhs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream->read(& rhsLtf, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} //end of oofem namespace
