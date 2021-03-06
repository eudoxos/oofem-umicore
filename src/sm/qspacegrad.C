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

#include "qspacegrad.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QSpaceGrad);

FEI3dHexaLin QSpaceGrad :: interpolation;

QSpaceGrad :: QSpaceGrad(int n, Domain *aDomain) :  QSpace(n, aDomain), GradDpElement()
    // Constructor.
{
    nPrimNodes = 8;
    nPrimVars = 2;
    nSecNodes = 4;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
}


IRResultType
QSpaceGrad :: initializeFrom(InputRecord *ir)
{
    IRResultType result = this->StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( ( numberOfGaussPoints != 8 ) && ( numberOfGaussPoints != 14 ) && ( numberOfGaussPoints != 27 ) && ( numberOfGaussPoints != 64 ) ) {
        numberOfGaussPoints = 27;
    }

    return IRRT_OK;
}


void
QSpaceGrad :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( inode <= nSecNodes ) {
        answer.setValues(4, D_u, D_v, D_w, G_0);
    } else {
        answer.setValues(3, D_u, D_v, D_w);
    }
}


void
QSpaceGrad :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 7);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
}



void
QSpaceGrad :: computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n(8);
    this->interpolation.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize(1, 8);
    answer.zero();

    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, i) = n.at(i);
    }
}

void
QSpaceGrad :: computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;
    IntArray a(8);
    for ( int i = 1; i < 9; i++ ) {
        a.at(i) = dofManArray.at(i);
    }

    answer.resize(3, 8);
    answer.zero();

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, i) = dnx.at(i, 1);
        answer.at(2, i) = dnx.at(i, 2);
        answer.at(3, i) = dnx.at(i, 3);
    }
}


void
QSpaceGrad :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int i)
// Returns the [60x60] nonlinear part of strain-displacement matrix {B} of the receiver,
// evaluated at gp

{
    FloatMatrix dnx;

    // compute the derivatives of shape functions
    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(60, 60);
    answer.zero();

    // put the products of derivatives of shape functions into the "nonlinear B matrix",
    // depending on parameter i, which is the number of the strain component
    if ( i <= 3 ) {
        for ( int k = 0; k < 20; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 60; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, i) * dnx.at( ( j - 1 ) / 3 + 1, i );
                }
            }
        }
    } else if ( i == 4 ) {
        for ( int k = 0; k < 20; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 60; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 2) * dnx.at( ( j - 1 ) / 3 + 1, 3 ) + dnx.at(k + 1, 3) * dnx.at( ( j - 1 ) / 3 + 1, 2 );
                }
            }
        }
    } else if ( i == 5 ) {
        for ( int k = 0; k < 20; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 60; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 1) * dnx.at( ( j - 1 ) / 3 + 1, 3 ) + dnx.at(k + 1, 3) * dnx.at( ( j - 1 ) / 3 + 1, 1 );
                }
            }
        }
    } else if ( i == 6 ) {
        for ( int k = 0; k < 20; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 60; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 1) * dnx.at( ( j - 1 ) / 3 + 1, 2 ) + dnx.at(k + 1, 2) * dnx.at( ( j - 1 ) / 3 + 1, 1 );
                }
            }
        }
    }

    return;
}
}
