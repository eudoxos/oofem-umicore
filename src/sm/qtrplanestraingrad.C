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

#include "qtrplanestraingrad.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "crosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
QTrPlaneStrainGrad :: QTrPlaneStrainGrad(int n, Domain *aDomain) : QTrPlaneStrain(n, aDomain), GradDpElement()
    // Constructor.
{
    nPrimNodes = 6;
    nPrimVars = 2;
    nSecNodes = 3;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
}


void
QTrPlaneStrainGrad :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( inode <= nSecNodes ) {
        answer.setValues(3, D_u, D_v, G_0);
    } else {
        answer.setValues(2, D_u, D_v);
    }
}

IRResultType
QTrPlaneStrainGrad :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                 // Required by IR_GIVE_FIELD macro
    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    //IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, _IFT_Element_nip);

    return IRRT_OK;
}

void
QTrPlaneStrainGrad :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
QTrPlaneStrainGrad :: computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    double l1, l2, l3;

    l1 = gp->giveCoordinate(1);
    l2 = gp->giveCoordinate(2);
    l3 = 1.0 - l1 - l2;

    answer.resize(1, 3);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 2) = l2;
    answer.at(1, 3) = l3;
}

void
QTrPlaneStrainGrad :: computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [2x6] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3, area;

    node1 = this->giveNode(1);
    node2 = this->giveNode(2);
    node3 = this->giveNode(3);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.resize(2, 3);

    answer.at(1, 1) = y2 - y3;
    answer.at(1, 2) = y3 - y1;
    answer.at(1, 3) = y1 - y2;
    answer.at(2, 1) = x3 - x2;
    answer.at(2, 2) = x1 - x3;
    answer.at(2, 3) = x2 - x1;


    answer.times( 1. / ( 2. * area ) );
}
}
