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

#include "planstrss.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "structuralms.h"
#include "strainvector.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "connectivitytable.h"
 #include "rcm2.h"
#endif

namespace oofem {
REGISTER_Element(PlaneStress2d);

FEI2dQuadLin PlaneStress2d :: interpolation(1, 2);

PlaneStress2d :: PlaneStress2d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(),
    DirectErrorIndicatorRCInterface(), EIPrimaryUnknownMapperInterface(),
    HuertaErrorEstimatorInterface(), HuertaRemeshingCriteriaInterface()
    // Constructor.
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 4;
}

PlaneStress2d :: ~PlaneStress2d()
// Destructor
{ }

void
PlaneStress2d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    int i;
    FloatMatrix dnx;
#ifdef  PlaneStress2d_reducedShearIntegration
    FloatArray coord;
#endif

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, 8);
    answer.zero();

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);
    }

#ifdef  PlaneStress2d_reducedShearIntegration
    coord.resize(2);
    coord.zero();
    this->interpolation.evaldNdx( dnx, coord, FEIElementGeometryWrapper(this) );
#endif

    for ( i = 1; i <= 4; i++ ) {
        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }
}


void
PlaneStress2d :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
//
// Returns the [4x8] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
    FloatMatrix dnx;
#ifdef  PlaneStress2d_reducedShearIntegration
    FloatArray coord;
#endif

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(4, 8);

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);     // dv/dy -2
    }

#ifdef  PlaneStress2d_reducedShearIntegration
    coord.resize(2);
    coord.zero();
    this->interpolation.evaldNdx( dnx, coord, FEIElementGeometryWrapper(this) );
#endif

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(3, 2 * i - 1) = dnx.at(i, 2);     // du/dy -6
        answer.at(4, 2 * i - 0) = dnx.at(i, 1);     // dv/dx -9
    }
}

void
PlaneStress2d :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
PlaneStress2d :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    /*
     *
     * computes interpolation matrix for element edge.
     * we assemble locally this matrix for only nonzero
     * shape functions.
     * (for example only two nonzero shape functions for 2 dofs are
     * necessary for linear plane stress tringle edge).
     * These nonzero shape functions are then mapped to
     * global element functions.
     *
     * Using mapping technique will allow to assemble shape functions
     * without regarding particular side
     */

    FloatArray n(2);
    this->interpolation.edgeEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 4);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 3) = n.at(2);
    answer.at(2, 2) = n.at(1);
    answer.at(2, 4) = n.at(2);
}


void
PlaneStress2d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(4);
    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 3;
        answer.at(2) = 4;
        answer.at(3) = 5;
        answer.at(4) = 6;
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer.at(1) = 5;
        answer.at(2) = 6;
        answer.at(3) = 7;
        answer.at(4) = 8;
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer.at(1) = 7;
        answer.at(2) = 8;
        answer.at(3) = 1;
        answer.at(4) = 2;
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}

double
PlaneStress2d :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    return result *gp->giveWeight();
}

void
PlaneStress2d :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}


int
PlaneStress2d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    double dx, dy, length;
    Node *nodeA, *nodeB;
    IntArray edgeNodes(2);

    answer.resize(2, 2);
    answer.zero();

    this->interpolation.computeEdgeMapping(edgeNodes, dofManArray, iEdge);

    // edge nodes are global numbers, not local element numbers
    nodeA   = domain->giveNode( edgeNodes.at(1) );
    nodeB   = domain->giveNode( edgeNodes.at(2) );

    dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);

    answer.at(1, 1) = dx / length;
    answer.at(1, 2) = -dy / length;
    answer.at(2, 1) = answer.at(1, 2);
    answer.at(2, 2) = dx / length;

    return 1;
}



double
PlaneStress2d :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );


    weight      = gp->giveWeight();
    thickness   = this->giveCrossSection()->give(CS_Thickness, gp);
    volume      = determinant * weight * thickness;

    return volume;
}

IRResultType
PlaneStress2d :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    IRResultType result = this->NLStructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING1("Number of Gauss points enforced to 4");
    }

    return IRRT_OK;
}


double
PlaneStress2d :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
//
// returns receivers characteristic length in gp (for some material models)
// for crack formed in plane with normal normalToCrackPlane.
//
{
    // return this -> giveLenghtInDir(normalToCrackPlane) / sqrt ((double) gp->giveIntegrationRule()->giveNumberOfIntegrationPoints());
    if ( normalToCrackPlane.at(3) < 0.999999 ) { //ensure that characteristic length is in the plane of element
        return this->giveLenghtInDir(normalToCrackPlane);
    } else { //otherwise compute out-of-plane characteristic length from element area
        return DirectErrorIndicatorRCI_giveCharacteristicSize();
    }
}

double
PlaneStress2d :: giveCharacteristicSize(GaussPoint *gp, FloatArray &normalToCrackPlane, ElementCharSizeMethod method)
//
// returns receiver's characteristic size at gp (for some material models)
// for crack formed in plane with normal normalToCrackPlane
// using the selected method
//
{
    if ( method == ECSM_SquareRootOfArea ) {
        // square root of element area
        return DirectErrorIndicatorRCI_giveCharacteristicSize();
    }

    if ( method == ECSM_Projection ) {
        // standard projection method
        return this->giveCharacteristicLenght(gp, normalToCrackPlane);
    }

    // evaluate average strain and its maximum principal direction
    FloatArray sumstrain, averageNormal;
    IntegrationRule *iRule = giveDefaultIntegrationRulePtr();
    int nGP = iRule->giveNumberOfIntegrationPoints();
    for ( int i = 0; i < nGP; i++ ) {
        GaussPoint *gpi = iRule->getIntegrationPoint(i);
        StructuralMaterialStatus *matstatus = dynamic_cast< StructuralMaterialStatus * >( gpi->giveMaterialStatus() );
        if ( matstatus ) {
            sumstrain.add( matstatus->giveTempStrainVector() );
        }
    }

    StrainVector sumstrainvec(sumstrain, _PlaneStress);
    sumstrainvec.computeMaxPrincipalDir(averageNormal);

    if ( method == ECSM_ProjectionCentered ) {
        // projection method based on principal direction of average strain
        normalToCrackPlane = averageNormal;
        return this->giveLenghtInDir(normalToCrackPlane);
    }

    if ( method == ECSM_Oliver1 || method == ECSM_Oliver1modified ) {
        // method based on derivative of auxiliary function phi at each Gauss point
        // in the maximum principal strain direction determined at
        // ECSM_Oliver1 ... at each Gauss point
        // ECSM_Oliver1modified ... at element center (from average strain)

        // coordinates of the element center
        FloatArray center(2);
        double cx = 0., cy = 0.;
        for ( int i = 1; i <= 4; i++ ) {
            cx += giveNode(i)->giveCoordinate(1);
            cy += giveNode(i)->giveCoordinate(2);
        }

        cx /= 4.;
        cy /= 4.;

        // nodal values of function phi (0 or 1)
        FloatArray phi(4);
        for ( int i = 1; i <= 4; i++ ) {
            if ( ( ( giveNode(i)->giveCoordinate(1) - cx ) * averageNormal.at(1) + ( giveNode(i)->giveCoordinate(2) - cy ) * averageNormal.at(2) ) > 0. ) {
                phi.at(i) = 1.;
            } else {
                phi.at(i) = 0.;
            }
        }

        // gradient of function phi at the current GP
        FloatMatrix dnx;
        this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        FloatArray gradPhi(2);
        gradPhi.zero();
        for ( int i = 1; i <= 4; i++ ) {
            gradPhi.at(1) += phi.at(i) * dnx.at(i, 1);
            gradPhi.at(2) += phi.at(i) * dnx.at(i, 2);
        }

        // scalar product of the gradient with crack normal (at GP)
        double dPhidN = 0.;
        if ( method == ECSM_Oliver1modified ) {
            normalToCrackPlane = averageNormal;
        }

        for ( int i = 1; i <= 2; i++ ) {
            dPhidN += gradPhi.at(i) * normalToCrackPlane.at(i);
        }

        if ( dPhidN == 0. ) {
            _error("Zero value of dPhidN in PlaneStress2d :: giveCharacteristicSize\n");
        }

        return 1. / fabs(dPhidN);
    }

    _error("PlaneStress2d :: giveCharacteristicSize: invalid method");
    return 0.;
}

void
PlaneStress2d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}


Interface *
PlaneStress2d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == DirectErrorIndicatorRCInterfaceType ) {
        return static_cast< DirectErrorIndicatorRCInterface * >(this);
    } else if ( interface == EIPrimaryUnknownMapperInterfaceType ) {
        return static_cast< EIPrimaryUnknownMapperInterface * >(this);
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return static_cast< HuertaErrorEstimatorInterface * >(this);
    } else if ( interface == HuertaRemeshingCriteriaInterfaceType ) {
        return static_cast< HuertaRemeshingCriteriaInterface * >(this);
    }

    return NULL;
}


void
PlaneStress2d :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode)
{
    Element *element = this->HuertaErrorEstimatorI_giveElement();
    int inode, nodes = 4, iside, sides = 4, nd1, nd2;
    FloatArray *corner [ 4 ], midSide [ 4 ], midNode, cor [ 4 ];
    double x = 0.0, y = 0.0;

    static int sideNode [ 4 ] [ 2 ] = { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 1 } };

    if ( sMode == HuertaErrorEstimatorInterface :: NodeMode ||
        ( sMode == HuertaErrorEstimatorInterface :: BCMode && aMode == HuertaErrorEstimator :: HEE_linear ) ) {
        for ( inode = 0; inode < nodes; inode++ ) {
            corner [ inode ] = element->giveNode(inode + 1)->giveCoordinates();
            if ( corner [ inode ]->giveSize() != 3 ) {
                cor [ inode ].resize(3);
                cor [ inode ].at(1) = corner [ inode ]->at(1);
                cor [ inode ].at(2) = corner [ inode ]->at(2);
                cor [ inode ].at(3) = 0.0;

                corner [ inode ] = & ( cor [ inode ] );
            }

            x += corner [ inode ]->at(1);
            y += corner [ inode ]->at(2);
        }

        for ( iside = 0; iside < sides; iside++ ) {
            midSide [ iside ].resize(3);

            nd1 = sideNode [ iside ] [ 0 ] - 1;
            nd2 = sideNode [ iside ] [ 1 ] - 1;

            midSide [ iside ].at(1) = ( corner [ nd1 ]->at(1) + corner [ nd2 ]->at(1) ) / 2.0;
            midSide [ iside ].at(2) = ( corner [ nd1 ]->at(2) + corner [ nd2 ]->at(2) ) / 2.0;
            midSide [ iside ].at(3) = 0.0;
        }

        midNode.resize(3);

        midNode.at(1) = x / nodes;
        midNode.at(2) = y / nodes;
        midNode.at(3) = 0.0;
    }

    this->setupRefinedElementProblem2D(element, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midSide, midNode,
                                       localNodeId, localElemId, localBcId,
                                       controlNode, controlDof, aMode, "PlaneStress2d");
}


#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void PlaneStress2d :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void PlaneStress2d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    WCRec p [ 4 ];
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}



void PlaneStress2d :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v [ 4 ];
    double s [ 4 ], defScale;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }

        indx = context.giveIntVarIndx();

        for ( i = 1; i <= 4; i++ ) {
            s [ i - 1 ] = v [ i - 1 ].at(indx);
        }

        if ( context.getScalarAlgo() == SA_ISO_SURF ) {
            for ( i = 0; i < 4; i++ ) {
                if ( context.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = context.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = 0.;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    p [ i ].z = 0.;
                }
            }

            //EASValsSetColor(gc.getYieldPlotColor(ratio));
            context.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);

            /*  } else if (context.getScalarAlgo() == SA_ISO_LINE) {
             *
             * EASValsSetColor(context.getActiveCrackColor());
             * EASValsSetLineWidth(OOFEG_ISO_LINE_WIDTH);
             *
             * for (i=0; i< 4; i++) {
             * if (context.getInternalVarsDefGeoFlag()) {
             * // use deformed geometry
             * defScale = context.getDefScale();
             * p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,defScale);
             * p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,defScale);
             * p[i].z = 0.;
             *
             * } else {
             * p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
             * p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
             * p[i].z = 0.;
             * }
             * }
             *
             * // isoline implementation
             * oofeg_drawIsoLinesOnQuad (p, s);
             */
        } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
            double landScale = context.getLandScale();

            for ( i = 0; i < 4; i++ ) {
                if ( context.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = context.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = s [ i ] * landScale;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    p [ i ].z = s [ i ] * landScale;
                }

                // this fixes a bug in ELIXIR
                if ( fabs(s [ i ]) < 1.0e-6 ) {
                    s [ i ] = 1.0e-6;
                }
            }

            if ( context.getScalarAlgo() == SA_ZPROFILE ) {
                EASValsSetColor( context.getDeformedElementColor() );
                EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
                tr =  CreateQuad3D(p);
                EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            } else {
                context.updateFringeTableMinMax(s, 4);
                tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                EGWithMaskChangeAttributes(LAYER_MASK, tr);
            }

            EMAddGraphicsToModel(ESIModel(), tr);
        }
    } else if ( context.giveIntVarMode() == ISM_local ) {
        if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() != 4 ) {
            return;
        }

        int ip;
        GaussPoint *gp;
        IntArray ind(4);
        FloatArray *gpCoords;
        WCRec pp [ 9 ];

        for ( i = 0; i < 4; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                pp [ i ].z = 0.;
            } else {
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                pp [ i ].z = 0.;
            }
        }

        for ( i = 0; i < 3; i++ ) {
            pp [ i + 4 ].x = 0.5 * ( pp [ i ].x + pp [ i + 1 ].x );
            pp [ i + 4 ].y = 0.5 * ( pp [ i ].y + pp [ i + 1 ].y );
            pp [ i + 4 ].z = 0.5 * ( pp [ i ].z + pp [ i + 1 ].z );
        }

        pp [ 7 ].x = 0.5 * ( pp [ 3 ].x + pp [ 0 ].x );
        pp [ 7 ].y = 0.5 * ( pp [ 3 ].y + pp [ 0 ].y );
        pp [ 7 ].z = 0.5 * ( pp [ 3 ].z + pp [ 0 ].z );

        pp [ 8 ].x = 0.25 * ( pp [ 0 ].x + pp [ 1 ].x + pp [ 2 ].x + pp [ 3 ].x );
        pp [ 8 ].y = 0.25 * ( pp [ 0 ].y + pp [ 1 ].y + pp [ 2 ].y + pp [ 3 ].y );
        pp [ 8 ].z = 0.25 * ( pp [ 0 ].z + pp [ 1 ].z + pp [ 2 ].z + pp [ 3 ].z );

        for ( ip = 1; ip <= integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip - 1);
            gpCoords = gp->giveCoordinates();
            if ( ( gpCoords->at(1) > 0. ) && ( gpCoords->at(2) > 0. ) ) {
                ind.at(1) = 0;
                ind.at(2) = 4;
                ind.at(3) = 8;
                ind.at(4) = 7;
            } else if ( ( gpCoords->at(1) < 0. ) && ( gpCoords->at(2) > 0. ) ) {
                ind.at(1) = 4;
                ind.at(2) = 1;
                ind.at(3) = 5;
                ind.at(4) = 8;
            } else if ( ( gpCoords->at(1) < 0. ) && ( gpCoords->at(2) < 0. ) ) {
                ind.at(1) = 5;
                ind.at(2) = 2;
                ind.at(3) = 6;
                ind.at(4) = 8;
            } else {
                ind.at(1) = 6;
                ind.at(2) = 3;
                ind.at(3) = 7;
                ind.at(4) = 8;
            }

            if ( giveIPValue(v [ 0 ], gp, context.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            indx = context.giveIntVarIndx();

            for ( i = 1; i <= 4; i++ ) {
                s [ i - 1 ] = v [ 0 ].at(indx);
            }

            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ ind.at(i + 1) ].x;
                p [ i ].y = pp [ ind.at(i + 1) ].y;
                p [ i ].z = pp [ ind.at(i + 1) ].z;
            }

            context.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}




void
PlaneStress2d :: drawSpecial(oofegGraphicContext &gc)
{
    int i;
    WCRec l [ 2 ];
    GraphicObj *tr;
    GaussPoint *gp;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarType() == IST_CrackState ) {
        // ask if any active crack exist
        int igp, crackStatus;
        double ax, ay, bx, by, norm, xc, yc, length;
        FloatArray crackDir;
        FloatArray gpglobalcoords;

        for ( igp = 1; igp <= integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints(); igp++ ) {
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(igp - 1);

            if ( this->giveIPValue(cf, gp, IST_CrackedFlag, tStep) == 0 ) {
                return;
            }

            if ( ( int ) cf.at(1) == 0 ) {
                return;
            }

            if ( this->giveIPValue(crackDir, gp, IST_CrackDirs, tStep) ) {
                this->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);
                for ( i = 1; i <= 3; i++ ) {
                    crackStatus = ( int ) crackStatuses.at(i);
                    if ( ( crackStatus != pscm_NONE ) && ( crackStatus != pscm_CLOSED ) ) {
                        // draw a crack
                        // this element is 2d element in x-y plane
                        //
                        // compute perpendicular line to normal in xy plane
                        ax = crackDir.at(i);
                        ay = crackDir.at(3 + i);
                        if ( fabs(ax) > 1.e-6 ) {
                            by = 1.;
                            bx = -ay * by / ax;
                            norm = sqrt(bx * bx + by * by);
                            bx = bx / norm; // normalize to obtain unit vector
                            by = by / norm;
                        } else {
                            by = 0.0;
                            bx = 1.0;
                        }

                        // obtain gp global coordinates
                        if ( gc.getInternalVarsDefGeoFlag() ) {
                            double ksi, eta, n1, n2, n3, n4;
                            ksi = gp->giveCoordinate(1);
                            eta = gp->giveCoordinate(2);

                            n1 = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
                            n2 = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
                            n3 = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
                            n4 = ( 1. + ksi ) * ( 1. - eta ) * 0.25;

                            gpglobalcoords.resize(2);

                            gpglobalcoords.at(1) = ( n1 * this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                    n2 * this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                    n3 * this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                    n4 * this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale) );
                            gpglobalcoords.at(2) = ( n1 * this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                    n2 * this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                    n3 * this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                    n4 * this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale) );
                        } else {
                            computeGlobalCoordinates( gpglobalcoords, * ( gp->giveCoordinates() ) );
                        }

                        xc = gpglobalcoords.at(1);
                        yc = gpglobalcoords.at(2);
                        length = sqrt( computeVolumeAround(gp) / this->giveCrossSection()->give(CS_Thickness, gp) ) / 2.;
                        l [ 0 ].x = ( FPNum ) xc + bx * length;
                        l [ 0 ].y = ( FPNum ) yc + by * length;
                        l [ 0 ].z = 0.;
                        l [ 1 ].x = ( FPNum ) xc - bx * length;
                        l [ 1 ].y = ( FPNum ) yc - by * length;
                        l [ 1 ].z = 0.;
                        EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
                        EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
                        if ( ( crackStatus == pscm_SOFTENING ) || ( crackStatus == pscm_OPEN ) ) {
                            EASValsSetColor( gc.getActiveCrackColor() );
                        } else {
                            EASValsSetColor( gc.getCrackPatternColor() );
                        }

                        tr = CreateLine3D(l);
                        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
                        EMAddGraphicsToModel(ESIModel(), tr);
                    }
                }
            }
        }
    }
}

#endif

void
PlaneStress2d :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
PlaneStress2d :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 5; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        _error("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: node unknown");
    }
}

int
PlaneStress2d :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}


SPRPatchType
PlaneStress2d :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


int
PlaneStress2d :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    int result;
    FloatArray lcoords;
    result = this->computeLocalCoordinates(lcoords, coords);

    return result;
}

double
PlaneStress2d :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(2), gcoords;
    double dist;
    int size, gsize;

    lcoords.at(1) = lcoords.at(2) = 0.0;
    this->computeGlobalCoordinates(gcoords, lcoords);

    if ( ( size = coords.giveSize() ) < ( gsize = gcoords.giveSize() ) ) {
        _error("SpatialLocalizerI_giveDistanceFromParametricCenter: coordinates size mismatch");
    }

    if ( size == gsize ) {
        dist = coords.distance(gcoords);
    } else {
        FloatArray helpCoords = coords;

        helpCoords.resizeWithValues(gsize);
        dist = helpCoords.distance(gcoords);
    }

    return dist;
}


double
PlaneStress2d :: DirectErrorIndicatorRCI_giveCharacteristicSize()
{
    IntegrationRule *iRule;
    GaussPoint *gp;
    double volume = 0.0;

    iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        volume += this->computeVolumeAround(gp) / this->giveCrossSection()->give(CS_Thickness, gp);
    }

    return sqrt(volume);
}


int
PlaneStress2d :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                  TimeStep *tStep, const FloatArray &coords,
                                                                  FloatArray &answer)
{
    FloatArray lcoords, u;
    FloatMatrix n;
    FloatArray ni(4);
    int result;

    result = this->computeLocalCoordinates(lcoords, coords);

    this->interpolation.evalN( ni, lcoords, FEIElementGeometryWrapper(this) );

    n.resize(2, 8);
    n.zero();

    n.at(1, 1) = n.at(2, 2) = ni.at(1);
    n.at(1, 3) = n.at(2, 4) = ni.at(2);
    n.at(1, 5) = n.at(2, 6) = ni.at(3);
    n.at(1, 7) = n.at(2, 8) = ni.at(4);

    this->computeVectorOf(EID_MomentumBalance, mode, tStep, u);
    answer.beProductOf(n, u);

    return result;
}


void
PlaneStress2d :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
}
} // end namespace oofem
