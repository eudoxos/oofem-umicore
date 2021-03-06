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

#include "tr1bubblestokes.h"
#include "node.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fei2dtrlin.h"
#include "masterdof.h"
#include "fluidcrosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Tr1BubbleStokes);

// Set up interpolation coordinates
FEI2dTrLin Tr1BubbleStokes :: interp(1, 2);
// Set up ordering vectors (for assembling)
IntArray Tr1BubbleStokes :: momentum_ordering;
IntArray Tr1BubbleStokes :: conservation_ordering;
IntArray Tr1BubbleStokes :: edge_ordering [ 3 ] = {
    IntArray(4), IntArray(4), IntArray(4)
};
bool Tr1BubbleStokes :: __initialized = Tr1BubbleStokes :: initOrdering();

Tr1BubbleStokes :: Tr1BubbleStokes(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 3;
    this->numberOfGaussPoints = 7;

    this->bubble = new ElementDofManager(1, aDomain, this);
    this->bubble->appendDof( new MasterDof(1, this->bubble, V_u) );
    this->bubble->appendDof( new MasterDof(2, this->bubble, V_v) );
}

Tr1BubbleStokes :: ~Tr1BubbleStokes()
{
    delete this->bubble;
}

void Tr1BubbleStokes :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

int Tr1BubbleStokes :: computeNumberOfDofs()
{
    return 11;
}

void Tr1BubbleStokes :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( ut == EID_MomentumBalance ) {
        answer.setValues(2, V_u, V_v);
    } else if ( ut == EID_ConservationEquation ) {
        answer.setValues(1, P_f);
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        answer.setValues(3, V_u, V_v, P_f);
    } else {
        answer.resize(0);
    }
}

void Tr1BubbleStokes :: giveInternalDofManDofIDMask(int i, EquationID eid, IntArray &answer) const
{
    if ( eid == EID_MomentumBalance_ConservationEquation || eid == EID_MomentumBalance ) {
        answer.setValues(2, V_u, V_v);
    } else {
        answer.resize(0);
    }
}

double Tr1BubbleStokes :: computeVolumeAround(GaussPoint *gp)
{
    double detJ = fabs( this->interp.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ *gp->giveWeight();
}

void Tr1BubbleStokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                                 TimeStep *tStep)
{
    // Compute characteristic vector for this element. I.e the load vector(s)
    if ( mtrx == ExternalForcesVector ) {
        this->computeExternalForcesVector(answer, tStep);
    } else if ( mtrx == InternalForcesVector ) {
        this->computeInternalForcesVector(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}

void Tr1BubbleStokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                 CharType mtrx, TimeStep *tStep)
{
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}

void Tr1BubbleStokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatArray a_pressure, a_velocity, devStress, epsp, N, dNv(8);
    double r_vol, pressure;
    FloatMatrix dN, B(3, 8);
    B.zero();

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a_velocity);
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, a_pressure);

    FloatArray momentum, conservation;

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs( this->interp.evaldNdx( dN, * lcoords, FEIElementGeometryWrapper(this) ) );
        this->interp.evalN( N, * lcoords, FEIElementGeometryWrapper(this) );
        double dA = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < 3; j++, k += 2 ) {
            dNv(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }

        // Bubble contribution;
        dNv(6) = B(0, 6) = B(2, 7) = 27. * ( dN(0, 0) * N(1) * N(2) + N(0) * dN(1, 0) * N(2) + N(0) * N(1) * dN(2, 0) );
        dNv(7) = B(1, 7) = B(2, 6) = 27. * ( dN(0, 1) * N(1) * N(2) + N(0) * dN(1, 1) * N(2) + N(0) * N(1) * dN(2, 1) );

        pressure = N.dotProduct(a_pressure);
        epsp.beProductOf(B, a_velocity);

        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);

        momentum.plusProduct(B, devStress, dA);
        momentum.add(-pressure * dA, dNv);
        conservation.add(r_vol * dA, N);
    }

    answer.resize(11);
    answer.zero();
    answer.assemble(momentum, this->momentum_ordering);
    answer.assemble(conservation, this->conservation_ordering);
}

void Tr1BubbleStokes :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    int load_number, load_id;
    Load *load;
    bcGeomType ltype;
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.resize(11);
    answer.zero();
    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        load_number = this->boundaryLoadArray.at(2 * i - 1);
        load_id = this->boundaryLoadArray.at(2 * i);
        load = this->domain->giveLoad(load_number);
        ltype = load->giveBCGeoType();

        if ( ltype == EdgeLoadBGT ) {
            this->computeBoundaryLoadVector(vec, static_cast< BoundaryLoad * >(load), load_id, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeLoadVector(vec, load, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }
}


void Tr1BubbleStokes :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatArray N, gVector, temparray(8);

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(k);
            FloatArray *lcoords = gp->giveCoordinates();

            double rho = mat->give('d', gp);
            double detJ = fabs( this->interp.giveTransformationJacobian( * lcoords, FEIElementGeometryWrapper(this) ) );
            double dA = detJ * gp->giveWeight();

            this->interp.evalN( N, * lcoords, FEIElementGeometryWrapper(this) );
            for ( int j = 0; j < 3; j++ ) {
                temparray(2 * j)     += N(j) * rho * gVector(0) * dA;
                temparray(2 * j + 1) += N(j) * rho * gVector(1) * dA;
            }

            temparray(6) += N(0) * N(1) * N(2) * rho * gVector(0) * dA;
            temparray(7) += N(0) * N(1) * N(2) * rho * gVector(1) * dA;
        }
    }

    answer.resize(11);
    answer.zero();
    answer.assemble(temparray, this->momentum_ordering);
}

void Tr1BubbleStokes :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int iEdge, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad = static_cast< BoundaryLoad * >(load);

        int numberOfEdgeIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 2. ) / 2. );

        GaussIntegrationRule iRule(1, this, 1, 1);
        FloatArray N, t, f(4);
        IntArray edge_mapping;

        f.zero();
        iRule.SetUpPointsOnLine(numberOfEdgeIPs, _Unknown);

        for ( int i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            this->interp.edgeEvalN( N, iEdge, * lcoords, FEIElementGeometryWrapper(this) );
            double detJ = fabs( this->interp.boundaryGiveTransformationJacobian( iEdge, * lcoords, FEIElementGeometryWrapper(this) ) );
            double dS = gp->giveWeight() * detJ;

            if ( boundaryLoad->giveFormulationType() == Load :: FT_Entity ) { // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interp.boundaryLocal2Global( gcoords, iEdge, * lcoords, FEIElementGeometryWrapper(this) );
                boundaryLoad->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < 2; j++ ) {
                f(2 * j)     += N(j) * t(0) * dS;
                f(2 * j + 1) += N(j) * t(1) * dS;
            }
        }

        answer.resize(11);
        answer.zero();
        answer.assemble(f, this->edge_ordering [ iEdge - 1 ]);
    } else {
        OOFEM_ERROR("Tr1BubbleStokes :: computeEdgeBCSubVectorAt - Strange boundary condition type");
    }
}

void Tr1BubbleStokes :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // Note: Working with the components; [K, G+Dp; G^T+Dv^T, C] . [v,p]
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatMatrix B(3, 8), EdB, K, G, Dp, DvT, C, Ed, dN;
    FloatArray dNv(8), N, Ep, Cd, tmpA, tmpB;
    double Cp;
    B.zero();

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        // Compute Gauss point and determinant at current element
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs( this->interp.evaldNdx( dN, * lcoords, FEIElementGeometryWrapper(this) ) );
        double dA = detJ * gp->giveWeight();
        this->interp.evalN( N, * lcoords, FEIElementGeometryWrapper(this) );
        for ( int j = 0, k = 0; j < 3; j++, k += 2 ) {
            dNv(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }

        // Bubble contribution;
        dNv(6) = B(0, 6) = B(2, 7) = 27. * ( dN(0, 0) * N(1) * N(2) + N(0) * dN(1, 0) * N(2) + N(0) * N(1) * dN(2, 0) );
        dNv(7) = B(1, 7) = B(2, 6) = 27. * ( dN(0, 1) * N(1) * N(2) + N(0) * dN(1, 1) * N(2) + N(0) * N(1) * dN(2, 1) );

        // Computing the internal forces should have been done first.
        mat->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep); // dsigma_dev/deps_dev
        mat->giveDeviatoricPressureStiffness(Ep, TangentStiffness, gp, tStep); // dsigma_dev/dp
        mat->giveVolumetricDeviatoricStiffness(Cd, TangentStiffness, gp, tStep); // deps_vol/deps_dev
        mat->giveVolumetricPressureStiffness(Cp, TangentStiffness, gp, tStep); // deps_vol/dp

        EdB.beProductOf(Ed, B);
        K.plusProductSymmUpper(B, EdB, dA);
        G.plusDyadUnsym(dNv, N, -dA);
        C.plusDyadSymmUpper(N, Cp * dA);

        tmpA.beTProductOf(B, Ep);
        Dp.plusDyadUnsym(tmpA, N, dA);

        tmpB.beTProductOf(B, Cd);
        DvT.plusDyadUnsym(N, tmpB, dA);
    }

    K.symmetrized();
    C.symmetrized();

    FloatMatrix GTDvT, GDp;
    GTDvT.beTranspositionOf(G);
    GTDvT.add(DvT);
    GDp = G;
    GDp.add(Dp);

    answer.resize(11, 11);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(C, this->conservation_ordering);
}

FEInterpolation *Tr1BubbleStokes :: giveInterpolation() const
{
    return & interp;
}

FEInterpolation *Tr1BubbleStokes :: giveInterpolation(DofIDItem id) const
{
    return & interp;
}

void Tr1BubbleStokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tr1BubbleStokes :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case ZZNodalRecoveryModelInterfaceType:
        return static_cast< ZZNodalRecoveryModelInterface * >(this);

    case SpatialLocalizerInterfaceType:
        return static_cast< SpatialLocalizerInterface * >(this);

    case EIPrimaryUnknownMapperInterfaceType:
        return static_cast< EIPrimaryUnknownMapperInterface * >(this);

    default:
        return FMElement :: giveInterface(it);
    }
}

int Tr1BubbleStokes :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

void Tr1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                              TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin;
    this->interp.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->interp.evalN( n_lin, lcoords, FEIElementGeometryWrapper(this) );
    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer(0) += n.at(i) * this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(mode, tStep);
        answer(1) += n.at(i) * this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(mode, tStep);
    }

    for ( int i = 1; i <= n_lin.giveSize(); i++ ) {
        answer(2) += n_lin.at(i) * this->giveNode(i)->giveDofWithID(P_f)->giveUnknown(mode, tStep);
    }
}

int Tr1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode, TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    bool ok;
    FloatArray lcoords, n, n_lin;
    ok = this->computeLocalCoordinates(lcoords, gcoords);
    if ( !ok ) {
        answer.resize(0);
        return false;
    }

    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void Tr1BubbleStokes :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    answer.setValues(3, V_u, V_v, P_f);
}

double Tr1BubbleStokes :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray center;
    FloatArray lcoords;
    lcoords.setValues(3, 0.333333, 0.333333, 0.333333);
    this->interp.local2global( center, lcoords, FEIElementGeometryWrapper(this) );
    return center.distance(coords);
}
} // end namespace oofem
