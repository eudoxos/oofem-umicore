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

#ifndef cct_h
#define cct_h

#include "mathfem.h"
#include "nlstructuralelement.h"
#include "layeredcrosssection.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "fei2dtrlin.h"
#include "zzerrorestimator.h"

#define _IFT_CCTPlate_Name "cctplate"

namespace oofem {
/**
 * This class implements an triangular three-node plate CCT finite element.
 * Each node has 3 degrees of freedom.
 *
 * Tasks:
 * - calculating its B,D,N matrices and dV.
 */
class CCTPlate : public NLStructuralElement,
public LayeredCrossSectionInterface, public ZZNodalRecoveryModelInterface,
public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface
{
protected:
    static FEI2dTrLin interp_lin;
    //static FEI2dTrRot interp_rot;

    double area;

public:
    CCTPlate(int n, Domain * d);
    virtual ~CCTPlate() { }

    virtual FEInterpolation *giveInterpolation() const { return & interp_lin; }
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    virtual MaterialMode giveMaterialMode()  { return _2dPlate; }
    virtual int giveApproxOrder() { return 2; }
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

protected:
    virtual void computeGaussPoints();
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);


    virtual void giveNodeCoordinates(double &x1, double &x2, double &x3,
                                     double &y1, double &y2, double &y3,
                                     double *z = NULL);


    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    //virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *gp) { answer.resize(0, 0); }
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    //virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const { answer.resize(0); }
    //virtual IntegrationRule *GetSurfaceIntegrationRule(int i) { return NULL; }
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    //virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf) { return 0.; }
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    //virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf) { answer.resize(0); }
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);

public:
    // definition & identification
    virtual const char *giveClassName() const { return "CCTPlate"; }
    virtual const char *giveInputRecordName() const { return _IFT_CCTPlate_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs() { return 9; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp);

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    virtual Interface *giveInterface(InterfaceType it);

    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                           InternalStateType type, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP() { return 1; }
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();
    // ZZErrorEstimatorInterface
    virtual Element *ZZErrorEstimatorI_giveElement() { return this; }

    // ZZRemeshingCriteriaInterface
    virtual double ZZRemeshingCriteriaI_giveCharacteristicSize();
    virtual int ZZRemeshingCriteriaI_givePolynOrder() { return 1; };

    // layered cross section support functions
    virtual void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                            GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &context);
#endif
};
} // end namespace oofem
#endif // cct_h
