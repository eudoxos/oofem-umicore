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

#ifndef tr21_2d_supg_h
#define tr21_2d_supg_h

#include "supgelement2.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "fei2dtrquad.h"
#include "fei2dtrlin.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "leplic.h"
#include "levelsetpcs.h"

#define _IFT_TR21_2D_SUPG_Name "tr21supg"

namespace oofem {
/**
 * Class representing 2d triangular element  with quadratic velocity
 * and linear pressure approximation for solving incompressible fluid problems
 * with SUPG solver.
 */
class TR21_2D_SUPG : public SUPGElement2, public LevelSetPCSElementInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dTrQuad velocityInterpolation;
    static FEI2dTrLin pressureInterpolation;
    IntArray pressureDofManArray;

public:
    TR21_2D_SUPG(int n, Domain * aDomain);
    virtual ~TR21_2D_SUPG();

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    // definition
    virtual const char *giveClassName() const { return "TR21_2D_SUPG"; }
    virtual const char *giveInputRecordName() const { return _IFT_TR21_2D_SUPG_Name; }
    virtual MaterialMode giveMaterialMode() { return _2dFlow; }

    virtual void giveElementDofIDMask(EquationID, IntArray &answer) const;
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int computeNumberOfDofs();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void updateYourself(TimeStep *tStep);
    /// Used to check consistency and initialize some element geometry data (area,b,c).
    virtual int checkConsistency();

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual double LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep);
    virtual void LS_PCS_computedN(FloatMatrix &answer);
    virtual double LS_PCS_computeVolume();
    virtual void LS_PCS_computeVolume(double &answer,  const FloatArray **coordinates);
    virtual double LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep);
    virtual void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi);


    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }


    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                           InternalStateType type, TimeStep *tStep);

    /// @name Helping functions for computing VOFFractions.
    //@{
    void computeIntersection(int iedge, FloatArray &intcoords, FloatArray &fi);
    void computeMiddlePointOnParabolicArc(FloatArray &answer, int iedge, FloatArray borderpoints);
    void computeCenterOf(FloatArray &C, FloatArray c, int dim);
    void computeQuadraticRoots(FloatArray Coeff, double &r1, double &r2);
    void computeCoordsOfEdge(FloatArray &answer, int iedge);
    void computeQuadraticFunct(FloatArray &answer, int iedge);
    void computeQuadraticFunct(FloatArray &answer, FloatArray line);
    //@{

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep);
    // Graphics output
    //void drawYourself(oofegGraphicContext&);
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawScalar(oofegGraphicContext &context);
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

    virtual double computeCriticalTimeStep(TimeStep *tStep);

    // three terms for computing their norms due to computing t_supg
    virtual void computeAdvectionTerm(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeAdvectionDeltaTerm(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassDeltaTerm(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeLSICTerm(FloatMatrix &answer, TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType);

protected:
    virtual void giveLocalVelocityDofMap(IntArray &map);
    virtual void giveLocalPressureDofMap(IntArray &map);

    virtual void computeGaussPoints();
    virtual void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);
    virtual void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp);
    virtual void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);
    virtual void computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);
    virtual int  giveNumberOfSpatialDimensions();
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void initGeometry();

    virtual void updateStabilizationCoeffs(TimeStep *tStep);
};
} // end namespace oofem
#endif // tr21_2d_supg_h
