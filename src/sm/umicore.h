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

#ifndef umicore_h
#define umicore_h


#include "idm1.h"

///@name Input fields for UmicoreMaterial
//@{
#define _IFT_UmicoreMaterial_Name "umicore"
#define _IFT_UmicoreMaterial_plastDmgCoeff "plastdmgcoeff"
#define _IFT_UmicoreMaterial_deactDmgInCompression "deactdmgincompression"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to UmicoreMaterial.
 * Stores the characteristic length of the element.
 */
class UmicoreMaterialStatus : public IsotropicDamageMaterial1Status
{
public:
    /// Constructor
    UmicoreMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~UmicoreMaterialStatus() { }

    FloatArray plasticStrain, tempPlasticStrain;
    double cumPlasticDeformation, tempCumPlasticDeformation;
    bool damageDeactivatedFlag, tempDamageDeactivatedFlag;
    bool elasticFlag, tempElasticFlag;


    // definition
    virtual const char *giveClassName() const { return "UmicoreMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);



};

/**
 * This class implements a simple local isotropic damage model for concrete in tension.
 * A model is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class UmicoreMaterial : public IsotropicDamageMaterial1
{
protected:
	/// interpolation coefficient between plasticity (value=1) and damage (value=0) model
	double plastDmgCoeff;

	/// if do or do not deactivate damage influence in compression
	bool deactDmgInCompression;


public:

    /// Constructor
    UmicoreMaterial(int n, Domain * d);
    /// Destructor
    virtual ~UmicoreMaterial();

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "UmicoreMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_UmicoreMaterial_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;



    virtual void giveRealStressVector(FloatArray &answer,  GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    void giveRealStressVectorDmg(FloatArray &answer,  GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    void giveRealStressVectorPlast(FloatArray &answer,  GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    double computeStressScalarForYieldFunction(const FloatArray& stress, GaussPoint *gp, TimeStep *tStep, const FloatMatrix& de, const FloatMatrix& ce);
    double computeStressScalarForYieldFunctionFromStrain(double lambda, const FloatArray& elasticStrain, const FloatArray& dStrain, const FloatMatrix& de, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    void stressReturn(FloatArray& stress, UmicoreMaterialStatus* status, double ft, GaussPoint *gp, TimeStep *tStep);

};
} // end namespace oofem
#endif // umicore_h
