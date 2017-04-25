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

/**
 * Select the mapping algorithm. The IDM_USE_MMAShapeFunctProjection does not work, since
 * this mapper does not preserve the max. property of damage and equivalent strain.
 */
#define IDM_USE_MMAContainingElementProjection

/*
 * Selects the use of mapped strain or projected strain from element.
 */
//#define IDM_USE_MAPPEDSTRAIN

#include "material.h"
#include "linearelasticmaterial.h"
#include "isodamagemodel.h"
#include "structuralms.h"
#include "randommaterialext.h"
#include "materialmapperinterface.h"

#ifdef IDM_USE_MMAClosestIPTransfer
 #include "mmaclosestiptransfer.h"
#endif

#ifdef IDM_USE_MMAContainingElementProjection
 #include "mmacontainingelementprojection.h"
#endif

#ifdef IDM_USE_MMAShapeFunctProjection
 #include "mmashapefunctprojection.h"
#endif

#ifdef IDM_USE_MMALeastSquareProjection
 #include "mmaleastsquareprojection.h"
#endif

///@name Input fields for UmicoreMaterial
//@{
#define _IFT_UmicoreMaterial_Name "umicore"
#define _IFT_UmicoreMaterial_e0 "e0"
#define _IFT_UmicoreMaterial_ef "ef"
#define _IFT_UmicoreMaterial_wf "wf"
#define _IFT_UmicoreMaterial_equivstraintype "equivstraintype"
#define _IFT_UmicoreMaterial_damageLaw "damlaw"
#define _IFT_UmicoreMaterial_k "k"
#define _IFT_UmicoreMaterial_md "md"
#define _IFT_UmicoreMaterial_ecsm "ecsm"
#define _IFT_UmicoreMaterial_At "at"
#define _IFT_UmicoreMaterial_Bt "bt"
#define _IFT_UmicoreMaterial_ft "ft"
#define _IFT_UmicoreMaterial_w1wf "w1wf"
#define _IFT_UmicoreMaterial_e1ef "e1ef"
#define _IFT_UmicoreMaterial_s1ft "s1ft"
#define _IFT_UmicoreMaterial_s1 "s1"
#define _IFT_UmicoreMaterial_w1 "w1"
#define _IFT_UmicoreMaterial_e1 "e1"
#define _IFT_UmicoreMaterial_ek "ek"
#define _IFT_UmicoreMaterial_gf "gf"
#define _IFT_UmicoreMaterial_gft "gft"
#define _IFT_UmicoreMaterial_ep "ep"
#define _IFT_UmicoreMaterial_e2 "e2"
#define _IFT_UmicoreMaterial_nd "nd"
#define _IFT_UmicoreMaterial_checkSnapBack "checksnapback"
#define _IFT_UmicoreMaterial_n "griff_n"
#define _IFT_UmicoreMaterial_plastDmgCoeff "plastdmgcoeff"
#define _IFT_UmicoreMaterial_deactDmgInCompression "deactdmgincompression"
//@}

namespace oofem {
#define IDM1_ITERATION_LIMIT 1.e-9

/**
 * This class implements associated Material Status to UmicoreMaterial.
 * Stores the characteristic length of the element.
 */
class UmicoreMaterialStatus : public IsotropicDamageMaterialStatus, public RandomMaterialStatusExtensionInterface
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

    virtual Interface *giveInterface(InterfaceType it);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};

/**
 * This class implements a simple local isotropic damage model for concrete in tension.
 * A model is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class UmicoreMaterial : public IsotropicDamageMaterial,
public RandomMaterialExtensionInterface,
public MaterialModelMapperInterface
{
protected:
	/// interpolation coefficient between plasticity (value=1) and damage (value=0) model
	double plastDmgCoeff;

	/// if do or do not deactivate damage influence in compression
	bool deactDmgInCompression;

    /// Equivalent strain at stress peak (or a similar parameter).
    double e0;
    /// Determines ductility -> corresponds to fracturing strain.
    double ef;
    /// Determines ductility -> corresponds to crack opening in the cohesive crack model.
    double wf;

    /**
     * Determines the softening -> corresponds to the initial fracture energy. For a linear law, it is the area
     * under the stress/strain curve. For an exponential law, it is the area bounded by the elastic range
     * and a tangent to the softening part of the curve at the peak stress. For a bilinear law,
     * gf corresponds to area bounded by elasticity and the first linear softening line projected to zero stress.
     */
    double gf;

    /// Determines the softening for the bilinear law -> corresponds to the total fracture energy.
    double gft;

    /// Determines the softening for the bilinear law -> corresponds to the strain at the knee point.
    double ek;

    /// Type characterizing the algorithm used to compute equivalent strain measure.
    enum EquivStrainType { EST_Unknown, EST_Mazars, EST_Rankine_Smooth, EST_Rankine_Standard, EST_ElasticEnergy, EST_ElasticEnergyPositiveStress, EST_ElasticEnergyPositiveStrain, EST_Mises, EST_Griffith };
    /// Parameter specifying the definition of equivalent strain.
    EquivStrainType equivStrainType;

    /// Parameter used in Mises definition of equivalent strain.
    double k;

    /// Parameter used in Griffith's criterion
    double griff_n;

    /// Remporary parameter reading type of softening law, used in other isotropic damage material models.
    int damageLaw;

    /** Type characterizing the formula for the damage law. For example, linear softening can be specified
     *   with fracturing strain or crack opening.
     */
    enum SofteningType { ST_Unknown, ST_Exponential, ST_Linear, ST_Mazars, ST_Smooth, ST_SmoothExtended, ST_Exponential_Cohesive_Crack, ST_Linear_Cohesive_Crack, ST_BiLinear_Cohesive_Crack, ST_Disable_Damage };

    /// Parameter specifying the type of softening (damage law).
    SofteningType softType;

    /// Parameters used in Mazars damage law.
    double At, Bt;
    /// Parameter used in "smooth damage law".
    double md;

    /// Parameters used if softType = 7 (extended smooth damage law)
    double e1, e2, s1, nd;
    /// Check possible snap back flag
    int checkSnapBack;

    /// Method used for evaluation of characteristic element size
    ElementCharSizeMethod ecsMethod;

#ifdef IDM_USE_MMAClosestIPTransfer
    /// Mapper used to map internal variables in adaptivity.
    static MMAClosestIPTransfer mapper;
#endif
#ifdef IDM_USE_MMAContainingElementProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMAContainingElementProjection mapper;
#endif
#ifdef IDM_USE_MMAShapeFunctProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMAShapeFunctProjection mapper;
#endif
#ifdef IDM_USE_MMALeastSquareProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMALeastSquareProjection mapper;
#endif

public:

    /// Constructor
    UmicoreMaterial(int n, Domain * d);
    /// Destructor
    virtual ~UmicoreMaterial();

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "UmicoreMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_UmicoreMaterial_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    /**
     * Computes invariants I1 and J2 of the strain tensor
     * from the strain components stored in a vector.
     * @param strainVector Input strain components.
     * @param[out] I1e Output value of strain invariant I1.
     * @param[out] J2e Output value of strain invariant J2.
     */
    static void computeStrainInvariants(const FloatArray &strainVector, double &I1e, double &J2e);

    bool isCrackBandApproachUsed() { return ( this->softType == ST_Exponential_Cohesive_Crack || this->softType == ST_Linear_Cohesive_Crack || this->gf != 0. ); }
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual void computeEta(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);
    /**
     * computes the value of damage parameter omega,
     * based on a given value of equivalent strain,
     * using iterations to achieve objectivity,
     * based on the crack band concept (effective element size used)
     * @param[out] omega Contains the resulting damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    void computeDamageParamForCohesiveCrack(double &omega, double kappa, GaussPoint *gp);
    /**
     * Returns the value of damage parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double damageFunction(double kappa, GaussPoint *gp);
    /**
     * Returns the value of compliance parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * The compliance parameter gamma is defined as
     * gamma = omega/(1-omega)
     * where omega is the damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    /**
     * Returns the value of derivative of damage function
     * wrt damage-driving variable kappa corresponding
     * to a given value of the  kappa, depending on
     * the type of selected damage law.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double damageFunctionPrime(double kappa, GaussPoint *gp);
    /**
     * Returns the value of compliance parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * The compliance parameter gamma is defined as
     * gamma = omega/(1-omega)
     * where omega is the damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double complianceFunction(double kappa, GaussPoint *gp);

    virtual Interface *giveInterface(InterfaceType it);

    virtual int MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep);
    virtual int MMI_update(GaussPoint *gp, TimeStep *tStep, FloatArray *estrain = NULL);
    virtual int MMI_finish(TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;

    virtual double give(int aProperty, GaussPoint *gp);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual void giveRealStressVector(FloatArray &answer,  GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    void giveRealStressVectorDmg(FloatArray &answer,  GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    void giveRealStressVectorPlast(FloatArray &answer,  GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
	 double computeStressScalarForYieldFunction(const FloatArray& stress, GaussPoint *gp, TimeStep *tStep, const FloatMatrix& de, const FloatMatrix& ce);
	double computeStressScalarForYieldFunctionFromStrain(double lambda, const FloatArray& elasticStrain, const FloatArray& dStrain, const FloatMatrix& de, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
	void stressReturn(FloatArray& stress, UmicoreMaterialStatus* status, double ft, GaussPoint *gp, TimeStep *tStep);

protected:
    /**
     * Performs initialization, when damage first appear. The characteristic length is
     * computed from the direction of largest positive principal strain and stored
     * in corresponding status.
     * @param kappa Scalar measure of strain level.
     * @param totalStrainVector Current total strain vector.
     * @param gp Integration point.
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp);
};
} // end namespace oofem
#endif // umicore_h
