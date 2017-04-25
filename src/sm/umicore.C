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

#include "umicore.h"
#include "gausspoint.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "classfactory.h"

#include<iostream>

namespace oofem {
REGISTER_Material(UmicoreMaterial);

UmicoreMaterial :: UmicoreMaterial(int n, Domain *d) : IsotropicDamageMaterial1(n, d)
    //
    // constructor
    //
{}


UmicoreMaterial :: ~UmicoreMaterial()
//
// destructor
//
{ }

IRResultType
UmicoreMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";     // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IsotropicDamageMaterial1 :: initializeFrom(ir);

    plastDmgCoeff = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, plastDmgCoeff, _IFT_UmicoreMaterial_plastDmgCoeff);
    if (plastDmgCoeff<0. || plastDmgCoeff>1.) {
      OOFEM_ERROR("Inadmissible value, it should be between 0 and 1");
    }
    
    deactDmgInCompression = false;
    IR_GIVE_OPTIONAL_FIELD(ir, deactDmgInCompression, _IFT_UmicoreMaterial_deactDmgInCompression);
    
    return IRRT_OK;
}


MaterialStatus *
UmicoreMaterial :: CreateStatus(GaussPoint *gp) const
{
    UmicoreMaterialStatus *answer = new UmicoreMaterialStatus(1, UmicoreMaterial :: domain, gp);
    return answer;
}



void
UmicoreMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *tStep)
{
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );

	// compute stresses from damage and plasticity model
	FloatArray stressDmg, stressPlast;
	giveRealStressVectorDmg(stressDmg,gp,totalStrain,tStep);
	giveRealStressVectorPlast(stressPlast,gp,totalStrain,tStep);

	// interpolate using plastDmgCoeff parameter
	stressPlast.times(plastDmgCoeff);
	stressDmg.times(1-plastDmgCoeff);
	answer = stressDmg;
	answer.add(stressPlast);

	// deactivate damage influence in compression
	double sigmaV = answer.at(1)+answer.at(2)+answer.at(3);
	status->tempDamageDeactivatedFlag = false;
	if (deactDmgInCompression && sigmaV < 0.) {
		double omega = status->giveTempDamage();
		answer.times(1/(1-omega*(1-plastDmgCoeff)));
		status->tempDamageDeactivatedFlag = true;
	}

	// update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}


void
UmicoreMaterial :: giveRealStressVectorDmg(FloatArray &answer, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
{

  IsotropicDamageMaterial1 :: giveRealStressVector(answer, gp, totalStrain, tStep);

}


void
UmicoreMaterial :: giveRealStressVectorPlast(FloatArray &answer, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *tStep)
{
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );
    FloatArray strainVector, reducedTotalStrainVector;
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);
	 //
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
	 FloatArray strain, stress;
    StructuralMaterial :: giveFullSymVectorForm( strain, reducedTotalStrainVector, gp->giveMaterialMode() );
    if ( gp->giveMaterialMode() == _PlaneStress ) {
        double nu = lmat->give(NYxz, gp);
        strain.at(3) = -nu * ( strain.at(1) + strain.at(2) ) / ( 1. - nu );
    } else if ( gp->giveMaterialMode() == _1dMat ) {
        double nu = lmat->give(NYxz, gp);
        strain.at(2) = -nu *strain.at(1);
        strain.at(3) = -nu *strain.at(1);
    }

	//
	 FloatArray elasticStrain = strain;
	 elasticStrain.subtract(status->plasticStrain);

    FloatMatrix de;
    lmat->give3dMaterialStiffnessMatrix(de, SecantStiffness, gp, tStep);
	 FloatMatrix ce; ce.beInverseOf(de);

	 stress.beProductOf(de,elasticStrain);
	 double omega = status->giveTempDamage();
	 double ft = lmat->give('E', gp)*std::max(e0,status->giveTempKappa())*(1-omega);
	 status->tempElasticFlag = true;
	 double y = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce) - ft;
	 if (y > 0) {
		 status->tempElasticFlag = false;
			stressReturn(stress,status,ft,gp,tStep);
	 }

    StructuralMaterial :: giveReducedSymVectorForm(answer, stress, gp->giveMaterialMode() );
	 FloatArray dPlasticStrain = status->tempPlasticStrain;
	 dPlasticStrain.subtract(status->plasticStrain);
	 status->tempCumPlasticDeformation = status->cumPlasticDeformation + dPlasticStrain.computeNorm();
}

void
UmicoreMaterial :: stressReturn(FloatArray& stress, UmicoreMaterialStatus* status, double ft, GaussPoint *gp, TimeStep *tStep) {
	 LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
	 FloatMatrix de;
	 lmat->give3dMaterialStiffnessMatrix(de, SecantStiffness, gp, tStep);
	 FloatMatrix ce; ce.beInverseOf(de);
	 StrainVector elasticStrain(_3dMat);
	 elasticStrain.beProductOf(ce,stress);
	 FloatArray princVals;
	 FloatMatrix princDirs;
	 elasticStrain.computePrincipalValDir(princVals,princDirs);
	 //
	 double lambda = .5;
	 double y0,y1,yC;
	 y0 = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce)-ft;
	 FloatArray e = elasticStrain;
	 e.times(-1);
	 e.add(elasticStrain);
	 stress.beProductOf(de,e);
	 y1 = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce)-ft;
	 if (y0<0) { _error("TODO y0"); }
	 if (y1>0) { _error("TODO y1"); }
	 double lambda0=0, lambda1=1;
	//
	 for (int i=0; i<30; i++) {
		 e = elasticStrain;
		 e.times(-lambda);
		 e.add(elasticStrain);
		 stress.beProductOf(de,e);
		 yC = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce)-ft;
		 if (yC > 0.) lambda0 = lambda;
		else lambda1 = lambda;
		lambda = .5*(lambda0+lambda1);
	 }
	 //
		FloatArray tempPlasticStrain = elasticStrain;
		tempPlasticStrain.times(lambda);
		tempPlasticStrain.add(status->plasticStrain);
		status->tempPlasticStrain = tempPlasticStrain;
	 return;
}

double
UmicoreMaterial :: computeStressScalarForYieldFunction(const FloatArray& stress, GaussPoint *gp, TimeStep *tStep, const FloatMatrix& de, const FloatMatrix& ce) {
	double ret = 0.;
    if ( this->equivStrainType == EST_Mazars ) {
		 LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
		 StrainVector strain(_3dMat);
		 strain.beProductOf(ce,stress);
		 FloatArray princVals;
		 strain.computePrincipalValues(princVals);
		 for (int i=1; i<=3; i++) {
			if ( princVals.at(i) > 0.0 ) {
				 ret += std::pow(princVals.at(i),2);
			 }
		 }
        double E = lmat->give('E', gp);
		 return E*std::sqrt(ret);
	 }
    if ( this->equivStrainType == EST_Rankine_Standard || this->equivStrainType == EST_Rankine_Smooth) {
		 StressVector stressV(stress,_3dMat);
		 FloatArray principalStress;
		 stressV.computePrincipalValues(principalStress);
       for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                if ( this->equivStrainType == EST_Rankine_Smooth ) {
                    ret += principalStress.at(i) * principalStress.at(i);
                } else if ( ret < principalStress.at(i) ) {
                    ret = principalStress.at(i);
                }
            } else if ( ret < principalStress.at(i) ) {
                ret = principalStress.at(i);
            }
        }

        if ( this->equivStrainType == EST_Rankine_Smooth ) {
            ret = sqrt(ret);
        }
		  return ret;
	 }
	  _error("computeStressScalarForYieldFunction: unknown EquivStrainType");
	return 0;
}

void
UmicoreMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
      tempDamage = status->giveTempDamage();
      tempDamage = min(tempDamage, maxOmega);
    }

    FloatMatrix dd,dp;
    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(dd, mode, gp, tStep);
    dp = dd;
    //
    if (!status->tempElasticFlag) dp.times(1-tempDamage);
    if (!status->tempDamageDeactivatedFlag) dd.times(1-tempDamage);
    //
    dp.times(plastDmgCoeff);
    dd.times(1-plastDmgCoeff);
    answer = dp;
    answer.add(dd);
}

int
UmicoreMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_CumPlasticStrain ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->cumPlasticDeformation;
        return 1;
    } else {
        return IsotropicDamageMaterial1 :: giveIPValue(answer, gp, type, tStep);
    }
    return 1;
}




UmicoreMaterialStatus :: UmicoreMaterialStatus(int n, Domain *d, GaussPoint *g) :
    IsotropicDamageMaterial1Status(n, d, g)
{

	 plasticStrain.resize(6); plasticStrain.zero();
	 tempPlasticStrain.resize(6); tempPlasticStrain.zero();
	 cumPlasticDeformation = tempCumPlasticDeformation = 0.;
}

void
UmicoreMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterialStatus :: initTempStatus();
    damageDeactivatedFlag = false;
}



void
UmicoreMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    cumPlasticDeformation = tempCumPlasticDeformation;
    damageDeactivatedFlag = tempDamageDeactivatedFlag;
    elasticFlag = tempElasticFlag;
}


}     // end namespace oofem
