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
#include "isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "strainvector.h"
#include "stressvector.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#include<iostream>

namespace oofem {
REGISTER_Material(UmicoreMaterial);

#ifdef IDM_USE_MMAClosestIPTransfer
MMAClosestIPTransfer UmicoreMaterial :: mapper;
#endif

#ifdef IDM_USE_MMAContainingElementProjection
MMAContainingElementProjection UmicoreMaterial :: mapper;
#endif

#ifdef IDM_USE_MMAShapeFunctProjection
MMAShapeFunctProjection UmicoreMaterial :: mapper;
#endif

#ifdef IDM_USE_MMALeastSquareProjection
MMALeastSquareProjection UmicoreMaterial :: mapper;
#endif

UmicoreMaterial :: UmicoreMaterial(int n, Domain *d) : IsotropicDamageMaterial(n, d),
    RandomMaterialExtensionInterface()
    //
    // constructor
    //
{
    // deleted by parent, where linearElasticMaterial instance declared
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    equivStrainType = EST_Unknown;
    softType = ST_Unknown;
    k = 0.;
    md = 1.;
    damageLaw = 0;
    e0 = 0.;
    ef = 0.;
    wf = 0.;
    gf = 0.;
    griff_n = 8.;
    ecsMethod = ECSM_Unknown;
}


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

    int equivStrainTypeRecord;
    IsotropicDamageMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    checkSnapBack = 1; //check by default
    IR_GIVE_OPTIONAL_FIELD(ir, checkSnapBack, _IFT_UmicoreMaterial_checkSnapBack);

    // specify the type of formula for equivalent strain
    equivStrainTypeRecord = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, equivStrainTypeRecord, _IFT_UmicoreMaterial_equivstraintype);
    if ( equivStrainTypeRecord == 0 ) {
        this->equivStrainType = EST_Mazars;
    } else if ( equivStrainTypeRecord == 1 ) {
        this->equivStrainType = EST_Rankine_Smooth;
    } else if ( equivStrainTypeRecord == 2 ) {
        this->equivStrainType = EST_ElasticEnergy;
    } else if ( equivStrainTypeRecord == 3 ) {
        this->equivStrainType = EST_Mises;
        IR_GIVE_FIELD(ir, k, _IFT_UmicoreMaterial_k);
    } else if ( equivStrainTypeRecord == 4 ) {
        this->equivStrainType = EST_Rankine_Standard;
    } else if ( equivStrainTypeRecord == 5 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStress;
    } else if ( equivStrainTypeRecord == 6 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStrain;
    } else if ( equivStrainTypeRecord == 7 ) {
        this->equivStrainType = EST_Griffith;
        IR_GIVE_OPTIONAL_FIELD(ir, griff_n, _IFT_UmicoreMaterial_n);
    } else {
        OOFEM_ERROR2("Unknown equivStrainType %d", equivStrainType);
    }

    // specify the type of formula for damage evolution law
    IR_GIVE_OPTIONAL_FIELD(ir, damageLaw, _IFT_UmicoreMaterial_damageLaw);
    if ( ( damageLaw != 6 ) && ( damageLaw != 7 ) ) {
        IR_GIVE_FIELD(ir, e0, _IFT_UmicoreMaterial_e0);
    }

    //applies only in this class
    switch ( damageLaw ) {
    case 0:     // exponential softening - default
        if ( ir->hasField(_IFT_UmicoreMaterial_wf) ) {
            this->softType = ST_Exponential_Cohesive_Crack;
            IR_GIVE_FIELD(ir, wf, _IFT_UmicoreMaterial_wf);
        } else if ( ir->hasField(_IFT_UmicoreMaterial_gf) ) {
            this->softType = ST_Exponential_Cohesive_Crack;
            IR_GIVE_FIELD(ir, gf, _IFT_UmicoreMaterial_gf);
        } else {
            this->softType = ST_Exponential;
            IR_GIVE_FIELD(ir, ef, _IFT_UmicoreMaterial_ef);
        }

        break;
    case 1:     // linear softening law
        if ( ir->hasField(_IFT_UmicoreMaterial_wf) ) {
            this->softType = ST_Linear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, wf, _IFT_UmicoreMaterial_wf);
        } else if ( ir->hasField(_IFT_UmicoreMaterial_gf) ) {
            this->softType = ST_Linear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, gf, _IFT_UmicoreMaterial_gf);
        } else {
            this->softType = ST_Linear;
            IR_GIVE_FIELD(ir, ef, _IFT_UmicoreMaterial_ef);
        }

        break;
    case 2:     // bilinear softening
        if ( ir->hasField(_IFT_UmicoreMaterial_gf) ) {
            this->softType = ST_BiLinear_Cohesive_Crack;
            IR_GIVE_FIELD(ir, gf, _IFT_UmicoreMaterial_gf);
            // ek is for the bilinear law, and corresponds to the strain at the knee point
            IR_GIVE_FIELD(ir, ek, _IFT_UmicoreMaterial_ek);
            // Gft is for the bilinear law, and corresponds to the total energy required to fail the specimen
            IR_GIVE_FIELD(ir, gft, _IFT_UmicoreMaterial_gft);
        } else {
            OOFEM_ERROR("Bilinear softening for wf and ef not implemented");
        }

        break;
    case 3: // reserved for Hordijk's law
        break;
    case 4:
        this->softType = ST_Mazars;
        IR_GIVE_FIELD(ir, At, _IFT_UmicoreMaterial_At);
        IR_GIVE_FIELD(ir, Bt, _IFT_UmicoreMaterial_Bt);
        break;
    case 5:
        this->softType = ST_Smooth;
        md = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, md, _IFT_UmicoreMaterial_md);
        break;
    case 6:
        this->softType = ST_Disable_Damage;
        break;
    case 7:
        this->softType = ST_SmoothExtended;
        double E, ep, ft; // auxiliary input variables
        IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);
        IR_GIVE_FIELD(ir, ep, _IFT_UmicoreMaterial_ep);
        IR_GIVE_FIELD(ir, ft, _IFT_UmicoreMaterial_ft);
        this->md = 1. / log(E * ep / ft);
        this->e0 = ep * pow(md, 1. / md);
        IR_GIVE_FIELD(ir, e1, _IFT_UmicoreMaterial_e1);
        this->s1 = e1 * exp( -pow(e1 / e0, md) );
        this->ef = -e1 / ( 1. - md * pow(e1 / e0, md) );
        IR_GIVE_FIELD(ir, e2, _IFT_UmicoreMaterial_e2);
        IR_GIVE_FIELD(ir, nd, _IFT_UmicoreMaterial_nd);
        break;
    default:
        OOFEM_ERROR2("Softening type number %d is unknown", damageLaw);
    }

    if ( ( softType == ST_Exponential_Cohesive_Crack ) || ( softType == ST_Linear_Cohesive_Crack ) || ( softType == ST_BiLinear_Cohesive_Crack ) ) {
        int ecsm = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, ecsm, _IFT_UmicoreMaterial_ecsm);
        switch ( ecsm ) {
        case 1: ecsMethod = ECSM_SquareRootOfArea;
            break;
        case 2: ecsMethod = ECSM_ProjectionCentered;
            break;
        case 3: ecsMethod = ECSM_Oliver1;
            break;
        case 4: ecsMethod = ECSM_Oliver1modified;
            break;
        default: ecsMethod = ECSM_Projection;
        }
    }

    this->mapper.initializeFrom(ir);

	plastDmgCoeff = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, plastDmgCoeff, _IFT_UmicoreMaterial_plastDmgCoeff);
	 if (plastDmgCoeff<0. || plastDmgCoeff>1.) {
		 OOFEM_ERROR("TODO");
	 }

	deactDmgInCompression = false;
    IR_GIVE_OPTIONAL_FIELD(ir, deactDmgInCompression, _IFT_UmicoreMaterial_deactDmgInCompression);

    return IRRT_OK;
}


void
UmicoreMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    // Done first so that record name is overwritten afterwards by femcomponent:
    this->linearElasticMaterial->giveInputRecord(input);
    this->mapper.giveInputRecord(input);
    IsotropicDamageMaterial :: giveInputRecord(input);
    RandomMaterialExtensionInterface :: giveInputRecord(input);

    input.setField(this->e0, _IFT_UmicoreMaterial_e0);
    input.setField(this->ef, _IFT_UmicoreMaterial_ef);
    input.setField(this->equivStrainType, _IFT_UmicoreMaterial_equivstraintype);

    input.setField(this->checkSnapBack, _IFT_UmicoreMaterial_checkSnapBack);

    input.setField(this->equivStrainType, _IFT_UmicoreMaterial_equivstraintype);
    if ( this->equivStrainType == EST_Mises ) {
        input.setField(this->k, _IFT_UmicoreMaterial_k);
    }

    input.setField(this->damageLaw, _IFT_UmicoreMaterial_damageLaw);
    if ( this->damageLaw != 7 ) {
        input.setField(this->e0, _IFT_UmicoreMaterial_e0);
    }

    switch ( damageLaw ) {
    case 0:
        if ( this->softType == ST_Exponential_Cohesive_Crack ) {
            input.setField(this->wf, _IFT_UmicoreMaterial_wf);
        } else if ( this->softType == ST_Exponential_Cohesive_Crack ) {
            input.setField(this->gf, _IFT_UmicoreMaterial_gf);
        } else {
            input.setField(this->ef, _IFT_UmicoreMaterial_ef);
        }

        break;
    case 1:     // linear softening law
        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            input.setField(this->wf, _IFT_UmicoreMaterial_wf);
        } else if ( this->softType == ST_Linear_Cohesive_Crack ) {
            input.setField(this->gf, _IFT_UmicoreMaterial_gf);
        } else {
            input.setField(this->ef, _IFT_UmicoreMaterial_ef);
        }

        break;
    case 2:     // bilinear softening
        input.setField(this->gf, _IFT_UmicoreMaterial_gf);
        input.setField(this->ek, _IFT_UmicoreMaterial_ek);
        input.setField(this->gft, _IFT_UmicoreMaterial_gft);

        break;
    case 4:
        input.setField(this->At, _IFT_UmicoreMaterial_At);
        input.setField(this->Bt, _IFT_UmicoreMaterial_Bt);
        break;
    case 5:
        input.setField(this->md, _IFT_UmicoreMaterial_md);
        break;
    case 7:
        input.setField(this->e1, _IFT_UmicoreMaterial_e1);
        input.setField(this->e2, _IFT_UmicoreMaterial_e2);
        input.setField(this->nd, _IFT_UmicoreMaterial_nd);
        break;
    }

    if ( softType == ST_Exponential_Cohesive_Crack || softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) {
        input.setField(this->ecsMethod, _IFT_UmicoreMaterial_ecsm);
    }
}


void
UmicoreMaterial :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray fullStrain;

    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }

    StructuralMaterial :: giveFullSymVectorForm( fullStrain, strain, gp->giveMaterialMode() );
    // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
    if ( gp->giveMaterialMode() == _PlaneStress ) {
        double nu = lmat->give(NYxz, gp);
        fullStrain.at(3) = -nu * ( fullStrain.at(1) + fullStrain.at(2) ) / ( 1. - nu );
    } else if ( gp->giveMaterialMode() == _1dMat ) {
        double nu = lmat->give(NYxz, gp);
        fullStrain.at(2) = -nu *fullStrain.at(1);
        fullStrain.at(3) = -nu *fullStrain.at(1);
    }

    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains;

        this->computePrincipalValues(principalStrains, fullStrain, principal_strain);

        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        kappa = sqrt(posNorm);
    } else if ( ( this->equivStrainType == EST_Rankine_Smooth ) || ( this->equivStrainType == EST_Rankine_Standard ) ) {
        // EST_Rankine equiv strain measure
        double sum = 0.;
        FloatArray stress, fullStress, principalStress;
        FloatMatrix de;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        stress.beProductOf(de, strain);
        StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                if ( this->equivStrainType == EST_Rankine_Smooth ) {
                    sum += principalStress.at(i) * principalStress.at(i);
                } else if ( sum < principalStress.at(i) ) {
                    sum = principalStress.at(i);
                }
            } else if ( sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
            }
        }

        if ( this->equivStrainType == EST_Rankine_Smooth ) {
            sum = sqrt(sum);
        }

        kappa = sum / lmat->give('E', gp);
    } else if ( ( this->equivStrainType == EST_ElasticEnergy ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStrain ) ) {
        // equivalent strain expressions based on elastic energy
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        if ( this->equivStrainType == EST_ElasticEnergy ) {
            // standard elastic energy
            stress.beProductOf(de, strain);
            sum = strain.dotProduct(stress);
        } else if ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) {
            // elastic energy corresponding to positive part of stress
            FloatArray fullStress, principalStress;
            StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
            this->computePrincipalValues(principalStress, fullStress, principal_stress);
            // TO BE FINISHED
            sum = 0.;
            OOFEM_ERROR("Elastic energy corresponding to positive part of stress not finished\n");
        } else {
            // elastic energy corresponding to positive part of strain
            // TO BE DONE
            sum = 0.;
            OOFEM_ERROR("Elastic energy corresponding to positive part of strain not finished\n");
        }

        kappa = sqrt( sum / lmat->give('E', gp) );
    } else if ( this->equivStrainType == EST_Mises ) {
        double nu = lmat->give(NYxz, NULL);
        FloatArray principalStrains;

        this->computePrincipalValues(principalStrains, fullStrain, principal_strain);
        double I1e, J2e;
        this->computeStrainInvariants(principalStrains, I1e, J2e);
        double a, b, c;
        a = ( k - 1 ) * I1e / ( 2 * k * ( 1 - 2 * nu ) );
        b = ( k - 1 ) * ( k - 1 ) * I1e * I1e / ( ( 1 - 2 * nu ) * ( 1 - 2 * nu ) );
        c = 12 * k * J2e / ( ( 1 + nu ) * ( 1 + nu ) );
        kappa = a + 1 / ( 2 * k ) * sqrt(b + c);
    } else if ( this->equivStrainType == EST_Griffith ) {
        kappa = 0.0;
        double kappa1 = 0.0, kappa2 = 0.0;
        FloatArray stress, fullStress, principalStress;
        FloatMatrix de;
        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        stress.beProductOf(de, strain);
        StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        double sum = 0., maxStress = 0.;
        //Compute equivalent strain for Rankine's criterion
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 && sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
            }
        }
        kappa1 = sum / lmat->give('E', gp);
        //Compute equivalent strain for Griffith's criterion
        //Check zero division first
        maxStress = max( fabs( principalStress.at(1) ), fabs( principalStress.at(3) ) );
        if ( maxStress == 0. || fabs( principalStress.at(3) ) < 1.e-6 * maxStress || fabs( principalStress.at(1) + principalStress.at(3) ) < 1.e-6 * maxStress ) {
            //Skip evaluation
        } else if ( principalStress.at(1) / principalStress.at(3) >= -0.33333 ) {
            kappa2 = -( principalStress.at(1) - principalStress.at(3) ) * ( principalStress.at(1) - principalStress.at(3) ) / this->griff_n / ( principalStress.at(1) + principalStress.at(3) ) / lmat->give('E', gp);
        }
        kappa = max(kappa1, 0.0);
        kappa = max(kappa, kappa2);
    } else {
        _error("computeEquivalentStrain: unknown EquivStrainType");
    }
}

//Computes derivative of the equivalent strain with regards to strain, used in tangent formulation
void
UmicoreMaterial :: computeEta(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();

    if ( strain.isEmpty() ) {
        answer.zero();
        return;
    }

    if ( this->equivStrainType == EST_Mazars ) {
        int dim = 0;
        double posNorm = 0.0;
        double nu = lmat->give(NYxz, gp);
        FloatArray principalStrains, fullStrain;
        FloatMatrix N, m;

        if ( gp->giveMaterialMode() == _1dMat ) {
            dim = 1;
            StrainVector fullStrain(strain, _1dMat);
            fullStrain.computePrincipalValDir(principalStrains, N);
            principalStrains.resizeWithValues(3);
            principalStrains.at(2) = -nu *principalStrains.at(1);
            principalStrains.at(3) = -nu *principalStrains.at(1);
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            dim = 2;
            StrainVector fullStrain(strain, _PlaneStress);
            fullStrain.computePrincipalValDir(principalStrains, N);
            principalStrains.resizeWithValues(3);
            principalStrains.at(3) = -nu * ( principalStrains.at(1) + principalStrains.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            dim = 3;
            StrainVector fullStrain(strain, _PlaneStrain);
            fullStrain.computePrincipalValDir(principalStrains, N);
        } else if ( gp->giveMaterialMode() == _3dMat ) {
            dim = 3;
            StrainVector fullStrain(strain, _3dMat);
            fullStrain.computePrincipalValDir(principalStrains, N);
        } else {
            dim = 0;
            OOFEM_ERROR("UmicoreMaterial :: computeEta - Unknown material mode.");
        }

        FloatArray n(dim);
        FloatMatrix Eta(dim, dim);
        Eta.zero();

        for ( int i = 1; i <= 3; i++ ) {
            if ( i <= dim ) {
                if ( principalStrains.at(i) > 0.0 ) {
                    for ( int j = 1; j < 3; j++ ) {
                        n.at(j) = N.at(j, i);
                    }

                    Eta.plusDyadSymmUpper( n, principalStrains.at(i) );
                }
            }

            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        Eta.symmetrized();

        double kappa = sqrt(posNorm);

        int numberOfEl = ( dim * ( dim - 1 ) / 2 + dim );
        answer.resize(numberOfEl);

        if ( kappa != 0 ) {
            Eta.times(1. / kappa);
        } else {
            answer.zero();
            return;
        }

        if ( gp->giveMaterialMode() == _1dMat ) {
            answer.at(1) = Eta.at(1, 1);
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            answer.at(1) = Eta.at(1, 1);
            answer.at(2) = Eta.at(2, 2);
            answer.at(3) = Eta.at(1, 2);
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            answer.resize(4);
            answer.at(1) = Eta.at(1, 1);
            answer.at(2) = Eta.at(2, 2);
            answer.at(3) = Eta.at(3, 3);
            answer.at(4) = Eta.at(1, 2);
        } else if ( gp->giveMaterialMode() == _3dMat ) {
            answer.at(1) = Eta.at(1, 1);
            answer.at(2) = Eta.at(2, 2);
            answer.at(3) = Eta.at(3, 3);
            answer.at(4) = Eta.at(2, 3);
            answer.at(5) = Eta.at(1, 3);
            answer.at(6) = Eta.at(1, 2);
        }
    } else if ( ( this->equivStrainType == EST_Rankine_Smooth ) || ( this->equivStrainType == EST_Rankine_Standard ) ) {
        int index = 0, dim = 0;
        double sum = 0.;
        FloatArray stress, principalStress, eta;
        FloatMatrix de, N, m, Eta;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        stress.beProductOf(de, strain);

        if ( gp->giveMaterialMode() == _1dMat ) {
            StressVector fullStress(stress, _1dMat);
            fullStress.computePrincipalValDir(principalStress, N);
            principalStress.resizeWithValues(3);
            dim = 1;
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            StressVector fullStress(stress, _PlaneStress);
            fullStress.computePrincipalValDir(principalStress, N);
            principalStress.resizeWithValues(3);
            dim = 2;
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            StressVector fullStress(stress, _PlaneStrain);
            fullStress.computePrincipalValDir(principalStress, N);
            dim = 3;
        } else if ( gp->giveMaterialMode() == _3dMat ) {
            StressVector fullStress(stress, _3dMat);
            fullStress.computePrincipalValDir(principalStress, N);
            dim = 3;
        } else {
            dim = 0;
            OOFEM_ERROR("UmicoreMaterial :: computeEta - Unknown material mode.");
        }

        FloatArray n(dim);
        n.zero();
        Eta.resize(dim, dim);
        Eta.zero();
        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                if ( this->equivStrainType == EST_Rankine_Smooth ) {
                    sum += principalStress.at(i) * principalStress.at(i);
                    if ( i <= dim ) {
                        for ( int j = 1; j <= dim; j++ ) {
                            n.at(j) = N.at(j, i);
                        }
                    }

                    Eta.plusDyadSymmUpper( n, principalStress.at(i) );
                } else if ( sum < principalStress.at(i) ) {
                    sum = principalStress.at(i);
                    index = i;
                }
            } else if ( sum < principalStress.at(i) ) {
                sum = principalStress.at(i);
                index = i;
            }
        }

        Eta.symmetrized();

        int numberOfEl = ( dim * ( dim - 1 ) / 2 + dim );
        eta.resize(numberOfEl);

        if ( this->equivStrainType == EST_Rankine_Smooth ) {
            sum = sqrt(sum);
            double kappa =  sum / lmat->give('E', gp);

            if ( kappa != 0 ) {
                Eta.times(1. / kappa);
            } else {
                answer.zero();
                return;
            }

            Eta.times( 1. / kappa / lmat->give('E', gp) / lmat->give('E', gp) );
        } else if ( this->equivStrainType == EST_Rankine_Standard ) {
            for ( int i = 1; i <= dim; i++ ) {
                n.at(i) = N.at(i, index);
            }

            Eta.beDyadicProductOf(n, n);
            Eta.times( 1. / lmat->give('E', gp) );
        }

        if ( gp->giveMaterialMode() == _1dMat ) {
            eta.at(1) = Eta.at(1, 1);
        } else if ( gp->giveMaterialMode() == _PlaneStress ) {
            eta.at(1) = Eta.at(1, 1);
            eta.at(2) = Eta.at(2, 2);
            eta.at(3) = Eta.at(1, 2);
        }  else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            eta.resize(4);
            eta.at(1) = Eta.at(1, 1);
            eta.at(2) = Eta.at(2, 2);
            eta.at(3) = Eta.at(3, 3);
            eta.at(4) = 2. * Eta.at(1, 2);
        } else if ( gp->giveMaterialMode() == _3dMat ) {
            eta.at(1) = Eta.at(1, 1);
            eta.at(2) = Eta.at(2, 2);
            eta.at(3) = Eta.at(3, 3);
            eta.at(4) = Eta.at(2, 3);
            eta.at(5) = Eta.at(1, 3);
            eta.at(6) = Eta.at(1, 2);
        }

        answer.beProductOf(de, eta);
    } else if ( this->equivStrainType == EST_ElasticEnergy ) {
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
        // standard elastic energy
        stress.beProductOf(de, strain);
        sum = strain.dotProduct(stress);

        double kappa = sqrt( sum / lmat->give('E', gp) );
        answer = stress;
        answer.times(1. / lmat->give('E', gp) / kappa);
    } else {
        _error("computeEta: unknown EquivStrainType");
    }
}

void
UmicoreMaterial :: computeStrainInvariants(const FloatArray &strainVector, double &I1e, double &J2e)
{
    I1e = strainVector.at(1) + strainVector.at(2) + strainVector.at(3);
    double s1 = strainVector.at(1) * strainVector.at(1);
    double s2 = strainVector.at(2) * strainVector.at(2);
    double s3 = strainVector.at(3) * strainVector.at(3);
    J2e = 1. / 2. * ( s1 + s2 + s3 ) - 1. / 6. * ( I1e * I1e );
}

void
UmicoreMaterial :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    if ( this->softType == ST_Disable_Damage ) { //dummy material with no damage
        omega = 0.;
    } else if ( isCrackBandApproachUsed() ) { // adjustment of softening law according to the element size, given crack opening or fracture energy
        computeDamageParamForCohesiveCrack(omega, kappa, gp);
    } else { // no adjustment according to element size, given fracturing strain
        omega = damageFunction(kappa, gp);
    }
}

void
UmicoreMaterial :: computeDamageParamForCohesiveCrack(double &omega, double kappa, GaussPoint *gp)
{
    const double e0 = this->give(e0_ID, gp);  // e0 is the strain at the peak stress
    const double E = this->giveLinearElasticMaterial()->give('E', gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);     // wf is the crack opening
    double Le;
    omega = 0.0;

    if ( kappa > e0 ) {
        if ( this->gf != 0. ) { //cohesive crack model
            if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
                wf = this->gf / E / e0; // wf is the crack opening
            } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi)linear softening law
                wf = 2. * gf / E / e0; // wf is the crack opening
            } else {
                OOFEM_ERROR2("Gf unsupported for softening type softType = %d", softType);
            }
        }

        UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );
        Le = status->giveLe();
        ef = wf / Le;    //ef is the fracturing strain
        if ( ef < e0 ) { //check that no snapback occurs
            double minGf = 0.;
            OOFEM_WARNING5("ef %e < e0 %e, this leads to material snapback in element %d, characteristic length %f", ef, e0, gp->giveElement()->giveNumber(), Le);
            if ( gf != 0. ) { //cohesive crack
                if ( softType == ST_Exponential_Cohesive_Crack ) { //exponential softening
                    minGf = E * e0 * e0 * Le;
                } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { //(bi)linear softening law
                    minGf = E * e0 * e0 * Le / 2.;
                } else {
                    OOFEM_WARNING2("Gf unsupported for softening type softType = %d", softType);
                }

                if ( checkSnapBack ) {
                    OOFEM_ERROR4("Material number %d, decrease e0, or increase Gf from %f to Gf=%f", this->giveNumber(), gf, minGf);
                }
            }

            if ( checkSnapBack ) { //given fracturing strain
                OOFEM_ERROR4("Material number %d, increase ef %f to minimum e0 %f", this->giveNumber(), ef, e0);
            }
        }

        if ( this->softType == ST_Linear_Cohesive_Crack ) {
            if ( kappa < ef ) {
                omega = ( ef / kappa ) * ( kappa - e0 ) / ( ef - e0 );
            } else {
                omega = 1.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
            }
        } else if (  this->softType == ST_BiLinear_Cohesive_Crack ) {
            double ek = this->give(ek_ID, gp);
            double gft = this->give(gft_ID, gp);
            double ef = 2 * gf / E / e0 / Le; //the first part corresponds to linear softening
            double sigmak = E * e0 * ( ef - ek ) / ( ef - e0 );
            double epsf = 2 * ( gft - gf ) / sigmak / Le + ef;
            if ( ( ek > ef ) || ( ek < e0 ) ) {
                OOFEM_WARNING4("ek %f is not between e0 %f and ef %f", ek, e0, ef);
            }

            if ( gft < gf ) {
                OOFEM_ERROR3("The total fracture energy gft %f must be greater than the initial fracture energy gf %f", gft, gf);
            }

            if ( kappa <= ek ) {
                omega = 1.0 - ( ( e0 / kappa ) * ( ek - kappa ) / ( ek - e0 ) + ( ( sigmak / ( E * kappa ) ) * ( kappa - e0 ) / ( ek - e0 ) ) );
            } else if ( kappa > ek && kappa <= epsf ) {
                omega = 1.0 - ( ( sigmak / ( E * kappa ) ) * ( epsf - kappa ) / ( epsf - ek ) );
            } else if ( kappa <= e0 ) {
                omega = 0.0;
            } else {
                omega = maxOmega;
            }
        } else if (  this->softType == ST_Exponential_Cohesive_Crack ) {
            // exponential cohesive crack - iteration needed
            double R, Lhs, help;
            int nite = 0;
            // iteration to achieve objectivity
            // we are looking for a state in which the elastic stress is equal to
            // the stress from crack-opening relation
            // ef has now the meaning of strain
            do {
                nite++;
                help = omega * kappa / ef;
                R = ( 1. - omega ) * kappa - e0 *exp(-help); //residuum
                Lhs = kappa - e0 *exp(-help) * kappa / ef; //- dR / (d omega)
                omega += R / Lhs;
                if ( nite > 40 ) {
                    _error("computeDamageParamForCohesiveCrack: algorithm not converging");
                }
            } while ( fabs(R) >= e0 * IDM1_ITERATION_LIMIT );
        } else {
            OOFEM_ERROR1("Unknown softening type for cohesive crack model.");
        }

        if ( omega > 1.0 ) {
            OOFEM_WARNING2("computeDamageParam: damage parameter is %f, which is greater than 1, snap-back problems", omega);
            omega = maxOmega;
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        }

        if ( omega < 0.0 ) {
            OOFEM_WARNING2("computeDamageParam: damage parameter is %f, which is smaller than 0, snap-back problems", omega);
            omega = 0.0;
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        }
    }
}
double
UmicoreMaterial :: damageFunction(double kappa, GaussPoint *gp)
{
    const double e0 = this->give(e0_ID, gp);
    double ef = 0.;
    if ( softType == ST_Linear || softType == ST_Exponential || softType == ST_SmoothExtended ) {
        ef = this->give(ef_ID, gp);         // ef is the fracturing strain
    }

    switch ( softType ) {
    case ST_Linear:
        if ( kappa <= e0 ) {
            return 0.0;
        } else if ( kappa < ef ) {
            return ( ef / kappa ) * ( kappa - e0 ) / ( ef - e0 );
        } else {
            return 1.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
        }

    case ST_Exponential:
        if ( kappa > e0 ) {
            return 1.0 - ( e0 / kappa ) * exp( -( kappa - e0 ) / ( ef - e0 ) );
        } else {
            return 0.0;
        }

    case ST_Mazars:
        return 1.0 - ( 1.0 - At ) * e0 / kappa - At *exp( -Bt * ( kappa - e0 ) );

    case ST_Smooth:
        return 1.0 - exp( -pow(kappa / e0, md) );

    case ST_SmoothExtended:
        // a special damage law used in Grassl and Jirasek (2010)
        if ( kappa <= e1 ) {
            return 1.0 - exp( -pow(kappa / e0, md) );
        } else {
            return 1.0 - s1 *exp( -( kappa - e1 ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) ) / kappa;
        }

    default:
        OOFEM_WARNING2("UmicoreMaterial::damageFunction ... undefined softening type %d\n", softType);
    }

    return 0.;         // to make the compiler happy
}

double
UmicoreMaterial :: damageFunctionPrime(double kappa, GaussPoint *gp)
{
    const double e0 = this->give(e0_ID, gp);
    double ef = 0.;
    const double E = this->giveLinearElasticMaterial()->give('E', gp);
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );
    const double Le = status->giveLe();

    if ( softType == ST_Linear || softType == ST_Exponential || softType == ST_SmoothExtended ) {
        ef = this->give(ef_ID, gp);         // ef is the fracturing strain
    }

    switch ( softType ) {
    case ST_Linear:
    {
        if ( kappa <= e0 ) {
            return 0.0;
        } else if ( kappa < ef ) {
            return ( ef * e0 ) / ( ef - e0 ) / ( kappa * kappa );
        } else {
            return 0.0; //maximum omega (maxOmega) is adjusted just for stiffness matrix in isodamagemodel.C
        }
    } break;
    case ST_Exponential:
    {
        if ( kappa > e0 ) {
            return ( e0 / ( kappa * kappa ) ) * exp( -( kappa - e0 ) / ( ef - e0 )  + e0 / ( kappa * ( ef - e0 ) ) ) * exp( -( kappa - e0 ) / ( ef - e0 ) );
        } else {
            return 0.0;
        }
    } break;
    case ST_Mazars:
    {
        return ( 1.0 - At ) * e0 / kappa / kappa - At *Bt *exp( -Bt * ( kappa - e0 ) );
    } break;
    case ST_Smooth:
    {
        return md / e0 *pow(kappa, md - 1.) * exp( -pow(kappa / e0, md) );
    } break;
    case ST_SmoothExtended:
    {
        // a special damage law used in Grassl and Jirasek (2010)
        if ( kappa <= 0 ) {
            return 0;
        } else if ( kappa <= e1 ) {
            return exp( -pow(kappa / e0, md) ) * md / pow(e0, md) * pow(kappa, md - 1.);
        } else {
            double a = ( ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) -  ef * nd * pow( ( kappa - e1 ) / e2, nd ) ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) );
            double answer =   s1 * exp( -( kappa - e1 ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) ) / kappa / kappa +  s1 *exp( -( kappa - e1 ) / ( ef * ( 1. + pow( ( kappa - e1 ) / e2, nd ) ) ) ) / kappa * a;
            return answer;
        }
    } break;
    case ST_Linear_Cohesive_Crack:
    {
        double wf = 2. * gf / E / e0; // wf is the crack opening
        if ( kappa <= e0 ) {
            return 0.0;
        } else if ( kappa < wf / Le ) {
            return ( e0 / ( kappa * kappa ) / ( 1. - Le * e0 / wf ) );
        } else {
            return 0.0;
        }
    } break;
    case ST_Exponential_Cohesive_Crack:
    {
        if ( kappa > e0 ) {
            double ef = gf / E / e0 / Le;
            double omega = status->giveTempDamage();
            double help = exp(omega * kappa / ef);
            double ret = -( ( omega * ef - ef ) * help - omega * e0 ) / ( ef * kappa * help - e0 * kappa );
            if ( isnan(ret) ) {
                return 0.;
            }
            return ret;
        } else {
            return 0.0;
        }
    } break;
    default:
        OOFEM_ERROR2("UmicoreMaterial::damageFunctionPrime ... undefined softening type %d\n", softType);
    }

    return 0.;         // to make the compiler happy
}

double
UmicoreMaterial :: complianceFunction(double kappa, GaussPoint *gp)
{
    double om = damageFunction(kappa, gp);
    return om / ( 1. - om );
}

void
UmicoreMaterial :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp)
{
    int i, indx = 1;
    double le = 0.;
    double E = this->giveLinearElasticMaterial()->give('E', gp);
    FloatArray principalStrains, crackPlaneNormal(3), fullStrain, crackVect(3);
    FloatMatrix principalDir(3, 3);
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );

    const double e0 = this->give(e0_ID, gp);
    const double ef = this->give(ef_ID, gp);
    const double gf = this->give(gf_ID, gp);
    double wf = this->give(wf_ID, gp);


    if ( softType == ST_Disable_Damage ) {
        return;
    }

    if ( gf != 0. ) { //cohesive crack model
        if ( softType == ST_Exponential_Cohesive_Crack ) { // exponential softening
            wf = gf / E / e0; // wf is the crack opening
        } else if ( softType == ST_Linear_Cohesive_Crack || softType == ST_BiLinear_Cohesive_Crack ) { // (bi) linear softening law
            wf = 2. * gf / E / e0; // wf is the crack opening
        } else {
            OOFEM_ERROR2("Gf unsupported for softening type softType = %d", softType);
        }
    }

    StructuralMaterial :: giveFullSymVectorForm( fullStrain, strainVector, gp->giveMaterialMode() );

    if ( ( kappa > e0 ) && ( status->giveDamage() == 0. ) ) {
        this->computePrincipalValDir(principalStrains, principalDir, fullStrain, principal_strain);
        // find index of max positive principal strain
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        // find index with minimal value but non-zero for plane-stress condition - this is the crack direction
        indx = 1;
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) < principalStrains.at(indx) && fabs( principalStrains.at(i) ) > 1.e-10 ) {
                indx = i;
            }
        }

        // Use orientation of the worst inclusion for Griffith criterion in compression.
        if ( this->equivStrainType == EST_Griffith ) {
            FloatArray stress, fullStress, principalStress, crackV(3), crackPlaneN(3);
            FloatMatrix de;
            LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
            lmat->giveStiffnessMatrix( de, SecantStiffness, gp, domain->giveEngngModel()->giveCurrentStep() );
            stress.beProductOf(de, strainVector);
            StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
            this->computePrincipalValDir(principalStress, principalDir, fullStress, principal_stress);
            if (  principalStress.at(1) <= 1.e-10 && principalStress.at(2) <= 1.e-10 && principalStress.at(3) <= 1.e-10 ) {
                int indexMax = principalStress.giveIndexMaxElem();
                int indexMin = principalStress.giveIndexMinElem();
                int indexMid = 0;
                if ( indexMin + indexMax == 3 ) {
                    indexMid = 3;
                } else if ( indexMin + indexMax == 4 ) {
                    indexMid = 2;
                } else if ( indexMin + indexMax == 5 ) {
                    indexMid = 1;
                }

                //inclination from maximum compressive stress (plane sig_1 and sig_3)
                double twoPsi = ( principalStress.at(indexMin) - principalStress.at(indexMax) ) / 2. / ( principalStress.at(indexMax) + principalStress.at(indexMin) );
                double psi = acos(twoPsi) / 2.;
                for ( int i = 1; i <= 3; i++ ) {
                    crackV.at(i) = principalDir.at(i, indexMin);
                    crackPlaneN = principalDir.at(i, indexMax);
                }

                //rotate around indexMid axis
                //see http://en.wikipedia.org/wiki/Rotation_matrix and Rodrigues' rotation formula
                FloatMatrix ux(3, 3), dyadU(3, 3), unitMtrx(3, 3), rotMtrx(3, 3);
                FloatArray u(3);
                ux.zero();
                ux.at(1, 2) = -principalDir.at(3, indexMid);
                ux.at(1, 3) = principalDir.at(2, indexMid);
                ux.at(2, 1) = principalDir.at(3, indexMid);
                ux.at(2, 3) = -principalDir.at(1, indexMid);
                ux.at(3, 1) = -principalDir.at(2, indexMid);
                ux.at(3, 2) = principalDir.at(1, indexMid);
                u.at(1) = principalDir.at(1, indexMid);
                u.at(2) = principalDir.at(2, indexMid);
                u.at(3) = principalDir.at(3, indexMid);
                dyadU.beDyadicProductOf(u, u);
                unitMtrx.beUnitMatrix();
                unitMtrx.times( cos(psi) );
                ux.times( sin(psi) );
                dyadU.times( 1. - cos(psi) );
                rotMtrx.zero();
                rotMtrx.add(unitMtrx);
                rotMtrx.add(ux);
                rotMtrx.add(dyadU);
                crackV.rotatedWith(rotMtrx, 'n');
                for ( int i = 1; i <= 3; i++ ) {
                    principalDir.at(i, indx) = crackV.at(i);
                }

                crackPlaneNormal.rotatedWith(rotMtrx, 'n');
            }
        }


        for ( i = 1; i <= 3; i++ ) {
            crackVect.at(i) = principalDir.at(i, indx);
        }

        status->setCrackVector(crackVect);

        if ( isCrackBandApproachUsed() ) { // le needed only if the crack band approach is used
            // old approach (default projection method)
            // le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
            // new approach, with choice of method
            le = gp->giveElement()->giveCharacteristicSize(gp, crackPlaneNormal, ecsMethod);
            // remember le in corresponding status
            status->setLe(le);
        }

        // compute and store the crack angle (just for postprocessing)
        double ca = 3.1415926 / 2.;
        if ( crackPlaneNormal.at(1) != 0.0 ) {
            ca = atan( crackPlaneNormal.at(2) / crackPlaneNormal.at(1) );
        }

        status->setCrackAngle(ca);

        if ( this->gf != 0. && e0 >= ( wf / le ) ) { // case for a given fracture energy
            OOFEM_WARNING6("Fracturing strain %e is lower than the elastic strain e0=%e, possible snap-back. Element number %d, wf %e, le %e", wf / le, e0, gp->giveElement()->giveLabel(), wf, le);
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        } else if ( wf == 0. && e0 >= ef ) {
            OOFEM_WARNING5( "Fracturing strain ef=%e is lower than the elastic strain e0=%f, possible snap-back. Increase fracturing strain to %f. Element number %d", ef, e0, e0, gp->giveElement()->giveLabel() );
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        } else if ( ef == 0. && e0 * le >= wf ) {
            OOFEM_WARNING5( "Crack opening at zero stress wf=%f is lower than the elastic displacement w0=%f, possible snap-back. Increase crack opening wf to %f. Element number %d", wf, e0 * le, e0 * le, gp->giveElement()->giveLabel() );
            if ( checkSnapBack ) {
                OOFEM_ERROR1("\n");
            }
        }
    }
}

double
UmicoreMaterial :: give(int aProperty, GaussPoint *gp)
{
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        return answer;
    } else if ( aProperty == e0_ID ) {
        return this->e0;
    } else if ( aProperty == ef_ID ) {
        return this->ef;
    } else if ( aProperty == ek_ID ) {
        return this->ek;
    } else if ( aProperty == gf_ID ) {
        return this->gf;
    } else if ( aProperty == wf_ID ) {
        return this->wf;
    } else if ( aProperty == gft_ID ) {
        return this->gft;
    } else {
        return IsotropicDamageMaterial :: give(aProperty, gp);
    }
}


Interface *
UmicoreMaterial :: giveInterface(InterfaceType type)
{
    if ( type == MaterialModelMapperInterfaceType ) {
        return static_cast< MaterialModelMapperInterface * >(this);
    } else {
        return NULL;
    }
}


MaterialStatus *
UmicoreMaterial :: CreateStatus(GaussPoint *gp) const
{
    UmicoreMaterialStatus *answer = new UmicoreMaterialStatus(1, UmicoreMaterial :: domain, gp);
    return answer;
}

MaterialStatus *
UmicoreMaterial :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus( status, this->giveNumber() );
            this->_generateStatusVariables(gp);
        }
    }

    return status;
}


int
UmicoreMaterial :: MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep)
{
    int result;
    FloatArray intVal, strainIncr(3);
    IntArray toMap(3);
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );


    toMap.at(1) = ( int ) IST_MaxEquivalentStrainLevel;
    toMap.at(2) = ( int ) IST_DamageTensor;
    toMap.at(3) = ( int ) IST_StrainTensor;
    this->mapper.init(oldd, toMap, gp, tStep);

    result = mapper.mapVariable(intVal, gp, IST_MaxEquivalentStrainLevel, tStep);
    if ( result ) {
        status->setTempKappa( intVal.at(1) );
    }

    result = mapper.mapVariable(intVal, gp, IST_DamageTensor, tStep);
    if ( result ) {
        status->setTempDamage( intVal.at(1) );
    }

#ifdef IDM_USE_MAPPEDSTRAIN
    result = mapper.mapVariable(intVal, gp, IST_StrainTensor, tStep);
    if ( result ) {
        status->letTempStrainVectorBe(intVal);
    }

#endif
    status->updateYourself(tStep);

    return result;
}




int
UmicoreMaterial :: MMI_update(GaussPoint *gp,  TimeStep *tStep, FloatArray *estrain)
{
    int result = 1;
    FloatArray intVal, strain;
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );

    // now update all internal vars accordingly
    strain = status->giveStrainVector();
#ifdef IDM_USE_MAPPEDSTRAIN
    this->giveRealStressVector(intVal, gp, strain, tStep);
#else
    this->giveRealStressVector(intVal, gp, * estrain, tStep);
#endif
    gp->updateYourself(tStep);
    return result;
}


int
UmicoreMaterial :: MMI_finish(TimeStep *tStep)
{
    this->mapper.finish(tStep);
    return 1;
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
	if (deactDmgInCompression && sigmaV<0.) {
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
// copied from isodamagemodel.C
{
    UmicoreMaterialStatus *status = static_cast< UmicoreMaterialStatus * >( this->giveStatus(gp) );
    //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray strainVector, reducedTotalStrainVector;
    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;

    this->initGpForNewStep(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);

    //crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, tStep);

    if ( llcriteria == idm_strainLevelCR ) {
        // compute value of loading function if strainLevel crit apply
        f = equivStrain - status->giveKappa();

        if ( f <= 0.0 ) {
            // damage does not grow
            tempKappa = status->giveKappa();
            omega     = status->giveDamage();
        } else {
            // damage grow
            tempKappa = equivStrain;
            this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
            // evaluate damage parameter
            this->computeDamageParam(omega, tempKappa, reducedTotalStrainVector, gp);
        }
    } else if ( llcriteria == idm_damageLevelCR ) {
        // evaluate damage parameter first
        tempKappa = equivStrain;
        this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
        this->computeDamageParam(omega, tempKappa, reducedTotalStrainVector, gp);
        if ( omega < status->giveDamage() ) {
            // unloading takes place
            omega = status->giveDamage();
            //printf (".");
        }
    } else {
        _error("giveRealStressVector: unsupported loading/uloading criteria");
    }


    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
		de.times(1-omega);

    answer.beProductOf(de, reducedTotalStrainVector);

    // update gp
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
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

	 // consider full strain vector
    StructuralMaterial :: giveFullSymVectorForm( strain, reducedTotalStrainVector, gp->giveMaterialMode() );
    if ( gp->giveMaterialMode() == _PlaneStress ) {
        double nu = lmat->give(NYxz, gp);
        strain.at(3) = -nu * ( strain.at(1) + strain.at(2) ) / ( 1. - nu );
    } else if ( gp->giveMaterialMode() == _1dMat ) {
        double nu = lmat->give(NYxz, gp);
        strain.at(2) = -nu *strain.at(1);
        strain.at(3) = -nu *strain.at(1);
    }

	// elastic strain = totalStrain - plasticStrain
	 FloatArray elasticStrain = strain;
	 elasticStrain.subtract(status->plasticStrain);

	 // stiffness and compliance matrix
    FloatMatrix de;
    lmat->give3dMaterialStiffnessMatrix(de, SecantStiffness, gp, tStep);
	 FloatMatrix ce; ce.beInverseOf(de);

	// trial stress
	 stress.beProductOf(de,elasticStrain);
	 double omega = status->giveTempDamage();
	 // current tensile strength
	 double ft = lmat->give('E', gp)*std::max(e0,status->giveTempKappa())*(1-omega);
	 status->tempElasticFlag = true;
	 // test yield function with trial stress ...
	 double y = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce) - ft;
	 if (y > 0) { // ... and perform stress return for not addmisible state
		 status->tempElasticFlag = false;
			stressReturn(stress,status,ft,gp,tStep);
	 }

	// update answer and gp
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
	 StrainVector elasticStrainTrial(_3dMat);
	 elasticStrainTrial.beProductOf(ce,stress);
	 FloatArray princVals;
	 //
	 // dPlasticStrain = lambda * elasticStrainTrial
	 // stress = Ce * elasticStrain = Ce * (elasticStrainTrial - dPlasticStrain) = Ce * elasticStrainTrial * (1-lambda)
	 // f(stress) = 0, solved by bisection method for lambda
	 //
	 double lambda = .5;
	 double y0,y1,yC;
	 y0 = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce)-ft;
	 FloatArray e = elasticStrainTrial;
	 e.times(0);
	 stress.beProductOf(de,e);
	 y1 = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce)-ft;
	 if (y0<0) { _error("TODO y0"); }
	 if (y1>0) { _error("TODO y1"); }
	 double lambda0=0, lambda1=1;
	//
	 for (int i=0; i<30; i++) {
		 e = elasticStrainTrial;
		 e.times(1-lambda);
		 stress.beProductOf(de,e);
		 yC = computeStressScalarForYieldFunction(stress,gp,tStep,de,ce)-ft;
		 if (yC > 0.) lambda0 = lambda;
		else lambda1 = lambda;
		lambda = .5*(lambda0+lambda1);
	 }
	 //
		 e = elasticStrainTrial;
		 e.times(1-lambda);
		 stress.beProductOf(de,e);
	 //

	 // update gp
		FloatArray tempPlasticStrain = elasticStrainTrial;
		tempPlasticStrain.times(lambda);
		tempPlasticStrain.add(status->plasticStrain);
		status->tempPlasticStrain = tempPlasticStrain;
	 return;
}

/**
stress dependend part of yield function
*/
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
		 // copied from idm1.C
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
        return IsotropicDamageMaterial :: giveIPValue(answer, gp, type, tStep);
    }
    return 1;
}




UmicoreMaterialStatus :: UmicoreMaterialStatus(int n, Domain *d, GaussPoint *g) :
    IsotropicDamageMaterialStatus(n, d, g), RandomMaterialStatusExtensionInterface()
{
    le = 0.0;

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

Interface *
UmicoreMaterialStatus :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return static_cast< RandomMaterialStatusExtensionInterface * >(this);
    } else {
        return NULL;
    }
}


contextIOResultType
UmicoreMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
UmicoreMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
}     // end namespace oofem
