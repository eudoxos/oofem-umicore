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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "abaqususermaterial.h"
#include "gausspnt.h"

#ifdef _WIN32 //_MSC_VER and __MINGW32__ included
 #include <Windows.h>
#else
 #include <dlfcn.h>
#endif

#include <cstring>

namespace oofem {
int AbaqusUserMaterial :: n = 1;

AbaqusUserMaterial :: ~AbaqusUserMaterial()
{
#ifdef _WIN32
    if ( this->umatobj ) {
        FreeLibrary( ( HMODULE ) this->umatobj );
    }
#else
    if ( this->umatobj ) {
        dlclose(this->umatobj);
    }

#endif
}

IRResultType AbaqusUserMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;
    std :: string umatname;
    std :: string umatfile;

    this->Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->numState, IFT_AbaqusUserMaterial_numState, "numstate");
    IR_GIVE_FIELD(ir, this->properties, IFT_AbaqusUserMaterial_properties, "properties");
    IR_GIVE_FIELD(ir, umatfile, IFT_AbaqusUserMaterial_userMaterial, "umat");
    umatname = "umat";
    IR_GIVE_OPTIONAL_FIELD(ir, umatname, IFT_AbaqusUserMaterial_userMaterial, "name");
    strncpy(this->cmname, umatname.c_str(), 80);

#ifdef _WIN32
    ///@todo Check all the windows support.
    this->umatobj = ( void * ) LoadLibrary( umatfile.c_str() );
    if ( !this->umatobj ) {
        OOFEM_ERROR2( "AbaqusUserMaterial :: initializeFrom - couldn't load \"%s\",\ndlerror: %s", umatfile.c_str() );
    }

//     * ( void ** )( & this->umat ) = GetProcAddress( ( HMODULE ) this->umatobj, "umat_" );
    *(FARPROC *)( & this->umat ) = GetProcAddress( ( HMODULE ) this->umatobj, "umat_" );//works for MinGW 32bit
    if ( !this->umat ) {
//         char *dlresult = GetLastError();
        DWORD dlresult = GetLastError();//works for MinGW 32bit
        OOFEM_ERROR2("AbaqusUserMaterial :: initializeFrom - couldn't load symbol umat,\nerror: %s\n", dlresult);
    }
    
#else
    this->umatobj = dlopen(umatfile.c_str(), RTLD_NOW);
    if ( !this->umatobj ) {
        OOFEM_ERROR3( "AbaqusUserMaterial :: initializeFrom - couldn't load \"%s\",\ndlerror: %s", umatfile.c_str(), dlerror() );
    }

    * ( void ** ) ( & this->umat ) = dlsym(this->umatobj, "umat_");
    char *dlresult = dlerror();
    if ( dlresult ) {
        OOFEM_ERROR2("AbaqusUserMaterial :: initializeFrom - couldn't load symbol umat,\ndlerror: %s\n", dlresult);
    }

#endif
    return IRRT_OK;
}

MaterialStatus *AbaqusUserMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new AbaqusUserMaterialStatus(n++, this->giveDomain(), gp, this->numState);
}

void AbaqusUserMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                    MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    AbaqusUserMaterialStatus *ms = dynamic_cast< AbaqusUserMaterialStatus * >( this->giveStatus(gp) );
    if ( !ms->hasTangent() ) { ///@todo Make this hack fit more nicely into OOFEM in general;
        // Evaluating the function once, so that the tangent can be obtained.
        MaterialMode mMode = gp->giveMaterialMode();
        int ncomp = 0;
        if ( mMode == _3dMat ) {
            ncomp = 6;
        } else if ( mMode == _PlaneStress ) {
            ncomp = 3;
        } else if ( mMode == _PlaneStrain ) {
            ncomp = 4;
        } /*else if ( mMode == _3dMat_F ) {
           * ncomp = 9;
           * } */
        else if ( mMode == _1dMat ) {
            ncomp = 1;
        }

        FloatArray stress(ncomp), strain(ncomp);
        strain.zero();
        this->giveRealStressVector(stress, form, gp, strain, tStep);
    }

    answer = ms->giveTempTangent();

#if 0
    double h = 1e-7;
    FloatArray strain, strainh, stress, stressh;
    strain = ( ( StructuralMaterialStatus * ) gp->giveMaterialStatus(AbaqusUserMaterialClass) )->giveTempStrainVector();
    stress = ( ( StructuralMaterialStatus * ) gp->giveMaterialStatus(AbaqusUserMaterialClass) )->giveTempStressVector();
    FloatMatrix En( strain.giveSize(), strain.giveSize() );
    for ( int i = 1; i <= strain.giveSize(); ++i ) {
        strainh = strain;
        strainh.at(i) += h;
        this->giveRealStressVector(stressh, form, gp, strainh, tStep);
        stressh.subtract(stress);
        stressh.times(1.0 / h);
        En.setColumn(stressh, i);
    }

    printf("En = ");
    En.printYourself();
    printf("Tangent = ");
    answer.printYourself();
#endif
}

void AbaqusUserMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                const FloatArray &reducedStrain, TimeStep *tStep)
{
    AbaqusUserMaterialStatus *ms = static_cast< AbaqusUserMaterialStatus * >( this->giveStatus(gp) );
    /* User-defined material name, left justified. Some internal material models are given names starting with
     * the “ABQ_” character string. To avoid conflict, you should not use “ABQ_” as the leading string for
     * CMNAME. */
    //char cmname[80];

    MaterialMode mMode = gp->giveMaterialMode();
    // Sizes of the tensors.
    int ndi;
    int nshr;
    ///@todo Check how to deal with large deformations.
    ///@todo Check order of entries in the Voigt notation (which order does Abaqus use? convert to that).
    if ( mMode == _3dMat ) {
        ndi = 3;
        nshr = 3;
    } else if ( mMode == _PlaneStress ) {
        ndi = 2;
        nshr = 1;
    } else if ( mMode == _PlaneStrain ) {
        ndi = 3;
        nshr = 1;
    } /*else if ( mMode == _3dMat_F ) {
       * ndi = 3;
       * nshr = 6;
       * } */
    else if ( mMode == _1dMat ) {
        ndi = 1;
        nshr = 0;
    } else {
        ndi = nshr = 0;
        OOFEM_ERROR2( "AbaqusUserMaterial :: giveRealStressVector : unsupported material mode (%s)", __MaterialModeToString(mMode) );
    }

    int ntens = ndi + nshr;
    FloatArray stress(ntens);
    FloatArray strain = ms->giveStrainVector();
    FloatArray strainIncrement;
    strainIncrement.beDifferenceOf(reducedStrain, strain);
    FloatArray state = ms->giveStateVector();
    FloatMatrix jacobian(ntens, ntens);
    int numProperties = this->properties.giveSize();

    // Temperature and increment
    double temp = 0.0, dtemp = 0.0;

    // Times and increment
    double dtime = tStep->giveTimeIncrement();
    ///@todo Check this. I'm just guessing. Maybe intrinsic time instead?
    double time [ 2 ] = {
        tStep->giveTargetTime() - dtime, tStep->giveTargetTime()
    };
    double pnewdt = 1.0; ///@todo Right default value? umat routines may change this (although we ignore it)

    /* Specific elastic strain energy, plastic dissipation, and “creep” dissipation, respectively. These are passed
     * in as the values at the start of the increment and should be updated to the corresponding specific energy
     * values at the end of the increment. They have no effect on the solution, except that they are used for
     * energy output. */
    double sse, spd, scd;

    // Outputs only in a fully coupled thermal-stress analysis:
    double rpl = 0.0; // Volumetric heat generation per unit time at the end of the increment caused by mechanical working of the material.
    FloatArray ddsddt(ntens); // Variation of the stress increments with respect to the temperature.
    FloatArray drplde(ntens); // Variation of RPL with respect to the strain increments.
    double drpldt = 0.0; // Variation of RPL with respect to the temperature.

    /* An array containing the coordinates of this point. These are the current coordinates if geometric
     * nonlinearity is accounted for during the step (see “Procedures: overview,” Section 6.1.1); otherwise,
     * the array contains the original coordinates of the point */
    FloatArray coords;
    gp->giveElement()->computeGlobalCoordinates( coords, * gp->giveCoordinates() ); ///@todo Large deformations?

    /* Rotation increment matrix. This matrix represents the increment of rigid body rotation of the basis
     * system in which the components of stress (STRESS) and strain (STRAN) are stored. It is provided so
     * that vector- or tensor-valued state variables can be rotated appropriately in this subroutine: stress and
     * strain components are already rotated by this amount before UMAT is called. This matrix is passed in
     * as a unit matrix for small-displacement analysis and for large-displacement analysis if the basis system
     * for the material point rotates with the material (as in a shell element or when a local orientation is used).*/
    FloatMatrix drot(3, 3);
    drot.beUnitMatrix();

    /* Characteristic element length, which is a typical length of a line across an element for a first-order
     * element; it is half of the same typical length for a second-order element. For beams and trusses it is a
     * characteristic length along the element axis. For membranes and shells it is a characteristic length in
     * the reference surface. For axisymmetric elements it is a characteristic length in the
     * plane only.
     * For cohesive elements it is equal to the constitutive thickness.*/
    double celent = 0.0; /// @todo Include the characteristic element length

    /* Array containing the deformation gradient at the beginning of the increment. See the discussion
     * regarding the availability of the deformation gradient for various element types. */
    FloatMatrix dfgrd0(3, 3);
    /* Array containing the deformation gradient at the end of the increment. The components of this array
     * are set to zero if nonlinear geometric effects are not included in the step definition associated with
     * this increment. See the discussion regarding the availability of the deformation gradient for various
     * element types. */
    FloatMatrix dfgrd1(3, 3);

    int noel = gp->giveElement()->giveNumber(); // Element number.
    int npt = 0; // Integration point number.

    int layer = 0; // Layer number (for composite shells and layered solids)..
    int kspt = 0; // Section point number within the current layer.
    int kstep = 0; // Step number.
    int kinc = 0; // Increment number.

    ///@todo No idea about these parameters
    double predef;
    double dpred;

    OOFEM_LOG_DEBUG("AbaqusUserMaterial :: giveRealStressVector - Calling subroutine");
    this->umat(stress.givePointer(), // STRESS
               state.givePointer(), // STATEV
               jacobian.givePointer(), // DDSDDE
               & sse, // SSE
               & spd, // SPD
               & scd, // SCD
               & rpl, // RPL
               ddsddt.givePointer(), // DDSDDT
               drplde.givePointer(), // DRPLDE
               & drpldt, // DRPLDT
               strain.givePointer(), // STRAN
               strainIncrement.givePointer(), // DSTRAN
               time, // TIME
               & dtime, // DTIME
               & temp, // TEMP
               & dtemp, // DTEMP
               & predef, // PREDEF
               & dpred, // DPRED
               this->cmname, // CMNAME
               & ndi, // NDI
               & nshr, // NSHR
               & ntens, // NTENS
               & numState, // NSTATV
               properties.givePointer(), // PROPS
               & numProperties, // NPROPS
               coords.givePointer(), // COORDS
               drot.givePointer(), // DROT
               & pnewdt, // PNEWDT
               & celent, // CELENT
               dfgrd0.givePointer(), // DFGRD0
               dfgrd1.givePointer(), // DFGRD1
               & noel, // NOEL
               & npt, // NPT
               & layer, // LAYER
               & kspt, // KSPT
               & kstep, // KSTEP
               & kinc); // KINC

    ms->letTempStrainVectorBe(reducedStrain);
    ms->letTempStressVectorBe(stress);
    ms->letTempStateVectorBe(state);
    ms->letTempTangentBe(jacobian);
    answer = stress;

    OOFEM_LOG_DEBUG("AbaqusUserMaterial :: giveRealStressVector - Calling subroutine was successful");
}

int AbaqusUserMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _1dMat;
}

void AbaqusUserMaterialStatus :: initTempStatus()
{
    stateVector.resize(numState);
    stateVector.zero();
    tempStateVector.resize(numState);
    stateVector.zero();
}

AbaqusUserMaterialStatus :: AbaqusUserMaterialStatus(int n, Domain *d, GaussPoint *gp, int numState) :
    StructuralMaterialStatus(n, d, gp),
    numState(numState), stateVector(numState), hasTangentFlag(false)
{
    this->initTempStatus();
}

void AbaqusUserMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    stateVector = tempStateVector;
}
} // end namespace oofem