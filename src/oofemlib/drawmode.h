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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//
// FILE: drawmode.h
//

#ifndef drawmode_h
#define drawmode_h

enum DrawMode {
    unknown,

    rawGeometry,
    deformedGeometry,
    eigenVectorGeometry,
    nodeAnnotation,
    appliedPrimaryBc,

    internalStateBegin,
    mxForce,
    myForce,
    mzForce,
    myzForce,
    mzxForce,
    mxyForce,
    //qxzForce,
    //qyzForce,
    sxForce,
    syForce,
    szForce,
    syzForce,
    szxForce,
    sxyForce,
    yieldState,
    crackedState,
    stressErrorState,
    requiredAdaptiveMeshSizeState,
    damageLevel,
    errorIndicatorLevel,
    relativeMeshSizeDensity,
    temperatureField,
    massConcentration1Field,
    velocityField,
    pressureField,
    vofField,
    densityField,

    hydrationDegreeState,
    humidityState,

    internalStateEnd
};

#endif // drawmode_h