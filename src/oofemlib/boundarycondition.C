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

#include "boundarycondition.h"
#include "timestep.h"
#include "loadtimefunction.h"
#include "verbose.h"
#include "classfactory.h"

namespace oofem {

REGISTER_BoundaryCondition( BoundaryCondition );

double BoundaryCondition :: give(Dof *dof, ValueModeType mode, TimeStep *stepN)
// Returns the value at stepN of the prescribed value of the kinematic
// unknown 'u'. Returns 0 if 'u' has no prescribed value.
{
    double factor;

    factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);
    return prescribedValue * factor;
}


IRResultType
BoundaryCondition :: initializeFrom(InputRecord *ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    GeneralBoundaryCondition :: initializeFrom(ir);

    if ( ir->hasField(_IFT_BoundaryCondition_PrescribedValue) ) {
        IR_GIVE_FIELD(ir, prescribedValue, _IFT_BoundaryCondition_PrescribedValue);
    } else {
        IR_GIVE_FIELD(ir, prescribedValue, _IFT_BoundaryCondition_PrescribedValue_d);
    }

    return IRRT_OK;
}



int
BoundaryCondition :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    GeneralBoundaryCondition :: giveInputRecordString(str, keyword);
    sprintf(buff, " prescribedvalue %e", this->prescribedValue);
    str += buff;

    return 1;
}
} // end namespace oofem