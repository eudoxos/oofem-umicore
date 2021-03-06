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

#include "constantpressureload.h"
#include "loadtimefunction.h"
#include "floatarray.h"
#include "timestep.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(ConstantPressureLoad);

ConstantPressureLoad :: ConstantPressureLoad(int i, Domain *d) : BoundaryLoad(i, d) {
    this->loadOffset = 0.0;
}

IRResultType
ConstantPressureLoad :: initializeFrom(InputRecord *ir)
{
    BoundaryLoad :: initializeFrom(ir);
    if ( componentArray.giveSize() != nDofs ) {
        _error("instanciateFrom: componentArray size mismatch");
    }

    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    IR_GIVE_OPTIONAL_FIELD(ir, this->loadOffset, _IFT_ConstantPressureLoad_LoadOffset);
    return IRRT_OK;
}

void
ConstantPressureLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode)
{
    // we overload general implementation on the boundary load level due
    // to implementation efficiency

    double factor;

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        _error("computeValueAt: mode not supported");
    }

    // ask time distribution

    /*
     * factor = this -> giveLoadTimeFunction() -> at(tStep->giveTime()) ;
     * if ((mode==VM_Incremental) && (!tStep->isTheFirstStep()))
     * //factor -= this->giveLoadTimeFunction()->at(tStep->givePreviousStep()->giveTime()) ;
     * factor -= this->giveLoadTimeFunction()->at(tStep->giveTime()-tStep->giveTimeIncrement()) ;
     */
    factor = this->giveLoadTimeFunction()->evaluate(tStep, mode);
    answer = componentArray;
    answer.times(factor);
}
} // end namespace oofem
