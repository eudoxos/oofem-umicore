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

#include "gausspoint.h"
#include "inputrecord.h"
#include "domain.h"
#include "material.h"
#include "randomfieldgenerator.h"
#include "randommaterialext.h"
#include "dynamicinputrecord.h"

namespace oofem {
bool
RandomMaterialStatusExtensionInterface :: _giveProperty(int key, double &value)
{
    if ( randProperties.includes(key) ) {
        value = randProperties.at(key);
        return true;
    } else {
        return false;
    }
}

void
RandomMaterialStatusExtensionInterface :: _setProperty(int key, double value)
{
    randProperties.at(key) = value;
}


IRResultType
RandomMaterialExtensionInterface :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    randVariables.resize(0);
    randomVariableGenerators.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, randVariables, _IFT_RandomMaterialExt_randVariables);
    IR_GIVE_OPTIONAL_FIELD(ir, randomVariableGenerators, _IFT_RandomMaterialExt_randGen);

    if ( randVariables.giveSize() != randomVariableGenerators.giveSize() ) {
        OOFEM_ERROR("RandomMaterialExtensionInterface::_initializeFrom: Incompatible size of randvars and randdist attrs");
    }

    return IRRT_OK;
}


void
RandomMaterialExtensionInterface :: giveInputRecord(DynamicInputRecord &ir)
{
    ir.setField(this->randVariables, _IFT_RandomMaterialExt_randVariables);
    ir.setField(this->randomVariableGenerators, _IFT_RandomMaterialExt_randGen);
}


bool
RandomMaterialExtensionInterface :: give(int key, GaussPoint *gp, double &value)
{
    ///@todo Have to wrap it through the material to ensure that it gets an actual material status (for now at least). Doesn't play nicely with layered cross-sections.
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterial()->giveStatus(gp) );
    RandomMaterialStatusExtensionInterface *interface = dynamic_cast< RandomMaterialStatusExtensionInterface * >
                                                        ( status->giveInterface(RandomMaterialStatusExtensionInterfaceType) );
    return interface->_giveProperty(key, value);
}

void
RandomMaterialExtensionInterface :: _generateStatusVariables(GaussPoint *gp) const
{
    // Have to wrap it through the material to ensure that it gets an actual material status (for now at least)
    int size = randVariables.giveSize();
    double value;
    MaterialStatus *matStat = static_cast< MaterialStatus * >( gp->giveMaterial()->giveStatus(gp) );
    RandomMaterialStatusExtensionInterface *status = static_cast< RandomMaterialStatusExtensionInterface * >
                                                     ( matStat->giveInterface(RandomMaterialStatusExtensionInterfaceType) );

    for ( int i = 1; i <= size; i++ ) {
        gp->giveElement()->giveDomain()->
        giveRandomFieldGenerator( randomVariableGenerators.at(i) )->generateRandomValueAt(value, gp);
        status->_setProperty(randVariables.at(i), value);
    }
}
} // end namespace oofem
