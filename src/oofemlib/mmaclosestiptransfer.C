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

#include "mmaclosestiptransfer.h"
#include "spatiallocalizer.h"
#include "domain.h"
#include "material.h"
#include "gausspoint.h"
#include "matstatmapperint.h"
#include "xfemelementinterface.h"
#include "structuralinterfacematerial.h"

namespace oofem {
MMAClosestIPTransfer :: MMAClosestIPTransfer() : MaterialMappingAlgorithm()
{ }

void
MMAClosestIPTransfer :: __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep, bool iCohesiveZoneGP)
{
    SpatialLocalizer *sl = dold->giveSpatialLocalizer();
    this->source = sl->giveClosestIP(coords, region, iCohesiveZoneGP);

    if ( iCohesiveZoneGP ) {
        XfemElementInterface *xFemEl = dynamic_cast< XfemElementInterface * >( source->giveElement() );

        if ( xFemEl == NULL ) {
            OOFEM_ERROR("In MMAClosestIPTransfer :: __init: xFemEl == NULL.\n");
        }

        mpMaterialStatus = xFemEl->mpCZMat->giveStatus(source);
    } else {
        mpMaterialStatus = source->giveMaterial()->giveStatus(source);
    }

    if ( !source ) {
        OOFEM_ERROR("MMAClosestIPTransfer::__init : no suitable source found");
    }
}

int
MMAClosestIPTransfer :: __mapVariable(FloatArray &answer, FloatArray &coords,
                                      InternalStateType type, TimeStep *tStep)
{
    if ( source ) {
        source->giveMaterial()->giveIPValue(answer, source, type, tStep);
        return 1;
    }

    return 0;
}

int
MMAClosestIPTransfer :: mapStatus(MaterialStatus &oStatus) const
{
    if ( mpMaterialStatus != NULL ) {
        MaterialStatusMapperInterface &interface = dynamic_cast< MaterialStatusMapperInterface & >(oStatus);
        interface.copyStateVariables(* mpMaterialStatus);

        return 1;
    } else {
        OOFEM_ERROR("Error in MMAClosestIPTransfer :: mapStatus(): source not set.\n");
    }

    return 0;
}
} // end namespace oofem
