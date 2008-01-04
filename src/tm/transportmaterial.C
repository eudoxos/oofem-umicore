/* $Header: */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

#include "transportmaterial.h"

void
TransportMaterial::updateInternalState (const FloatArray& vec, GaussPoint* gp, TimeStep*)
{
 TransportMaterialStatus *ms = (TransportMaterialStatus*) this->giveStatus(gp);
 if (ms) ms->letTempStateVectorBe (vec);
}


TransportMaterialStatus::TransportMaterialStatus (int n,Domain* d, GaussPoint* g) : 
MaterialStatus (n,d,g),stateVector(),tempStateVector()
{}

void  TransportMaterialStatus :: printOutputAt (FILE * File, TimeStep* tNow)
// Prints the strains and stresses on the data file.
{
 FloatArray helpVec;
  int  i,n ;

  MaterialStatus::printOutputAt (File,tNow);
  
  fprintf (File,"  state vector ") ;
 //((StructuralCrossSection*)
 // gp->giveCrossSection())->giveFullCharacteristicVector(helpVec, gp, strainVector);
  //n = helpVec.giveSize() ;
  //for (i=1 ; i<=n ; i++)
 //  fprintf (File," % .4e",helpVec.at(i)) ;
  n = stateVector.giveSize();
 for (i=1 ; i<=n ; i++)
  fprintf (File," % .4e", stateVector.at(i)) ;
 
  fprintf (File,"\n") ;
  
}

void  TransportMaterialStatus :: updateYourself (TimeStep* tStep)
// Performs end-of-step updates.
{
  MaterialStatus::updateYourself (tStep);
 stateVector = tempStateVector;
}


void
TransportMaterialStatus :: initTempStatus () 
//
// initialize record at the begining of new load step
//
{
  MaterialStatus::initTempStatus ();
 
 tempStateVector = stateVector;
}


contextIOResultType
TransportMaterialStatus :: saveContext (DataStream* stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
// saving the data in  TDictionary is left to material (yield crit. level).
{
 contextIOResultType iores;
 if (stream == NULL) _error ("saveContex : can't write into NULL stream");

 if ((iores = MaterialStatus::saveContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = stateVector.storeYourself(stream, mode)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}


contextIOResultType
TransportMaterialStatus :: restoreContext (DataStream* stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
  // FloatArray *s;
 contextIOResultType iores;
 if (stream == NULL) _error ("saveContex : can't write into NULL stream");

 if ((iores = MaterialStatus::restoreContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = stateVector.restoreYourself(stream, mode)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}



int
TransportMaterial::giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime)
// IST_Humidity must be overriden!
{
 if ((type == IST_Temperature) || (type == IST_MassConcentration_1) || (type == IST_Humidity) ) {
  FloatArray vec = ( (TransportMaterialStatus*) this->giveStatus(aGaussPoint) ) -> giveStateVector();
  answer.resize(1);
  answer.at(1) = vec.at((type == IST_Temperature)?1:2);
  return 1;
 } else return Material::giveIPValue (answer, aGaussPoint, type, atTime);
}

InternalStateValueType
TransportMaterial::giveIPValueType (InternalStateType type)
{
 if ((type == IST_Temperature) || (type == IST_MassConcentration_1) || (type == IST_Humidity) ) return ISVT_SCALAR;
 else return Material::giveIPValueType (type);
}

int
TransportMaterial::giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode)
{
 if ((type == IST_Temperature) || (type == IST_MassConcentration_1) || (type == IST_Humidity) ) {
  answer.resize(1);
  answer.at(1) = 1;
  return 1;
 } else
  return Material::giveIntVarCompFullIndx (answer, type, mmode);
}

int
TransportMaterial::giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint)
{
 if ((type == IST_Temperature) || (type == IST_MassConcentration_1) || (type == IST_Humidity) ) {
  return 1;
 } else
  return Material::giveIPValueSize (type, aGaussPoint);
}