/* $Header: /home/cvs/bp/oofem/oofemlib/src/primaryfield.C,v 1.2.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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

#include "primaryfield.h"
#include "spatiallocalizer.h"
#include "dofmanager.h"
#include "dof.h"
#include "element.h"
#include "timestep.h"
#include "datastream.h"

PrimaryField::PrimaryField (EngngModel* a, int idomain, 
                            FieldBaseID ft, EquationID ut, int nHist) : Field (ft), solutionVectors(nHist+1), solStepList(nHist+1)
{
 FloatArray *sv;

 this->actualStepNumber = -999;
 this->actualStepIndx = 0;
 this->nHistVectors = nHist;
 this->ut = ut;
 //atomicFields = atomFieds;

 emodel = a;
 domainIndx = idomain;

 for (int i=0; i<=nHist; i++) {
  sv = new FloatArray();
  solutionVectors.put(i+1, sv);
 }
}

PrimaryField :: ~PrimaryField()
{
}

double
PrimaryField :: giveUnknownValue (Dof* dof, ValueModeType mode, TimeStep* atTime) 
{
 int eq = dof->giveEquationNumber();
 if (eq == 0) _error ("giveUnknownValue: invalid equation number");
 
 if (mode == VM_Total) 
  return this->giveSolutionVector (atTime) -> at(eq);
 else if (mode == VM_Incremental) {
  int indxm1 = this->resolveIndx (atTime, -1);
  return (this->giveSolutionVector (atTime) -> at(eq) - this->giveSolutionVector (indxm1) -> at(eq));
 } else _error ("giveUnknownValue: unsupported mode");
 return 0.0;
}

int
PrimaryField :: evaluateAt (FloatArray& answer, FloatArray& coords, IntArray& dofId, 
              ValueModeType mode, TimeStep* atTime) 
{
 Element* bgelem;
 Domain* domain = emodel->giveDomain(domainIndx);
 SpatialLocalizer* sl = domain->giveSpatialLocalizer();
 // locate background element
 if ((bgelem = sl->giveElementContainingPoint (coords)) == NULL) {
   //_error ("PrimaryField::evaluateAt: point not found in domain\n");
   return 1;
 }
 
 EIPrimaryFieldInterface* interface = (EIPrimaryFieldInterface*) (bgelem->giveInterface(EIPrimaryFieldInterfaceType));
 if (interface) {
  interface -> EIPrimaryFieldI_evaluateFieldVectorAt (answer, *this, coords, dofId, mode, atTime);
 } else {
  _error ("ScalarPrimaryField::operator(): background element does not support EIPrimaryFiledInterface\n");
 }
 return 0; // ok
}


FloatArray* 
PrimaryField::giveSolutionVector (TimeStep* atTime)
{
 return this->giveSolutionVector (resolveIndx(atTime,0));
}

FloatArray* 
PrimaryField::giveSolutionVector (int i)
{
 FloatArray* answer = NULL;
 if ((i>=1) && (i<=(nHistVectors+1)))
  answer = solutionVectors.at(i); // alist 1-based access 
 else
  _error ("giveSolutionVector: index out of range");
 return answer; 
}



int 
PrimaryField::resolveIndx (TimeStep* atTime, int shift) 
{
 int stepNo = atTime->giveNumber();
 int relPos = actualStepNumber - stepNo - shift;
 if ((relPos>=0)&&(relPos<=nHistVectors)) {
  return (actualStepIndx+relPos)%(nHistVectors+1) + 1;
 } else {
   _error3 ("resolveIndx: History not available for relative step no. %d to step no. %d", shift, stepNo);
 }
 return 0;
}


void
PrimaryField::advanceSolution (TimeStep* atTime)
{
 TimeStep* newts;
 if ((actualStepNumber >= 0) && (actualStepNumber+1 != atTime->giveNumber()))
  _error ("advanceSolution: can not advance due to steps skipped");
 actualStepIndx = (actualStepIndx>0)?actualStepIndx-1:nHistVectors;
 actualStepNumber = atTime->giveNumber();
 if ((newts = solStepList.at(actualStepIndx+1))) {
  *newts = *atTime;
 } else {
  solStepList.put (actualStepIndx+1, new TimeStep(*atTime));
 }
}


contextIOResultType    
PrimaryField::saveContext (DataStream* stream, ContextMode mode)
{
  int i, type_id = PrimaryFieldClass;
  contextIOResultType iores;
  // write class header
  if (!stream->write(&type_id,1)) return CIO_IOERR;

  if (!stream->write(&actualStepNumber,1)) THROW_CIOERR(CIO_IOERR);
  if (!stream->write(&actualStepIndx,1)) THROW_CIOERR(CIO_IOERR);

  for (i=0; i<=nHistVectors; i++) 
    if ((iores = solutionVectors.at(i+1)->storeYourself(stream,mode))!= CIO_OK) THROW_CIOERR(iores);
 
  TimeStep *iStep;
  int flag;
  for (i=0; i<=nHistVectors; i++) {
    if ((iStep = solStepList.at(i+1))) flag = 1; else flag = 0;
    if (!stream->write(&flag,1)) THROW_CIOERR(CIO_IOERR);
    if (flag) 
      if ((iores = solStepList.at(i+1)->saveContext(stream,mode))!= CIO_OK) THROW_CIOERR(iores);
  }
  return CIO_OK;
}

contextIOResultType    
PrimaryField::restoreContext(DataStream* stream, ContextMode mode)
{
  int i, class_id;
 contextIOResultType iores;
  // read class header
  if (!stream->read(&class_id,1)) return CIO_IOERR;
  if (class_id != PrimaryFieldClass) return CIO_BADVERSION;

 if (!stream->read (&actualStepNumber,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->read (&actualStepIndx,1)) THROW_CIOERR(CIO_IOERR);

 for (i=0; i<=nHistVectors; i++) 
  if ((iores = solutionVectors.at(i+1)->restoreYourself(stream,mode))!= CIO_OK) THROW_CIOERR(iores);
 
 int flag;
 TimeStep *iStep;
 for (i=0; i<=nHistVectors; i++) {
  if (!stream->read (&flag,1)) THROW_CIOERR(CIO_IOERR);
  if (flag) {
   iStep = new TimeStep(emodel);
   if ((iores = iStep->restoreContext(stream,mode))!= CIO_OK) THROW_CIOERR(iores);
  } else {
    iStep = NULL;
  }
  solStepList.put (i+1, iStep);
 }
 return CIO_OK;
}
