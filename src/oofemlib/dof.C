/* $Header: /home/cvs/bp/oofem/oofemlib/src/dof.C,v 1.10.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   file DOF.CC

#include "dof.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"

#include "flotarry.h"
#include "dictionr.h"

#include "debug.h"
#include "cltypes.h"
#include "logger.h"
#include "datastream.h"

#ifndef __MAKEDEPEND
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#endif



Dof :: Dof (int i, DofManager* aNode, DofID id)
   // Constructor. Creates a new d.o.f., with number i, belonging
   // to aNode 
{
  number         = i ;
  dofManager     = aNode ;
  dofID          = (DofIDItem) id  ;;
}


int
Dof :: giveDofManNumber () const {return this->dofManager->giveNumber();} // termitovo


void  Dof :: printSingleOutputAt (FILE *File, TimeStep* stepN, char ch,
                                  EquationID type, ValueModeType mode, double scale) 
   // Prints in the data file the unknown 'u' (for example, the displacement
   // 'd') of the receiver, at stepN.
{
   double x ;
   x    = scale*this -> giveUnknown(type, mode, stepN) ;
   fprintf (File,"  dof %d   %c % .8e\n",number,ch,x) ;
}



void  Dof :: printMultipleOutputAt (FILE *File, TimeStep* stepN, char* ch,
                                    EquationID type, ValueModeType* mode, int nite) 
   // Prints in the data file the unknown 'u' (for example, the displacement
   // 'd') of the receiver, at stepN.
{
 int i;
 double x ;

 fprintf (File,"  dof %d", number);
 for (i=1; i<=nite; i++) {
  x    = this -> giveUnknown(type, mode[i-1], stepN) ;
  fprintf (File,"   %c % .8e",ch[i-1],x) ;
 }
 fprintf (File,"\n");
}


void  Dof :: printYourself () 
 // Prints the receiver on screen.
{
 printf ("dof %d  of %s %d :\n",number,dofManager->giveClassName(), dofManager->giveNumber()) ;
 // printf ("equation %d    bc %d \n",equationNumber,bc) ;

 // printOutputAt (node->giveDomain()->giveEngngModel()->giveCurrentStep());
}


void Dof :: error (const char* file, int line, const char *format, ...) 
{
  char buffer[MAX_ERROR_MSG_LENGTH];
	va_list args;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);

  __OOFEM_ERROR5 (file, line, "Class: %s, number %d, of DofManager %d\n%s", 
                 this->giveClassName(), number, dofManager->giveNumber(), buffer);
}

/*
UnknownType
Dof :: giveUnknownType () 
{
// returns CharType type corresponding to receiver's DofID value
 switch (this->giveDofID()) {
 case D_u:
 case D_v:
 case D_w:
 case R_u:
 case R_v:
 case R_w:
  return DisplacementVector;

 case T_f:
 case P_f:
  return FluxVector;

 case G_0:
 case G_1:
  return GeneralizedDisplacementVector;

 default:
  _error ("giveUnknownType:  undefined DofID value\n") ;
  exit (1);
 }
 return UnknownType_Unknown;
}
*/

char* 
Dof :: giveDofIDName  (char* s) 
{
// returns character string representing receiver's dofID
 switch (this->giveDofID()) {
 case D_u: 
  return strcpy (s,"D_u");

 case D_v:
  return strcpy (s,"D_v");

 case D_w:
  return strcpy (s,"D_w");

 case R_u:
  return strcpy (s,"R_u");

 case R_v:
  return strcpy (s,"R_v");

 case R_w:
  return strcpy (s,"R_w");

 case T_f:
  return strcpy (s,"T_f");

 case G_0:
  return strcpy (s,"G_0");
  
 case G_1:
  return strcpy (s,"G_1");

 default:
  sprintf (s,"%3d",this->number);
  return s;
 }
 // return s;
}

double  
Dof :: giveBcValue (ValueModeType mode, TimeStep* tStep) 
{
 double rel = 0.0;
 if (mode == VM_Incremental && tStep->isTheFirstStep() && hasIcOn(VM_Total)) 
  rel = giveIc() -> give(VM_Total);

 return this->giveBc()->give (this, mode, tStep) - rel;
}

contextIOResultType  
Dof::saveContext (DataStream* stream, ContextMode mode, void *obj)  
{
 if (stream == NULL) THROW_CIOERR(CIO_IOERR);

 if (mode & CM_Definition ) {
   int _val = dofID;

   if (!stream->write(&number,1)) THROW_CIOERR(CIO_IOERR);
   if (!stream->write(&_val,1)) THROW_CIOERR(CIO_IOERR);
 }
 return CIO_OK;
}


contextIOResultType  
Dof::restoreContext(DataStream* stream, ContextMode mode, void *obj) 
{
  if (mode & CM_Definition ) {
    int _val;

    if (!stream->read(&number,1)) THROW_CIOERR(CIO_IOERR);
    if (!stream->read(&_val,1)) THROW_CIOERR(CIO_IOERR);
    dofID = (DofIDItem) _val;
  }
  return CIO_OK;
}



  
#ifdef __PARALLEL_MODE
#endif