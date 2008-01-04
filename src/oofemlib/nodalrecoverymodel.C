/* $Header: /home/cvs/bp/oofem/oofemlib/src/nodalrecoverymodel.C,v 1.9.4.1 2004/04/05 15:19:43 bp Exp $ */
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

//
// file nodalrecoverymodel.C
//

#include "nodalrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "elementside.h"
#include "crosssection.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

NodalRecoveryModel::~NodalRecoveryModel()
{
 //printf ("NodalRecoveryModel::~NodalRecoveryModel()\n");
}


int
NodalRecoveryModel :: clear () 
{
 int i, nnodes = this->nodalValList.giveSize();
 FloatArray* valArray;

 if (nnodes) {
  TDictionaryIterator<int, FloatArray> iterator (this->nodalValList.at(1));
  for (i=1; i<= nnodes; i++) {
   // this->nodalValList->at(i)->clear(); 
   iterator.initialize (this->nodalValList.at(i));
   while ((valArray = iterator.next())) {
    // test if INCLUDED
    if (valArray) valArray->resize(0);
   }
  }
 }
 return 1;
} 
 
int
NodalRecoveryModel :: giveNodalVector (const FloatArray*& answer, int node, int region)
{
 if (this->includes (node, region)) {
  answer = this->nodalValList.at(node)->at(region);
  if (answer->giveSize()) return 1;
 } else {
  answer = NULL;
 }
 return 0;
}

FloatArray* 
NodalRecoveryModel :: giveNodalVectorPtr (int node, int region)
{
 if (this->includes (node, region))
  return this->nodalValList.at(node)->at(region);
 else {
   OOFEM_ERROR3 ("NodalRecoveryModel::giveNodalVectorPtr: No such record exist (dofMan %d, region %d)\n",node, region);
 }
 return NULL;
}

int 
NodalRecoveryModel :: includes (int node, int region)
{
 vectorDictType* dict;
 if (this->nodalValList.includes(node)) {
  dict = this->nodalValList.at(node);
  if (dict->includes(region)) return 1;
 }
 return 0;
}

int
NodalRecoveryModel :: init ()
{
 /* Initializes receiver for given domain */
 int i, nnodes = domain->giveNumberOfDofManagers();
 vectorDictType* dict;

 // allocate array of nodal dictionaries, containing nodal values for each region
 if (this->nodalValList.giveSize() != nnodes) {
  // printf("<NodalRecoveryModel:clearing>\n");
  this->nodalValList.clear ();
  this->nodalValList.growTo(nnodes);
 }
 
 for (i=1; i<= nnodes; i++) {
  if (this->nodalValList.at(i) == NULL) { 
   // printf("d+");
   dict = new vectorDictType ();
   this->nodalValList.put (i,dict);
  }
 }

 this->stateCounter = 0;
 //
 return 1;
}

//int NodalRecoveryModel :: giveElementRegion (Element* element) {return element->giveCrossSection()->giveNumber();} 
//int NodalRecoveryModel :: giveNumberOfRegions () {return domain->giveNumberOfCrossSectionModels();}


int
NodalRecoveryModel :: updateRegionRecoveredValues (const int ireg, const IntArray& regionNodalNumbers, 
                           int regionValSize, const FloatArray& rhs)
{
 int i, node, nnodes = domain->giveNumberOfDofManagers();
 FloatArray* nodalVal;
 
 // update recovered values 
 for (node=1; node<= nnodes; node++) {
  // find nodes in region
  if (regionNodalNumbers.at(node)) {
   if (this->includes (node, ireg)) {
    nodalVal = this->giveNodalVectorPtr (node, ireg);
    if (nodalVal) {
     // override record
     nodalVal->resize (regionValSize);
    } else {
     // create new entry 
     // printf ("a+");
     nodalVal = new FloatArray (regionValSize);
     this->nodalValList.at(node)->add (ireg, nodalVal);
    }
   } else {
    // create new entry 
    // printf ("a+");
    nodalVal = new FloatArray (regionValSize);
    this->nodalValList.at(node)->add (ireg, nodalVal);
   }
   for (i = 1; i<= regionValSize; i++) nodalVal->at(i) = rhs.at((regionNodalNumbers.at(node)-1)*regionValSize + i);
  }
 } // end update recovered values 
 
 return 1;
}

int
NodalRecoveryModel :: initRegionNodeNumbering (IntArray& regionNodalNumbers, int& regionDofMans, int reg)
{
 int ielem, nelem = domain->giveNumberOfElements();
 int nnodes = domain->giveNumberOfDofManagers();
 int elemNodes;
 int elementNode, node;
 Element* element;
 
 regionNodalNumbers.resize(nnodes);
 regionNodalNumbers.zero();
 regionDofMans = 0;

 for (ielem = 1; ielem<= nelem; ielem++) {
  element = domain->giveElement (ielem);
  if (element->giveRegionNumber() != reg) continue;
  elemNodes = element->giveNumberOfDofManagers ();
  
  // determine local region node numbering
  for (elementNode = 1; elementNode<= elemNodes; elementNode++) {
   node = element->giveDofManager(elementNode)->giveNumber();
   if (regionNodalNumbers.at(node) == 0) // assign new number
    regionNodalNumbers.at(node) = ++regionDofMans;
  }
 }
 return 1;
}

int
NodalRecoveryModel::giveRegionRecordSize (int reg, InternalStateType type)
{
 int ielem, nelem = domain->giveNumberOfElements();
 //Element* element;

 for (ielem = 1; ielem<= nelem; ielem++) {
  if (domain->giveElement (ielem)->giveRegionNumber() == reg) {
   return domain->giveElement (ielem)->giveIPValueType(type);
  }
 }
 OOFEM_WARNING2 ("NodalRecoveryModel::giveRegionRecordSize: bad region number (%d) or no element in given region\n",reg);
 return 0;
}



void
NodalRecoveryModel::giveRegionRecordMap (IntArray& answer, int reg, InternalStateType type)
{
 int ielem, nelem = domain->giveNumberOfElements();
 //Element* element;

 for (ielem = 1; ielem<= nelem; ielem++) {
  if (domain->giveElement (ielem)->giveRegionNumber() == reg) {
   domain->giveElement (ielem)->giveIntVarCompFullIndx(answer, type);
   return;
  }
 }
 OOFEM_WARNING2 ("NodalRecoveryModel::giveRegionRecordMap: bad region number (%d) or no element in given region\n",reg);
 answer.resize(0);
}