/* $Header: /home/cvs/bp/oofem/oofemlib/src/modulemanager.h,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   *********************************
//   *** CLASS ModuleManager ***
//   *********************************


#ifndef modulemanager_h
#define modulemanager_h

#include "alist.h"
#include "datareader.h"
#include "oofem_limits.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <time.h>
#endif

namespace oofem {

class EngngModel;

/**
 * Class representing and implementing ModuleManager. It is attribute of EngngModel.
 * It manages the modules of given type.
 */
template <class M>
class ModuleManager
{
protected:
    /// module list
    AList< M > *moduleList;
    /// number of modules
    int numberOfModules;
    /// Associated Engineering model.
    EngngModel *emodel;

public:
    ModuleManager(EngngModel *emodel) {
      this->emodel = emodel;
      
      moduleList = new AList< M >(0);
      numberOfModules = 0;
    }      
    virtual ~ModuleManager() {
      delete moduleList;
    }
    /** Creates new instance of module of given name, belonging to given EngngModel */
    virtual M* CreateModuleOfType (char *name, EngngModel *emodel) = 0;
    /**
     * Reads receiver description from input stream and creates corresponding modules components accordingly.
     * It scans input file, each line is assumed to be single record describing particular module.
     * The record line is converted to lowercase letters.
     * Corresponding component is created using ofType function.
     * After new output module object is created, its initializeForm member function is
     * called with its record as parameter.
     * @param inputStream input stream with domain description
     * @initString the e-model record containing export module manager record
     * @return nonzero if o.k.
     */
    virtual int                instanciateYourself(DataReader *dr, InputRecord *ir) {
      const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
      IRResultType result;                     // Required by IR_GIVE_FIELD macro
      
      int i;
      char name [ MAX_NAME_LENGTH ];
      M *module;
      InputRecord *mir;
      
      // read modules
      moduleList->growTo(numberOfModules);
      for ( i = 0; i < numberOfModules; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_expModuleRec, i + 1);
        result = mir->giveRecordKeywordField(name, MAX_NAME_LENGTH);
        if ( result != IRRT_OK ) {
	  IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }
	
        // read type of module
        module = this->CreateModuleOfType(name, emodel);
        if ( module == NULL ) {
            OOFEM_ERROR2("InitModuleManager::instanciateYourself: unknown module (%s)", name);
        }
	
        module->initializeFrom(mir);
        moduleList->put(i + 1, module);
      }
      
#  ifdef VERBOSE
      VERBOSE_PRINT0("Instanciated output modules ", numberOfModules)
#  endif
	return 1;
    }
    
    /**
     * Instanciates the receiver from input record. Called from instanciateYourself to initialize yourself
     * from corresponding record. Should be caled before instanciateYourself.
     */
    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    /** returns classname of the receiver */
    virtual const char *giveClassName() const { return "ModuleManager"; }
protected:

    /**
     * Returns the required module.
     * @param num module number
     */
    M *giveModule(int num) {
      M *elem = NULL;
      
      if ( moduleList->includes(num) ) {
        elem = moduleList->at(num);
      } else {
        OOFEM_ERROR2("ModuleManager::giveOuputModule: No module no. %d defined", num);
      }
      
      return elem;
    }
};

} // end namespace oofem
#endif // modulemanager_h