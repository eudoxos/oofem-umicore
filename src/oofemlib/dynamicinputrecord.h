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

#ifndef dynamicinputrecord_h
#define dynamicinputrecord_h

#include "inputrecord.h"

#include <map>
#include <set>

namespace oofem {
/**
 * Class representing the a dynamic Input Record.
 * The input record is represented as a list of fields.
 * This is intended for internal usage, where new input records and such are created dynamically.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT DynamicInputRecord : public InputRecord
{
protected:
    std :: string recordKeyword;
    int recordNumber;

    // Record representation.
    std :: set< std :: string >emptyRecord; ///< Fields without values
    std :: map< std :: string, int >intRecord;
    std :: map< std :: string, double >doubleRecord;
    std :: map< std :: string, bool >boolRecord;
    std :: map< std :: string, std :: string >stringRecord;
    std :: map< std :: string, FloatArray >floatArrayRecord;
    std :: map< std :: string, IntArray >intArrayRecord;
    std :: map< std :: string, FloatMatrix >matrixRecord;
    std :: map< std :: string, std :: vector< std :: string > >stringListRecord;
    std :: map< std :: string, Dictionary >dictionaryRecord;
    std :: map< std :: string, std :: list< Range > >rangeRecord;

public:
    /// Constructor. Creates an empty input record.
    DynamicInputRecord();
    /// Copy constructor.
    DynamicInputRecord(const DynamicInputRecord &);
    /// Destructor.
    virtual ~DynamicInputRecord();
    /// Assignment operator.
    DynamicInputRecord &operator = ( const DynamicInputRecord & );

    virtual InputRecord *GiveCopy() { return new DynamicInputRecord(*this); }
    virtual void finish(bool wrn = true);

    virtual IRResultType giveRecordKeywordField(std :: string &answer, int &value);
    virtual IRResultType giveRecordKeywordField(std :: string &answer);
    virtual IRResultType giveField(int &answer, InputFieldType id);
    virtual IRResultType giveField(double &answer, InputFieldType id);
    virtual IRResultType giveField(bool &answer, InputFieldType id);
    virtual IRResultType giveField(std :: string &answer, InputFieldType id);
    virtual IRResultType giveField(FloatArray &answer, InputFieldType id);
    virtual IRResultType giveField(IntArray &answer, InputFieldType id);
    virtual IRResultType giveField(FloatMatrix &answer, InputFieldType id);
    virtual IRResultType giveField(std :: vector< std :: string > &answer, InputFieldType id);
    virtual IRResultType giveField(Dictionary &answer, InputFieldType id);
    virtual IRResultType giveField(std :: list< Range > &answer, InputFieldType id);

    virtual bool hasField(InputFieldType id);
    virtual void printYourself();
    // Setters, unique for the dynamic input record
    virtual void setRecordKeywordField(const std :: string &keyword, int number);
    virtual void setRecordKeywordNumber(int number);
    virtual void setField(int item, InputFieldType id);
    virtual void setField(double item, InputFieldType id);
    virtual void setField(bool item, InputFieldType id);
    virtual void setField(const std :: string &item, InputFieldType id);
    virtual void setField(const FloatArray &item, InputFieldType id);
    virtual void setField(const IntArray &item, InputFieldType id);
    virtual void setField(const FloatMatrix &item, InputFieldType id);
    virtual void setField(const std :: vector< std :: string > &item, InputFieldType id);
    virtual void setField(const Dictionary &item, InputFieldType id);
    virtual void setField(const std :: list< Range > &item, InputFieldType id);
    /// Sets an empty field with given id.
    virtual void setField(InputFieldType id);
    /// Removes given field from record.
    virtual void unsetField(InputFieldType id);

    virtual void report_error(const char *_class, const char *proc, InputFieldType id,
                              IRResultType result, const char *file, int line);

    /// Returns record as string.
    std :: string giveRecordAsString() const;
};
} // end namespace oofem
#endif // dynamicinputrecord_h
