#pragma once
#include "../fields/COutFieldCollection.hpp"

class COutFieldManager {

protected:
  mutable COutFieldCollection fields;
  /*! \brief Vector to cache the positions of the field in the data array */
  mutable std::vector<short>                            cacheIndexVector;
  /*! \brief Current value of the cache index */
  mutable unsigned short                                cacheIndex;
  /*! \brief Boolean to store whether the field index cache should be build. */
  mutable bool                                          buildIndexCache;

  mutable bool cacheEnabled = false;
public:
  using FieldRefVector = COutFieldCollection::InsertionVector;
  using FieldRef       = COutFieldCollection::Map::iterator;

  COutFieldManager() = default;

  inline void SetFieldValue(const std::string& refName, su2double value){
    if (cacheEnabled){
      if (buildIndexCache){
        const int index = fields.GetIndex(refName);
        cacheIndexVector.push_back(index);
        fields.SetValueByIndex(index, value);
      } else {
        const int index = cacheIndexVector[cacheIndex++];
        if (cacheIndex == cacheIndexVector.size()) cacheIndex = 0;
        fields.SetValueByIndex(index, value);
      }
    } else {
      fields.SetValueByKey(refName, value);
    }
  }

  inline su2double GetFieldValue(const std::string& refName) const{
    if (cacheEnabled){
      if (buildIndexCache){
        const int index = fields.GetIndex(refName);
        cacheIndexVector.push_back(index);
        return fields.GetItemByIndex(index).value;
      } else {
        const int index = cacheIndexVector[cacheIndex++];
        if (cacheIndex == cacheIndexVector.size()) cacheIndex = 0;
        return fields.GetItemByIndex(index).value;
      }
    }
    return fields.GetValueByKey(refName);
  }

  inline const COutFieldCollection& GetCollection() const { return fields; }

  void SetCaching(bool cacheEnable){
    cacheEnabled = cacheEnable;
    cacheIndex = 0;
    cacheIndexVector.clear();
  }

  void StartCaching(){
    buildIndexCache = cacheIndexVector.empty();
  }

};

class CHistoryOutFieldManager : public COutFieldManager {

public:

  CHistoryOutFieldManager () = default;

  inline FieldRef AddField(const std::string& refName, const std::string& headerName,
                           ScreenOutputFormat format, const std::string& groupName,
                           const std::string& description, FieldType type){
    return fields.AddItem(refName, COutputField(headerName, format, groupName, description, type));
  }
};

class CVolumeOutFieldManager : public COutFieldManager {

public:

  CVolumeOutFieldManager () = default;

  inline FieldRef AddField(const std::string& refName, const std::string& headerName,
                           const std::string& groupName, const std::string& description,
                           FieldType type){

    return fields.AddItem(refName, COutputField(headerName, groupName, description, type));
  }

  inline FieldRef AddField(const std::string&& refName, const std::string& headerName,
                           const std::string& groupName, const std::string& description,
                           FieldType type){
    /*--- Throw an error if the refName not smaller than 15 characters.
     * Otherwise, small string optimization will not be used, which can lead to some loss in performance,
     * when accessing field by a string literal.
     * The reasoning behind having it in the function that takes an rvalue reference is that fields
     *  defined using a rvalue reference are likely also to be accessed using an rvalue reference (i.e. string literal). ---*/
    if(refName.size() > 15){
      SU2_MPI::Error("Field names must be smaller than 15 characters: " + refName, CURRENT_FUNCTION);
    }
    return fields.AddItem(refName, COutputField(headerName, groupName, description, type));
  }
};
