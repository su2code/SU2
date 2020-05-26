#pragma once

#include "../../../Common/include/toolboxes/CIndexedMap.hpp"
#include "../../../Common/include/datatype_structure.hpp"
#include "../../../Common/include/mpi_structure.hpp"
#include "COutField.hpp"

class COutFieldCollection : public CIndexedMap<std::string, COutputField> {

public:
  using CurrentType     = CIndexedMap<std::string, COutputField>;
  using InsertionVector = typename CurrentType::InsertionVector;
  using KeyVector       = typename CurrentType::KeyVector;
  using KeyType         = typename CurrentType::KeyType;
  using ItemType        = typename CurrentType::ItemType;
  using Map             = typename CurrentType::Map;

protected:
  using CurrentType::insertionVector;
  using CurrentType::map;

  /*! \brief Vector to cache the positions of the field in the data array */
  std::vector<short>                            cacheIndexVector;
  /*! \brief Current value of the cache index */
  unsigned short                                cacheIndex;
  /*! \brief Boolean to store whether the field index cache should be build. */
  bool                                          buildIndexCache;

  bool cacheEnabled = false;

private:

  /*!
   * \brief Check if an iterator position matches a field type
   * \param type - The field type to check
   * \param mapIndex - The iterator position
   * \return <TRUE> if iterator position is of specified type
   */
  static bool findFieldWithType(const FieldType& type, const typename Map::iterator& mapIndex){
    return (mapIndex->second.fieldType == type);
  }

  /*!
   * \brief Check if an iterator position matches a group
   * \param group - The name of the group
   * \param mapIndex - The iterator position
   * \return <TRUE> if iterator position is in specified group
   */
  static bool findFieldWithGroup(const std::string& group, const typename Map::iterator& mapIndex){
    return (mapIndex->second.outputGroup == group);
  }

  /*!
   * \brief Check if an iterator position matches a certain key
   * \param name - The name of the key
   * \param mapIndex - The iterator position
   * \return <TRUE> if iterator position has specified name
   */
  static bool findFieldWithName(const std::string& key, const typename Map::iterator& mapIndex){
    return (mapIndex->first == key);
  }

public:

  /*!
   * \brief Get fields by using a list of keys
   * \param[in] fieldKey   - A vector of keys
   * \param[out] notFound  - Vector containing a list of keys that have not been found
   * \param[in] refVector  - Vector containing references to entries of an indexed map
   * \return Vector containing references to the entries with the specified keys and/or group names
   */
  static InsertionVector GetFieldsByKey(const KeyVector& fieldKey, KeyVector& notFound,
                                        const InsertionVector& refVector) {
    return GetReferences(fieldKey, notFound, refVector, findFieldWithName);
  }

  /*!
   * \brief Get fields by using a list of keys
   * \param[in] fieldKey   - A vector of keys
   * \param[in] refVector  - Vector containing references to entries of an indexed map
   * \return Vector containing references to the entries with the specified keys and/or group names
   */
  static InsertionVector GetFieldsByKey(const KeyVector& fieldKey, const InsertionVector& refVector) {
    KeyVector dummy;
    return GetFieldsByKey(fieldKey, dummy, refVector);
  }

  /*!
   * \brief Get fields by using a list of keys
   * \param[in] fieldKey   - A vector of keys
   * \param[out] notFound  - Vector containing a list of keys that have not been found
   * \return Vector containing references to the entries with the specified keys and/or group names
   */
  InsertionVector GetFieldsByKey(const KeyVector& fieldKey, KeyVector& notFound) const {
    return CurrentType::GetReferences(fieldKey, notFound, insertionVector, findFieldWithName);
  }

  /*!
   * \brief Get fields by using a list of keys
   * \param[in] fieldKey   - A vector of keys
   * \param[in] searchGroup - Boolean indicating whether also group names should be also included in the search
   * \return Vector containing references to the entries with the specified keys and/or group names
   */
  InsertionVector GetFieldsByKey(const KeyVector& fieldKey) const {
    KeyVector dummy;
    return GetFieldsByKey(fieldKey, dummy);
  }

  /*!
   * \brief Get fields by using a list of group names
   * \param[in] groupList  - A vector of group names
   * \param[out] notFound  - Vector containing a list of group names that have not been found
   * \param[in] refVector  - Vector containing references to entries of an indexed map
   * \return Vector containing references to the entries with the specified group names
   */
  static InsertionVector GetFieldsByGroup(const std::vector<std::string>& groupList,
                                          std::vector<std::string>& notFound, const InsertionVector& refVector) {
    return CurrentType::GetReferences(groupList, notFound, refVector, findFieldWithGroup);
  }

  /*!
   * \brief Get fields by using a list of group names
   * \param[in] groupList  - A vector of group names
   * \param[out] notFound  - Vector containing a list of group names that have not been found
   * \return Vector containing references to the entries with the specified group names
   */
  InsertionVector GetFieldsByGroup(const std::vector<std::string>& groupList, std::vector<std::string>& notFound) const {
    return CurrentType::GetReferences(groupList, notFound, insertionVector, findFieldWithGroup);
  }

  /*!
   * \brief Get fields by using a list of group names
   * \param[in] groupList  - A vector of group names
   * \return Vector containing references to the entries with the specified group names
   */
  InsertionVector GetFieldsByGroup(const std::vector<std::string>& groupList) const {
    std::vector<std::string> dummy;
    return GetFieldsByGroup(groupList, dummy);
  }

  /*!
   * \brief Get fields by using a list of field types
   * \param[in] groupList  - A vector of field types
   * \param[in] refVector  - Vector containing references to entries of an indexed map
   * \return Vector containing references to the entries with the specified field types
   */
  static InsertionVector GetFieldsByType(const std::vector<FieldType>& type, const InsertionVector& refVector) {
    std::vector<FieldType> dummy;
    return CurrentType::GetReferences(type, dummy, refVector, findFieldWithType);
  }

  /*!
   * \brief Get fields by using a list of field types
   * \param[in] groupList  - A vector of field types
   * \return Vector containing references to the entries with the specified field types
   */
  InsertionVector GetFieldsByType(const std::vector<FieldType>& type) const {
    std::vector<FieldType> dummy;
    return CurrentType::GetReferences(type, dummy, insertionVector, findFieldWithType);
  }

  void SetCaching(bool cacheEnable){
    cacheEnabled = cacheEnable;
    cacheIndex = 0;
    cacheIndexVector.clear();
  }

  void StartCaching(){
    buildIndexCache = cacheIndexVector.empty();
  }

  ItemType& GetItemByKey(const KeyType& key){
    if (cacheEnabled){
      if (buildIndexCache){
        const int index = CurrentType::GetIndex(key);
        cacheIndexVector.push_back(index);
        return CurrentType::GetItemByIndex(index);
      } else {
        const int index = cacheIndexVector[cacheIndex++];
        if (cacheIndex == cacheIndexVector.size()) cacheIndex = 0;
        return CurrentType::GetItemByIndex(index);
      }
    } else {
      return CurrentType::GetItemByKey(key);
    }
  }

  su2double GetValueByKey(const KeyType& key){
    return GetItemByKey(key).value;
  }

  /*!
   * \brief Set the value of specific field by using its key
   * \param[in] key   - The key of the field
   * \param[in] value - The new value for this field
   */
  void SetValueByKey(const KeyType& key, su2double value){
    GetItemByKey(key).value = value;
  }

  /*!
   * \brief Set the value of specific field by using its index
   * \param[in] i     - The index of the field
   * \param[in] value - The new value for this field
   */
  void SetValueByIndex(const int i, su2double value){
    insertionVector[i]->second.value = value;
  }
};

/*!
 * \brief Typedef for a collection of history fields
 */
//typedef COutFieldCollection<COutputField> COutFieldCollection;

/*!
 * \brief Typedef for a collection of volume fields
 */
//typedef COutFieldCollection<COutputField> COutFieldCollection;


