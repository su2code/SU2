#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include "../mpi_structure.hpp"

/*! The CIndexedMap combines a associative map and
 * and a sequentially indexed container. This allows access to items using a
 * key (e.g. a string) and at the same time the possibiilty to sequentially loop
 * over all key/item pairs in a fixed order, based on their insertion.
 *
 */
template <class Key, class Item>
class CIndexedMap {

public:
  using KeyType          = Key;
  using ItemType         = Item;
  using Map              = std::unordered_map<Key, Item>;
  using KeyVector        = std::vector<Key>;
  using InsertionVector  = std::vector<typename Map::iterator>;

protected:

  Map map;
  InsertionVector insertionVector;

public:

  /*!
   * \brief Get the size of the container
   * \return The current number of keys/item pairs.
   */
  int GetSize() const { return insertionVector.size(); }

  /*!
   * \brief Add a new item to the container using a key
   * \param[in] key - The key associated to the item
   * \param[in] item - The item that should be stored
   */
  void AddItem(const Key& key, const Item& item){
    bool success;
    typename Map::iterator iter;
    std::tie(iter, success) = map.insert({key, item});
    insertionVector.emplace_back(map.insert({key, item}).first);
    if (!success){
      SU2_MPI::Error(std::string("Key with name ") + key + std::string(" already exists"), CURRENT_FUNCTION);
    }
  }

  /*!
   * \brief Check if a specific key exists in the container
   * \param[in] key - The key that should be checked
   * \return <TRUE> then when the key could be found, otherwise <FALSE>
   */
  bool CheckKey(const Key& key) const {
    return (map.count(key) > 0);
  }

  /*!
   * \brief Get the key at a specific position
   * \param[in] i - The position of the key
   * \return The key at the requested position
   */
  const Key& GetKey(const int i) const {
    return insertionVector[i]->first;
  }

  /*!
   * \brief Get the item using the associated key
   * \param[in] key - The key value
   * \return The item associated with the specified key
   */
  const Item& GetItemByKey(const Key& key) const {
    return map.at(key);
  }

  /*!
   * \brief Get the item at a specific position
   * \param[in] i - The position of the item
   * \return The item at the specified position
   */
  const Item& GetItemByIndex(const int i) const {
    return insertionVector[i]->second;
  }

  /*!
   * \brief Get the item using the associated key using the bracket operator
   * \param[in] key - The key value
   * \return The item associated with the specified key
   */
  const Item& operator[](const Key& key) const {
    return GetItemByKey(key);
  }

  /*!
   * \brief Get references to all key/item pairs
   * \return vector containing iterator references to all key/item pairs
   */
  const InsertionVector& GetReferencesAll() const {
    return insertionVector;
  }

  /*!
   * \brief Get the index of a certain key
   * \param[in] key - The key
   * \return The index of the item with a certain key
   */
  int GetIndex(const Key& key){
    int index;
    for (index = 0; index < insertionVector.size(); index++){
      if (insertionVector[index]->first == key){
        return index;
      }
    }
    return -1;
  }

  /*!
   * \brief Get references to specific key/item pairs that have a certain property
   * \param[in] propertyList - Vector containing a number of properties
   * \param[out] notFound - Vector containing properties that could not be found
   * \param[in] f - Function that takes a key and a iterator position to check for a certain condition
   * \return vector containing iterator references to matched key/item pairs
   */
  template<typename P, typename F>
  static InsertionVector GetReferences(const std::vector<P> &propertyList,
                                       std::vector<P>& notFound,
                                       const InsertionVector& refVector,
                                       const F &f) {
    InsertionVector references;
    notFound.clear();
    std::vector<bool> found(propertyList.size(), false);
    for (const auto& mapIndex : refVector){
      int i = 0;
      for (const auto& prop : propertyList){
        if (f(prop, mapIndex)){
          references.push_back(mapIndex);
          found[i] = true;
        }
        i++;
      }
    }
    for (unsigned int i = 0; i < found.size(); i++){
      if (!found[i]){
        notFound.push_back(propertyList[i]);
      }
    }
    return references;
  }

  /*!
   * \brief Extract the keys from an insertion vector
   * \param[in] inVector - The insertion vector
   * \return Vector containing the keys
   */
  static KeyVector GetKeys(const InsertionVector& inVector){
    KeyVector keyList;
    for (auto mapIndex : inVector){
      keyList.push_back(mapIndex->first);
    }
    return keyList;
  }

  /*!
   * \brief Combine two insertionVectors into one and remove duplicates.
   * \param[in] vecOne - The first vector for the combination
   * \param[in] vecTwo - The second vector for the combination
   * \return The combined vector
   */
  static InsertionVector Combine(const InsertionVector& vecOne, const InsertionVector& vecTwo){
    InsertionVector newVector;
    for (const auto& field : vecOne){
      newVector.push_back(field);
    }
    for (const auto& field: vecTwo){
      if (std::find(newVector.begin(), newVector.end(), field) == newVector.end()){
        newVector.push_back(field);
      }
    }
    return newVector;
  }
};