#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>

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
    insertionVector.emplace_back(map.insert({key, item}).first);
  }

  /*!
   * \brief Check if a specific key exists in the container
   * \param[in] key - The key that should be checked
   * \return <TRUE> then when the key could be found, otherwise <FALSE>
   */
  bool FindKey(const Key& key) const {
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
  const Item& GetItem(const Key& key) const {
    return map.at(key);
  }

  /*!
   * \brief Get the item at a specific position
   * \param[in] i - The position of the item
   * \return The item at the specified position
   */
  const Item& GetItem(const int i) const {
    return insertionVector[i]->second;
  }

  /*!
   * \brief Get the item using the associated key using the bracket operator
   * \param[in] key - The key value
   * \return The item associated with the specified key
   */
  const Item& operator[](const Key& key) const {
    return GetItem(key);
  }

  /*!
   * \brief Get the key at a specific position using the bracket operator
   * \param[in] i - The position of the key
   * \return The key at the requested position
   */
  const Item& operator[](const int i) const {
    return GetItem(i);
  }

  /*!
   * \brief Get references to all key/item pairs
   * \return vector containing iterator references to all key/item pairs
   */
  InsertionVector GetReferencesAll() const {
    return insertionVector;
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
    bool found = false;
    for (auto prop : propertyList){
      found = false;
      for (auto mapIndex : refVector){
        if (f(prop, mapIndex)){
          references.push_back(mapIndex);
          found = true;
        }
      }
      if (!found) notFound.push_back(prop);
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
};