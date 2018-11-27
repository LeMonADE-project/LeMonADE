#ifndef LATTICE_PREDICATES_H
#define LATTICE_PREDICATES_H
/******************************************************************************/
/**
 *@author Toni
 *@date 2018/02/26
 * 
 * @brief provides functors what is written on the lattice 
/******************************************************************************/

/**
 * @struct Bool
 *
 * @brief Functor returning always \a True.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 **/

struct Bool
{
  template <class MoleculesType>
  bool operator() (const MoleculesType& m, int i){return true;}
};

/**
 * @struct Attribute
 *
 * @brief Functor returning  attribute tag.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 **/

struct Attribute
{
  template <class MoleculesType>
  uint32_t operator() (const MoleculesType& m, int i){ return m[i].getAttributeTag();}
};

/**
 * @struct MonomerID
 *
 * @brief Functor returning monomer id .
 *
 * @deprecated
 *
 *
 **/

struct MonomerID
{
  template <class MoleculesType>
  uint32_t operator() (const MoleculesType& m, int i){ return i+1;}
};

/**
 * @struct AttributeANDMonomerID
 *
 * @brief Functor returning attribute and monomer id .
 *
 * @details This functor writes on the first two bits the attribute and on the remaining bits the monomer id.
 * The id's start at 1. To get the monomer id use: LatticeEntry>>6 and to get the attribute use LatticeEntry&63.
 * Up to know there can be only attributes stored from 0 to 63 and monomer id's up to 67.108.864
 * @deprecated
 *
 *
 **/

struct AttributeANDMonomerID
{
  template <class MoleculesType>
  uint32_t operator() (const MoleculesType& m, int i){ 
    uint32_t LatticeEntry(i+1);
    uint32_t Attribute(m[i].getAttributeTag());
    LatticeEntry=(LatticeEntry <<6) + Attribute;
    return LatticeEntry;
  }
};
#endif
