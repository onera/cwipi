#ifndef __FACTORY_H__
#define __FACTORY_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2026  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "singleton.hpp"

#include <exception>
#include <map>

/**
 * \cond
 */

namespace cwipi {

  /**
   * \defgroup FactoryErrorPoliciesGroup Factory Error Policies
   * \ingroup  FactoryGroup
   * \brief    Manages the "Unknown Type" error in an object factory.
   *
   * \class    DefaultFactoryError factory.hpp "factory.hpp"
   * \ingroup  FactoryErrorPoliciesGroup
   * \brief    Default policy that throws an exception.
   *
   */

  template <typename IdentifierType, class AbstractProduct>
  struct DefaultFactoryError {

    /**
     * \struct  Exception
     * \ingroup FactoryErrorPoliciesGroup
     * \brief   Throws an exception.
     *
     */
    struct Exception : public std::exception
    {

      /**
       * \brief what function.
       *
       * \return "Unknown Type"
       */
      const char* what() const throw() { return "Unknown Type"; }

    };


    /**
     * \brief Throws an exception.
     *
     * \return AbstractProduct*
     */
    static AbstractProduct* ofUnknownType(IdentifierType)
    {
      throw Exception();
    }

  };


  /**
   * \class   Factory factory.hpp "factory.hpp"
   * \ingroup FactoryGroup
   * \brief   Generic object factory.
   *
   * Implements a generic object factory
   *
   */

  template <
    class                           AbstractProduct,
    typename                        IdentifierType,
    typename                        ProductCreator = AbstractProduct* (*)(),
    template<typename, class> class FactoryErrorPolicy = DefaultFactoryError
    >
  class Factory :
    public Singleton < Factory <
                                AbstractProduct,
                                IdentifierType,
                                ProductCreator,
                                FactoryErrorPolicy > >,
    public FactoryErrorPolicy < IdentifierType, AbstractProduct >
  {

  public:

    /**
     * \brief Add a new couple identifier / concrete product.
     *
     * \tparam     ConcreteProduct Concrete product
     * \param [in] id              Identifier
     *
     */
    template < class ConcreteProduct >
    bool add(const IdentifierType& id)
    {
      return _associations.insert(
        typename IdToProductMap::value_type(id,
                        newCreator<ConcreteProduct>)).second != 0;
    }


    /**
     * \brief Remove a concrete product from its identifier.
     *
     * \param [in]  id Identifier
     *
     */
    bool remove(const IdentifierType& id)
    {
      return _associations.erase(id) != 0;
    }


    /**
     * \brief Create an object of Concrete product from its identifier.
     *
     * \param [in]  id Identifier
     *
     */
    AbstractProduct* objectCreate(const IdentifierType& id)
    {
      typename IdToProductMap::iterator i = _associations.find(id);
      if (i != _associations.end())
        {
          return (i->second)();
        }
      return this->ofUnknownType(id);
    }

  private:

    typedef std::map<IdentifierType, ProductCreator> IdToProductMap;   /*!< Container type (map) */
    IdToProductMap _associations;                                     /*!< Container */

    /**
     * \brief Create an object of Concrete product from its identifier.
     *
     * \tparam ConcreteProduct Concrete product
     * \return A new concrete product with the abstract product interface
     *
     */
    template <class ConcreteProduct>
    static AbstractProduct * newCreator() {
      return new ConcreteProduct;
    }

  };

}

/**
 * \endcond
 */

#endif // __FACTORY_H__