#ifndef __SINGLETON_H__
#define __SINGLETON_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-2026  ONERA

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

/**
 * \cond
 */

namespace cwipi {

  template <typename T>

  /**
   * \class Singleton singleton.hpp "singleton.hpp"
   * \brief Singleton pattern implementation.
   *
   * Template singleton class ensuring a single instance per type
   *
   */
  class Singleton {

  public:

    /**
     * \brief Singleton instance.
     */
    Singleton() {}


    /**
     * \brief Get the singleton instance.
     *
     * Creates the instance on first call and returns a reference to it
     *
     * \return Reference to the singleton instance
     */
    static T &getInstance()
    {
      if (_singleton == nullptr)
        _singleton = new T;

      return *(static_cast<T*>(_singleton));
    }


    /**
     * \brief Destroy the singleton instance.
     *
     * Deletes the singleton instance if it exists
     */
    static void kill()
    {
      if (_singleton != nullptr) {
        delete _singleton;
        _singleton = nullptr;
      }
    }

  protected:

    /**
     * \brief Destructor.
     */
    virtual ~Singleton() {}

  private:

    static T* _singleton; /*!< Pointer to the singleton instance */

  };

  template <typename T>
  T* Singleton<T>::_singleton = nullptr;

}

/**
 * \endcond
 */

#endif // __SINGLETON_H__