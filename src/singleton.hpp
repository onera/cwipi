#ifndef __SINGLETON_H__
#define __SINGLETON_H__

#include <iostream>
#include <string>
#include <bft_printf.h>

namespace couplings {
  template <typename T>
  class Singleton
  {
  protected:
    Singleton () { }
    virtual ~Singleton () { bft_printf( "destroying singleton.\n" ); }
    Singleton (const Singleton & other) { }
    Singleton & operator=(const Singleton & other) { }

  public:
    static T &getInstance ()
    {
      if (NULL == _singleton)
        _singleton = new T;

      return *(static_cast<T*> (_singleton));
    }

    static void kill ()
    {
      if (NULL != _singleton) {
        delete _singleton;
        _singleton = NULL;
      }
    }

  private:
    static T *_singleton;
  };

  template <typename T>
  T *Singleton<T>::_singleton = NULL;
} // namespace couplings

#endif
