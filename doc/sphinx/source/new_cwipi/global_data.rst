.. _global data:

Exchange of Global Data
~~~~~~~~~~~~~~~~~~~~~~~

A *GlobalData* object is to be used to send/recv (non-blocking) interlaced strided arrays of any
type through a Coupling object. Arrays to exchange are characterized the number of elements,
the stride, *i.e.* the number of values by element, and the size of one value in bytes.

Each Coupling object can handle an arbitrary number of GlobalData objects.
