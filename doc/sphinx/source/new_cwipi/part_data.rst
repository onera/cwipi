.. _part data:

Exchange of Partitioned Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A *PartData* object is to be used to send/recv (non-blocking) interlaced strided arrays of any
type through a Coupling object. Contrary to the GlobalData object, this object accounts for the
global numbering of entities to exchange data between entities having the same global number. Note that
when an entity is defined on multiple ranks (eg vertices) only the value on the lowest MPI rank is kept.
Multiple exchanges, identified by an integer ``exch_id``, can be performed simultaneously. Arrays to exchange
are characterized by the stride, *i.e.* the number of values by element, and the size of one value in bytes.

Each Coupling object can handle an arbitrary number of PartData objects.
