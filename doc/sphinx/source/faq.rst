.. _faq:

CWIPI's FAQ
===========

This is a list of Frequently Asked Questions about CWIPI.  Feel free to
suggest new entries!

Why...
------

... does CWIPI fail with the given data entry?
   Before accusing CWIPI have you verifyed if your data entry complies with the CWIPI
   data entry written in the documentation?

... does CWIPI fail at certain runs?
   This is due to the asynchronous exchanges. It is not possible to get a data dat has
   no been set by the other code. A possible fix is to do the get in a while loop to wait
   until the data is available.



