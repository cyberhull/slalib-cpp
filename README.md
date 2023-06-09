
C++ Port of the SLALIB Library
==============================

This project is a C++ port of the FORTRAN library `SLALIB`. This is what the original
library's documentation says about its purpose:

> SLALIB is a library used by writers of positional-astronomy applications.
> Most of the 188 routines are concerned with astronomical position and time, but
> a number have wider trigonometrical, numerical or general applications.

The original (FORTRAN) library comes with extensive documentation, available
from many sites; for instance,
[here](http://star-www.rl.ac.uk/docs/sun67.htx/sun67.html). Not only has it
quite extensive coverage of every function, but is also contains large section
explaining [math behind library
routines](http://star-www.rl.ac.uk/docs/sun67.htx/sun67se4.html#x193-5840004).
That said, in our opinion, the very best source of information on positional
astronomy and related topics is the *Fundamental Astronomy* book by *H.
Karttunen et. al.* [available from
Springer](https://www.springer.com/gp/book/9783662530443). Highly recommended.

The FORTRAN source code used for the port is `SLABIB` version `2.5-7` bundled
with the `starlink` software by P.T. Wallace that can be [found on
GitHub](https://github.com/Starlink/starlink/tree/master/libraries/sla).

Hope you will find this C++ library useful.

The CyberHULL Team.
