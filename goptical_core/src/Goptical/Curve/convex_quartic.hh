/*

      This file is part of the Goptical Core library.
  
      The Goptical library is free software; you can redistribute it
      and/or modify it under the terms of the GNU General Public
      License as published by the Free Software Foundation; either
      version 3 of the License, or (at your option) any later version.
  
      The Goptical library is distributed in the hope that it will be
      useful, but WITHOUT ANY WARRANTY; without even the implied
      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
      See the GNU General Public License for more details.
  
      You should have received a copy of the GNU General Public
      License along with the Goptical library; if not, write to the
      Free Software Foundation, Inc., 59 Temple Place, Suite 330,
      Boston, MA 02111-1307 USA
  
      Copyright (C) 2010-2018 Free Software Foundation, Inc
      Author: Arnaubec Aurelien

*/


#ifndef GOPTICAL_CURVE_CONVEX_QUARTIC_HH_
#define GOPTICAL_CURVE_CONVEX_QUARTIC_HH_

#include "Goptical/common.hh"

#include "Goptical/Curve/conic_base.hh"

namespace _Goptical {

  namespace Curve {

    /**
       @short Convex quartic curve model (a4*r^4 + a2*r^2 with a4 & a2 of the same sign for convexity)
       @header Goptical/Curve/ConvexQuartic
       @module {Core}
       @main

       This class provides an efficient convex quartic curve implementation.
     */
    class ConvexQuartic : public Rotational
    {
    public:
      /** Creates a ConvexQuartic curve with given a2 and a4 */
      ConvexQuartic(double a2, double a4);

      bool intersect(Math::Vector3 &point, const Math::VectorPair3 &ray) const;

      double sagitta(double r) const;
      double derivative(double r) const;
    private:
      double _a2,_a4;
    };

  }
}

#endif

