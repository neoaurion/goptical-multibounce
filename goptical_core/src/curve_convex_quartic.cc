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
      Author: Alexandre Becoulet

*/


#include <Goptical/Curve/ConvexQuartic>
#include <Goptical/Math/Vector>
#include <Goptical/Math/VectorPair>
#include <Goptical/Math/VectorPair>
#include <limits>
#include <complex>

namespace _Goptical {

namespace Curve {

ConvexQuartic::ConvexQuartic(double a2, double a4)
{
    _a2=a2;
    _a4=a4;
    assert(a4!=0);
    assert(a2*a4>=0);
}

double ConvexQuartic::sagitta(double r) const
{
    return _a2*Math::square(r)+_a4*Math::square(Math::square(r));
}

double ConvexQuartic::derivative(double r) const
{
    return 2*_a2*r+4*_a4*r*r*r;
}

bool ConvexQuartic::intersect(Math::Vector3 &point, const Math::VectorPair3 &ray) const
{

    const double eps=1e-6;

    const double      px = ray.origin().x();
    const double      py = ray.origin().y();
    const double      pz = ray.origin().z();
    const double      dx = ray.direction().x();
    const double      dy = ray.direction().y();
    const double      dz = ray.direction().z();

    /*
        find intersection point between quartic and line, it is a specific case of 4th order polynomial root finding (see wikipedia)
      */

    // define 4th order polynomial problem
    std::complex<double> a(_a4*dx*dx*dx*dx + 2*_a4*dx*dx*dy*dy + _a4*dy*dy*dy*dy);
    std::complex<double> b(4*_a4*px*dx*dx*dx + 4*_a4*py*dx*dx*dy + 4*_a4*px*dx*dy*dy + 4*_a4*py*dy*dy*dy);
    std::complex<double> c(6*_a4*dx*dx*px*px + 2*_a4*dx*dx*py*py + _a2*dx*dx
                          + 8*_a4*dx*dy*px*py + 2*_a4*dy*dy*px*px + 6*_a4*dy*dy*py*py + _a2*dy*dy);
    std::complex<double> d(4*_a4*dx*px*px*px + 4*_a4*dy*px*px*py + 4*_a4*dx*px*py*py
                          + 2*_a2*dx*px + 4*_a4*dy*py*py*py + 2*_a2*dy*py - dz);
    std::complex<double> e(_a4*px*px*px*px + 2*_a4*px*px*py*py + _a2*px*px + _a4*py*py*py*py + _a2*py*py - pz);

    std::complex<double> delta0=c*c -3.0*b*d +12.0*a*e;
    std::complex<double> delta1=2.0*c*c*c -9.0*b*c*d +27.0*b*b*e +27.0*a*d*d -72.0*a*c*e;

    std::complex<double> p=(8.0*a*c-3.0*b*b)/(8.0*a*a);
    std::complex<double> q=(b*b*b-4.0*a*b*c+8.0*a*a*d)/(8.0*a*a*a);

    std::complex<double> Q=pow( (delta1+sqrt(delta1*delta1-4.0*delta0*delta0*delta0))/2.0 , 1.0/3.0);
    std::complex<double> S=0.5*sqrt(-2.0*p/3.0+(Q+delta0/Q)/(3.0*a));

    // 4 roots of the general problem
    std::vector< std::complex<double> > t_temp;
    t_temp.push_back( -b/(4.0*a)-S+0.5*sqrt(-4.0*S*S-2.0*p+q/S) );
    t_temp.push_back( -b/(4.0*a)-S-0.5*sqrt(-4.0*S*S-2.0*p+q/S) );
    t_temp.push_back( -b/(4.0*a)+S+0.5*sqrt(-4.0*S*S-2.0*p-q/S) );
    t_temp.push_back( -b/(4.0*a)+S-0.5*sqrt(-4.0*S*S-2.0*p-q/S) );

    // select only real positive roots
    std::vector<double> t;
    for (unsigned int i=0; i<t_temp.size(); i++)
    {
        if (fabs(std::imag(t_temp[i]))<eps){
            double t_temp_r=std::real(t_temp[i]);
            if(t_temp_r > eps)
            {
                t.push_back(t_temp_r);
            }
        }
    }

    if (t.size()==0)
        return false;

    double t_min = std::numeric_limits<double>::max();
    for (unsigned int i=0; i<t.size(); i++)
    {
        if(t[i]<t_min)
            t_min = t[i];
    }

    point = ray.origin() + ray.direction() * t_min;

    return true;
}

}

}

