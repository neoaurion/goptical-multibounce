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
  
      Copyright (C) 2010-2011 Free Software Foundation, Inc
      Author: Alexandre Becoulet

*/


#ifndef GOPTICAL_TRACE_RESULT_HH_
#define GOPTICAL_TRACE_RESULT_HH_

#include <set>
#include <deque>

#include "Goptical/common.hh"

#include "Goptical/Sys/element.hh"
#include "Goptical/Sys/surface.hh"
#include "Goptical/Trace/ray.hh"

namespace _Goptical {

  namespace Trace {

    /**
       @short Store light propagation result
       @header Goptical/Trace/Result
       @module {Core}
       @main

       This class encapsulates rays data propagation result.

       It must be properly configured before light propagation as
       needed by the analysis currently being performed. All
       requested light propagation informations will be store for
       further processing.

       All @ref Ray object are allocated by this class. It is able
       to remember which element intercepted and generated each ray.
    */
    class Result
    {
      friend class Tracer;

    public:
      typedef std::vector<const Sys::Source *> sources_t;

      /** Crate a new empty result object */
      Result();

      ~Result();

      /** Get the list of rays striking a given surface */
      inline const rays_queue_t & get_intercepted(const Sys::Surface &s) const;
      /** Get the list of rays generated by a given element */
      inline const rays_queue_t & get_generated(const Sys::Element &s) const;

      /** Get list of sources used for ray tracing */
      inline const Trace::Result::sources_t & get_source_list() const;

      /** Get window which include all ray intercepted on a surface */
      Math::VectorPair3 get_intercepted_window(const Sys::Surface &s) const;
      /** Get center of window */
      Math::Vector3 get_intercepted_center(const Sys::Surface &s) const;
      /** Get centroid of all ray intercepted on a surface */
      Math::Vector3 get_intercepted_centroid(const Sys::Surface &s) const;

      /** Clear all result data */
      void clear();

      /** List of rays striking this surface must be saved when tracing rays */
      void set_intercepted_save_state(const Sys::Element &e, bool enabled = true);
      /** Return true if generated rays must be saved for this element */
      bool get_intercepted_save_state(const Sys::Element &e);

      /** List of rays generated by this element must be saved when tracing rays */
      void set_generated_save_state(const Sys::Element &e, bool enabled = true);
      /** Return true if generated rays must be saved for this element */
      bool get_generated_save_state(const Sys::Element &e);

      /** Set all save states to false */
      void clear_save_states();

      /** Get maximum intensity for a single ray FIXME */
      double get_max_ray_intensity() const;

      /* Get raytracing mode used FIXME */
      //  inline Tracer::Mode get_mode() const;

      /** Allocate a new Trace::Ray object from result */
      inline Ray & new_ray();

      /** Allocate a new Trace::Ray object from result */
      inline Ray & new_ray(const Light::Ray &r);

      /** Declare a new ray interception */
      inline void add_intercepted(const Sys::Surface &s, Ray &ray);
      /** Declare a new ray generation */
      inline void add_generated(const Sys::Element &s, Ray &ray);

      /** Declare ray wavelen used for tracing */
      inline void add_ray_wavelen(double wavelen);

      /** Get ray wavelen in use set */
      inline const std::set<double> & get_ray_wavelen_set() const;

      /** Get reference to tracer parameters used */
      inline const Params & get_params() const;

      /** Draw all tangential rays using specified renderer. Only rays
          which end up hitting the image plane are drawn when @tt
          hit_image is set. */
      void draw_2d(Io::Renderer &r, bool hit_image = false,
                   const Sys::Element *ref = 0) const;
      /** Draw all rays using specified renderer. Only rays
          which end up hitting the image plane are drawn when @tt
          hit_image is set.*/
      void draw_3d(Io::Renderer &r, bool hit_image = false,
                   const Sys::Element *ref = 0) const;

    private:
      void init(const Sys::System &system);
      void init(const Sys::Element &element);

      void prepare();

      struct element_result_s
      {
        rays_queue_t *_intercepted; // list of rays for each intercepted surfaces
        rays_queue_t *_generated; // list of rays for each generator surfaces
        bool _save_intercepted_list;
        bool _save_generated_list;
      };

      inline struct element_result_s & get_element_result(const Sys::Element &e);
      inline const struct element_result_s & get_element_result(const Sys::Element &e) const;

      vector_pool<Ray, 256> _rays; // rays allocation pool
      std::vector<struct element_result_s> _elements;
      std::set<double>          _wavelengths;
      rays_queue_t              *_generated_queue;
      Trace::Result::sources_t  _sources;
      unsigned int              _bounce_limit_count;
      const Sys::System         *_system;
      const Trace::Params       *_params;
      //  Tracer::Mode          _mode;
    };
  }
}

#endif

