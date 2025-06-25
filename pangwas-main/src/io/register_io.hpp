#ifndef PANGWAS_IO_REGISTER_IO_HPP_INCLUDED
#define PANGWAS_IO_REGISTER_IO_HPP_INCLUDED

/**
 * \file register_io.hpp
 * Includes a function to call to register IO handlers for vg types.
 * Copied from vg
 */

namespace pangwas {

namespace io {

using namespace std;

/**
 * Register libvg types with libvgio.
 * Must be called by library users before doing IO.
 * Does not magically statically call itself.
 * Returns true on success.
 */
bool register_libvg_io();

}

}

#endif
