#-------------------------------------------------------------------------------
#   This file is part of the "Base Functions and Types" library, intended to
#   simplify and enhance portability, memory and I/O use for scientific codes.
#
#   Copyright (C) 2011       ONERA
#   Copyright (C) 2004-2008  EDF
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Recommended for VPATH support in autoconf documentation ; does not
# seem necessary with current Makefile functionality using
# autoconf 2.59 and automake 1.8.
# srcdir = @srcdir@
# VPATH = @srcdir@

# To disable dependency tracking
# AUTOMAKE_OPTIONS = no-dependencies

#include directories

#INCLUDES = -I$(top_srcdir)/include

# Add any debug flags

CC = $(CXX)

if DEBUG
  DEBUG_CPPFLAGS = -DDEBUG
else
  DEBUG_CPPFLAGS = -DNDEBUG
endif

# Main part

AM_CPPFLAGS = $(DEBUG_CPPFLAGS)
AM_CFLAGS   =
AM_LDFLAGS  =

# Public header files (to be installed)

bftincludedir = $(includedir)
bftinclude_HEADERS = \
bft_backtrace.h \
bft_error.h \
bft_file.h \
bft_fp_trap.h \
bft_mem.h \
bft_mem_usage.h \
bft_printf.h \
bft_sys_info.h \
bft_timer.h \
bft_version.h

noinst_HEADERS = \
bft_config_defs.h

# Library source files

lib_LTLIBRARIES = libbfto.la
libbfto_la_SOURCES   = \
bft_backtrace.cpp \
bft_error.cpp \
bft_file.cpp \
bft_fp_trap.cpp \
bft_mem.cpp \
bft_mem_usage.cpp \
bft_printf.cpp \
bft_sys_info.cpp \
bft_timer.cpp \
bft_version.cpp

libbfto_la_LDFLAGS = -no-undefined -version-info 1:0:0

