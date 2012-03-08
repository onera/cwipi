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
bftinclude_HEADERS =
#bftinclude_HEADERS = \
#bftc_backtrace.h \
#bftc_error.h \
#bftc_file.h \
#bftc_fp_trap.h \
#bftc_mem.h \
#bftc_mem_usage.h \
#bftc_printf.h \
#bftc_sys_info.h \
#bftc_timer.h \
#bftc_version.h

noinst_HEADERS = \
bftc_config_defs.h \
bftc_backtrace.h \
bftc_error.h \
bftc_file.h \
bftc_fp_trap.h \
bftc_mem.h \
bftc_mem_usage.h \
bftc_printf.h \
bftc_sys_info.h \
bftc_timer.h \
bftc_version.h

# Library source files

lib_LTLIBRARIES = libbftc.la
libbftc_la_SOURCES   = \
bftc_backtrace.c \
bftc_error.c \
bftc_file.c \
bftc_fp_trap.c \
bftc_mem.c \
bftc_mem_usage.c \
bftc_printf.c \
bftc_sys_info.c \
bftc_timer.c \
bftc_version.c

libbftc_la_LDFLAGS = -no-undefined -version-info 1:0:0

