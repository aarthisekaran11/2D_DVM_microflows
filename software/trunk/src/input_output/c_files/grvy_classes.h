// -*-c++-*-
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// libGRVY - a utility library for scientific computing.
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, 
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// grvy_classes.h: Basic class definitions
//
// $Id: grvy_classes.h 13660 2010-09-19 16:05:05Z karl $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRVY_CLASSES_H_
#define GRVY_CLASSES_H_

#include<grvy.h>

#include<cstdarg>
#include<limits>
#include<map>
#include<vector>
#include<string>
#include<stack>
#include<config_grvy.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

namespace GRVY {

  //--------------------------
  // Logging Class
  //--------------------------

  class GRVY_Log_Class {
  private:
    int log_level;			               // Current log level priority
    std::map<int,std::string> LogMask;                      // String masks for log messages
  public:
    GRVY_Log_Class();		

    bool isLog(int priority) {                                      // inlined log priority test
      return(priority <= log_level); }                              
    void msg       (int priority, std::string msg);                      // post new log message with a priority
    int  msg_printf(int priority, const char *format,va_list argp); // post printf style log message
    void change_priority(int priority);                             // change current log level priority

  };

//--------------------------
// HDF5 Utility Class
//--------------------------

  class GRVY_HDF5_Class {
  private: 

#ifdef HAVE_HDF5
    hid_t m_fileId;		        // hdf5 file handle
    H5E_auto2_t error_orig_func;	// error-handle func
    void       *error_orig_data;        // error-handle stack data
#endif

  public:
    GRVY_HDF5_Class      ();
#ifdef HAVE_HDF5
    int  file_create(const char *filename, bool overwrite_existing);
    int  close();
    void silence_hdf_error_handler();
    void restore_hdf_error_handler();
#endif
  };

}

#endif
