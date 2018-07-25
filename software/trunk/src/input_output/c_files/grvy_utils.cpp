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
// grvy_utils: Internal Utility Functions
//
// $Id: grvy_utils.cpp 13640 2010-09-18 16:49:30Z karl $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// Notes: these are primarily for internal use by GRVY functions.  Hint: these
// may not be the external droids you are looking for.


#include<errno.h>
#include<ftw.h>
#include<libgen.h>
#include<stdio.h>
#include<string.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<unistd.h>
#include<stack>

#include<grvy_classes.h>
#include<grvy_int.h>
#include<grvy.h>

#include"fortran_string_order.h"

using namespace std;

namespace GRVY {

  //-------------------------
  // Miscellaneous utilities
  //-------------------------

  extern "C" int grvy_check_file_path(const char *pathname)
  {
    const int MAX_DEPTH = 50;

    char *pathlocal;
    char *parents;
    char *dirstring;
    char *token;
    int depth = 0;

    // Save a copy of pathname and look for the parent directories.

    pathlocal = strdup(pathname);
    dirstring = strdup(pathname);
    parents   = dirname(pathlocal);

    if(strcmp(parents,".") == 0)
      {
	free(pathlocal);
	free(dirstring);
	return 0;
      }

    if( (token = strtok(parents,"/")) != NULL )
      {
	if ( _GRVY_CheckDir(token) )
	  {
	    free(pathlocal);
	    free(dirstring);
	    return -1;
	  }

	// Now, search for any remaining parent directories.

	sprintf(dirstring,"%s",token);

	while ( (token = strtok(0,"/")) && (depth < MAX_DEPTH) )
	  {
	    dirstring = strcat(dirstring,"/");
	    if(_GRVY_CheckDir(strcat(dirstring,token)))
	      {
		free(pathlocal);
		free(dirstring);
		return -1;
	      }
	    depth++;
	  };

	if(depth >= MAX_DEPTH )
	  { 
	    free(pathlocal);
	    free(dirstring);
	    return -1;
	  }
      }

    // Clean Up
    free(pathlocal);
    free(dirstring);

    return 0;
  }

  int _GRVY_CheckDir(const char *dirname)
  {
    struct stat st;

    if(stat(dirname,&st) != 0)
      {
	if( mkdir(dirname,0700) != 0 )
	  {
	    return -1;
	  }
      }
    else if (!S_ISDIR(st.st_mode))
      {
	return -1;
      }

    return 0;
  }

  extern "C" int grvy_create_unique_dir(char *name_template)
  {

    int status = verify_string_ends_with_6_Xs(name_template);
    if (status) return status;

    if (mkdtemp(name_template) == NULL)  /* mkdtemp sets errno */
      return -1;
    else
      return 0;
  }

  // Utility to provide a unique scratch directory for the user - note
  // that this utility has a very specific requirement on the input
  // string char template to include 6 trailing X's. We check to make
  // sure this is the case as some OS's (eg. OSX) will create the
  // directory without the trailing X's, but it will not have any
  // derived "uniqueness" at runtime.  Linux does not do this, so we
  // check explicitly to hopefully avoid confusion for the user.

#ifdef TLS_TODO /* Use thread local storage for the stack, if possible */
  static TLS std::stack<char *> _grvy_create_scratch_dir_paths;
#else
  static     std::stack<char *> _grvy_create_scratch_dir_paths;
#endif

  extern "C" int grvy_create_scratch_dir(char *name_template)
  {

    int status = verify_string_ends_with_6_Xs(name_template);
    if (status) return status;
    
    status = grvy_create_unique_dir(name_template);
    if (status) return status;
    
    char * name_copy = strdup(name_template);
    if (name_copy == NULL) {
      errno = ENOMEM;
      return -1;
    }
    
    if (_grvy_create_scratch_dir_paths.empty())
      {
	if (status = atexit(_GRVY_create_scratch_dir_atexit_handler))
	  {
	    free(name_copy);
	    return status;
	  }
      }

    _grvy_create_scratch_dir_paths.push(name_copy);

    return status;
  }

  void _GRVY_create_scratch_dir_atexit_handler()
  {
    while(!_grvy_create_scratch_dir_paths.empty())
    {
      char * path = _grvy_create_scratch_dir_paths.top();
      if (_GRVY_RemoveAll(path)) perror(path);
      free(path);
      _grvy_create_scratch_dir_paths.pop();
    }
  }

  int _GRVY_RemoveAll(const char *path)
  {
    struct stat st;
    int status;

    if (0 != lstat(path,&st) )
      {
	/* Succeed silently on "removal" of non-existent path */
	status = (errno == ENOENT) ? 0 : -1; /* errno set by lstat */
      }
    else if (S_ISREG(st.st_mode) || S_ISLNK(st.st_mode))
      {
	status = unlink(path);               /* errno set by unlink */
      }
    else if (S_ISDIR(st.st_mode))
      {
	status = nftw(path, _GRVY_RemoveAll_nftw_helper, _POSIX_OPEN_MAX,
	              FTW_DEPTH | FTW_PHYS); /* errno set by nftw */
      }
    else
      {
	/* unable to handle st.st_mode */
        status = -1;
	errno = EINVAL;
      }

    return status;
  }

  int _GRVY_RemoveAll_nftw_helper(const char *path,
                                  const struct stat * st,
				  int flag,
				  struct FTW *f)
  {
    int status;

    switch (flag) {
    case FTW_F:   /* File */
    case FTW_SL:  /* Symlink */
    case FTW_SLN: /* Symlink (dangling) */
      status = unlink(path);                 /* errno set by unlink */
      if (status) perror(path);
      break;
    case FTW_D:   /* Directory */
    case FTW_DP:  /* Directory that used to have subdirectories */
      status = rmdir(path);                  /* errno set by rmdir */
      if (status) perror(path);
      break;
    case FTW_DNR: /* Unreadable directory */
    case FTW_NS:  /* Permissions failed */
    default:      /* Unknown flag from nftw */
      /* Presumably related errors will arise on later rmdir attempts */
      status = 0;
    }

    return status;
  }

  //-----------------------------------------------------------------
  //                     Fortran Interfaces
  //-----------------------------------------------------------------

#ifdef _GRVY_FORTRAN_STRING_ORDER1
extern "C" void grvy_check_file_path_(char *pathname,int *flag,int _namelen)
#else
extern "C" void grvy_check_file_path_(char *pathname,int _namelen,int *flag)
#endif
  {
    char *name = grvy_f2c_char(pathname,_namelen);

    *flag = grvy_check_file_path(name);

    delete[] name;
    return;
  }

#ifdef _GRVY_FORTRAN_STRING_ORDER1
extern "C" void grvy_create_unique_dir_(char *name_template,int *flag,int _namelen)
#else
extern "C" void grvy_create_unique_dir_(char *name_template,int _namelen,int *flag)
#endif
  {
    char *c_name_template = grvy_f2c_char(name_template,_namelen);

    *flag = grvy_create_unique_dir(c_name_template);
    strncpy(name_template,c_name_template,_namelen);

    delete[] c_name_template;
  }

#ifdef _GRVY_FORTRAN_STRING_ORDER1
extern "C" void grvy_create_scratch_dir_(char *name_template,int *flag,int _namelen)
#else
extern "C" void grvy_create_scratch_dir_(char *name_template,int _namelen,int *flag)
#endif
  {
    char *c_name_template = grvy_f2c_char(name_template,_namelen);

    *flag = grvy_create_scratch_dir(c_name_template);
    strncpy(name_template,c_name_template,_namelen);

    delete[] c_name_template;
  }

  // ----------------------------------------------------------------
  // -------------------- Convenience Functions ---------------------
  // ----------------------------------------------------------------

  // grvy_f2c_char(): Convert evil Fortran character strings to C

  char *grvy_f2c_char(char*input,int len)
  {
    char* output = new char[len+1];
    strncpy(output,input,len);
    output[len]='\0';
    
    return(output);
  }

  
  // verify_string_ends_with_6_Xs(): 
  //
  // Verify that the template satisfies the necessary requirements to
  // end with 6 X's for use with mkdtemp() system function.

  int verify_string_ends_with_6_Xs(char *name_template)
  {
    const int NUM_REQUIRED_TRAILING_X = 6;
    int length_template;
    
    if(name_template == NULL)
      {
	return 1;
      }

    length_template = strlen(name_template);
    
    for(int i=length_template-NUM_REQUIRED_TRAILING_X;i<length_template;i++)
      {
	if(name_template[i] != 'X')
	  {
	    return 1;
	  }
      }
    
    return 0;
  }


}
