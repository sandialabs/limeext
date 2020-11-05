// ********************************************************************************************
// LIME 1.0:  Lightweight Integrating Multiphysics Environment for coupling codes, Version 1.0
// Copyright (c) 2012, Sandia National Laboratories
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted 
// provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright notice, this list of 
//      conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright notice, this list of 
//      conditions and the following disclaimer in the documentation and/or other materials  
//      provided with the distribution.
//    * Neither the name of Sandia National Laboratories nor the names of its contributors may 
//      be used to endorse or promote products derived from this software without specific prior 
//      written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// (Note: This is a BSD 3-Clause License.  For info. see www.opensource.org/licenses/BSD-3-Clause )
// --------------------------------------------------------------------------------------------
//
// LIME_Broadcast_File.cpp
//
// Description: Broadcast the contents of the named file to all processes in this
// MPI communicator.  The root MPI rank and the file arguments are checked for errors.
// Two broadcasts occur, the first sends the size of the file in bytes so that receiving
// processes can allocate buffers big enough to hold the file contents.  The second
// broadcast sends the file contents.  This routine won't work if the file to be read
// is order 2 Gb or more in size.
//
// 
// ********************************************************************************************

//-----------------------------------------------------------------------------
// Include Directives

#include <stdexcept>
#include <string>
#include <vector>

#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "LIME_Broadcast_File.hpp"

using std::runtime_error;
using std::string;
using std::vector;

//-----------------------------------------------------------------------------
// 

vector<char> 
LIME::broadcast_file(string filename, MPI_Comm comm, unsigned int root)
{
  // buffer to hold file (don't allocate yet)
  vector<char> buffer(0);

#ifdef HAVE_MPI
  // check passed arguments
  int size = 0;
  MPI_Comm_size(comm, &size);

  if (static_cast<unsigned int>(size) <= root) {
    string s = "Error: broadcast_file passed invalid processor rank";
    throw runtime_error(s);
  }

  // get my rank
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  // on root processor
  if (root == static_cast<unsigned int>(rank)) {
    // current working directory
    char* cwd = getcwd(0, 0);

    // form fqdn (fully qualified domain name)
    string setup_file = cwd;
    setup_file += "/" + filename;

    // release getcwd allocated memory
    free(cwd);

    // throw an execption if the file doesn't exist
    struct stat sb;
    if (stat(setup_file.c_str(), &sb)) {
      string s = strerror(errno);
      throw runtime_error(s.c_str());
    }

    // open file to be broadcast
    int fd = open(setup_file.c_str(), O_RDONLY, S_IRWXU | S_IRWXG | S_IRWXO);
    if (fd < 0) {
      string s = "Error: unable to open file " + setup_file + " for reading";
      throw runtime_error(s.c_str());
    }

    // allocate enough room for file
    buffer.resize(sb.st_size);

    // read file into buffer
    ssize_t bytes_read = read(fd, &buffer[0], sb.st_size);

    // check if problems reading data, or not all data read in
    if (bytes_read < 0 || bytes_read != sb.st_size) {
      string s = "Error: unable to read file " + setup_file;
      throw runtime_error(s.c_str());
    }

    // done with file
    close(fd);
  }

  // broadcast how big file is, and the buffer
  size = buffer.size();
  MPI_Bcast(&size, 1, MPI_INT, root, comm);

  // non-root processors allocate buffer here
  if (root != static_cast<unsigned int>(rank)) {
     buffer.resize(size);
  }

  // root pass copy of input file to all other processors
  MPI_Bcast(&buffer[0], size, MPI_CHAR, root, comm);
#endif

  return buffer;
}
