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
// LIME_Strings.hpp
//
// Description: Certain types of problems deal with string data and, occasionally, this
// string data needs to be known to multiple model evaluators.  For these kinds of problems
// I create this file to act as the repository of common variable names needed by those
// shared physics.  All strings are std::strings and are declared extern here in the header
// so we don't get multiply defined symbol names from the compiler or linker.  Note that
// we use the _str suffix to help readers recognize text string data from other kinds of data.
// 
// ********************************************************************************************

#ifndef lime_strings_hpp
#define lime_strings_hpp

#include <string>

//-----------------------------------------------------------------------------

// other frequently used strings

extern std::string us;

// physics related strings

extern std::string constants_str;
extern std::string physics_str;
extern std::string name_str;
extern std::string value_str;
extern std::string instance_str;
extern std::string input_file_str;
extern std::string input_mesh_file_str;
extern std::string output_mesh_file_str;
extern std::string fluid_str;
extern std::string solid_str;
extern std::string delta_str;
extern std::string sfactor_str;
extern std::string sample_size_str;
extern std::string resample_str;
extern std::string distribution_str;
extern std::string lower_str;
extern std::string upper_str;
extern std::string unit_str;
extern std::string mean_str;
extern std::string stdev_str;
extern std::string lambda_str;

// transfer related strings

extern std::string transfer_str;
extern std::string source_str;
extern std::string target_str;
extern std::string type_str;
extern std::string pre_str;
extern std::string post_str;
extern std::string mixed_str;

#endif
