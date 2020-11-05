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
// LIME_Strings.cpp
//
// Description: Here we initialize the global string data shared by multiple model evaluators.
// 
// ********************************************************************************************

#include <LIME_Strings.hpp>

#include <string>

//-----------------------------------------------------------------------------

using std::string;

string us = "_";

string constants_str = "Constants";
string physics_str = "Physics";  // also used by transfer operators
string name_str = "Name";
string value_str = "value";
string instance_str = "Instance";
string input_file_str = "Input File";
string input_mesh_file_str = "Input Mesh File";
string output_mesh_file_str = "Output Mesh File";
string fluid_str = "Fluid";
string solid_str = "Solid";
string delta_str = "Delta";
string sfactor_str = "Sfactor";
string sample_size_str = "Sample Size";
string resample_str = "Resample";
string distribution_str = "distribution";
string lower_str = "lower";
string upper_str = "upper";
string unit_str = "unit";
string mean_str = "mean";
string stdev_str = "stdev";
string lambda_str = "lambda";

string transfer_str = "Transfer";
string source_str = "Source";
string target_str = "Target";
string type_str = "Type";
string pre_str = "Pre";
string post_str = "Post";
