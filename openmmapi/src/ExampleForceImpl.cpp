/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "internal/ExampleForceImpl.h"
#include "ExampleKernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <julia.h>
#include <iostream>

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

ExampleForceImpl::ExampleForceImpl(const ExampleForce& owner) : owner(owner) {
    ip_path = owner.ace_path;
}

ExampleForceImpl::~ExampleForceImpl() {
}

jl_function_t* _atoms_from_c; 
jl_value_t* _energyfcn;
jl_value_t* _forcefcn;
jl_value_t* _stressfcn;

void ExampleForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcExampleForceKernel::Name(), context);
    kernel.getAs<CalcExampleForceKernel>().initialize(context.getSystem(), owner);
    jl_init();  // start up Julia process and load packages
    jl_eval_string("try using ACE1 catch; using ACE end;");
    jl_eval_string("using JuLIP");
    // define function pointers for creating Julia Atoms objects, energy, forces and stress
    _atoms_from_c = jl_eval_string("(X, Z, cell, bc) -> Atoms(X = X, Z = Z, cell=cell, pbc = Bool.(bc))");
    _energyfcn = (jl_value_t*)jl_get_function(jl_main_module, "energy");
    _forcefcn = (jl_value_t*)jl_eval_string("(calc, at) -> mat(forces(calc, at))[:]");
    _stressfcn = (jl_value_t*)jl_eval_string("(calc, at) -> vcat(stress(calc, at)...)");

    const char* comm1 = "IP = read_dict( load_dict(\"";
    const char* comm2 = "\")[\"IP\"])";
    char* read_ip;
    read_ip = (char*)calloc(strlen(comm1) + strlen(this->ip_path.c_str()) + strlen(comm2) + 1, sizeof(char));
    strcpy(read_ip, comm1);
    strcat(read_ip, this->ip_path.c_str());
    strcat(read_ip, comm2);
    jl_eval_string(read_ip);
    free(read_ip);
    //jl_eval_string("IP = read_dict( load_dict(\"/home/cdt1906/Documents/phd/ACE_dev/interfaces/test_openmm/CH_ace_test.json\")[\"IP\"])");
    
    // check if there were any Julia exceptions and print them
    if (jl_exception_occurred()){
        printf("ERROR DURING INITIALIZATION: %s \n", jl_typeof_str(jl_exception_occurred()));
        throw OpenMMException("Julia error during initialization of ACE potential");
    }
}

double ExampleForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        // Call the ACE energy and force calculation, pass the necessary function pointers
        return kernel.getAs<CalcExampleForceKernel>().execute(context, includeForces, includeEnergy, _atoms_from_c, _energyfcn, _forcefcn, _stressfcn);
    return 0.0;
}

std::vector<std::string> ExampleForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcExampleForceKernel::Name());
    return names;
}

vector<pair<int, int> > ExampleForceImpl::getBondedParticles() const {
    throw OpenMMException("getBondedParticles: Not defined bonded particles in ACE");
    // int numBonds = owner.getNumBonds();
    // vector<pair<int, int> > bonds(numBonds);
    // for (int i = 0; i < numBonds; i++) {
    //     double length, k;
    //     owner.getBondParameters(i, bonds[i].first, bonds[i].second, length, k);
    // }
    // return bonds;
}

void ExampleForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcExampleForceKernel>().copyParametersToContext(context, owner);
}
