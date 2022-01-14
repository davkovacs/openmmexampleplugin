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

#include "ReferenceExampleKernels.h"
#include "ExampleForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include <iostream>
#include <julia.h>


using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

double unbox_float64(jl_value_t* jlval) {
   if jl_typeis(jlval, jl_float64_type) {
      return jl_unbox_float64(jlval);
   } 
   jl_errorf("Can't identify the return type!!! \n");
   return 0.0;
}

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

void ReferenceCalcExampleForceKernel::initialize(const System& system, const ExampleForce& force) {
    std::map<int, int> mass_to_Z {{1, 1}, {12, 6}, {14, 7}, {16, 8}, {19, 9}, {31, 15}, {32, 16}, {35, 17}, {80, 35}, {127, 53}};

    // Initialize bond parameters.
    
    // int numBonds = force.getNumBonds();
    // particle1.resize(numBonds);
    // particle2.resize(numBonds);
    // length.resize(numBonds);
    // k.resize(numBonds);
    // for (int i = 0; i < numBonds; i++)
    //     force.getBondParameters(i, particle1[i], particle2[i], length[i], k[i]);
}

double ReferenceCalcExampleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, jl_function_t*& _atoms_from_c, jl_value_t*& _energyfcn, jl_value_t*& _forcefn, jl_value_t*& _stressfcn) {
    vector<Vec3>& posData = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();  
    vector<int> Z;
    for (int i = 0; i < numParticles; i++)
    {
        double mass = context.getSystem().getParticleMass(i); 
        // cout << std::round(mass) << "\n";
        Z.push_back(std::round(mass));
        // force[i] += RealVec(2.0, 2.0, 2.0); 
    }

    jl_value_t* jl_float64 = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
    jl_value_t* jl_int32 = jl_apply_array_type((jl_value_t*)jl_int32_type, 1);
    jl_array_t* _X = jl_alloc_array_1d(jl_float64, 3 * numParticles);
    jl_array_t* _cell = jl_alloc_array_1d(jl_float64, 9);
    jl_array_t* _bc = jl_alloc_array_1d(jl_int32, 3);
    jl_array_t* _Z = jl_alloc_array_1d(jl_int32, numParticles);

    double *XData = (double*)jl_array_data(_X);
    for (int i = 0; i < numParticles; i++){
        for (int j = 0; j < 3; j++)
        {
            XData[3*i+j] = posData[i][j];
        }
    }

    double cell [9] = {50.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 50.0};
    double *cellData = (double*)jl_array_data(_cell);
    for (int i = 0; i < 9; i++) cellData[i] = cell[i];

    int bc [3] = {0, 0, 0};
    int32_t *bcData = (int32_t*)jl_array_data(_bc);
    for (int i = 0; i < 3; i++) bcData[i] = bc[i];

    int32_t *ZData = (int32_t*)jl_array_data(_Z);
    for (int i = 0; i < numParticles; i++) ZData[i] = Z[i]; 

    jl_value_t* at_args[] = {(jl_value_t*)_X, (jl_value_t*)_Z, (jl_value_t*)_cell, (jl_value_t*)_bc};   
    jl_value_t* at = jl_call(_atoms_from_c, at_args, 4);

    jl_value_t** args; 
    JL_GC_PUSHARGS(args, 2);

    jl_value_t* calc = jl_eval_string("IP");
    args[0] = calc; 
    args[1] = at;

    jl_set_global(jl_main_module, jl_symbol("at"), at);

    double E = 0.0; 
    jl_value_t* jlE;
    cout << "Calling the Juila energy function" << "\n";

    if (jl_exception_occurred())
        printf("Exception before calling energy function : %s \n", 
                jl_typeof_str(jl_exception_occurred()));

    jlE = jl_call2(_energyfcn, calc, at);   
    
    if (jl_exception_occurred()) {
        printf("Exception at jl_call2(_energyfcn, calc, at) : %s \n", 
                jl_typeof_str(jl_exception_occurred()));
        jl_errorf("Something when wrong calling the energy\n");
    } else {
        E = unbox_float64(jlE); 
    }
    
    JL_GC_POP();

    return E;

    // vector<RealVec>& pos = extractPositions(context);
    // vector<RealVec>& force = extractForces(context);
    // int numBonds = particle1.size();
    // double energy = 0;
    
    // // Compute the interactions.
    
    // for (int i = 0; i < numBonds; i++) {
    //     int p1 = particle1[i];
    //     int p2 = particle2[i];
    //     RealVec delta = pos[p1]-pos[p2];
    //     RealOpenMM r2 = delta.dot(delta);
    //     RealOpenMM r = sqrt(r2);
    //     RealOpenMM dr = (r-length[i]);
    //     RealOpenMM dr2 = dr*dr;
    //     energy += k[i]*dr2*dr2;
    //     RealOpenMM dEdR = 4*k[i]*dr2*dr;
    //     dEdR = (r > 0) ? (dEdR/r) : 0;
    //     force[p1] -= delta*dEdR;
    //     force[p2] += delta*dEdR;
    // }
    // return energy;
}

void ReferenceCalcExampleForceKernel::copyParametersToContext(ContextImpl& context, const ExampleForce& force) {
    if (force.getNumBonds() != particle1.size())
        throw OpenMMException("updateParametersInContext: The number of Example bonds has changed");
    for (int i = 0; i < force.getNumBonds(); i++) {
        int p1, p2;
        force.getBondParameters(i, p1, p2, length[i], k[i]);
        if (p1 != particle1[i] || p2 != particle2[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
}
