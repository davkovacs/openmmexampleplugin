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

// turn Float64 returned by Julia into C Float number
double unbox_float64(jl_value_t* jlval) {
   if jl_typeis(jlval, jl_float64_type) {
      return jl_unbox_float64(jlval);
   } 
   jl_errorf("Can't identify the return type!!! \n");
   throw OpenMMException("Cannot identify return type of energy");;
}

// map masses to atomic numbers (necessary for ACE)
std::map<int, int> mass_to_Z {{1, 1}, {12, 6}, {14, 7}, {16, 8}, {19, 9}, {31, 15}, {32, 16}, {35, 17}, {80, 35}, {127, 53}};

// unit conversion factor between eV and kJ/mol
const double eV_to_kJ_per_mol = 96.48533212331002; // unit conversion factor from WolframAlpha

// get the poistions of the atoms
static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

// get the forces on the atoms
static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

void ReferenceCalcExampleForceKernel::initialize(const System& system, const ExampleForce& force) {
}

double ReferenceCalcExampleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, jl_function_t*& _atoms_from_c, jl_value_t*& _energyfcn, jl_value_t*& _forcefcn, jl_value_t*& _stressfcn) {
    vector<Vec3>& posData = extractPositions(context);  // get atomic positions
    vector<RealVec>& force = extractForces(context);  // get forces on the atoms
    int numParticles = context.getSystem().getNumParticles();  // get number of particles
    int Z[numParticles];  // create vector to store atomic numbers
    // loop over the particles
    for (int i = 0; i < numParticles; i++)
    {
        double mass = context.getSystem().getParticleMass(i);  // get the mass of the particle
        Z[i] = mass_to_Z[std::round(mass)];  // add atomic number to array
    }

    // create objects to store necessary Julia arrays with correct types 
    jl_value_t* jl_float64 = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
    jl_value_t* jl_int32 = jl_apply_array_type((jl_value_t*)jl_int32_type, 1);
    jl_array_t* _X = jl_alloc_array_1d(jl_float64, 3 * numParticles);  // positions
    jl_array_t* _cell = jl_alloc_array_1d(jl_float64, 9);  // cell
    jl_array_t* _bc = jl_alloc_array_1d(jl_int32, 3);  // boundary conditions
    jl_array_t* _Z = jl_alloc_array_1d(jl_int32, numParticles);  // atomic numbers

    // Pass positions to Julia
    double *XData = (double*)jl_array_data(_X);
    for (int i = 0; i < numParticles; i++){
        for (int j = 0; j < 3; j++)
        {
            XData[3*i+j] = 10 * posData[i][j];  // nm to A
        }
    }
    //cout << "############## NO PROBLEM HERE ##################" << "\n";

    // Set the default cell size for non-periodic simulations to 1.0 (doesn't matter what it is)
    double cell [9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    bool uses_pbc  = context.getSystem().usesPeriodicBoundaryConditions();
    int bc [3] = {0, 0, 0};
    Vec3 box[3];
    if (uses_pbc == true)  // WARNING! PBC assumed in ALL direction or NON
    {
        context.getPeriodicBoxVectors(box[0], box[1], box[2]);  // get periodic box vectors in nm
        for (int i = 0; i < 3; i++)
        {
            cell[3*i] = box[i][0] * 10.0;  // nm to A
            cell[3*i+1] = box[i][1] * 10.0;
            cell[3*i+2] = box[i][2] * 10.0;
            bc[i] = 1;
        }
    } 
    
    // Pass the cell to Julia
    double *cellData = (double*)jl_array_data(_cell);
    for (int i = 0; i < 9; i++) cellData[i] = cell[i];

    // Pass boundary conditions to Julia
    int32_t *bcData = (int32_t*)jl_array_data(_bc);
    for (int i = 0; i < 3; i++) bcData[i] = bc[i];

    // Pass atomic numbers to Julia
    int32_t *ZData = (int32_t*)jl_array_data(_Z);
    for (int i = 0; i < numParticles; i++) ZData[i] = Z[i]; 

    // create atoms object in Julia
    jl_value_t* at_args[] = {(jl_value_t*)_X, (jl_value_t*)_Z, (jl_value_t*)_cell, (jl_value_t*)_bc};   
    jl_value_t* at = jl_call(_atoms_from_c, at_args, 4);

    // ensure args doesn't get garbage collected
    jl_value_t** args; 
    JL_GC_PUSHARGS(args, 2);

    jl_value_t* calc = jl_eval_string("IP");
    args[0] = calc; 
    args[1] = at;

    // create a global variable for the atoms object at
    jl_set_global(jl_main_module, jl_symbol("at"), at);
    //jl_eval_string("@show at");  // can be used to print the Atoms object from Julia

    // compute the energy
    double E = 0.0;  
    jl_value_t* jlE;

    jlE = jl_call2(_energyfcn, calc, at);  // calling the energy function energy(IP, at) 

    if (jl_exception_occurred()) {
        printf("Exception at jl_call2(_energyfcn, calc, at) : %s \n", 
                jl_typeof_str(jl_exception_occurred()));
        jl_errorf("Something when wrong calling the energy\n");
        throw OpenMMException("Something went wrong calling the energy");;
    } else {
        E = unbox_float64(jlE);  // unbox energy from Julia Float to C++
    }

    // compute the forces
    jl_array_t* _F = (jl_array_t*)jl_call2(_forcefcn, calc, at);
    double *Fdata = (double*)jl_array_data(_F); 

    // write forces to the force array. (unit conversion done here)
    for (int i = 0; i < numParticles; i++)
    {
        force[i] += RealVec(Fdata[3*i] * eV_to_kJ_per_mol * 10, Fdata[3*i+1] * eV_to_kJ_per_mol * 10, Fdata[3*i+2] * eV_to_kJ_per_mol * 10); 
    }
    
    JL_GC_POP();

    return E * eV_to_kJ_per_mol;  // Return energy in kJ/mol
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
