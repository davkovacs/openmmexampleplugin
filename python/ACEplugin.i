%module ACEplugin

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include <std_string.i>

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "ACEForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import openmm as mm
import simtk.unit as unit
%}

/*
 * Add units to function outputs.
*/
%pythonappend ACEPlugin::ACEForce::getBondParameters(int index, int& particle1, int& particle2,
                                                             double& length, double& k) const %{
    val[2] = unit.Quantity(val[2], unit.nanometer)
    val[3] = unit.Quantity(val[3], unit.kilojoule_per_mole/unit.nanometer**4)
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


namespace ACEPlugin {

class ACEForce : public OpenMM::Force {
public:
    ACEForce(const char* IP_path);

    const std::string& get_IP_path() const;

    void setUsesPeriodicBoundaryConditions(bool periodic);

    void setAtomicNumbers(std::vector<int> atomic_numbs);

    void setAtomInds(std::vector<int> at_inds);

    bool usesPeriodicBoundaryConditions() const;

    void updateParametersInContext(OpenMM::Context& context);

    /*
     * Add methods for casting a Force to an ACEForce.
    */
    %extend {
        static ACEPlugin::ACEForce& cast(OpenMM::Force& force) {
            return dynamic_cast<ACEPlugin::ACEForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<ACEPlugin::ACEForce*>(&force) != NULL);
        }
    }
};

}
