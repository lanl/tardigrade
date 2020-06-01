/*!
 * Output the results of the simulation in XDMF format
 */

#ifndef XDMF_HEADER
#define XDMF_HEADER

//MOOSE includes
#include "AdvancedOutput.h"

//XDMF includes
//#include "XdmfDomain.hpp"
//#include "XdmfWriter.hpp"

//Forward declarations
class Xdmf;

template <>
InputParameters validParams<Xdmf>();

/**
 * Class for output data to the XDMF format
 */

class Xdmf : public AdvancedOutput
{
    public:
        static InputParameters validParams();

        enum class OutputDimension : int
        {
            DEFAULT,
            ONE,
            TWO,
            THREE,
            PROBLEM_DIMENSION
        };

        /*!
         * Class constructor
         */

        Xdmf( const InputParameters & parameters );

        /*!
         * Initialization method
         */

        void initialSetup();

        /*!
         * Class destructor
         */

        ~Xdmf(){};
};

#endif
