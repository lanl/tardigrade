
#ifndef MULTIAPPOVERLAPCOUPLINGTRANSFER_H
#define MULTIAPPOVERLAPCOUPLINGTRANSFER_H

    #include "MultiAppTransfer.h"
    #include "OverlapCoupling.h"

    class MultiAppOverlapCouplingTransfer;

    template <>
    InputParameters validParams< MultiAppOverlapCouplingTransfer >( );

    /*!
     * Copies attributes from an overlap coupling user object in the main domain
     * to an overlap coupling user object in the sub domain. 
     */

    class MultiAppOverlapCouplingTransfer : public MultiAppTransfer{

        public:

            MultiAppOverlapCouplingTransfer( const InputParameters & parameters );

            virtual void execute( ) override;

        protected:
            const std::string &_source_user_object_name;
            const std::string &_target_user_object_name;

    };

#endif
