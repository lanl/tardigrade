/* ---------------------------------------------------------------------
 *       _                           _
 *   ___| |__   __ _ _ __ ___   ___ (_)___
 *  / __| '_ \ / _` | '_ ` _ \ / _ \| / __|
 * | (__| | | | (_| | | | | | | (_) | \__ \
 *  \___|_| |_|\__,_|_| |_| |_|\___/|_|___/
 *
 * Chamois - a MOOSE interface to constitutive models developed at the
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of chamois.
 * ---------------------------------------------------------------------
 */

#include "GradientEnhancedDamagedMicromorphicContinuumAction.h"
#include <string>
#include <vector>
#include "FEProblem.h"
#include "Factory.h"

registerMooseAction( "tardigradeApp", GradientEnhancedDamagedMicromorphicContinuumAction, "add_kernel" );

registerMooseAction( "tardigradeApp", GradientEnhancedDamagedMicromorphicContinuumAction, "add_material" );

const std::vector< std::string > GradientEnhancedDamagedMicromorphicContinuumAction::excludedParameters = {
    "model_name",
    "material_fparameters",
    /* "save_in_disp_x", */
    /* "save_in_disp_y", */
    /* "save_in_disp_z", */
};

InputParameters
GradientEnhancedDamagedMicromorphicContinuumAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addRequiredCoupledVar( "displacements",
                                "The components of the displacement vector field" );
  params.addRequiredCoupledVar(
      "micro_displacement_gradient", "The components of the micro displacement gradient" );

  /* params.addRequiredCoupledVar( "nonlocal_damage", "The nonlocal damage field" ); */

  params.addParam< std::string >( "base_name", "Material property base name" );
  /* params.addParam< std::vector< AuxVariableName > >( "save_in_disp_x", */
  /*                                                    "Store displacement residuals" ); */
  /* params.addParam< std::vector< AuxVariableName > >( "save_in_disp_y", */
  /*                                                    "Store displacement residuals" ); */
  /* params.addParam< std::vector< AuxVariableName > >( "save_in_disp_z", */
  /*                                                    "Store displacement residuals" ); */
  params.addParam< std::vector< SubdomainName > >( "block",
                                                   "The list of ids of the blocks (subdomain) "
                                                   "that the kernels will be "
                                                   "applied to" );

    // Vectors of material properties
    params.addRequiredParam<std::vector<Real>>(
        "material_fparameters", "The vector of floating point material parameters required for the stiffness matrices");

    params.addRequiredParam<std::string>(
        "model_name", "The material model name");

    // The state variable array
    params.addParam< int >(
        "number_SDVS", 0, "The number of solution-dependent state variables" );

  return params;
}

GradientEnhancedDamagedMicromorphicContinuumAction::GradientEnhancedDamagedMicromorphicContinuumAction(
    const InputParameters & parameters )
  : Action( parameters ),
    _ndisp( getParam< std::vector< VariableName > >( "displacements" ).size() ),
    _nmicro_disp_grad ( getParam< std::vector< VariableName > >( "micro_displacement_gradient" ).size() )
{
  if ( _ndisp != 3 || _nmicro_disp_grad  != 9 )
    mooseError( "Gradient-enhanced micromorphic kernels are implemented only for 3D!" );

  /* if ( parameters.isParamSetByUser( "use_displaced_mesh" ) ) */
  /* { */
  /*   bool use_displaced_mesh_param = getParam< bool >( "use_displaced_mesh" ); */
  /*   if ( use_displaced_mesh_param ) */
  /*     mooseError( "use_displaced_mesh must be set to false" ); */
  /* } */
}

void
GradientEnhancedDamagedMicromorphicContinuumAction::act()
{
  if ( _current_task == "add_kernel" )
    addKernels();
  else if ( _current_task == "add_material" )
    addMaterial();
}

void
GradientEnhancedDamagedMicromorphicContinuumAction::addKernels()
{

    std::cout << "add kernel" << std::endl;
  /* std::string pki_stress_kernel( "GradientEnhancedDamagedMicromorphicPKIDivergence" ); */

  /* for ( unsigned int i = 0; i < _ndisp; ++i ) */
  /* { */
  /*   const std::string kernel_name = name() + "_div_pki_stress_" + Moose::stringify( i ); */

  /*   InputParameters pki_stress_kernel_params = _factory.getValidParams( pki_stress_kernel ); */
  /*   pki_stress_kernel_params.applyParameters( parameters(), excludedParameters ); */

  /*   pki_stress_kernel_params.set< unsigned int >( "component" ) = i; */
  /*   pki_stress_kernel_params.set< NonlinearVariableName >( "variable" ) = */
  /*       getParam< std::vector< VariableName > >( "displacements" )[i]; */
  /*   pki_stress_kernel_params.set< std::string >( "tensor" ) = "pk_i_stress"; */

  /*   if ( i == 0 && isParamValid( "save_in_disp_x" ) ) */
  /*     pki_stress_kernel_params.set< std::vector< AuxVariableName > >( "save_in" ) = */
  /*         getParam< std::vector< AuxVariableName > >( "save_in_disp_x" ); */

  /*   if ( i == 1 && isParamValid( "save_in_disp_y" ) ) */
  /*     pki_stress_kernel_params.set< std::vector< AuxVariableName > >( "save_in" ) = */
  /*         getParam< std::vector< AuxVariableName > >( "save_in_disp_y" ); */

  /*   if ( i == 2 && isParamValid( "save_in_disp_z" ) ) */
  /*     pki_stress_kernel_params.set< std::vector< AuxVariableName > >( "save_in" ) = */
  /*         getParam< std::vector< AuxVariableName > >( "save_in_disp_z" ); */

  /*   _problem->addKernel( pki_stress_kernel, kernel_name, pki_stress_kernel_params ); */
  /* } */

  std::string internal_force_kernel_name( "InternalForce" );
  for ( unsigned int i = 0; i < 3 ; ++i )
  /* for ( unsigned int j = 0; j < 3 ; ++j ) */
  {
    const std::string the_kernel_name = name() + "_internal_force_" + Moose::stringify( i  );

    InputParameters internal_force_kernel_params = _factory.getValidParams( internal_force_kernel_name );
    internal_force_kernel_params.applyParameters( parameters(), excludedParameters );

    internal_force_kernel_params.set< int >( "component" ) = i;
    internal_force_kernel_params.set< int >( "dof_num" ) = i;

    internal_force_kernel_params.set< NonlinearVariableName >( "variable" ) = getParam< std::vector< VariableName > >( "displacements" )[i];

    internal_force_kernel_params.set< std::vector< VariableName> >( "u1" ) = {getParam< std::vector< VariableName > >( "displacements" )[0]};
    internal_force_kernel_params.set< std::vector< VariableName> >( "u2" ) = {getParam< std::vector< VariableName > >( "displacements" )[1]};
    internal_force_kernel_params.set< std::vector< VariableName> >( "u3" ) = {getParam< std::vector< VariableName > >( "displacements" )[2]};

    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_11" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[0] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_12" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[1] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_13" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[2] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_21" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[3] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_22" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[4] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_23" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[5] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_31" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[6] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_32" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[7] } ;
    internal_force_kernel_params.set< std::vector< VariableName> >( "phi_33" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[8] } ;

    _problem->addKernel( internal_force_kernel_name, the_kernel_name, internal_force_kernel_params );

  }

  std::string internal_couple_kernel_name( "InternalCouple" );
  for ( unsigned int i = 0; i < 3 ; ++i )
  for ( unsigned int j = 0; j < 3 ; ++j )
  {
    const std::string the_kernel_name = name() + "_internal_couple_" + Moose::stringify( i * 3 + j );

    InputParameters internal_couple_kernel_params = _factory.getValidParams( internal_couple_kernel_name );
    internal_couple_kernel_params.applyParameters( parameters(), excludedParameters );

    internal_couple_kernel_params.set< int >( "component_i" ) = i;
    internal_couple_kernel_params.set< int >( "component_j" ) = j;
    internal_couple_kernel_params.set< int >( "dof_num" ) = 3 + 3*i + j;

    internal_couple_kernel_params.set< NonlinearVariableName >( "variable" ) = getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[i*3 + j];

    internal_couple_kernel_params.set< std::vector< VariableName> >( "u1" ) = {getParam< std::vector< VariableName > >( "displacements" )[0]};
    internal_couple_kernel_params.set< std::vector< VariableName> >( "u2" ) = {getParam< std::vector< VariableName > >( "displacements" )[1]};
    internal_couple_kernel_params.set< std::vector< VariableName> >( "u3" ) = {getParam< std::vector< VariableName > >( "displacements" )[2]};

    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_11" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[0] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_12" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[1] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_13" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[2] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_21" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[3] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_22" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[4] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_23" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[5] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_31" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[6] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_32" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[7] } ;
    internal_couple_kernel_params.set< std::vector< VariableName> >( "phi_33" ) = { getParam< std::vector< VariableName > >( "micro_displacement_gradient" )[8] } ;

    _problem->addKernel( internal_couple_kernel_name, the_kernel_name, internal_couple_kernel_params );

  }

  /* std::string nonlocal_damage_kernel( "GradientEnhancedDamagedMicromorphicDamage" ); */
  /* InputParameters nonlocal_damage_kernel_params = _factory.getValidParams( nonlocal_damage_kernel ); */

  /* nonlocal_damage_kernel_params.applyParameters( parameters(), excludedParameters ); */

  /* const std::string kernel_name = name() + "_nonlocal_damage"; */

  /* nonlocal_damage_kernel_params.set< NonlinearVariableName >( "variable" ) = */
  /*     getParam< std::vector< VariableName > >( "nonlocal_damage" )[0]; */

  /* _problem->addKernel( nonlocal_damage_kernel, kernel_name, nonlocal_damage_kernel_params ); */
}

void
GradientEnhancedDamagedMicromorphicContinuumAction::addMaterial()
{
  std::string materialType = "GradientEnhancedDamagedMicromorphicMaterial";

  auto materialParameters = _factory.getValidParams( materialType );
  materialParameters.applyParameters( parameters() );

  materialParameters.set< std::string >( "model_name" ) =
      getParam< std::string >( "model_name" );
  materialParameters.set< std::vector< Real > >( "material_fparameters" ) =
      getParam< std::vector< Real > >( "material_fparameters" );

  _problem->addMaterial( materialType, name() + "_material", materialParameters );
}
