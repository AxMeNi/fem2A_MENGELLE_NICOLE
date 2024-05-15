#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quad = false;
    const bool t_elmap = false;
    const bool t_shpf = false;
    const bool t_asmblmat = false;
    const bool t_lcltgbl = false;
    const bool t_cdtdirich = false;
    const bool t_assemble_elementary_vector = false;
    const bool t_local_to_global_vector = false;
    const bool t_assemble_neumann = false;
    const bool t_poisson = true;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if( t_quad ) Tests::test_quadrature(2, false);
    if( t_elmap) Tests::test_ElementMapping();
    if( t_shpf) Tests::test_ShapeFunctions();
    if( t_asmblmat) Tests::test_assemble_elementary_matrix();
    if( t_lcltgbl) Tests::test_local_to_global_matrix() ;
    if( t_cdtdirich) Tests::test_apply_dirichlet_boundary_conditions();
    if( t_assemble_elementary_vector) Tests::test_assemble_elementary_vector();
    if( t_local_to_global_vector ) Tests::test_local_to_global_vector();
    if( t_assemble_neumann ) Tests::test_assemble_elementary_neumann_vector();
    if( t_poisson) Tests::test_poisson_problem("data/square_fine.mesh");
}

void run_simu()
{

    const bool simu_pure_dirichlet = false;
    const bool simu_dirichlet_source_term = false;
    const bool simu_sinus_bump = false;
    const bool simu_analytic_sinus_bumb = false;
    const bool simu_diff_sin = false;
    const bool simu_neumann_pb = false;
    const bool simu_mug_pb = false;
    //const bool simu_diff_poisson_orders = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) 
    {
        Simu::pure_dirichlet_pb("data/square.mesh", verbose);
    }
    
    if (simu_dirichlet_source_term ) 
    {
    	Simu::dirichlet_pb_source_term("data/square.mesh", verbose );
    }
    
    if (simu_sinus_bump ) 
    {
    	Simu::dirichlet_sinus_bump_pb("data/square_fine.mesh", verbose );
    }
    
    if (simu_analytic_sinus_bumb)
    {
    	Simu::sinus_bump_pb_analytic("data/square_fine.mesh", verbose);
    }
    
    if (simu_diff_sin)
    {
    	Simu::diff_pb_sin("data/square_fine.mesh", verbose);
    }
        
    if (simu_neumann_pb )
    {
    	Simu::neumann_pb("data/square.mesh", verbose );
    }
    
    if (simu_mug_pb )
    {
    	Simu::mug_pb("data/mug_0_5.mesh", verbose);
    }
    
    //if (simu_diff_poisson_orders )
    //{
    //	Simu::diff_poisson_orders("data/square_fine.mesh", verbose);
    //}
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }
    
    /// C'EST ICI QUE LES COMMANDES SONT EXECUTEES
    

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }
    return 0;
}





















