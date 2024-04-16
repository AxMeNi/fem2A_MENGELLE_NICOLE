#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################
        
        double unit_fct(FEM2A::vertex vertex) {
    		return 1.0; // Toujours renvoyer 1 pour n'importe quelle valeur de vertex
	}

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            	std::cout << "Solving a pure Dirichlet problem" << std::endl;
            	if ( verbose ) {
                	std::cout << " with lots of printed details..." << std::endl;
            	}
            	std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
        	Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping elt_mapping(mesh, false, 4);
		    
		ShapeFunctions reference_functions(2, 1);
		    
		Quadrature quad;
		Quadrature quadrature = quad.get_quadrature(2, false);
		    
		DenseMatrix Ke_in;
		    
		assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke_in);
	   	
		SparseMatrix K (mesh.nb_vertices());
        	local_to_global_matrix(mesh, 4, Ke_in, K);
        	
        	
        	
        	// Apply dirichlet condition
        	std::vector< bool >& attribute_is_dirichlet (mesh.nb_triangles()); 
        	for (int i = 0; i < mesh.nb_triangles() ; i++) {
        		mesh.set_attribute(unit_fct, 1, false)
        		attribute_is_dirichlet[i] = true;        		
        	}
        	td::vector< double > values (mesh.nb_vertices);
        	for (int i = 0; i < mesh.nb_vertices() ; i++) {
        		values[i] = mesh.get_vertex( i ).x + mesh.get_vertex( i ).y
        	}
        	
        	
        	apply_dirichlet_boundary_conditions(mesh, unit_fct, attribute_is_dirichlet, values)
        }

    }

}
