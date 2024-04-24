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

	double sinus_bump_coefficient(FEM2A::vertex vertex) {
        	double coeff = 2*pow(M_PI,2)*sin(M_PI*vertex.x)*sin(M_PI*vertex.y);
    		return coeff; // Toujours renvoyer 1 pour n'importe quelle valeur de vertex
	}
	
        //#################################
        //  Simulations
        //#################################
        




        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            	std::cout << "Solving a pure Dirichlet problem" << std::endl;
            	if ( verbose ) {
                	std::cout << " with lots of printed details..." << std::endl;
            	}
            	std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
        	Mesh mesh;
		mesh.load(mesh_filename);
		mesh.set_attribute(unit_fct, 1, true);      
		
		
		//Initiatlisation de F et de K
		// K
		SparseMatrix K (mesh.nb_vertices());
		
		// F
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Création de K et F
		for (int i = 0; i < mesh.nb_triangles(); i++) {
			ElementMapping elt_mapping( mesh, false, i);
			ShapeFunctions reference_functions(2, 1);
			Quadrature quadrature = Quadrature::get_quadrature(2);
			// K
			DenseMatrix Ke_in;
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, unit_fct, Ke_in);
			local_to_global_matrix(mesh, i, Ke_in, K);
	   	}       	
        	
        	// Apply dirichlet condition
        	// cf ligne 414 du fem.cpp
        	std::vector< bool > attribute_is_dirichlet (2); 
        	attribute_is_dirichlet[0] = false; 
        	attribute_is_dirichlet[1] = true; 
        	
        	std::vector< double > values (mesh.nb_vertices());
        	for (int i = 0; i < mesh.nb_vertices() ; i++) {
        		values[i] = xy_fct(mesh.get_vertex(i));
        	}
        	
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	std::vector<double> x (mesh.nb_vertices());
        	FEM2A::solve(K,F,x);
        	mesh.save("solutions/pure_dirich_fine_square.mesh");
        	save_solution( x, "solutions/pure_dirich_fine_square.bb" );
        }
	
	void dirichlet_pb_source_term( const std::string& mesh_filename, bool verbose )
	{
            	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) {
                	std::cout << " with lots of printed details..." << std::endl;
            	}
            	std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
        	Mesh mesh;
		mesh.load(mesh_filename);
		mesh.set_attribute(unit_fct, 1, true);      
		
		
		//Initiatlisation de F et de K
		// K
		SparseMatrix K (mesh.nb_vertices());
		
		// F
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Création F et Fe
		for (int i = 0; i < mesh.nb_triangles(); i++) {
			ElementMapping elt_mapping( mesh, false, i);
			ShapeFunctions reference_functions(2, 1);
			Quadrature quadrature = Quadrature::get_quadrature(2);
			// F
			std::vector <double> Fe_in;
			assemble_elementary_vector (elt_mapping, reference_functions, quadrature, sinus_bump_coefficient, Fe_in);
			local_to_global_vector(mesh, false, i, Fe_in, F );
			// K
			DenseMatrix Ke_in;
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, sinus_bump_coefficient, Ke_in);
			local_to_global_matrix(mesh, i, Ke_in, K);
			
	   	}       	
        	
        	// Apply dirichlet condition
        	// cf ligne 414 du fem.cpp
        	std::vector< bool > attribute_is_dirichlet (2); 
        	attribute_is_dirichlet[0] = false; 
        	attribute_is_dirichlet[1] = true; 
        	
        	std::vector< double > values (mesh.nb_vertices());
        	for (int i = 0; i < mesh.nb_vertices() ; i++) {
        		values[i] = 0;
        	}
        	
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	std::vector<double> x (mesh.nb_vertices());
        	FEM2A::solve(K,F,x);
        	mesh.save("solutions/dirich_source_term_sinusbump_fine_square.mesh");
        	save_solution( x, "solutions/dirich_source_term_sinusbump_fine_square.bb" );
        }
        
        void neumann_condition_problem(const std::string& mesh_filename, bool verbose )
	{
            	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) {
                	std::cout << " with lots of printed details..." << std::endl;
            	}
            	std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
        	Mesh mesh;
		mesh.load(mesh_filename);
		mesh.set_attribute(unit_fct, 1, true); 
		
		//Initialisation de F et de K
		// K
		SparseMatrix K (mesh.nb_vertices());
		
		// F
		std::vector< double > F (mesh.nb_vertices(), 0.);   
		
		///Création F et Fe
		for (int i = 0; i < mesh.nb_triangles(); i++) {
			ElementMapping elt_mapping( mesh, false, i);
			ShapeFunctions reference_functions(2, 1);
			Quadrature quadrature = Quadrature::get_quadrature(2);
			// F
			std::vector <double> Fe_in;
			assemble_elementary_vector (elt_mapping, reference_functions, quadrature, sinus_bump_coefficient, Fe_in);
			local_to_global_vector(mesh, false, i, Fe_in, F );
			// K
			DenseMatrix Ke_in;
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, sinus_bump_coefficient, Ke_in);
			local_to_global_matrix(mesh, i, Ke_in, K);
			
	   	}       
		// Apply dirichlet condition
		// A MODIFIER POUR ADAPTER A NEUMANN NUL
        	std::vector< bool > attribute_is_dirichlet (2); 
        	attribute_is_dirichlet[0] = false; 
        	attribute_is_dirichlet[1] = true; 
        	
        	std::vector< double > values (mesh.nb_vertices());
        	for (int i = 0; i < mesh.nb_vertices() ; i++) {
        		values[i] = 0;
        	}
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Apply neumann condition
        	//NEUMANN NUL = YA RIEN A FAIRE, ne pas mettre quoi que ce soit pour neumann num
        	std::vector< bool > attribute_is_neumann
	}
    }

}










































