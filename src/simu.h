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
    		return coeff; 
	}
	
	double region_top (vertex v)
	{
		if (v.y == 1) {return 1;}
		else {return 0;};
	}
	
	double region_bottom (vertex v)
	{
		if (v.y = 0) {return 1;}
		else {return 0;};
	}
	
	double region_right (vertex v)
	{
		if (v.x == 1) {return 1;}
		else {return 0;};
	}
	
	double region_left (vertex v)
	{
		if (v.x == 0) {return 1;}
		else {return 0;};
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
        	Mesh mesh;
		mesh.load(mesh_filename);
		
		//Choix des attributs
		mesh.set_attribute(unit_fct, 1, true);      
		
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		DenseMatrix Ke_in;
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Création de K et F
		for (int i = 0; i < mesh.nb_triangles(); i++) {
			ElementMapping elt_mapping( mesh, false, i);
			ShapeFunctions reference_functions(2, 1);
			Quadrature quadrature = Quadrature::get_quadrature(2,false);
			// K
			
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
        	
        	// Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices());
        	FEM2A::solve(K,F,x);
        	
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/pure_dirich" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        }
        
        
        
        
        
	
	void dirichlet_pb_source_term( const std::string& mesh_filename, bool verbose )
	{
            	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
        	Mesh mesh;
		mesh.load(mesh_filename);
		
		//Choix des attributs
		mesh.set_attribute(unit_fct, 1, true);      
		
		//Initialisation des shape functions
		ShapeFunctions reference_functions(2, 1);
		
		//Initialisation de la quadrature
		Quadrature quadrature = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Création F Fe K Ke
		for (int i = 0; i < mesh.nb_triangles(); i++) 
		{
			ElementMapping elt_mapping( mesh, false, i);		
				
			// K
			DenseMatrix Ke_in;
			Ke_in.set_size(3,3);
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, unit_fct, Ke_in);
			local_to_global_matrix(mesh, i, Ke_in, K);
			
			// F
			std::vector <double> Fe_in;
			assemble_elementary_vector (elt_mapping, reference_functions, quadrature, unit_fct, Fe_in);
			local_to_global_vector(mesh, false, i, Fe_in, F );
			
	   	}       	
        	
        	//Application des conditions de Dirichlet
        	std::vector< bool > attribute_is_dirichlet (2); 
        	attribute_is_dirichlet[0] = false; 
        	attribute_is_dirichlet[1] = true; 
        	
        	std::vector< double > values (mesh.nb_vertices(),0);
 
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices(),0);
        	
        	FEM2A::solve(K,F,x);
        	
        	//Création du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/dirich_source_term_unit_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        }
        
        
        
        void dirichlet_sinus_bump_pb( const std::string& mesh_filename, bool verbose )
	{
            	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
        	Mesh mesh;
		mesh.load(mesh_filename);
		
		//Choix des attributs
		mesh.set_attribute(unit_fct, 1, true);      
		
		//Initialisation des shape functions
		ShapeFunctions reference_functions(2, 1);
		
		//Initialisation de la quadrature
		Quadrature quadrature = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Création F et Fe
		for (int i = 0; i < mesh.nb_triangles(); i++) 
		{
			ElementMapping elt_mapping( mesh, false, i);		
				
			// K
			DenseMatrix Ke_in;
			Ke_in.set_size(3,3);
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, unit_fct, Ke_in);
			local_to_global_matrix(mesh, i, Ke_in, K);
			
			// F
			std::vector <double> Fe_in;
			assemble_elementary_vector (elt_mapping, reference_functions, quadrature, sinus_bump_coefficient, Fe_in);
			local_to_global_vector(mesh, false, i, Fe_in, F );
			
	   	}       	
        	
        	//Application des conditions de Dirichlet
        	std::vector< bool > attribute_is_dirichlet (2); 
        	attribute_is_dirichlet[0] = false; 
        	attribute_is_dirichlet[1] = true; 
        	
        	std::vector< double > values (mesh.nb_vertices(),0);
 
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Recherche de la Solution
        	std::vector<double> x (mesh.nb_vertices(),0);
        	
        	FEM2A::solve(K,F,x);
        	
        	//Création du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/dirich_source_term_sinus_bump_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        	
        	//for(int vertice ; vertice<mesh.nb_vertices();++vertice)
            	//{
                //	x[i] -= sin(M_PI*mesh.get_vertex(vertice).x)*sin(M_PI*mesh.get_vertex(vertice).y);
            	//}
            	//save_solution( x, sol_path + "bb");
            		
        }
        
        
        void neumann_pb( const std::string& mesh_filename, bool verbose )
	{
            	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
        	Mesh mesh;
		mesh.load(mesh_filename);
		
		//Choix des attributs
		mesh.set_attribute(region_right, 1, true);  	//Cond de Dirichlet
		mesh.set_attribute(region_left, 2, true); 	//Cond de Neumann
		mesh.set_attribute(region_top, 3, true);	//Cond de Neumann nulle
		mesh.set_attribute(region_bottom, 3, true);	//Cond de Neumann nulle
		
		//Initialisation des shape functions
		ShapeFunctions reference_functions(2, 1);
		
		//Initialisation de la quadrature
		Quadrature quadrature = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Création F et Fe
		for (int i = 0; i < mesh.nb_triangles(); i++) 
		{
			ElementMapping elt_mapping( mesh, false, i);		
				
			// K
			DenseMatrix Ke_in;
			Ke_in.set_size(3,3);
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, unit_fct, Ke_in);
			local_to_global_matrix(mesh, i, Ke_in, K);
			
			// F
			std::vector <double> Fe_in;
			assemble_elementary_vector (elt_mapping, reference_functions, quadrature, sinus_bump_coefficient, Fe_in);
			local_to_global_vector(mesh, false, i, Fe_in, F );
			
	   	}       	
        	
        	//Application des conditions de Dirichlet
        	std::vector< bool > attribute_is_dirichlet (2); 
        	attribute_is_dirichlet[0] = false; 
        	attribute_is_dirichlet[1] = true;
        	attribute_is_dirichlet[2] = false;
        	attribute_is_dirichlet[3] = false;
        	
        	std::vector< double > values (mesh.nb_vertices(),0);

        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Application des conditions de Neumann
        	
        	
        	//Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices(),0);
        	
        	FEM2A::solve(K,F,x);
        	
        	//Création du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/neumann_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        	
        	//for(int vertice ; vertice<mesh.nb_vertices();++vertice)
            	//{
                //	x[i] -= sin(M_PI*mesh.get_vertex(vertice).x)*sin(M_PI*mesh.get_vertex(vertice).y);
            	//}
            	//save_solution( x, sol_path + "bb");
            		
        }
        
        
        
        
        
    }

}










































