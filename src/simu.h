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
	
	//Regions pour le pb du carré avec les conditions de Neumann
	double region_top (vertex v)
	{
		if (v.y == 1) {return 1.;}
		else {return -1.;};
	}
	
	double region_bottom (vertex v)
	{
		if (v.y == 0) {return 1.;}
		else {return -1.;};
	}
	
	double region_right (vertex v)
	{
		if (v.x == 1) {return 1.;}
		else {return -1.;};
	}
	
	double region_left (vertex v)
	{
		if (v.x == 0) {return 1.;}
		else {return -1.;};
	}
	
	//Fonction de Neumann pour le pb du carré
	double neumann_fct (vertex v)
	{
		double nf = sin(M_PI*v.y);
		return nf;
	}
	
	//Regions du mug
	double region_hot_liquid (vertex v) //Region chaude car en contact avec le liquide
	{
		if (((v.y == 1) && (v.x >= 1) && (v.x <= 20)) 	//Bas
		 ||((v.x == 1) && (v.y >= 1) && (v.y <= 10)) 	//Cote gauche
		 ||((v.x == 20) && (v.y >= 1) && (v.y <= 10))) 	//Cote droit
		 	{return 1.;}
		else {return -1.;}
		
	}
	
	double region_free_air (vertex v) //Region a lair libre
	{
		if (((v.y == 1) && (v.x >= 1) && (v.x <= 20)) 	//Bas
		 ||((v.x == 1) && (v.y >= 1) && (v.y <= 10)) 	//Cote gauche
		 ||((v.x == 20) && (v.y >= 1) && (v.y <= 10))) 	//Cote droit
		 	{return -1.;}
		else {return 1.;}
		
	}
	
	double constant_flux (vertex v) //Flux constant neumann pour le pb du mug
	{
		return -0.1;
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
        	
        	//Initialisation du vecteur de values pour Dirichlet
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
        	
        	//Initialisation du vecteur de values pour Dirichlet
        	std::vector< double > values (mesh.nb_vertices(),0.);
 
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices(),0.);
        	
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
        	
        	//Initialisation du vecteur de values pour Dirichlet
        	std::vector< double > values (mesh.nb_vertices(),0.);
 
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Recherche de la Solution
        	std::vector<double> x (mesh.nb_vertices(),0.);
        	
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
		
		//Initialisation des shape functions
		ShapeFunctions reference_functions(2, 1); // Pour les triangles
		ShapeFunctions reference_functions_1D(1, 1); // Pour les edges
		
		//Initialisation de la quadrature
		Quadrature quadrature = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		//Création F et Fe
		for (int tri = 0; tri < mesh.nb_triangles(); tri++) 
		{
			ElementMapping elt_mapping( mesh, false, tri);		
				
			// K
			DenseMatrix Ke_in;
			Ke_in.set_size(3,3);
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, unit_fct, Ke_in);
			local_to_global_matrix(mesh, tri, Ke_in, K);
			
			// F
			std::vector <double> Fe_in;
			assemble_elementary_vector (elt_mapping, reference_functions, quadrature, unit_fct, Fe_in);
			local_to_global_vector(mesh, false, tri, Fe_in, F );
			
	   	}       	
        	
        	//Initialisation du vecteur de values pour Dirichlet
        	std::vector< double > values (mesh.nb_vertices(),0);
        	
        	//Choix des attributs
		mesh.set_attribute(region_hot_liquid, 0, true);  	//Cond de Dirichlet
		mesh.set_attribute(region_free_air, 1, true); 	//Cond de Neumann
        	
        	//Création des booléens qui indiquent quelle action sera à effectuer
        	std::vector< bool > attribute_is_dirichlet (3); 
        	attribute_is_dirichlet[0] = true; 
        	attribute_is_dirichlet[1] = false;
        	attribute_is_dirichlet[2] = false;
        	
        	std::vector< bool > attribute_is_neumann (3);
        	attribute_is_neumann[0] = false; 
        	attribute_is_neumann[1] = true;
        	attribute_is_neumann[2] = false;
        	
        	std::vector< bool > attribute_is_neumann_null (3);
        	attribute_is_neumann_null[0] = false; 
        	attribute_is_neumann_null[1] = false;
        	attribute_is_neumann_null[2] = true;
        	
		
		//Application des conditions de Dirichlet sur la frontière droite
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Application des conditions de Neumann sur les autres forntières
        	for (int edge_i = 0; edge_i < mesh.nb_edges(); ++ edge_i) //On itère sur tous les bords du maillage
        	{ 
        		ElementMapping elt_mapping_1D(mesh, true, edge_i); // On travaille sur les borders avec Neumann
        		
        		//On applique d'abord la conditon de Neumann non nulle sur la frontière gauche
        		if (attribute_is_neumann[mesh.get_edge_attribute(edge_i)]) 
        		{
        			std::vector <double> Fe_in;
				assemble_elementary_neumann_vector(elt_mapping_1D, reference_functions_1D, quadrature, neumann_fct, Fe_in);
				local_to_global_vector(mesh, true, edge_i, Fe_in, F );
        		}
        		
        		//On applique ensuite la condition de Neumann nulle sur les frontières en haut et en bas
        		if (attribute_is_neumann_null[mesh.get_edge_attribute(edge_i)])
        		{
        			std::vector <double> Fe_in;
        			assemble_elementary_neumann_vector(elt_mapping_1D, reference_functions_1D, quadrature, zero_fct, Fe_in);
				local_to_global_vector(mesh, true, edge_i, Fe_in, F );
        		}
        	
        	}
        	
        	//Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices(),0.);
        	FEM2A::solve(K,F,x);
        	
        	//Création du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/neumann_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        	

        }
        
        
        
        
        void mug_pb( const std::string& mesh_filename, bool verbose )
	{
            	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
        	Mesh mesh;
		mesh.load(mesh_filename);
		
		//Initialisation des shape functions
		ShapeFunctions reference_functions(2, 1); // Pour les triangles
		ShapeFunctions reference_functions_1D(1, 1); // Pour les edges
		
		//Initialisation de la quadrature
		Quadrature quadrature = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		//Création F et Fe
		for (int tri = 0; tri < mesh.nb_triangles(); tri++) 
		{
			ElementMapping elt_mapping( mesh, false, tri);		
				
			// K
			DenseMatrix Ke_in;
			Ke_in.set_size(3,3);
			assemble_elementary_matrix (elt_mapping, reference_functions, quadrature, unit_fct, Ke_in);
			local_to_global_matrix(mesh, tri, Ke_in, K);
			
			//On touche pas a F car pas de terme source ici (cf le premier pb)
			
	   	}       	
        	
        	//Initialisation du vecteur de values pour Dirichlet
        	std::vector< double > values (mesh.nb_vertices(),100.);
        	
        	//Choix des attributs
		mesh.set_attribute(region_hot_liquid, 0, true);  	//Cond de Dirichlet
		mesh.set_attribute(region_free_air, 1, true); 	//Cond de Neumann
        	
        	//Création des booléens qui indiquent quelle action sera à effectuer
        	std::vector< bool > attribute_is_dirichlet (2); 
        	attribute_is_dirichlet[0] = true; 
        	attribute_is_dirichlet[1] = false;
        	
        	std::vector< bool > attribute_is_neumann (2);
        	attribute_is_neumann[0] = false; 
        	attribute_is_neumann[1] = true;
		
		//Application des conditions de Dirichlet sur la frontière droite
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Application des conditions de Neumann sur les autres forntières
        	for (int edge_i = 0; edge_i < mesh.nb_edges(); ++ edge_i) //On itère sur tous les bords du maillage
        	{ 
        		ElementMapping elt_mapping_1D(mesh, true, edge_i); // On travaille sur les borders avec Neumann
        		
        		//On applique d'abord la conditon de Neumann non nulle sur la frontière gauche
        		if (attribute_is_neumann[mesh.get_edge_attribute(edge_i)]) 
        		{
        			std::vector <double> Fe_in;
				assemble_elementary_neumann_vector(elt_mapping_1D, reference_functions_1D, quadrature, constant_flux, Fe_in);
				local_to_global_vector(mesh, true, edge_i, Fe_in, F );
        		}
        	
        	}
        	
        	//Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices(),0);
        	FEM2A::solve(K,F,x);
        	
        	//Création du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/mug_pb_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        	

        }
    }

}










































