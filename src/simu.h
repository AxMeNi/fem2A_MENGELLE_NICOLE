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
		if (v.y >= 0.999999999) {return 1.;}
		else {return -1.;};
	}
	
	double region_bottom (vertex v)
	{
		if (v.y <= 0.000000001) {return 1.;}
		else {return -1.;};
	}
	
	double region_right (vertex v)
	{
		if (std::abs(v.x-1.) <= 0.000000001) {return 1.;}
		else {return -1.;};
	}
	
	double region_left (vertex v)
	{
		if (std::abs(v.x) <= 0.000000001) {return 1.;}
		else {return -1.;};
	}
	
	//Fonction de Neumann pour le pb du carre
	double neumann_fct_square (vertex v)
	{
		if (region_left(v) == 1)	//Neumann non nul pour le bord gauche
			{return sin(M_PI*v.y);} 
		else {return 0;} 		//Neumann nul pour le bord droit
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
	
	//Fonction de neumann pour le pb du mug
	
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
		DenseMatrix Ke;
		
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Creation de K et F
		for (int i = 0; i < mesh.nb_triangles(); i++) {
			ElementMapping elt_mapping( mesh, false, i);
			ShapeFunctions shp_fcts(2, 1);
			Quadrature quad = Quadrature::get_quadrature(2,false);
			// K
			
			assemble_elementary_matrix (elt_mapping, shp_fcts, quad, unit_fct, Ke);
			local_to_global_matrix(mesh, i, Ke, K);
	   	}       	
        	
        	// Apply dirichlet condition
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
        	std::string sol_path = "solutions/pure_dirich_" + solution_name;
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

		//Initialisation des shape functions
		ShapeFunctions shp_fcts(2, 1);
		
		//Initialisation de la quadrature
		Quadrature quad = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		DenseMatrix Ke;
		SparseMatrix K (mesh.nb_vertices());
		
		std::vector <double> Fe(3,0);
		std::vector< double > F (mesh.nb_vertices(), 0.);
		
		///Creation F Fe K Ke
		for (int i = 0; i < mesh.nb_triangles(); i++) 
		{
			ElementMapping elt_mapping( mesh, false, i);		
				
			// K
			assemble_elementary_matrix (elt_mapping, shp_fcts, quad, unit_fct, Ke);
			local_to_global_matrix(mesh, i, Ke, K);
			
			// F
			assemble_elementary_vector (elt_mapping, shp_fcts, quad, unit_fct, Fe);
			local_to_global_vector(mesh, false, i, Fe, F );
			
	   	}       	
	   	
	   	//Choix des attributs
		mesh.set_attribute(unit_fct, 1, true);      
        	
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
        	
        	//Creation du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/unit_src_trm_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        }
        
        
        
        
       	
        void dirichlet_sinus_bump_pb( const std::string& mesh_filename, bool verbose )
	{
            	std::cout << "Solving a Dirichlet problem with sinus source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
        	Mesh mesh;
		mesh.load(mesh_filename);   
		
		//Initialisation des shape functions
		ShapeFunctions shp_fcts(2, 1);
		
		//Initialisation de la quadrature
		Quadrature quad = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		DenseMatrix Ke;
		
		std::vector< double > F (mesh.nb_vertices(), 0.);
		std::vector <double> Fe (3,0);
		
		///Creation F et Fe
		for (int i = 0; i < mesh.nb_triangles(); i++) 
		{
			ElementMapping elt_mapping( mesh, false, i);		
				
			// K
			
			Ke.set_size(3,3);
			assemble_elementary_matrix (elt_mapping, shp_fcts, quad, unit_fct, Ke);
			local_to_global_matrix(mesh, i, Ke, K);
			
			// F
			
			assemble_elementary_vector (elt_mapping, shp_fcts, quad, sinus_bump_coefficient, Fe);
			local_to_global_vector(mesh, false, i, Fe, F );
			
	   	}     
	   	
	   	//Choix des attributs
		mesh.set_attribute(unit_fct, 1, true);     	
        	
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
        	
        	//Creation du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/sinus_src_trm_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        	
            		
        }
        
        void sinus_bump_pb_analytic(const std::string& mesh_filename, bool verbose) {
                std::cout << "Exact sinus source term solution" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
       		Mesh mesh;
            	mesh.load(mesh_filename);
            	
            	std::vector<double> x(mesh.nb_vertices(),0);
            	
            	for (int i=0 ; i<mesh.nb_vertices() ; ++i)
            	{
            	    x[i] = std::sin(M_PI*mesh.get_vertex(i).x)*std::sin(M_PI*mesh.get_vertex(i).y);
            	}
            	
                std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/analytic_sinus_src_trm_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
            }
        
        
        
        
        void diff_pb_sin ( const std::string& mesh_filename, bool verbose )
        {
        	std::cout << "Solving a Dirichlet problem with sinus source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
        	Mesh mesh;
		mesh.load(mesh_filename);   
		
		//Initialisation des shape functions
		ShapeFunctions shp_fcts(2, 1);
		
		//Initialisation de la quadrature
		Quadrature quad = Quadrature::get_quadrature(2, false);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		DenseMatrix Ke;
		
		std::vector< double > F (mesh.nb_vertices(), 0.);
		std::vector <double> Fe (3,0);
		
		for (int i = 0; i < mesh.nb_triangles(); i++) 
		{
			ElementMapping elt_mapping( mesh, false, i);		
				
			// K
			
			Ke.set_size(3,3);
			assemble_elementary_matrix (elt_mapping, shp_fcts, quad, unit_fct, Ke);
			local_to_global_matrix(mesh, i, Ke, K);
			
			// F
			
			assemble_elementary_vector (elt_mapping, shp_fcts, quad, sinus_bump_coefficient, Fe);
			local_to_global_vector(mesh, false, i, Fe, F );
			
	   	}     
	   	
	   	//Choix des attributs
		mesh.set_attribute(unit_fct, 1, true);     	
        	
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
        	
        	//La solution analytique
		std::vector<double> x_anal(mesh.nb_vertices(),0);
            	
            	for (int i=0 ; i<mesh.nb_vertices() ; ++i)
            	{
            	    	x_anal[i] = std::sin(M_PI * mesh.get_vertex(i).x) * std::sin(M_PI * mesh.get_vertex(i).y);
            	}
		
		//Calcul de la diffence
		std::vector<double> x_diff(mesh.nb_vertices(),0);
		
        	for (int i = 0 ; i < mesh.nb_vertices(); i++)
        	{
        		x_diff[i] = std::abs(x_anal[i] - x[i]);
        	}
        	
        	//Creation du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/diff_sinus_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x_diff, sol_path + "bb");
        	
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
		ShapeFunctions shp_fcts_2D(2, 1); // Pour les triangles
		ShapeFunctions shp_fcts_1D(1, 1); // Pour les edges
		
		//Initialisation de la quadrature
		Quadrature quad_2D = Quadrature::get_quadrature(2, false);
		Quadrature quad_1D = Quadrature::get_quadrature(2, true);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		DenseMatrix Ke;

		std::vector< double > F (mesh.nb_vertices(), 0.);
		std::vector <double> Fe_2D (3, 0.);
		std::vector <double> Fe_1D (2, 0.);
		
		for (int tri = 0; tri < mesh.nb_triangles(); tri++) 
		{
			ElementMapping elt_mapping_2D( mesh, false, tri);		
				
			// K
			assemble_elementary_matrix (elt_mapping_2D, shp_fcts_2D, quad_2D, unit_fct, Ke);
			local_to_global_matrix(mesh, tri, Ke, K);
			
			// F
			assemble_elementary_vector (elt_mapping_2D, shp_fcts_2D, quad_2D, unit_fct, Fe_2D);
			local_to_global_vector(mesh, false, tri, Fe_2D, F );
			
	   	}       	
        	
        	//Initialisation du vecteur de values pour Dirichlet
        	std::vector< double > values (mesh.nb_vertices(),0);
        	
        	//Choix des attributs
		mesh.set_attribute(region_right, 1, true);  	//Cond de Dirichlet
		mesh.set_attribute(region_left, 2, true); 	//Cond de Neumann - non nulle
		mesh.set_attribute(region_top, 2, true); 	//Cond de Neumann - nulle
		mesh.set_attribute(region_bottom, 2, true); 	//Cond de Neumann - nulle

        	//Creation des booléens qui indiquent quelle action sera a effectuer
        	std::vector< bool > attribute_is_dirichlet (3, false); 
        	attribute_is_dirichlet[1] = true;
        	
        	std::vector< bool > attribute_is_neumann (3, false);
        	attribute_is_neumann[2] = true;
		
		//Application des conditions de Dirichlet sur la frontière droite
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Application des conditions de Neumann sur les autres forntières
        	for (int edge_i = 0; edge_i < mesh.nb_edges(); ++ edge_i) //On itere sur tous les bords du maillage
        	{ 
        		
        		//On applique dabord la conditon de Neumann non nulle sur la frontière gauche
        		if (attribute_is_neumann[mesh.get_edge_attribute(edge_i)]) 
        		{
        			ElementMapping elt_mapping_1D(mesh, true, edge_i); // On travaille sur les borders avec Neumann
				assemble_elementary_neumann_vector(elt_mapping_1D, shp_fcts_1D, quad_1D, neumann_fct_square, Fe_1D);
				local_to_global_vector(mesh, true, edge_i, Fe_1D, F );
        		}
        	}
        	
        	//Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices(),0.);
        	FEM2A::solve(K,F,x);
        	
        	//Creation du fichier de solution
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
		
		//Initialisation des shape fcts
		ShapeFunctions shp_fcts_2D(2, 1); // Pour les triangles
		ShapeFunctions shp_fcts_1D(1, 1); // Pour les edges
		
		//Initialisation de la quadrature
		Quadrature quad_2D = Quadrature::get_quadrature(2, false);
		Quadrature quad_1D = Quadrature::get_quadrature(2, true);
		
		//Initiatlisation de F et de K
		SparseMatrix K (mesh.nb_vertices());
		DenseMatrix Ke;
		
		std::vector< double > F (mesh.nb_vertices(), 0.);
		std::vector <double> Fe_1D(2,0);
		
		for (int tri = 0; tri < mesh.nb_triangles(); tri++) 
		{
			ElementMapping elt_mapping_2D( mesh, false, tri);		
				
			// K
			
			Ke.set_size(3,3);
			assemble_elementary_matrix (elt_mapping_2D, shp_fcts_2D, quad_2D, unit_fct, Ke);
			local_to_global_matrix(mesh, tri, Ke, K);
			
			//On touche pas a F car pas de terme source ici (cf le premier pb)
			
	   	}       	
        	
        	//Initialisation du vecteur de values pour Dirichlet
        	std::vector< double > values (mesh.nb_vertices(),100.);
        	
        	//Choix des attributs
		mesh.set_attribute(region_hot_liquid, 1, true);  	//Cond de Dirichlet
		mesh.set_attribute(region_free_air, 2, true); 	//Cond de Neumann
        	
        	//Creation des booléens qui indiquent quelle action sera a effectuer
        	std::vector< bool > attribute_is_dirichlet (3, false); 
        	attribute_is_dirichlet[1] = true;
        	
        	std::vector< bool > attribute_is_neumann (3, false);
        	attribute_is_neumann[2] = true;
        	
		
		//Application des conditions de Dirichlet sur la frontiere droite
        	apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
        	
        	//Application des conditions de Neumann sur les autres forntieres
        	for (int edge_i = 0; edge_i < mesh.nb_edges(); ++ edge_i) //On itere sur tous les bords du maillage
        	{ 
        		ElementMapping elt_mapping_1D(mesh, true, edge_i); // On travaille sur les borders avec Neumann
        		
        		//On applique la conditon de Neumann 
        		if (attribute_is_neumann[mesh.get_edge_attribute(edge_i)]) 
        		{
				assemble_elementary_neumann_vector(elt_mapping_1D, shp_fcts_1D, quad_1D, constant_flux, Fe_1D);
				local_to_global_vector(mesh, true, edge_i, Fe_1D, F );
        		}
        	}
        	
        	//Recherche de la solution
        	std::vector<double> x (mesh.nb_vertices(),0);
        	FEM2A::solve(K,F,x);
        	
        	//Creation du fichier de solution
        	std::string solution_name;
        	solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end()-4);
        	std::string sol_path = "solutions/mug_pb_" + solution_name;
        	mesh.save(sol_path + "mesh");
        	save_solution( x, sol_path + "bb");
        	

        }
    }

}










































