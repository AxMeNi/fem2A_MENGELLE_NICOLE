#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }


        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        
        bool test_quadrature(int order, bool border = false) 
        {
		std::cout << "[TEST DE QUADRATURE] \n";
        	Quadrature quad;
        	quad = quad.get_quadrature(order, border); // Je crée une quadrature
        	std::cout << "Valeur de la somme des poids : ";
        	double sum_w = 0;
        	for (int i = 0; i < quad.nb_points(); i++) { // Les i sont les index des points quad.point(i).x et quad.point(i).y sont les coordonnées du point i
        		sum_w += quad.weight(i); // Je somme les poids de la quadrature
        	}
        	std::cout << sum_w;
        	if (sum_w == 0.5) { //Si la somme vaut 0.5 ça fonctionne
        		return true;
        	}	
        	return false;
        }
        
        
        bool test_ElementMapping() {
		std::cout << "[TEST DE ELEMENT MAPPING] \n";
        	Mesh mesh;
        	mesh.load("data/square.mesh");
        	// Test constructeur
		std::cout << ">> TEST DU COSNTRUCTEUR \n";
        	ElementMapping el_triangle( mesh, false, 4 );
        	ElementMapping el_border( mesh, true, 4 );
        	for (int i = 0; i < 3 ;++i) {
        		std::cout << "v" << i << " : ";
        		std::cout << el_triangle.get_vertices()[i].x ;
        		std::cout << ", ";
        		std::cout << el_triangle.get_vertices()[i].y ;
        		std::cout << "\n";
        	}
        	for (int i = 0; i < 2 ;++i) {
        		std::cout << "v" << i << " : ";
        		std::cout << el_border.get_vertices()[i].x ;
        		std::cout << ", ";
        		std::cout << el_border.get_vertices()[i].y ;
        		std::cout << "\n";
        	}
 
        	// Test transform
		std::cout << ">> TEST DE TRANSFORM \n";
        	vertex vert1;
        	vert1.x = 0.2;
        	vert1.y = 0.4;
        	vertex out = el_triangle.transform(vert1);
        	std::cout << "Dans l'espace global >> ";
        	std::cout << "x : " << out.x << " y :" << out.y ;
        	std::cout << "\n";
        	
        	// Test jacobian matrix
		std::cout << ">> TEST DE JACOBIAN MATRIX \n";
        	DenseMatrix jacob_matrix = el_triangle.jacobian_matrix(vert1);
        	jacob_matrix.print();
        	
        	// Test det of jacobian matrix
		std::cout << ">> TEST DE JACOBIAN \n";
        	double det1 = el_triangle.jacobian(vert1);
        	std::cout << "determinant with the triangle :" << det1 << "\n";
        	double det2 = el_border.jacobian(vert1);
        	std::cout << "determinant with the border :" << det2 << "\n";
        	return true;
        }
        
        
        bool test_ShapeFunctions() {
		std::cout << "[TEST DE SHAPE FUNCTIONS] \n";
        	vertex vert1;
        	vert1.x = 0.2;
        	vert1.y = 0.4;
        	ShapeFunctions shpf_ref_tri (2, 1); // dim, 
        	ShapeFunctions shpf_ref_seg (1, 1); // dim
        	
        	int nb_functions_test =  shpf_ref_tri.nb_functions();
        	int nb_functions_test2 =  shpf_ref_seg.nb_functions();
        	std::cout << "nb_functions = " << nb_functions_test << "\n";
        	std::cout << "nb_functions = " << nb_functions_test2 << "\n";
        	
        	double evaluation_tri = shpf_ref_tri.evaluate(0, vert1);
        	double evaluation_seg = shpf_ref_seg.evaluate(0, vert1);
        	std::cout << "evaluation_tri = " << evaluation_tri << "\n";
        	std::cout << "evaluation_seg = " << evaluation_seg << "\n";
        	
        	vec2 evaluation_grad_tri = shpf_ref_tri.evaluate_grad(0, vert1);
        	vec2 evaluation_grad_seg = shpf_ref_seg.evaluate_grad(1, vert1);
        	std::cout << "evaluation_grad_tri = " << evaluation_grad_tri.x << "," << evaluation_grad_tri.y << "\n";
        	std::cout << "evaluation_grad_seg = " << evaluation_grad_seg.x << "," << evaluation_grad_tri.y << "\n";
        	return true;
        }
        
        
        // Définition de la fonction de pointeur renvoyant toujours 1
	double constant_coefficient(FEM2A::vertex vertex) {
    		return 1.0; // Toujours renvoyer 1 pour n'importe quelle valeur de vertex
	}
	
	
	
        bool test_assemble_elementary_matrix() {
		std::cout << "[TEST DE ASSEMBLE ELEMENTARY MATRIX] \n";
		Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping elt_mapping(mesh, false, 4);
		
		ShapeFunctions reference_functions(2, 1);
		
		Quadrature quad;
		Quadrature quadrature = quad.get_quadrature(2, false);
		
		DenseMatrix Ke;
		
		// Assemblage de la matrice élémentaire en utilisant la fonction de pointeur constant_coefficient
		assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, constant_coefficient, Ke);
		Ke.print()
		return true;
	
        }
        
        
        bool test_local_to_global_matrix() {
		std::cout << "[TEST DE LOCAL TO GLOBAL MATRIX] \n";
		Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping elt_mapping(mesh, false, 4);
		    
		ShapeFunctions reference_functions(2, 1);
		    
		Quadrature quad;
		Quadrature quadrature = quad.get_quadrature(2, false);
		    
		DenseMatrix Ke;
		    
	   	assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, constant_coefficient, Ke);
	   	
	   	SparseMatrix K (mesh.nb_vertices());
        	local_to_global_matrix(mesh, 4, Ke, K);
        	K.print()
        	return true;
        }
        
        
        bool test_apply_dirichlet_boundary_conditions()
        {
		std::cout << "[TEST DE APPLY DIRICHLET BOUNDARY CONDITIONS] \n";
		std::cout << "Remarque : ce test affiche x, le vecteur des valeurs finales après l'application de la conditon de Dirichlet. Pour effectuer un simulation et produire un fichier de simulation, passer par la commande -m.\n";
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		mesh.set_attribute( constant_coefficient, 0, true );
		
		ShapeFunctions reference_functions(2,1);
		Quadrature quadrat = Quadrature::get_quadrature(0,false);
		DenseMatrix Ke ;
		SparseMatrix K(mesh.nb_vertices() * quadrat.nb_points());
		
		
		for (int i=0 ; i<mesh.nb_triangles(); ++i)
		{
		ElementMapping element(mesh, false, i);
		assemble_elementary_matrix(element, reference_functions, quadrat, constant_coefficient, Ke);
		local_to_global_matrix(mesh, i, Ke, K);
		}
		
		std::vector< double > F(mesh.nb_vertices(), 1);
		
		std::vector< bool > attribute_bool(1, true);
		std::vector< double > values(mesh.nb_vertices());
		
		for (int i =0 ; i<mesh.nb_vertices(); ++i)
		{
		values[i] = mesh.get_vertex(i).x + mesh.get_vertex(i).y; 
		}
		
		apply_dirichlet_boundary_conditions(mesh, attribute_bool, values, K, F);
		
		std::vector< double > x(mesh.nb_vertices(), 0);
		solve(K,F, x);
		
		for (double xi :x)
		{
		std::cout << xi << '\t';
		}
		
		return true;
        }
        
        bool test_assemble_elementary_vector() 
        {
		std::cout << "[TEST DE ASSEMBLE ELEMENTARY MATRIX] \n";
		Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping elt_mapping(mesh, false, 4);
		
		ShapeFunctions reference_functions(2, 1);
		
		Quadrature quad;
		Quadrature quadrature = quad.get_quadrature(2, false);
		
		std::vector <double> Fe (3,0);
		
		// Assemblage de la matrice élémentaire en utilisant la fonction de pointeur constant_coefficient
		assemble_elementary_vector(elt_mapping, reference_functions, quadrature, constant_coefficient, Fe ) ;
		for (int i = 0 ; i < 3 ; i++){
		
		std::cout << Fe[i] << "\n";
		}
		
		return true;
	
        }
        
        bool test_local_to_global_vector() 
        {
		std::cout << "[TEST DE LOCAL TO GLOBAL VECTOR] \n";
        	Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping elt_mapping(mesh, false, 4);
		    
		ShapeFunctions reference_functions(2, 1);
		    
		Quadrature quad;
		Quadrature quadrature = quad.get_quadrature(2, false);
		    
		std::vector <double> Fe;
		    
	   	assemble_elementary_vector(elt_mapping, reference_functions, quadrature, constant_coefficient, Fe);
	   	
	   	std::vector <double> F (mesh.nb_vertices());
        	local_to_global_vector(mesh, false, 4, Fe, F );
        	for (int i = 0; i < mesh.nb_vertices() ; ++i) {
        		std::cout << F[i] << "\t"; 
        	}
        	for (int i = 0; i < 3 ; ++i) {
        		std::cout << mesh.get_triangle_vertex_index(4,i) << "\n"; 
        	}
        		
        	return true;
        }
        
        bool test_assemble_elementary_neumann_vector() 
        {
		std::cout << "[TEST DE ASSEMBLE ELEMENTARY NEUMANN VECTOR] \n";
       		Mesh mesh;
        	mesh.load("data/square.mesh");
       		ElementMapping elt_mapping(mesh,true,4);

       		ShapeFunctions reference_functions(1,1);

            	Quadrature quad;
            	Quadrature quadrature = quad.get_quadrature(2,true);

            	std::vector<double> F(mesh.nb_vertices(),0);
            	std::vector<double> Fe(2,0);

		assemble_elementary_neumann_vector(elt_mapping,reference_functions,quadrature,constant_coefficient,Fe);
            	local_to_global_vector(mesh,true,4,Fe,F);
            	
            	for (double i : Fe)
            	{
            		std::cout << i << "\t";
            	}
        	
        	return true;
        }
        
        
        ////////////////////////////////////////////////////////////////////////////
        ////// ALL THE FUNCTIONS ABOVE ARE FOR SOLVING THE POISSON PROBLEM//////////
        ////////////////////////////////////////////////////////////////////////////
        
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
	
	double unit_fct( vertex v )
	{
    		return 1.;
	}
	
	double zero_fct( vertex v )
	{
    		return 0.;
	}
	
	double neumann_fct_square (vertex v)
	{
		if (region_left(v) == 1)	//Neumann non nul pour le bord gauche
			{return sin(M_PI*v.y);} 
		else {return 0;} 		//Neumann nul pour le bord droit
	}
        
        bool test_poisson_problem( const std::string& mesh_filename)
        {
		std::cout << "[TEST DE POISSON PROBELM] \n";
		Mesh M;
		M.load(mesh_filename);
		
		//Definition des attributs
		M.set_attribute(region_right, 1, true);  	//Cond de Dirichlet
		M.set_attribute(region_left, 2, true); 		//Cond de Neumann - non nulle
		M.set_attribute(region_top, 2, true); 		//Cond de Neumann - nulle
		M.set_attribute(region_bottom, 2, true); 	//Cond de Neumann - nulle
		
		//Recherche de la solution
        	std::vector<double> solution (M.nb_vertices(),0);
        	
		solve_poisson_problem(M, unit_fct, unit_fct, zero_fct, neumann_fct_square, solution, false);
            
		//Creation du fichier de solution
		std::string solution_name;
		solution_name.assign(mesh_filename.begin() + 5, mesh_filename.end() - 4);
		std::string sol_path = "solutions/poisson_" + solution_name;
		M.save(sol_path + "mesh");
		save_solution(solution, sol_path + "bb");
		std::cout << "Successfully saved the poisson file in data/output/ \n" ;
		return true;
        }
        

    }
}














