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
        	Mesh mesh;
        	mesh.load("data/square.mesh");
        	// Test constructeur
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
        	vertex vert1;
        	vert1.x = 0.2;
        	vert1.y = 0.4;
        	vertex out = el_triangle.transform(vert1);
        	std::cout << "Dans l'espace global >> ";
        	std::cout << "x : " << out.x << " y :" << out.y ;
        	std::cout << "\n";
        	
        	// Test jacobian matrix
        	DenseMatrix jacob_matrix = el_triangle.jacobian_matrix(vert1);
        	jacob_matrix.print();
        	
        	// Test det of jacobian matrix
        	double det1 = el_triangle.jacobian(vert1);
        	std::cout << "determinant with the triangle :" << det1 << "\n";
        	double det2 = el_border.jacobian(vert1);
        	std::cout << "determinant with the border :" << det2 << "\n";
        	return true;
        }
        
        
        bool test_ShapeFunctions() {
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
	    Mesh mesh;
	    mesh.load("data/square.mesh");
	    ElementMapping elt_mapping(mesh, false, 4);
	    
	    ShapeFunctions reference_functions(2, 1);
	    
	    Quadrature quad;
	    Quadrature quadrature = quad.get_quadrature(2, false);
	    
	    DenseMatrix Ke_in;
	    
	    // Assemblage de la matrice élémentaire en utilisant la fonction de pointeur constant_coefficient
	    assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, constant_coefficient, Ke_in);
	    
	    return true;
	
        }
        
        
        bool test_local_to_global_matrix() {
		Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping elt_mapping(mesh, false, 4);
		    
		ShapeFunctions reference_functions(2, 1);
		    
		Quadrature quad;
		Quadrature quadrature = quad.get_quadrature(2, false);
		    
		DenseMatrix Ke_in;
		    
	   	assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, constant_coefficient, Ke_in);
	   	
	   	SparseMatrix K (mesh.nb_vertices());
        	local_to_global_matrix(mesh, 4, Ke_in, K);
        	
        	return true;
        }
        
        
        bool test_apply_dirichlet_boundary_conditions()
        {
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
        
        bool test_assemble_elementary_vector() {
       	    Mesh mesh;
	    mesh.load("data/square.mesh");
	    ElementMapping elt_mapping(mesh, false, 4);
	    
	    ShapeFunctions reference_functions(2, 1);
	    
	    Quadrature quad;
	    Quadrature quadrature = quad.get_quadrature(2, false);
	    
	    std::vector <double> Fe_in;
	    
	    // Assemblage de la matrice élémentaire en utilisant la fonction de pointeur constant_coefficient
	    assemble_elementary_vector(elt_mapping, reference_functions, quadrature, constant_coefficient, Fe_in ) ;
	    for (int i = 0 ; i < 3 ; i++){
	    	
	    	std::cout << Fe_in[i] << "\n";
	    }
	    
	    return true;
	
        }
        
        bool test_local_to_global_vector() {
        	Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping elt_mapping(mesh, false, 4);
		    
		ShapeFunctions reference_functions(2, 1);
		    
		Quadrature quad;
		Quadrature quadrature = quad.get_quadrature(2, false);
		    
		std::vector <double> Fe_in;
		    
	   	assemble_elementary_vector(elt_mapping, reference_functions, quadrature, constant_coefficient, Fe_in);
	   	
	   	std::vector <double> F (mesh.nb_vertices());
        	local_to_global_vector(mesh, false, 4, Fe_in, F );
        	for (int i = 0; i < mesh.nb_vertices() ; ++i) {
        		std::cout << F[i] << "\t"; 
        	}
        	for (int i = 0; i < 3 ; ++i) {
        		std::cout << mesh.get_triangle_vertex_index(4,i) << "\n"; 
        	}
        		
        	return true;
        }
        
        bool test_assemble_elementary_neumann_vector() {
        	// A IMPLEMENTER
        	return true;
        }
    }
}














