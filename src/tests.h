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
        	return true;
        }

    }
}
