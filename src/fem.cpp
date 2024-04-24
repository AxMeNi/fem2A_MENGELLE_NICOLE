#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        //std::cout << "[ElementMapping] constructor for element " << i << " ";
        //if ( border ) std::cout << "(border)";
        std::cout << '\n';
        std::vector< vertex > vertices ;
        //Objectif : récupérer les vertices de l'élément de mapping  et les stocker dans l'attribut vertex
        if ( border ) { // Si l'élément est un edge
        	for (int index_local_vertex = 0; index_local_vertex < 2 ; index_local_vertex++) {
        	vertices.push_back(M.get_edge_vertex(i, index_local_vertex));
        	}
        }
        else { //Si l'élément est un triangle
        	for (int index_local_vertex = 0; index_local_vertex < 3 ; index_local_vertex++) {
        	vertices.push_back(M.get_triangle_vertex(i, index_local_vertex));
        	}
        }
        
        vertices_ = vertices;
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] transform reference to world space" << '\n';
        //TODO
        vertex r ;
        if (border_) {
        	r.x = vertices_[0].x * (1 - x_r.x) + vertices_[1].x * x_r.x;
        	r.y = vertices_[0].y * (1 - x_r.x) + vertices_[1].y * x_r.x;
        }
        else {
        	r.x = (1 - x_r.x - x_r.y) * vertices_[0].x + vertices_[1].x * x_r.x + vertices_[2].x * x_r.y;
        	r.y = (1 - x_r.x - x_r.y) * vertices_[0].y + vertices_[1].y * x_r.x + vertices_[2].y * x_r.y;
        }
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        // TODO
        DenseMatrix J ;
        if (border_) {
        	J.set_size(2, 1) ;
        	J.set(0, 0, -1 * vertices_[0].x + 1 * vertices_[1].x);
        	J.set(1, 0, -1 * vertices_[0].y + 1 * vertices_[1].y);
        }
        else {
        	J.set_size(2, 2) ;
        	J.set(0, 0, -1 * vertices_[0].x + 1 * vertices_[1].x);
        	J.set(0, 1, -1 * vertices_[0].x + 1 * vertices_[2].x);
        	J.set(1, 0, -1 * vertices_[0].y + 1 * vertices_[1].y);
        	J.set(1, 1, -1 * vertices_[0].y + 1 * vertices_[2].y);
         }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        // TODO
        double det;
        DenseMatrix jacob_matrix = jacobian_matrix( x_r );
        if (border_) {
        	det = pow( (pow( jacob_matrix.get(0,0), 2.) + pow (jacob_matrix.get(1,0), 2.)),0.5);
        }
        else {
        	det = jacob_matrix.det_2x2();
        }
        return det ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        //std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
        // TODO
        assert ((dim_ == 1) || (dim == 2));
        assert (order == 1);
    }

    int ShapeFunctions::nb_functions() const
    {
        //std::cout << "[ShapeFunctions] number of functions" << '\n';
        if ( dim_ == 1 ) {
        	return 2;
        }
        else {
        	return 3;        
        }
        return 0 ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
       	if (dim_ == 1) {
       		if (i == 0) {
       			double value = 1 - x_r.x;
       			return value;
       		}
       		else if (i == 1) {
       			return x_r.x;
       		}
       	}
       	else if (dim_ == 2) {
       		if (i == 0) {
       			double value = 1 - x_r.x - x_r.y;
       			return value;
       		}
       		else if (i == 1) {
       			return x_r.x;
       		}
       		else if (i == 2) {
       			return x_r.y;
       		}
       	}
        return 0. ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        // TODO
        vec2 g ;
        if (dim_ == 1) {
       		if (i == 0) {
       			g.x = -1;
       			g.y = 0;
       		}
       		else if (i == 1) {
       			g.x = 1;
       			g.y = 0;
       		}
       	}
       	else if (dim_ == 2) {
       		if (i == 0) {
       			g.x = -1;
       			g.y = -1;
       		}
       		else if (i == 1) {
       			g.x = 1;
       			g.y = 0;
       		}
       		else if (i == 2) {
       			g.x = 0;
       			g.y = 1;
       		}
       	}
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
    	double wq;
    	vertex Me;
    	double coeff;
    	DenseMatrix Je;
    	double absJac;
    	vec2 grad_phi_i;
    	vec2 grad_phi_j;
    	double Ke_ij = 0;
    	DenseMatrix trans_inv_Je;
    	trans_inv_Je.set_size(2,2);
    	int max_i = reference_functions.nb_functions();
    	int max_j = reference_functions.nb_functions();
    	vec2 matricegauche;
    	vec2 matricedroite;
    	Ke.set_size(max_i, max_j);
    	
     	//std::cout << "compute elementary matrix" << '\n';
     	for (int i = 0; i < max_i; i++) {
     	
     		for (int j = 0; j < max_j; j++) {
     			
			for (int q = 0; q < quadrature.nb_points(); q++) {
				
				wq = quadrature.weight(q);
				Me = elt_mapping.transform(quadrature.point(q));
				coeff = coefficient(Me);
				Je = elt_mapping.jacobian_matrix(quadrature.point(q));
				trans_inv_Je = Je.transpose().invert_2x2();
				absJac = abs(elt_mapping.jacobian(quadrature.point(q)));
				grad_phi_i = reference_functions.evaluate_grad(i, quadrature.point(q));
				grad_phi_j = reference_functions.evaluate_grad(j, quadrature.point(q));
				//On multiplie les membre de droite et de gauche par les coefficients
				grad_phi_i.x *= wq * coeff;
				grad_phi_i.y *= wq * coeff;
				grad_phi_j.x *= absJac;
				grad_phi_j.y *= absJac;
				// Traitement de la matrice de gauche
				matricegauche = trans_inv_Je.mult_2x2_2(grad_phi_i);
				// Traitement de la matrice de droite
				matricedroite = trans_inv_Je.mult_2x2_2(grad_phi_j);
				Ke_ij += dot(matricegauche, matricedroite);
			}
			Ke.set(i, j, Ke_ij);
			Ke_ij = 0;
		}
	
        } 
    	// Ke.print();     
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        //std::cout << "Ke -> K" << '\n';
        //TODO
        // Ke et K sont symétriques
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			int glob_id1 = M.get_triangle_vertex_index(t, i);
			int glob_id2 = M.get_triangle_vertex_index(t, j);
			K.add(glob_id1, glob_id2, Ke.get(i,j));
		}
	}
	// K.print();
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
        int max = reference_functions.nb_functions();
        
        for (int i = 0; i < max; ++i){
            double Fe_i = 0.;
            for (int q = 0; q < quadrature.nb_points(); ++q){
                vertex pt_dintegration = quadrature.point(q);
                double w = quadrature.weight(q);
                double shape_i = reference_functions.evaluate(i, pt_dintegration);
                double det = elt_mapping.jacobian(pt_dintegration);
               
                Fe_i += w*shape_i*source(pt_dintegration)*det;
            }
            Fe.push_back(Fe_i);
        }
    }


    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        int max = reference_functions_1D.nb_functions();
        
        for (int i = 0; i < max; ++i){
            double Fe_i = 0.;
            for (int q = 0; q < quadrature_1D.nb_points(); ++q){
                vertex pt_dintegration = quadrature_1D.point(q);
                double w = quadrature_1D.weight(q);
                double shape_i = reference_functions.evaluate(i, pt_dintegration);
                double det = elt_mapping_1D.jacobian(pt_dintegration);
               
                Fe_i += w*shape_i*source(pt_dintegration)*det;
            }
            Fe.push_back(Fe_i);
        }
    }
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::cout << "Fe -> F" << '\n';
        if (border ) {
        	for (int j = 0; j < Fe.size(); ++j) {
        		F[M.get_edge_vertex_index(i, j)] += Fe[j];
        	}
        }
        else {
        	for (int j = 0; j < Fe.size(); ++j) {
        		F[M.get_triangle_vertex_index(i, j)] += Fe[j];
        	}
        }
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: must be M.nb_vertices() */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';
        // TODO prendre p =10 000
        std::vector< bool > processed_vertices(values.size(), false);
        double penalty_coefficient = 10000.;
        for ( int edge = 0; edge < M.nb_edges(); edge++ ) {
        	int edge_attribute = M.get_edge_attribute(edge);
        	if ( attribute_is_dirichlet[edge_attribute] ) {
        	// Si edge_attribute vaut 0 alors c'est false et si edge_attribute vaut 1 alors c'est true avec la fonction unit_fct
        		for (int v = 0; v < 2; v++) {
        			int vertex_index = M.get_edge_vertex_index(edge, v);
        			if ( !processed_vertices[vertex_index] ) {
        				processed_vertices[vertex_index] = true;
        				K.add(vertex_index, vertex_index, penalty_coefficient);
        				F[vertex_index] += penalty_coefficient * values[vertex_index];
        			}
        		}
        	}
        }
        
    }





    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        
        
    }

}

























