///only move sematic, not redundant
//https://docs.microsoft.com/en-gb/cpp/cpp/explicitly-defaulted-and-deleted-functions
//http://blog.matejzavrsnik.com/using_smart_pointers_in_building_trees_in_which_child_nodes_need_an_access_to_the_parent.html
//http://www.opengl-tutorial.org/miscellaneous/clicking-on-objects/picking-with-custom-ray-obb-function/

//#include "stdafx.h"
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/median.h>
#include <igl/igl_inline.h>
#include <igl/matrix_to_list.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>

#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/is_irregular_vertex.h>

#include "tutorial_shared_path.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <float.h>
#include <string>
#include "fsdf.h"

Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
	//Read data
	MatrixXd V;
	MatrixXi F;
	MatrixXd N;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", V, F);
	igl::per_vertex_normals(V, F, N);

	//An affine transformation
	Eigen::Affine3d T = Eigen::Affine3d::Identity();
	T.rotate(Eigen::AngleAxisd(0.5, Eigen::Vector3d(-1, 1, 3)));
	T.translate(Eigen::Vector3d(0, 0.125*cos(0.5), 0));
	V = V * T.matrix().block(0, 0, 3, 3).transpose();
	Eigen::RowVector3d trans = T.matrix().block(0, 3, 3, 1).transpose();
	V = (V.rowwise() + trans).eval();

	cout << "Welcome to a prototype of the Offset Surface Diameter Function" << endl;
	cout << "Author: Adrian Szatmari" << endl;
	cout << "Date: 25 March 2018" << endl << endl;

	cout << "Note: \t this prototype can be easily greatly optimised." << endl;
	cout << "\t some conceptual shortcomings in the original paper transpire here" << endl;

	cout << "Starting program..." << endl;
	cout << "The mesh has " << V.rows() << " vertices." << endl;
	
	OBBNode root(V, N);
	VectorXd area;
	igl::doublearea(V, F, area);
	double param_max_side = area.maxCoeff();
	int max_level = 15;

	root.init(max_level, param_max_side);
	cout << "Starting to build an OBB Tree with max depth : " << max_level << endl;

	root.build_tree();
	
	cout << "Finished building OBB Tree..." << endl;
	cout << "Starting to querying every vertex in the mesh..." << endl;
	
	VectorXd queries = root.query();

	cout << "Finished querying every vertex in the mesh..." << endl;
	cout << "The patches in blue are queries 'closer' queries." << endl << endl;

	MatrixXd C(F.rows(), 3);

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);

	// Compute per-vertex colors
	igl::jet(queries, true, C);
	// Add per-vertex colors
	viewer.data().set_colors(C);
	// Launch the viewer
	viewer.launch();	
}