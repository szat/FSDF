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

#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/floor.h>

#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/is_irregular_vertex.h>

#include "tutorial_shared_path.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <float.h>
#include <ctime>
#include <string>
#include "fsdf.h"

Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

/*
int SplitOBB(const Ref<const MatrixXd> PCL, const OBB & box, vector<int> & idx_in, vector<int> & idx_L, vector<int> & idx_R){
	//side2, side1, and side0 should be orthonormal
	Vector3d side0 = box.side0;
	Vector3d side1 = box.side1;
	Vector3d side2 = box.side2;
	double volume = side0.norm() * side1.norm() * side2.norm();
	double area = side1.norm() * side2.norm();

	//return the index of the side if there is a split
	//return -1 if there is no split 

	//I know this is hugely unoptimal, but it is the simplest
	MatrixXd V(idx_in.size(), 3);
	for (size_t i = 0; i < idx_in.size(); ++i) {
		V.row(i) = PCL.row(idx_in.at(i));
	}

	Vector3d center = V.colwise().mean(); //side2 is the longest side
	double t_center = (center - box.corner).dot(box.side2)/box.side2.norm();
	double t_pt;
	int count1 = 0;
	int count2 = 0;
	for (size_t i = 0; i < V.rows(); ++i) {
		t_pt = (V.row(i) - box.corner.transpose()).dot(box.side2) / box.side2.norm();
		if (t_pt < t_center) {
			count1++;
			idx_L.push_back(idx_in.at(i));
		}
		if (t_pt >= t_center) {
			count2++;
			idx_R.push_back(idx_in.at(i));
		}
	}

	return 2;
}
*/

int main(int argc, char *argv[])
{
	//Read data
	MatrixXd VA;
	MatrixXi FA;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", VA, FA);
	MatrixXd CA(FA.rows(), 3);
	CA << RowVector3d(0.2, 0.3, 0.8).replicate(FA.rows(), 1);

	//An affine transformation
	Eigen::Affine3d T = Eigen::Affine3d::Identity();
	T.rotate(Eigen::AngleAxisd(0.5, Eigen::Vector3d(-1, 1, 3)));
	T.translate(Eigen::Vector3d(0, 0.125*cos(0.5), 0));
	VA = VA * T.matrix().block(0, 0, 3, 3).transpose();
	Eigen::RowVector3d trans = T.matrix().block(0, 3, 3, 1).transpose();
	VA = (VA.rowwise() + trans).eval();

	vector<int> idx;
	for (size_t i = 0; i < VA.rows(); ++i) 
		idx.push_back(i);

 	//box = ComputeOBB(VA,idx);
	MatrixXd NA;
	igl::per_vertex_normals(VA, FA, NA);
	OBBNode root(VA, NA);
	VectorXd area;
	igl::doublearea(VA, FA, area);
	double param_max_side = area.maxCoeff();
	root.init(15, param_max_side);

	//root.compute_obb();
	root.build_tree();
	vector<int> truth = root.validate();
	cout << "build done " << endl;
	cout << "sum of idxs " << truth.size() << endl;
	cout << "size of initial vertices " << VA.rows() << endl;
	Vector3d source = VA.row(0);
	Vector3d dir = NA.row(0);
	cout << "test in box " << root.in_box(source) << endl;
	cout << "test in box left " << root.left->in_box(source) << endl;
	cout << "test in box right " << root.right->in_box(source) << endl;
	cout << "test intersect left " << root.left->intersect_box(source, dir) << endl;
	cout << "test intersect right " << root.right->intersect_box(source, dir) << endl;

	vector<int> temp = root.ray_intersect(source, dir);
	cout << "vector<int> temp length = " << temp.size() << endl;
	//root.compute_obb();
	//MatrixXd box = root.get_obb();
	
	cin.ignore();
	
	/*
	CB << RowVector3d(0.2, 0.3, 0.8).replicate(FB.rows(), 1);

	MatrixXd VU(VA.rows() + VB.rows(), VA.cols());
	VU << VA,
		VB;
	MatrixXi FU(FA.rows() + FB.rows(), VA.cols());
	FU << FA,
		(FB.array() + VA.rows());
	MatrixXd CU(FU.rows(), 3);
	CU << CA,
		CB;

		*/

	// Plot the mesh
	//igl::opengl::glfw::Viewer viewer;
	//viewer.data().set_mesh(VU, FU);

	// Use the z coordinate as a scalar field over the surface
	//Eigen::VectorXd Z = VV.col(2);

	// Compute per-vertex colors
	//igl::jet(Z, true, C);
	// Add per-vertex colors
	//viewer.data().set_colors(CA);
	// Launch the viewer
	//viewer.launch();	
}