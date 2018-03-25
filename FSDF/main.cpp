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

	root.compute_obb();
	root.build_tree();
	//root.compute_obb();
	MatrixXd box = root.get_obb();
	cout << "gh " << root.get_obb().size() << endl;

	cin.ignore();
	
	//if (side2.norm() > side0.norm() && side2.norm() > side1.norm())
	//	cout << "side 2 is the longest" << endl;
	//else
	//	cout << "side 2 is NOT the longest" << endl;

	//vector<int> idx; 
	//for (size_t i = 0; i < VA.rows(); ++i) {
	//	idx.push_back((int)i);
	//}

	//vector<int> idx1;
	//vector<int> idx2;

	//Vector3d center = VA.colwise().mean(); //side2 is the longest side
	//double t_center = (center - corner).dot(side2)/side2.norm();
	//double t_pt;
	//int count1 = 0;
	//int count2 = 0;
	//for (size_t i = 0; i < VA.rows(); ++i) {
	//	t_pt = (VA.row(i) - corner.transpose()).dot(side2) / side2.norm();
	//	if (t_pt < t_center) {
	//		count1++;
	//		idx1.push_back(idx.at(i));
	//	}
	//	if (t_pt >= t_center) {
	//		count2++;
	//		idx2.push_back(idx.at(i));
	//	}
	//}

	//ComputeOBB(VA, corner, side0, side1, side2, evals); left
	//ComputeOBB(VA, corner, side0, side1, side2, evals); right

	//cout << "number of pts = " << VA.rows() << endl;
	//cout << "pts in group 1 = " << count1 << endl;
	//cout << "pts in group 2 = " << count2 << endl;*/

	//for any pt in pcl 
	//	if t_pt < t_center
	//		go into group 1
	//	else 
	//		go into group 2

	/*
	MatrixXd VB(8,3);
	MatrixXi FB(12,3);
	MatrixXd CB(12,3);

	VB.row(0) = corner.transpose();
	VB.row(1) = corner.transpose() + side0.transpose();
	VB.row(2) = corner.transpose() + side1.transpose();
	VB.row(3) = corner.transpose() + side2.transpose();
	VB.row(4) = corner.transpose() + side0.transpose() + side1.transpose();
	VB.row(5) = corner.transpose() + side0.transpose() + side2.transpose();
	VB.row(6) = corner.transpose() + side1.transpose() + side2.transpose();
	VB.row(7) = corner.transpose() + side0.transpose() + side1.transpose() + side2.transpose();

	FB << 0, 2, 1,
		2, 4, 1,
		0, 3, 6,
		2, 0, 6,
		3, 5, 7,
		6, 3, 7,
		5, 1, 4,
		7, 5, 4,
		0, 5, 3,
		0, 1, 5,
		7, 4, 2,
		6, 7, 2;

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
	//viewer.data().set_mesh(VA, FA);

	// Use the z coordinate as a scalar field over the surface
	//Eigen::VectorXd Z = VV.col(2);

	// Compute per-vertex colors
	//igl::jet(Z, true, C);
	// Add per-vertex colors
	//viewer.data().set_colors(CA);
	// Launch the viewer
	//viewer.launch();	
}