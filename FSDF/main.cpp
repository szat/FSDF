//#include "stdafx.h"
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <igl/median.h>
#include <igl/igl_inline.h>
#include <igl/matrix_to_list.h>

#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/floor.h>

#include "tutorial_shared_path.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <float.h>
#include <ctime>
#include <string>

Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

struct OBB {
	Eigen::Vector3d corner;
	Eigen::Vector3d side0;
	Eigen::Vector3d side1;
	Eigen::Vector3d side2;
	Eigen::Vector3d evals;
};

OBB * ComputeOBB(const Ref<const MatrixXd> PCL, vector<int> & idx)
{
	OBB * box = new OBB(); //on the heap

	//I know this is hugely unoptimal, but it is the simplest
	MatrixXd V(idx.size(), 3);
	for (size_t i = 0; i < idx.size(); ++i) {
		V.row(i) = PCL.row(idx.at(i));
	}

	Vector3d center = V.colwise().mean();
	MatrixXd centered = V.rowwise() - center.transpose();
	MatrixXd cov = (centered.adjoint() * centered) / double(V.rows() - 1);

	// Get eigenvectors from point cloud
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	box->evals = eig.eigenvalues();
	Vector3d ev0 = eig.eigenvectors().col(0); //min
	Vector3d ev1 = eig.eigenvectors().col(1); //mid
	Vector3d ev2 = eig.eigenvectors().col(2); //max

	//cout << "Inside ComputeOBB" << endl;
	//cout << "||ev0|| = " << ev0.norm() << endl;
	//cout << "||ev1|| = " << ev1.norm() << endl;
	//cout << "||ev2|| = " << ev2.norm() << endl;
	//cout << "evals = " << evals << endl;

	Vector3d tMin = DBL_MAX * Vector3d::Ones(3);
	Vector3d tMax = -DBL_MAX * Vector3d::Ones(3);

	// Create oriented bounding box by projecting points onto eigenvectors.
	for (size_t i = 0; i < centered.rows(); ++i) { //Still in col major, potential optimiation
		Vector3d pt = centered.row(i);

		// a + ab * dot(ap,ab)/dot(ab,ab), a = means, p = pt, b = amax (or amid or amin)
		double t;
		t = pt.dot(ev0)/ev0.norm();
		if (t < tMin[0])
			tMin[0] = t;
		if (t > tMax[0])
			tMax[0] = t;

		t = pt.dot(ev1)/ev1.norm();
		if (t < tMin[1])
			tMin[1] = t;
		if (t > tMax[1])
			tMax[1] = t;
		
		t = pt.dot(ev2)/ev2.norm();
		if (t < tMin[2])
			tMin[2] = t;
		if (t > tMax[2])
			tMax[2] = t;
	}

	box->corner = center + tMin[0] * ev0 + tMin[1] * ev1 + tMin[2] * ev2;
	box->side0 = (tMax[0] - tMin[0]) * ev0;
	box->side1 = (tMax[1] - tMin[1]) * ev1;
	box->side2 = (tMax[2] - tMin[2]) * ev2;
	return box;
}

void SplitPCL(const Ref<const MatrixXd> V, const Ref<const Vector3d> corner, const Ref<const Vector3d> side0, const Ref<const Vector3d> side1, const Ref<const Vector3d> side2, const Ref<const Vector3d> evals) {

}

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

	Vector3d corner, side0, side1, side2, evals;
	OBB *box = new OBB();

	vector<int> idx;
	for (size_t i = 0; i < VA.rows(); ++i) {
		idx.push_back(i);
	}

 	box = ComputeOBB(VA,idx);

	corner = box->corner;
	side0 = box->side0;
	side1 = box->side1;
	side2 = box->side2;
	evals = box->evals;

	
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

/*
	cout << "number of pts = " << VA.rows() << endl;
	cout << "pts in group 1 = " << count1 << endl;
	cout << "pts in group 2 = " << count2 << endl;*/

	//for any pt in pcl 
	//	if t_pt < t_center
	//		go into group 1
	//	else 
	//		go into group 2

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

	//VectorXi vec(11);
	//vec << -3, 4, 5, 2, 2, 7,1,2,3,4,9;
	//// 100 random indices into rows of F
	//int I;
	//// median
	//igl::median(vec, I);
	//std::cout << "vec " << vec << endl;
	//std::cout << "I " << I << endl;

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(VU, FU);

	// Use the z coordinate as a scalar field over the surface
	//Eigen::VectorXd Z = VV.col(2);

	// Compute per-vertex colors
	//igl::jet(Z, true, C);
	// Add per-vertex colors
	viewer.data().set_colors(CU);
	// Launch the viewer
	viewer.launch();	
}