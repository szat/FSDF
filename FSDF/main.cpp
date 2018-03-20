//#include "stdafx.h"
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include "tutorial_shared_path.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <float.h>
#include <ctime>
//#include <igl/copyleft/cgal/mesh_boolean.h>

Eigen::MatrixXd VA;
Eigen::MatrixXi FA;
Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

void ComputeOBB(const Ref<const MatrixXd> V, Ref<Vector3d> corner, Ref<Vector3d> side0, Ref<Vector3d> side1, Ref<Vector3d> side2, Ref<Vector3d> evals)
{
	Vector3d center = V.colwise().mean();
	MatrixXd centered = V.rowwise() - center.transpose();
	MatrixXd cov = (centered.adjoint() * centered) / double(V.rows() - 1);

	// Get eigenvectors from point cloud
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	evals = eig.eigenvalues();
	Vector3d ev0 = eig.eigenvectors().col(0);
	Vector3d ev1 = eig.eigenvectors().col(1);
	Vector3d ev2 = eig.eigenvectors().col(2);

	cout << "Inside ComputeOBB" << endl;
	cout << "||ev0|| = " << ev0.norm() << endl;
	cout << "||ev1|| = " << ev1.norm() << endl;
	cout << "||ev2|| = " << ev2.norm() << endl;
	cout << "evals = " << evals << endl;

	Vector3d tMin = DBL_MAX * Vector3d::Ones(3);
	Vector3d tMax = -DBL_MAX * Vector3d::Ones(3);

	// Create oriented bounding box by projecting points onto eigenvectors.
	for (size_t i = 0; i < centered.rows(); ++i) { //Still in col major, potential optimiation
		Vector3d pt = centered.row(i);

		// a + ab * dot(ap,ab)/dot(ab,ab), a = means, p = pt, b = amax (or amid or amin)
		double t;
		t = pt.dot(ev0)/(ev0).dot(ev0);
		if (t < tMin[0])
			tMin[0] = t;
		if (t > tMax[0])
			tMax[0] = t;

		t = (pt).dot(ev1)/(ev1).dot(ev1);
		if (t < tMin[1])
			tMin[1] = t;
		if (t > tMax[1])
			tMax[1] = t;
		
		t = pt.dot(ev2)/(ev2).dot(ev2);
		if (t < tMin[2])
			tMin[2] = t;
		if (t > tMax[2])
			tMax[2] = t;
	}

	corner = center + tMin[0] * ev0 + tMin[1] * ev1 + tMin[2] * ev2;
	side0 = (tMax[0] - tMin[0]) * ev0;
	side1 = (tMax[1] - tMin[1]) * ev1;
	side2 = (tMax[2] - tMin[2]) * ev2;
}


int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	igl::readOFF(TUTORIAL_SHARED_PATH "/screwdriver.off", VA, FA);
	
	MatrixXd VS;
	MatrixXi FS;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/sphere.obj", VS, FS);
	VS = VS * 0.02;

	Vector3d means = VA.colwise().mean();

	VS = VS + means.transpose();

	cout << "nb rows " << VA.rows() << endl;
	cout << "nb cols " << VA.cols() << endl;
	cout << "mean " << VA.mean() << endl;
	cout << "mean cols" << VA.colwise().mean()  << endl;

	Vector3d corner, side0, side1, side2, evals;
	ComputeOBB(VA, corner, side0, side1, side2, evals);

	/*
	MatrixXd VB(8,3);
	MatrixXi FB(3,3);

	VB.row(0) = corner.transpose();
	VB.row(1) = corner.transpose() + max.transpose();
	VB.row(2) = corner.transpose() + mid.transpose();
	VB.row(3) = corner.transpose() + min.transpose();

	FB.row(0) << 0, 1, 2;
	FB.row(1) << 0, 2, 3;
	FB.row(2) << 0, 1, 3;
	// Concatenate (VA,FA) and (VB,FB) into (V,F)
	Eigen::MatrixXd V(VA.rows() + VS.rows(), 3);
	V << VA, 
		VS;

	Eigen::MatrixXi F(FA.rows() + FS.rows(), 3);
	F << FA, 
		(FS.array() + VA.rows());

	// blue color for faces of first mesh, orange for second
	Eigen::MatrixXd C(F.rows(), 3);
	C <<
		Eigen::RowVector3d(0.2, 0.3, 0.8).replicate(FA.rows(), 1),
		Eigen::RowVector3d(1.0, 0.7, 0.2).replicate(FS.rows(), 1);



	//cout << "corner " << corner << endl;
	//cout << Vbox << endl;

	//Viz by concatenation and by opening the half of the cube
	//Viz in Wireframe

	//cin.ignore();

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);

	// Use the z coordinate as a scalar field over the surface
	//Eigen::VectorXd Z = VV.col(2);

	// Compute per-vertex colors
	//igl::jet(Z, true, C);

	// Add per-vertex colors
	//viewer.data().set_colors(C);

	// Launch the viewer
	viewer.launch();
	*/
}