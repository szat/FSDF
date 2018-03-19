//#include "stdafx.h"
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include "tutorial_shared_path.h"
#include <iostream>
#include <Eigen/Dense>


Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	igl::readOFF(TUTORIAL_SHARED_PATH "/screwdriver.off", V, F);
	
	cout << "nb rows " << V.rows() << endl;
	cout << "nb cols " << V.cols() << endl;
	cout << "mean " << V.mean() << endl;
	cout << "mean cols" << V.colwise().mean()  << endl;
	Vector3d means = V.colwise().mean();
	MatrixXd centered = V.rowwise() - means.transpose();

	cout << "centered rows" << centered.rows() << endl;
	cout << "centered cols" << centered.cols() << endl;
	cout << "new mean " << centered.colwise().mean() << endl;

	MatrixXd cov = (centered.adjoint() * centered) / double(V.rows() - 1);

	cin.ignore();

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);

	// Use the z coordinate as a scalar field over the surface
	Eigen::VectorXd Z = V.col(2);

	// Compute per-vertex colors
	igl::jet(Z, true, C);

	// Add per-vertex colors
	viewer.data().set_colors(C);

	// Launch the viewer
	viewer.launch();
}
