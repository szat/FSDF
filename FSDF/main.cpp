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

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

void ComputeOBB(const Ref<const MatrixXd> V, Ref<Vector3d> corner, Ref<Vector3d> max, Ref<Vector3d> mid, Ref<Vector3d> min, Ref<Vector3d> size)
{
	Vector3d means = V.colwise().mean();
	MatrixXd centered = V.rowwise() - means.transpose();
	MatrixXd cov = (centered.adjoint() * centered) / double(V.rows() - 1);

	// Get eigenvectors from point cloud
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	Vector3d evals = eig.eigenvalues();
	MatrixXd evecs = eig.eigenvectors();

	Vector3d vmax = evecs.col(0);
	Vector3d vmid = evecs.col(1);
	Vector3d vmin = evecs.col(2);

	Vector3d amax = means + vmax;
	Vector3d amid = means + vmid;
	Vector3d amin = means + vmin;

	Vector3d tMin = DBL_MAX * Vector3d::Ones(3);
	Vector3d tMax = -DBL_MAX * Vector3d::Ones(3);

	// Create oriented bounding box by projecting points onto eigenvectors.
	for (size_t i = 0; i < V.rows(); ++i) { //Still in col major, potential optimiation
		Vector3d pt = V.row(i);

		// a + ab * dot(ap,ab)/dot(ab,ab), a = means, p = pt, b = amax (or amid or amin)
		double t;
		t = ((pt - means).dot(amax - means) / (amax - means).dot(amax - means));
		if (t < tMin[0])
			tMin[0] = t;
		if (t > tMax[0])
			tMax[0] = t;
		t = ((pt - means).dot(amid - means) / (amid - means).dot(amid - means));
		if (t < tMin[1])
			tMin[1] = t;
		if (t > tMax[1])
			tMax[1] = t;
		t = ((pt - means).dot(amin - means) / (amin - means).dot(amin - means));
		if (t < tMin[2])
			tMin[2] = t;
		if (t > tMax[2])
			tMax[2] = t;
	}

	corner = means + tMin[0] * vmax + tMin[1] * vmid + tMin[2] * vmin;
	max = (tMax[0] - tMin[0]) * vmax;
	mid = (tMax[1] - tMin[1]) * vmid;
	min = (tMax[2] - tMin[2]) * vmin;
}


int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	igl::readOFF(TUTORIAL_SHARED_PATH "/screwdriver.off", V, F);
	
	cout << "nb rows " << V.rows() << endl;
	cout << "nb cols " << V.cols() << endl;
	cout << "mean " << V.mean() << endl;
	cout << "mean cols" << V.colwise().mean()  << endl;

	Vector3d corner, max, mid, min, size;
	ComputeOBB(V, corner, max, mid, min, size);

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