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

struct OBJ {
	Eigen::MatrixXd V;
	Eigen::MatrixXd C;
	Eigen::MatrixXi F;
};

Eigen::MatrixXd C;

using namespace Eigen;
using namespace std;

OBJ Union(const Ref<const MatrixXd> VA, const Ref<const MatrixXi> FA, const Ref<const MatrixXd> CA,
	const Ref<const MatrixXd> VB, const Ref<const MatrixXi> FB, const Ref<const MatrixXd> CB) {
	
	MatrixXd VU(VA.rows() + VB.rows(), VA.cols());
	VU << VA,
		VB;

	MatrixXi FU(FA.rows() + FA.rows(), VA.cols());
	FU << FA,
		(FB.array() + VA.rows());

	MatrixXd CU(FU.rows(), 3);
	CU << CA,
		CB;

	OBJ out;
	out.V = VU;
	out.F = FU;
	out.C = CU;

	return out;
}

void ComputeOBB(const Ref<const MatrixXd> V, Ref<Vector3d> corner, Ref<Vector3d> side0, Ref<Vector3d> side1, Ref<Vector3d> side2, Ref<Vector3d> evals)
{
	Vector3d center = V.colwise().mean();
	MatrixXd centered = V.rowwise() - center.transpose();
	MatrixXd cov = (centered.adjoint() * centered) / double(V.rows() - 1);

	// Get eigenvectors from point cloud
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	evals = eig.eigenvalues();
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
	MatrixXd VA;
	MatrixXi FA;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", VA, FA);
	MatrixXd CA(FA.rows(), 3);
	CA << RowVector3d(0.2, 0.3, 0.8).replicate(FA.rows(), 1);
	
	MatrixXd VS1;
	MatrixXi FS1;
	igl::readOBJ(TUTORIAL_SHARED_PATH "/sphere.obj", VS1, FS1);
	VS1 = VS1 * 0.05;
	MatrixXd CS1(FS1.rows(), 3);
	CS1 << RowVector3d(0.2, 0.8, 0.3).replicate(FS1.rows(), 1);

	MatrixXd VS2, VS3, VS4, VS5, VS6, VS7, VS8 = VS1;
	MatrixXi FS2, FS3, FS4, FS5, FS6, FS7, FS8 = FS1;
	MatrixXd CS2, CS3, CS4, CS5, CS6, CS7, CS8 = CS1;

	Vector3d corner, side0, side1, side2, evals;
	ComputeOBB(VA, corner, side0, side1, side2, evals);

	cout << "corner " << corner << endl;
	cout << "side0 " << side0 << endl;

	VS1 = VS1.rowwise() + (corner.transpose() +side0.transpose());// +side0.transpose();

	OBJ out = Union(VA, FA, CA, VS1, FS1, CS1);

	cout << "out size " << out.V.size() << endl;

	//Union
	MatrixXd VU(VA.rows() + VS1.rows(), VA.cols());
	VU << VA,
		VS1;

	MatrixXi FU(FA.rows() + FS1.rows(), VA.cols());
	FU << FA,
		(FS1.array() + VA.rows());

	MatrixXd CU(FU.rows(), 3);
	CU << CA,
		CS1;

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(out.V, out.F);

	// Use the z coordinate as a scalar field over the surface
	//Eigen::VectorXd Z = VV.col(2);

	// Compute per-vertex colors
	//igl::jet(Z, true, C);
	// Add per-vertex colors
	viewer.data().set_colors(out.C);
	// Launch the viewer
	viewer.launch();	
}