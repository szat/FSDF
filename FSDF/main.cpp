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


//MatrixXd is dynamic, Vector3d is not
//void ComputeOBB(MatrixXd & V, Vector3d & corner, Vector3d & max, Vector3d & mid, Vector3d & min, Vector3d & size)
//{
//	int i;
//	int numPts, pointId;
//	Vector3d x, mean, xp, *v[3], v0, v1, v2;
//	double *a[3], a0[3], a1[3], a2[3];
//	Vector3d tMin, tMax, closest;
//	double t;
//
//	Vector3d mean = V.colwise().mean();
//
//	MatrixXd centered = V.rowwise() - mean;
//	
//	MatrixXd cov = (centered.adjoint() * centered) / double(mat.rows() - 1);
//
//	
//	 Compute covariance matrix
//	
//	a[0] = a0; a[1] = a1; a[2] = a2;
//	for (i = 0; i < 3; i++)
//	{
//		a0[i] = a1[i] = a2[i] = 0.0;
//	}
//
//	for (pointId = 0; pointId < numPts; pointId++)
//	{
//		pts->GetPoint(pointId, x);
//		xp[0] = x[0] - mean[0]; xp[1] = x[1] - mean[1]; xp[2] = x[2] - mean[2];
//		for (i = 0; i < 3; i++)
//		{
//			a0[i] += xp[0] * xp[i];
//			a1[i] += xp[1] * xp[i];
//			a2[i] += xp[2] * xp[i];
//		}
//	}//for all points
//
//	for (i = 0; i < 3; i++)
//	{
//		a0[i] /= numPts;
//		a1[i] /= numPts;
//		a2[i] /= numPts;
//	}
//
//	
//	 Extract axes (i.e., eigenvectors) from covariance matrix.
//	
//	v[0] = v0; v[1] = v1; v[2] = v2;
//	vtkMath::Jacobi(a, size, v);
//	max[0] = v[0][0]; max[1] = v[1][0]; max[2] = v[2][0];
//	mid[0] = v[0][1]; mid[1] = v[1][1]; mid[2] = v[2][1];
//	min[0] = v[0][2]; min[1] = v[1][2]; min[2] = v[2][2];
//
//	for (i = 0; i < 3; i++)
//	{
//		a[0][i] = mean[i] + max[i];
//		a[1][i] = mean[i] + mid[i];
//		a[2][i] = mean[i] + min[i];
//	}
//
//	
//	 Create oriented bounding box by projecting points onto eigenvectors.
//	
//	tMin[0] = tMin[1] = tMin[2] = VTK_DOUBLE_MAX;
//	tMax[0] = tMax[1] = tMax[2] = -VTK_DOUBLE_MAX;
//
//	for (pointId = 0; pointId < numPts; pointId++)
//	{
//		pts->GetPoint(pointId, x);

		//do it for amax, amid, amin
//		for (i = 0; i < 3; i++)
//		{
//			vtkLine::DistanceToLine(x, mean, a[i], t, closest);
//			if (t < tMin[i])
//			{
//				tMin[i] = t;
//			}
//			if (t > tMax[i])
//			{
//				tMax[i] = t;
//			}
//		}
//	}//for all points
//
//	for (i = 0; i < 3; i++)
//	{
//		corner[i] = mean[i] + tMin[0] * max[i] + tMin[1] * mid[i] + tMin[2] * min[i];
//
//		max[i] = (tMax[0] - tMin[0]) * max[i];
//		mid[i] = (tMax[1] - tMin[1]) * mid[i];
//		min[i] = (tMax[2] - tMin[2]) * min[i];
//	}
//}


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
	cout << "cov rows " << cov.rows() << endl;
	cout << "cov cols " << cov.cols() << endl;

	//	 Extract axes (i.e., eigenvectors) from covariance matrix.

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	
	Vector3d evals = eig.eigenvalues();
	MatrixXd evecs = eig.eigenvectors();
	//Eigen::MatrixXd firstTwo = evecs.rightCols(2);

	cout << "evals " << evals << endl;
	cout << "evecs " << evecs << endl;

	// Get the two major eigenvectors and omit the others.
	//Eigen::MatrixXf evecs = eig.eigenvectors();
	
	Vector3d vmax = evecs.col(0);
	Vector3d vmid = evecs.col(1);
	Vector3d vmin = evecs.col(2);

	Vector3d amax = means + vmax;
	Vector3d amid = means + vmid;
	Vector3d amin = means + vmin;

	//	 Create oriented bounding box by projecting points onto eigenvectors.

	Vector3d tMin = DBL_MAX * Vector3d::Ones(3);
	Vector3d tMax = -DBL_MAX * Vector3d::Ones(3);

	cout << "tmin " << tMin << endl;
	cout << "tMax " << tMax << endl;

	for (size_t i = 0; i < V.rows(); ++i) {
		Vector3d pt = V.row(i);
		// Like in a distance to line computation
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

	cout << "tmin " << tMin << endl;
	cout << "tMax " << tMax << endl;
	
	Vector3d corner;
	corner = means + tMin[0] * vmax + tMin[1] * vmid + tMin[2] * vmin;
	vmax = (tMax[0] - tMin[0]) * vmax;
	vmid = (tMax[1] - tMin[1]) * vmid;
	vmin = (tMax[2] - tMin[2]) * vmin;

	//
		//	for (pointId = 0; pointId < numPts; pointId++)
		//	{
		//		pts->GetPoint(pointId, x);
		//		for (i = 0; i < 3; i++)
		//		{
		//			vtkLine::DistanceToLine(x, mean, a[i], t, closest); //t is not updated
		//			if (t < tMin[i])
		//			{
		//				tMin[i] = t;
		//			}
		//			if (t > tMax[i])
		//			{
		//				tMax[i] = t;
		//			}
		//		}
		//	}//for all points
		//
		//	for (i = 0; i < 3; i++)
		//	{
		//		corner[i] = mean[i] + tMin[0] * max[i] + tMin[1] * mid[i] + tMin[2] * min[i];
		//
		//		max[i] = (tMax[0] - tMin[0]) * max[i];
		//		mid[i] = (tMax[1] - tMin[1]) * mid[i];
		//		min[i] = (tMax[2] - tMin[2]) * min[i];
		//	}
		//}


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