#include "fsdf.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Node
//////////////////////////////////////////////////////////////////////////////////////////////////////////

OBBNode::OBBNode(MatrixXd & vertices, MatrixXd & normals) : vertices(vertices), normals(normals) {
	this->left = nullptr;
	this->right = nullptr;
	this->nb_pts = 0;
	this->depth = 0;
	this->idx.clear();
}

OBBNode::init() {

}

int OBBNode::get_depth() { 
	return this->depth;
}

void OBBNode::set_idx(vector<int> idx) { //want copy
	this->idx = idx;
}

void OBBNode::set_obb(MatrixXd box) {
	this->box = box;
}

MatrixXd OBBNode::get_obb() const {
	return this->box;
};

void OBBNode::compute_obb() {
	
	// Slicing 
	MatrixXd V(this->idx.size(), 3);
	for (size_t i = 0; i < this->idx.size(); ++i) {
		V.row(i) = this->vertices.row(this->idx.at(i));
	}

	Vector3d center = V.colwise().mean();
	MatrixXd centered = V.rowwise() - center.transpose();
	MatrixXd cov = (centered.adjoint() * centered) / double(V.rows() - 1);

	// Get eigenvectors from point cloud
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	Vector3d ev0 = eig.eigenvectors().col(0); //min
	Vector3d ev1 = eig.eigenvectors().col(1); //mid
	Vector3d ev2 = eig.eigenvectors().col(2); //max

											  // Vernier
	Vector3d tMin = DBL_MAX * Vector3d::Ones(3);
	Vector3d tMax = -DBL_MAX * Vector3d::Ones(3);

	// Create oriented bounding box by projecting points onto eigenvectors.
	for (size_t i = 0; i < centered.rows(); ++i) { //Still in col major, potential optimiation
		Vector3d pt = centered.row(i);

		// a + ab * dot(ap,ab)/dot(ab,ab), a = means, p = pt, b = amax (or amid or amin)
		double t;
		t = pt.dot(ev0) / ev0.norm();
		if (t < tMin[0])
			tMin[0] = t;
		if (t > tMax[0])
			tMax[0] = t;

		t = pt.dot(ev1) / ev1.norm();
		if (t < tMin[1])
			tMin[1] = t;
		if (t > tMax[1])
			tMax[1] = t;

		t = pt.dot(ev2) / ev2.norm();
		if (t < tMin[2])
			tMin[2] = t;
		if (t > tMax[2])
			tMax[2] = t;
	}

	this->box.row(0) = (tMax[0] - tMin[0]) * ev0; //side0
	this->box.row(1) = (tMax[1] - tMin[1]) * ev1; //side1
	this->box.row(2) = (tMax[2] - tMin[2]) * ev2; //side2
	this->box.row(3) = center + tMin[0] * ev0 + tMin[1] * ev1 + tMin[2] * ev2; //corner
	this->box.row(4) = eig.eigenvalues();
};

void OBBNode::build_tree() {
	//// First level
	//MatrixXd compute_obb(vector<int> & idx);

	////one iteration
	//root->left = unique_ptr<OBBNode>(new OBBNode(this->pcl));
	//root->right = unique_ptr<OBBNode>(new OBBNode(this->pcl));
	////stopping condition: max_level or obb statistics
	//return;
	return;
}

vector<int> OBBNode::ray_intersect(const Vector3d & source, const Vector3d & dir) const {
	vector<int> pt_list;
	return pt_list;
}

