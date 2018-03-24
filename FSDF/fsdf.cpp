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

void OBBNode::init(int param_max_level, double param_max_side) {
	this->param_max_level = param_max_level;
	this->param_max_side = param_max_side;

	this->idx.clear();
	for (size_t i = 0; i < this->vertices.rows(); ++i) {
		idx.push_back((int)i);
	}
	this->nb_pts = idx.size();
	this->depth = 0;
}

void OBBNode::build_tree() {
	cout << this->depth << endl;
	cout << "nb pts " << this->nb_pts << endl;
	this->compute_obb();

	Vector3d side0 = this->box.row(0);  //side0
	Vector3d side1 = this->box.row(1);  //side1
	Vector3d side2 = this->box.row(2);  //side2
	Vector3d corner = this->box.row(3); //side3

	// Stopping conditions
	if (this->depth >= this->param_max_level) {
		cout << "return 1" << endl;
		return;
	}
	if (this->nb_pts <= 3) {
		cout << "return 2" << endl;
		return;
	}
	if (side2.norm() < 3* this->param_max_side) { //the times 3 is rather arbitrary
		cout << "return 3" << endl;
		return;
	}

	// Otherwise split along the longuest side
	if (side2.norm() < side1.norm() || side1.norm() < side0.norm() || side2.norm() < side0.norm()) {
		throw std::invalid_argument("build_tree(): should be ||side2|| >= ||side1|| >= ||side0||.");
	}
	vector<int> idx_L;
	vector<int> idx_R;

	MatrixXd pcl_temp(this->idx.size(), 3); //there are better ways to slice
	for (size_t i = 0; i < idx.size(); ++i) {
		pcl_temp.row(i) = this->vertices.row(idx.at(i));
	}

	Vector3d center = pcl_temp.colwise().mean(); //side2 is the longest side
	double t_center = (center - corner).dot(side2)/side2.norm();
	double t_pt;
	for (size_t i = 0; i < idx.size(); ++i) {
		t_pt = (pcl_temp.row(i) - corner.transpose()).dot(side2) / side2.norm(); //project each point on side2
		if (t_pt < t_center) {
			idx_L.push_back(idx.at(i));
		}
		if (t_pt >= t_center) {
			idx_R.push_back(idx.at(i));
		}
	}

	if (static_cast<long int>(idx_L.size()) + static_cast<long int>(idx_R.size()) != static_cast<long int>(this->idx.size())) {
		throw std::invalid_argument("build_tree(): sum of splits is not equal to init size!");
	}
	if (idx_L.size() == 0 || idx_R.size() == 0) {
		cout << "return 4" << endl;
		return; // Still no split!
	}


	cout << "step one" << endl;
	// Split will happen here
	this->idx.clear(); //free memory

	this->left = unique_ptr<OBBNode>(new OBBNode(this->vertices, this->normals));
	this->left->param_max_level = this->param_max_level;
	this->left->param_max_side = this->param_max_side;
	this->left->set_idx(idx_L); 
	this->left->depth = this->depth + 1;
	this->left->nb_pts = idx_L.size();
	this->left->build_tree(); //<========= Recursion

	this->right = unique_ptr<OBBNode>(new OBBNode(this->vertices, this->normals));
	this->right->param_max_level = this->param_max_level;
	this->right->param_max_side = this->param_max_side;
	this->right->set_idx(idx_R);
	this->right->depth = this->depth + 1;
	this->right->nb_pts = idx_R.size();
	this->right->build_tree(); //<========= Recursion

	return;
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

	MatrixXd temp(5, 3);
	temp.row(0) = (tMax[0] - tMin[0]) * ev0; //side0
	temp.row(1) = (tMax[1] - tMin[1]) * ev1; //side1
	temp.row(2) = (tMax[2] - tMin[2]) * ev2; //side2
	temp.row(3) = center + tMin[0] * ev0 + tMin[1] * ev1 + tMin[2] * ev2; //corner
	temp.row(4) = eig.eigenvalues();

	this->box = temp;

	//cout << "box " << this->box << endl;
};


vector<int> OBBNode::ray_intersect(const Vector3d & source, const Vector3d & dir) const {
	vector<int> pt_list;
	return pt_list;
}

