#include "fsdf.h"

OBBNode::OBBNode(MatrixXd & vertices, MatrixXd & normals) : vertices(vertices), normals(normals) {
	this->left = nullptr;
	this->right = nullptr;
	this->nb_pts = 0;
	this->depth = 0;
	this->idx.clear();
	this->longest_side = -1;
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
	cout << "build_tree(): depth = " << this->depth << ", nb_pts = " << this->nb_pts << ", idx length = " << this->idx.size() << endl;
	this->compute_obb(); //all nodes have computed obb

	MatrixXd sides(3, 3);
	sides.row(0) = this->box.row(0);  //side0
	sides.row(1) = this->box.row(1);  //side1
	sides.row(2) = this->box.row(2);  //side2
	Vector3d corner = this->box.row(3); //side3

	int longest = -1;
	double measure = -1;
	for (int i = 0; i < 3; ++i) {
		if (sides.row(i).norm() > measure) {
			longest = i;
			measure = sides.row(i).norm();
		}
	}
	this->longest_side = longest;

	//cout << "build_tree 2" << endl;

	// Stopping conditions
	if (this->depth >= this->param_max_level) {
		//cout << "return 1" << endl;
		return;
	}
	if (this->nb_pts <= 3) {
		//cout << "return 2" << endl;
		return;
	}
	if (sides.row(longest).norm() < 3* this->param_max_side) { //the times 3 is rather arbitrary
		//cout << "return 3" << endl;
		return;
	}

	vector<int> idx_L;
	vector<int> idx_R;

	MatrixXd pcl_temp(this->idx.size(), 3); //there are better ways to slice
	for (size_t i = 0; i < idx.size(); ++i) {
		pcl_temp.row(i) = this->vertices.row(idx.at(i));
	}

	//cout << "build_tree 5" << endl;
	Vector3d center = pcl_temp.colwise().mean(); //side2 is the longest side
	double t_center = (center - corner).dot(sides.row(longest))/ sides.row(longest).norm();
	double t_pt;
	for (size_t i = 0; i < idx.size(); ++i) {
		t_pt = (pcl_temp.row(i) - corner.transpose()).dot(sides.row(longest)) / sides.row(longest).norm(); //project each point on side2
		if (t_pt < t_center) {
			idx_L.push_back(idx.at(i));
		}
		if (t_pt >= t_center) {
			idx_R.push_back(idx.at(i));
		}
	}

	//cout << "build_tree 6" << endl;
	if (static_cast<long int>(idx_L.size()) + static_cast<long int>(idx_R.size()) != static_cast<long int>(this->idx.size())) {
		throw std::invalid_argument("build_tree(): sum of splits is not equal to init size!");
	}
	if (idx_L.size() == 0 || idx_R.size() == 0) {
		cout << "return 4" << endl;
		return; // Still no split!
	}

	//cout << "step one" << endl;
	// Split will happen here, along the longest side
	this->idx.clear(); //free memory
	//cout << "after idx.clear()" << endl;

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
};

bool OBBNode::is_leaf() const {
	if (this->left == nullptr &&  this->right == nullptr || this->depth == this->param_max_level) {
		return true;
	}
	else {
		return false;
	}
}

bool OBBNode::in_box(const Vector3d & source) const
{
	bool in_box = true;
	Vector3d corner = this->box.row(3);
	for (int i = 0; i < 3; ++i) { //the three sides
		double t = (source - corner).dot(this->box.row(i)) / this->box.row(i).norm();
		if (t < 0 || t > 1) in_box = false;
	}
	return in_box;
}

bool OBBNode::intersect_box(const Vector3d & source, const Vector3d & dir) const 
{
	// From Geometric tools for computer graphics, p. 632.
	// Can be modified for the corner/side representation, but I am running out of time.
	// Create the same representation as used in GTCG (book above)
	Vector3d center = this->box.row(4) + (this->box.row(0) + this->box.row(1) + this->box.row(2)) / 2;
	Vector3d halves;
	halves << this->box.row(0).norm() / 2, this->box.row(1).norm() / 2, this->box.row(2).norm() / 2;
	MatrixXd axis(3, 3);
	axis.row(0) = this->box.row(0) / this->box.row(0).norm();
	axis.row(1) = this->box.row(1) / this->box.row(1).norm();
	axis.row(2) = this->box.row(2) / this->box.row(2).norm();
	Vector3d dirN = dir / dir.norm();

	double tNear = -DBL_MAX;
	double tFar = DBL_MAX; 
	double r,s, t0, t1;
	for (int i = 0; i < 3; ++i) { // Check for ray parallel to planes
		if (abs(dirN.dot(axis.row(i)) < numeric_limits<double>::epsilon())) {
			// Ray parallel to planes
			r = axis.row(i).dot(center - source);
			if (-r - halves[i] > 0 || -r + halves[i] > 0) {
				// No intersection
				return false;
			}
		}
		r = axis.row(i).dot(center - source);
		s = axis.row(i).dot(dirN);
		// Ray not parallel to planes, so find parameters of intersections
		t0 = (r + halves[i]) / s;
		t1 = (r - halves[i]) / s;
		// Check ordering
		if (t0 > t1) {
			// Swap them
			double tmp = t0;
			t0 = t1;
			t1 = tmp;
		}
		// Compare with current values
		if (t0 > tNear) {
			tNear = t0;
		}
		if (t1 < tFar) {
			tFar = t1;
		}
		// Check if ray misses entirely
		if (tNear > tFar) {
			return false;
		}
		if (tFar < 0) {
			return false;
		}
	}
	// Box definitely intersected
	if (tNear > 0) {
		double tIntersect = tNear;
	}
	else {
		double tIntersect = tFar;
	}
	return true;
}


vector<int> OBBNode::ray_intersect(const Vector3d & source, const Vector3d & dir) const {
	vector<int> cat_out; //concatenate by recurrence
	//if (this->is_leaf()) {
	//	if (this->in_box(source)) { //if the source ppt is in the obb, then we are in the "mother" obb
	//		vector<int> empty;
	//		return empty;
	//	}
	//	else {
	//		//check wether ray intersects the obb
	//		vector<int> sphere_list;
	//		//return a sphere (list)
	//		return sphere_list;
	//	}
	//}
	//else {
	//	//if source is in left
	//	this->left->ray_intersect()
	//	//if source is in right
	//	this->right->ray_intersect()
	//}
	////if line intersects
	//this->box;
	////return sphere index list
	////

	////do not return intersections with spheres in its own obb
	//
	vector<int> pt_list;
	return pt_list;
}
