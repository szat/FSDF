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
//	cout << "build_tree(): depth = " << this->depth << ", nb_pts = " << this->nb_pts << ", idx length = " << this->idx.size() << endl;
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
	if (this->nb_pts <= 8) {
		//cout << "return 2" << endl;
		return;
	}
	if (sides.row(longest).norm() < 3 * this->param_max_side) { //the times 3 is rather arbitrary
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
	double t_center = (center - corner).dot(sides.row(longest)) / sides.row(longest).norm();
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
	this->left->build_tree(); //<========= Recursion Left

	this->right = unique_ptr<OBBNode>(new OBBNode(this->vertices, this->normals));
	this->right->param_max_level = this->param_max_level;
	this->right->param_max_side = this->param_max_side;
	this->right->set_idx(idx_R);
	this->right->depth = this->depth + 1;
	this->right->nb_pts = idx_R.size();
	this->right->build_tree(); //<========= Recursion Right

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

		// a + ab * dot(ap,ab)/norm(ab,ab), a = means, p = pt, b = amax (or amid or amin)
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
	double ux = (this->box.row(0)).dot(source);
	double vx = (this->box.row(1)).dot(source);
	double wx = (this->box.row(2)).dot(source);
	double up1 = (this->box.row(0)).dot(this->box.row(3));
	double up2 = (this->box.row(0)).dot(this->box.row(3) + this->box.row(0));
	double vp1 = (this->box.row(1)).dot(this->box.row(3));
	double vp2 = (this->box.row(1)).dot(this->box.row(3) + this->box.row(1));
	double wp1 = (this->box.row(2)).dot(this->box.row(3));
	double wp2 = (this->box.row(2)).dot(this->box.row(3) + this->box.row(2));

	if (up1 <= ux && ux <= up2 && vp1 <= vx && vx <= vp2 && wp1 <= wx && wx <= wp2) {
		return true;
	}
	else {
		return false;
	}
}

bool OBBNode::intersect_box(const Vector3d & source, const Vector3d & dir) const {
	MatrixXd cornerz(8, 3);
	cornerz.row(0) = this->box.row(3);
	cornerz.row(1) = this->box.row(3) + this->box.row(0);
	cornerz.row(2) = this->box.row(3) + this->box.row(1);
	cornerz.row(3) = this->box.row(3) + this->box.row(2);
	cornerz.row(4) = this->box.row(3) + this->box.row(0) + this->box.row(1);
	cornerz.row(5) = this->box.row(3) + this->box.row(0) + this->box.row(2);
	cornerz.row(6) = this->box.row(3) + this->box.row(1) + this->box.row(2);
	cornerz.row(7) = this->box.row(3) + this->box.row(0) + this->box.row(1) + this->box.row(2);

	for (int i = 0; i < 8; ++i) {
		double t = (cornerz.row(i) - source.transpose()).dot(dir) / dir.norm();
		if (this->in_box(source + t * dir)) {
			return true;
		}
	}
	return false;
}

vector<int> OBBNode::validate() const {
	if (this->is_leaf()) {
		return this->idx;
	}
	else {
		vector<int> out;
		vector<int> temp;
		temp = this->left->validate();  //<========= Recursion Left
		out.insert(std::end(out), std::begin(temp), std::end(temp));
		temp = this->right->validate(); //<========= Recursion Right
		out.insert(std::end(out), std::begin(temp), std::end(temp));
		return out;
	}
}

vector<int> OBBNode::ray_intersect(const Vector3d & source, const Vector3d & dir) const {
	if (this->is_leaf()) {
		return this->idx;
	}
	else {
		vector<int> out; 
		if (this->left->intersect_box(source, dir)) { //<========= Recursion Left
			vector<int> temp = this->left->ray_intersect(source, dir);
			out.insert(end(out), begin(temp), end(temp));
		}
		if (this->right->intersect_box(source, dir)) { //<========= Recursion Right
			vector<int> temp = this->right->ray_intersect(source, dir);
			out.insert(end(out), begin(temp), end(temp));
		}
		return out;
	}
}