#include "stdafx.h"
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

class OBBNode {
private:
	// "Global" members
	MatrixXd & vertices; 
	MatrixXd & normals;
	int param_max_level;
	double param_max_side;

	// "Local" members
	vector<int> idx;
	int nb_pts;
	int depth;
	int longest_side; //along which side it was split
	MatrixXd box; //5 x 3


public:
	unique_ptr<OBBNode> left;
	unique_ptr<OBBNode> right;

	// Special member functions, move semantic only
	OBBNode(MatrixXd & vertices, MatrixXd & normals);
	OBBNode(const OBBNode& other) = delete;
	OBBNode& operator=(const OBBNode& rhs) = delete;

	// Tree functions
	void init(int param_max_level, double param_max_side);
	void build_tree();

	void compute_obb();
	void set_idx(vector<int> idx);
	void set_obb(MatrixXd box);
	MatrixXd get_obb() const;
	int get_depth();
	bool is_leaf() const;
	bool in_box(const Vector3d & source) const;
	vector<int> ray_intersect(const Vector3d & source, const Vector3d & dir) const;
	bool intersect_box(const Vector3d & source, const Vector3d & dir) const;
	vector<int> validate() const;
	vector<int> obb_nbh(const Vector3d & source) const;
	MatrixXd query() const;
};
