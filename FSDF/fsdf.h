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
	MatrixXd & vertices; //stored in the tree
	MatrixXd & normals;
	int max_level;

	// "Local" members
	vector<int> idx;
	int nb_pts;
	int depth;
	MatrixXd box; //5 x 3, r(0) = corner, r(1) = side1, r(2) = side2, r(3) = side3, r(4) = evals; 
	
	unique_ptr<OBBNode> left;
	unique_ptr<OBBNode> right;
public:

	// Special member functions, move semantic only
	OBBNode(MatrixXd & vertices, MatrixXd & normals);
	OBBNode(const OBBNode& other) = delete;
	OBBNode& operator=(const OBBNode& rhs) = delete;

	// Tree functions
	void compute_obb();
	void build_tree();
	void set_idx(vector<int> idx);
	void set_obb(MatrixXd box);
	MatrixXd get_obb() const;
	int get_depth();
	vector<int> ray_intersect(const Vector3d & source, const Vector3d & dir) const;
};
