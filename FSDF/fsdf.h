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
	// Members
	MatrixXd & vertices; //stored in the tree
	vector<int> idx;
	int nb_pts;
	int depth;
	MatrixXd box; //5 x 3, r(0) = corner, r(1) = side1, r(2) = side2, r(3) = side3, r(4) = evals; 
public:
	// Members
	unique_ptr<OBBNode> left;
	unique_ptr<OBBNode> right;

	// Special member functions, move semantic only
	OBBNode(MatrixXd & vertices);
	OBBNode(const OBBNode& other) = delete;
	OBBNode& operator=(const OBBNode& rhs) = delete;

	// Tree functions
	void compute_obb();
	void build_tree();
	void set_idx(vector<int> idx);
	void set_obb(MatrixXd box);
	MatrixXd get_obb() const;
};

class OBBTree {
private:
	// Members
	int max_depth; //of the tree
	int nb_leaves;
	double max_dist; //between two vertices linked by an edge
	MatrixXd pcl;
	MatrixXd normals;
	unique_ptr<OBBNode> root;

	// Inner functions
	vector<int> ray_intersect(const Vector3d & source, const Vector3d & dir) const;
public:
	// Special member functions, move semantic only
	OBBTree(MatrixXd pcl, MatrixXd normals); //Want copy of structures
	OBBTree(const OBBTree& other) = delete;
	OBBTree& operator=(const OBBTree& rhs) = delete;

	// User functions
	void build();
	VectorXd query() const;
};
