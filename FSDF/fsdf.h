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
typedef Matrix<double, 5, 3> Box;

private:
	// Members
	vector<int> idx;
	int nb_pts;
	int depth;
	Box box;
	//box.row(0) = corner;
	//box.row(1) = side1; 
	//box.row(2) = side2; 
	//box.row(3) = side3; 
	//box.row(4) = evals; 

public:
	// Members
	unique_ptr<OBBNode> left;
	unique_ptr<OBBNode> right;

	// Special member functions, move semantic only
	OBBNode();
	OBBNode(const OBBNode& other) = delete;
	OBBNode& operator=(const OBBNode& rhs) = delete;

	// User functions
	void compute_obb(); 
	Box get_obb() const;
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
	vector<int> ray_intersect(const Ref<const Vector3d> source, const Ref<const Vector3d> dir) const;
public:
	// Special member functions, move semantic only
	OBBTree(MatrixXd pcl, MatrixXd normals);
	OBBTree(const OBBTree& other) = delete;
	OBBTree& operator=(const OBBTree& rhs) = delete;

	// User functions
	void build_tree();	
	VectorXd query() const;
};
