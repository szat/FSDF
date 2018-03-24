#include "fsdf.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Node
//////////////////////////////////////////////////////////////////////////////////////////////////////////

OBBNode::OBBNode() {
	this->left = nullptr;
	this->right = nullptr;
	this->nb_pts = 0;
	this->depth = 0;
	this->idx.clear();
}

void OBBNode::compute_obb() {
	this->box.row(0) << 0, 0, 0;
	this->box.row(1) << 0, 0, 0;
	this->box.row(2) << 0, 0, 0;
	this->box.row(3) << 0, 0, 0;
	this->box.row(4) << 0, 0, 0;
}

OBBNode::Box OBBNode::get_obb() const {
	return this->box;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tree
//////////////////////////////////////////////////////////////////////////////////////////////////////////

OBBTree::OBBTree(MatrixXd pcl, MatrixXd normals) {
	this->pcl = pcl; //I actually want a new copy of pcl
	this->normals = normals;
	this->max_depth = 4;
	this->max_dist = 0;
	this->nb_leaves = 0;
	this->root = unique_ptr<OBBNode>(new OBBNode()); //owning
}

void OBBTree::build_tree() {
	//one iteration
	root->left = unique_ptr<OBBNode>(new OBBNode());
	root->right = unique_ptr<OBBNode>(new OBBNode());
	//stopping condition: max_level or obb statistics
	return;
}

vector<int> OBBTree::ray_intersect(const Ref<const Vector3d> source, const Ref<const Vector3d> dir) const {
	vector<int> pt_list;
	return pt_list;
}

VectorXd OBBTree::query() const {
	VectorXd out;
	return out;
}