//#pragma once
//
//#include <Eigen/Dense>
//
//class OBBNode { //;prevent man page generation
//public:
//	OBBNode();
//	~OBBNode();
//
//	double Corner[3]; 
//	double Axes[3][3]; //the axes defining the OBB - ordered from long->short
//	OBBNode *Parent; //parent node; nullptr if root
//	OBBNode **Kids; //two children of this node; nullptr if leaf
//
//private:
//	OBBNode(const OBBNode& other) = delete;
//	OBBNode& operator=(const OBBNode& rhs) = delete;
//};
//
//class OBBTree {
//public:
//	OBBTree();
//	~OBBTree();
//
//	void ComputeOBB(Eigen::MatrixXd & V, double corner[3], double max[3], double mid[3], double min[3], double size[3]);
//
//private:
//	OBBNode *Tree;
//	void BuildTree(vtkIdList *cells, OBBNode *parent, int level);
//	Eigen::MatrixXd V;
//	int *InsertedPoints;
//	int OBBCount;
//
//	void DeleteTree(OBBNode *OBBptr);
//
//private:
//	OBBTree(const OBBTree&) = delete;
//	void operator=(const OBBTree&) = delete;
//};
