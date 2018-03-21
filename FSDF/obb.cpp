//#include "obb.h"
//
//OBBNode::OBBNode() {
//	this->Parent = nullptr;
//	this->Kids = nullptr;
//}
//
//OBBNode::~OBBNode() {
//	delete[] this->Kids;
//}
//
//OBBTree::OBBTree()
//{
//	this->DataSet = nullptr;
//	this->Level = 4;
//	this->MaxLevel = 12;
//	this->Automatic = 1;
//	this->Tolerance = 0.01;
//	this->Tree = nullptr;
//	this->PointsList = nullptr;
//	this->InsertedPoints = nullptr;
//	this->OBBCount = this->Level = 0;
//}
//
//OBBTree::~OBBTree()
//{
//	if (this->Tree)
//	{
//		this->DeleteTree(this->Tree);
//		delete this->Tree;
//		this->Tree = nullptr;
//	}
//}
//
//void OBBTree::DeleteTree(OBBNode *OBBptr)
//{
//	if (OBBptr->Kids != nullptr)
//	{
//		this->DeleteTree(OBBptr->Kids[0]);
//		this->DeleteTree(OBBptr->Kids[1]);
//		delete OBBptr->Kids[0];
//		delete OBBptr->Kids[1];
//	}
//}
//
//void ComputeOBB(Eigen::MatrixXd & V, double corner[3], double max[3], double mid[3], double min[3], double size[3]) {
//	1 + 1;
//}