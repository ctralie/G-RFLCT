from Primitives3D import *
import math

class EMMaterial(object):
	def EMMaterial(object):
		self.R = 1.0;
		self.T = 0.0;

class EMNode(object):
	def EMNode(self):
		self.parent = None
		self.mesh = None
		self.children = []
		self.transformation = Matrix4()#Make identity matrix by default
		self.material = EMMaterial()

class EMScene(object):
	def EMScene(self):
		


class AcousticNode { 
public:
	AcousticNode();
	AcousticNode* parent;
	vector<AcousticNode*> children;
	R3Shape* shape;
	//NOTE: mesh and intersectTriangle are only used if it's a mesh
	R3Mesh* mesh;
	R3MeshFace* intersectedTriangle;
	R4Matrix transformation;
	R3Material* material;
	AcousticMaterial* acousticMaterial;
	//AcousticMaterial acousticMaterial;
	R3Box bbox;
	R3Node* r3node;//Used for visual rendering
};

class AcousticScene {
public:
	AcousticScene();
	int Read(const char *filename, AcousticNode* root = NULL);
	void loadLights();
	void drawSource();
	void drawListener();
	void rayIntersectScene(R3Ray* ray, double* t, R3Point* position, R3Vector* normal, AcousticNode** intersectNode, R4Matrix* transform);
	void rayIntersectNode(R3Ray* ray, AcousticNode* node, R4Matrix transforms, double* t, R3Point* position,
		R3Vector* normal, AcousticNode** intersectNode, R4Matrix* transform);
	void transformMeshes(AcousticNode* node, R4Matrix transforms);

	vector<R3AcousticLight*> lights;
	AcousticNode* root;
	AcousticListener* listener;
	AcousticSource* source;
	R3Box bbox;
	RNRgb ambient;

	ImpulseResponse leftearResponse;
	ImpulseResponse rightearResponse;

    vector<R3Ray> debugRays;
    vector<R3MeshEdge*> diffEdges;
};
