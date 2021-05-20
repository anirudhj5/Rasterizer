#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <algorithm>
#include <vector>
#include <cmath>

#define NORMALS
#define M_PI 3.14159265358979323846

//My movie can be accessed here https://drive.google.com/file/d/1u9tM6Wfkhm57TYX4YcKnmsaUWhJ5RbJ5/view?usp=sharing

using std::cerr;
using std::endl;

double ceil__441(double f)
{
    return ceil(f-0.00001);
}

double floor__441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int height, int width)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Matrix
{
public:
	double          A[4][4];  // A[i][j] means row i, column j

	void            TransformPoint(const double* ptIn, double* ptOut);
	static Matrix   ComposeMatrices(const Matrix&, const Matrix&);
	void            Print(ostream& o);
};

void
Matrix::Print(ostream& o)
{
	for (int i = 0; i < 4; i++)
	{
		char str[256];
		sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
		o << str;
	}
}

void
Matrix::TransformPoint(const double* ptIn, double* ptOut)
{
	ptOut[0] = ptIn[0] * A[0][0]
		+ ptIn[1] * A[1][0]
		+ ptIn[2] * A[2][0]
		+ ptIn[3] * A[3][0];
	ptOut[1] = ptIn[0] * A[0][1]
		+ ptIn[1] * A[1][1]
		+ ptIn[2] * A[2][1]
		+ ptIn[3] * A[3][1];
	ptOut[2] = ptIn[0] * A[0][2]
		+ ptIn[1] * A[1][2]
		+ ptIn[2] * A[2][2]
		+ ptIn[3] * A[3][2];
	ptOut[3] = ptIn[0] * A[0][3]
		+ ptIn[1] * A[1][3]
		+ ptIn[2] * A[2][3]
		+ ptIn[3] * A[3][3];
}

Matrix
Matrix::ComposeMatrices(const Matrix& M1, const Matrix& M2)
{
	Matrix rv;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			rv.A[i][j] = 0;
			for (int k = 0; k < 4; k++)
				rv.A[i][j] += M1.A[i][k] * M2.A[k][j];
		}

	return rv;
}

double dot_product(std::vector<double> v1, std::vector<double> v2) {
	if (v1.size() != v2.size()) {
		cout << "ERROR INCOMPATIBLE VECTORS" << endl;
		cout << v1.size() << "     :      " << v2.size() << endl;
		return -1.0;
	}
	else {
		double val = 0;
		for (int i = 0; i < v1.size(); i++) {
			val += v1[i] * v2[i];
		}
		return val;
	}
}

class Camera
{
public:
	double          near, far;
	double          angle;
	double          position[3];
	double          focus[3];
	double          up[3];

	Matrix          ViewTransform(double alpha, double n, double f);
	Matrix          CameraTransform(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> t);
	Matrix          DeviceTransform(int n, int m);
};

Matrix Camera::ViewTransform(double alpha, double n, double f) {
	Matrix view;
	view.A[0][0] = 1 / tan(alpha/2);
	view.A[0][1] = view.A[0][2] = view.A[0][3] = view.A[1][0] = view.A[1][2] = view.A[1][3] = view.A[2][0] = view.A[2][1] = view.A[3][0] = view.A[3][1] = view.A[3][3] = 0;
	view.A[1][1] = 1 / tan(alpha/2); 
	view.A[2][2] = (f + n) / (f - n);
	view.A[2][3] = -1;
	view.A[3][2] = 2 * f * n / (f - n);
	return view;
}

Matrix Camera::CameraTransform(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> t) {
	Matrix camera;
	for (int i = 0; i < 3; i++) {
		camera.A[i][0] = v1[i];
		camera.A[i][1] = v2[i];
		camera.A[i][2] = v3[i];
		camera.A[i][3] = 0;
	}
	camera.A[3][0] = dot_product(v1, t);
	camera.A[3][1] = dot_product(v2, t);
	camera.A[3][2] = dot_product(v3, t);
	camera.A[3][3] = 1;
	return camera;
}

Matrix Camera::DeviceTransform(int n, int m) {
	Matrix device;
	device.A[0][0] = n / 2;
	device.A[0][1] = device.A[0][2] = device.A[0][3] = device.A[1][0] = device.A[1][2] = device.A[1][3] = device.A[2][0] = device.A[2][1] =  device.A[2][3] = device.A[3][2] = 0;
	device.A[1][1] = m / 2;
	device.A[2][2] = device.A[3][3] = 1;
	device.A[3][0] = n / 2;
	device.A[3][1] = m / 2;
	return device;
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
	int nNonRamp = nFrames - 2 * ramp;
	double height = 1. / (nNonRamp + 4 * ramp / M_PI);
	if (curFrame < ramp)
	{
		double factor = 2 * height * ramp / M_PI;
		double eval = cos(M_PI / 2 * ((double)curFrame) / ramp);
		return (1. - eval) * factor;
	}
	else if (curFrame > nFrames - ramp)
	{
		int amount_left = nFrames - curFrame;
		double factor = 2 * height * ramp / M_PI;
		double eval = cos(M_PI / 2 * ((double)amount_left / ramp));
		return 1. - (1 - eval) * factor;
	}
	double amount_in_quad = ((double)curFrame - ramp);
	double quad_part = amount_in_quad * height;
	double curve_part = height * (2 * ramp) / M_PI;
	return quad_part + curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
	double t = SineParameterize(frame, nframes, nframes / 10);
	Camera c;
	c.near = 5;
	c.far = 200;
	c.angle = M_PI / 6;
	c.position[0] = 40 * sin(2 * M_PI * t);
	c.position[1] = 40 * cos(2 * M_PI * t);
	c.position[2] = 40;
	c.focus[0] = 0;
	c.focus[1] = 0;
	c.focus[2] = 0;
	c.up[0] = 0;
	c.up[1] = 1;
	c.up[2] = 0;
	return c;
}

struct LightingParameters
{
	LightingParameters(void)
	{
		lightDir[0] = -0.6;
		lightDir[1] = 0;
		lightDir[2] = -0.8;
		Ka = 0.3;
		Kd = 0.7;
		Ks = 2.8;
		alpha = 50.5;
	};

	double lightDir[3]; // The direction of the light source
	double Ka;          // The coefficient for ambient lighting
	double Kd;          // The coefficient for diffuse lighting
	double Ks;          // The coefficient for specular lighting
	double alpha;       // The exponent term for specular lighting
}; LightingParameters lp;

LightingParameters
GetLighting(Camera c)
{
	LightingParameters lp;
	lp.lightDir[0] = c.position[0] - c.focus[0];
	lp.lightDir[1] = c.position[1] - c.focus[1];
	lp.lightDir[2] = c.position[2] - c.focus[2];
	double mag = sqrt(lp.lightDir[0] * lp.lightDir[0]
		+ lp.lightDir[1] * lp.lightDir[1]
		+ lp.lightDir[2] * lp.lightDir[2]);
	if (mag > 0)
	{
		lp.lightDir[0] /= mag;
		lp.lightDir[1] /= mag;
		lp.lightDir[2] /= mag;
	}

	return lp;
}

std::vector<double> cross_product(std::vector<double> v1, std::vector<double> v2) {
	std::vector<double> v3(3);
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

	return v3;
}

std::vector<double> normalize(std::vector<double> v1) {
	double c = 0;
	for (int i = 0; i < v1.size(); i++) {
		c += v1[i] * v1[i];
	}
	double v = sqrt(c);
	std::vector<double> normal(v1.size());
	for (int j = 0; j < normal.size(); j++) {
		normal[j] = (v1[j] / v);
	}
	return normal;
}

double calc_shading(LightingParameters lp, std::vector<double> viewDir, double *normals) {
	std::vector<double> v1(3);
	std::vector<double> v2(3);
	std::vector<double> v4(3);
	for (int i = 0; i < 3; i++) {
		v1[i] = lp.lightDir[i];
		v2[i] = normals[i];
	}
	double diff = dot_product(v1, v2);
	double tot_spec;
	for (int j = 0; j < 3; j++) {
		v4[j] = 2 * diff * v2[j] - v1[j];
	}
	double spec = dot_product(v4, viewDir);
	if (isnan(lp.Ks * pow(spec, lp.alpha))) {
		tot_spec = 0;
	}else{
		tot_spec = lp.Ks * pow(spec, lp.alpha);
	}
	if (diff < 0) {
		diff = 0;
	}
	return lp.Ka + lp.Kd*diff + tot_spec;
}

double lin_interp(double fA, double fB, double A, double B, double x) {
	double t = (x - A) / (B - A);
	return fA + t * (fB - fA);
}

class Triangle
{
public:
	double         X[3];
	double         Y[3];
	double         Z[3];
	double colors[3][3];
	double normals[3][3];
	double shading[3];

	int isRight() {
		int left = 0;
		double val = X[0];
		for (int i = 1; i < 3; i++) {
			if (X[i] < val) {
				left = i;
				val = X[i];
			}
		}

		std::vector<int> v;
		for (int i = 0; i < 3; i++) {
			if (i != left) {
				v.push_back(i);
			}
		}
		
		for (int i = 0; i < 2; i++) {
			if (X[left] == X[v[i]]) {
				return 1;
			}
		}
		return 0;

	}

	int isArb() {
		double val = X[0];
		for (int j = 1; j < 3; j++) {
			if (X[j] == val) {
				return 0;
			}

		}
		if (X[1] == X[2]) {
			return 0;
		}
		return 1;
	}

	int find_mid() {
		int left = 0;
		double val1 = X[0];
		for (int j = 1; j < 3; j++) {
			if (X[j] < val1) {
				left = j;
				val1 = X[j];
			}
		}
	
		int right = 0;
		double val2 = X[0];
		for (int j = 1; j < 3; j++) {
			if (X[j] > val2) {
				right = j;
				val2 = X[j];
			}
		}
		return 3 - left - right;
	}

	int get_point(int is_Right) {
		int point = 0;
		double val = X[0];
		for (int i = 1; i < 3; i++) {
			if (is_Right == 1) {
				if (X[i] > val) {
					point = i;
					val = X[i];
				}
			}
			else {
				if (X[i] < val) {
					point = i;
					val = X[i];
				}
			}
		}
		return point;
	}

	int get_bot(int direction) {
		int bot = 0;
		double val2 = 2000;
		for (int j = 0; j < 3; j++) {
			if (j != get_point(direction)) {
				if (Y[j] < val2) {
					bot = j;
					val2 = Y[j];
				}
			}
		}
		return bot;
	}

	int get_top(int direction) {
		int top = 0;
		double val2 = 0;
		for (int j = 0; j < 3; j++) {
			if (j != get_point(direction)) {
				if (Y[j] > val2) {
					top = j;
					val2 = Y[j];
				}
			}
		}
		return top;
	}

	double slope(double x, int p1, int p2) {
		double m = (Y[p2] - Y[p1]) / (X[p2] - X[p1]);
		double b = Y[p2] - m * X[p2];
		return m * x + b;
	}
};

class Screen
{
public:
	unsigned char* buffer;
	int width, height;
	// would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
	vtkPolyDataReader* rdr = vtkPolyDataReader::New();
	rdr->SetFileName("proj1f_geometry.vtk");
	cerr << "Reading" << endl;
	rdr->Update();
	cerr << "Done reading" << endl;
	if (rdr->GetOutput()->GetNumberOfCells() == 0)
	{
		cerr << "Unable to open file!!" << endl;
		exit(EXIT_FAILURE);
	}
	vtkPolyData* pd = rdr->GetOutput();

	int numTris = pd->GetNumberOfCells();
	vtkPoints* pts = pd->GetPoints();
	vtkCellArray* cells = pd->GetPolys();
	vtkDoubleArray* var = (vtkDoubleArray*)pd->GetPointData()->GetArray("hardyglobal");
	double* color_ptr = var->GetPointer(0);
	//vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
	//float *color_ptr = var->GetPointer(0);
	vtkFloatArray* n = (vtkFloatArray*)pd->GetPointData()->GetNormals();
	float* normals = n->GetPointer(0);
	std::vector<Triangle> tris(numTris);
	vtkIdType npts;
	vtkIdType* ptIds;
	int idx;
	for (idx = 0, cells->InitTraversal(); cells->GetNextCell(npts, ptIds); idx++)
	{
		if (npts != 3)
		{
			cerr << "Non-triangles!! ???" << endl;
			exit(EXIT_FAILURE);
		}
		double* pt = NULL;
		pt = pts->GetPoint(ptIds[0]);
		tris[idx].X[0] = pt[0];
		tris[idx].Y[0] = pt[1];
		tris[idx].Z[0] = pt[2];
#ifdef NORMALS
		tris[idx].normals[0][0] = normals[3 * ptIds[0] + 0];
		tris[idx].normals[0][1] = normals[3 * ptIds[0] + 1];
		tris[idx].normals[0][2] = normals[3 * ptIds[0] + 2];
#endif
		pt = pts->GetPoint(ptIds[1]);
		tris[idx].X[1] = pt[0];
		tris[idx].Y[1] = pt[1];
		tris[idx].Z[1] = pt[2];
#ifdef NORMALS
		tris[idx].normals[1][0] = normals[3 * ptIds[1] + 0];
		tris[idx].normals[1][1] = normals[3 * ptIds[1] + 1];
		tris[idx].normals[1][2] = normals[3 * ptIds[1] + 2];
#endif
		pt = pts->GetPoint(ptIds[2]);
		tris[idx].X[2] = pt[0];
		tris[idx].Y[2] = pt[1];
		tris[idx].Z[2] = pt[2];
#ifdef NORMALS
		tris[idx].normals[2][0] = normals[3 * ptIds[2] + 0];
		tris[idx].normals[2][1] = normals[3 * ptIds[2] + 1];
		tris[idx].normals[2][2] = normals[3 * ptIds[2] + 2];
#endif

		// 1->2 interpolate between light blue, dark blue
		// 2->2.5 interpolate between dark blue, cyan
		// 2.5->3 interpolate between cyan, green
		// 3->3.5 interpolate between green, yellow
		// 3.5->4 interpolate between yellow, orange
		// 4->5 interpolate between orange, brick
		// 5->6 interpolate between brick, salmon
		double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
		double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
		unsigned char RGB[8][3] = { { 71, 71, 219 },
									{ 0, 0, 91 },
									{ 0, 255, 255 },
									{ 0, 128, 0 },
									{ 255, 255, 0 },
									{ 255, 96, 0 },
									{ 107, 0, 0 },
									{ 224, 76, 76 }
		};
		for (int j = 0; j < 3; j++)
		{
			float val = color_ptr[ptIds[j]];
			int r;
			for (r = 0; r < 7; r++)
			{
				if (mins[r] <= val && val < maxs[r])
					break;
			}
			if (r == 7)
			{
				cerr << "Could not interpolate color for " << val << endl;
				exit(EXIT_FAILURE);
			}
			double proportion = (val - mins[r]) / (maxs[r] - mins[r]);
			tris[idx].colors[j][0] = (RGB[r][0] + proportion * (RGB[r + 1][0] - RGB[r][0])) / 255.0;
			tris[idx].colors[j][1] = (RGB[r][1] + proportion * (RGB[r + 1][1] - RGB[r][1])) / 255.0;
			tris[idx].colors[j][2] = (RGB[r][2] + proportion * (RGB[r + 1][2] - RGB[r][2])) / 255.0;
		}
	}

	return tris;
}


void fill_right(Triangle t1, unsigned char* buffer, double* zbuffer) {
	int dir = 1;
	double tri_left = t1.X[t1.get_top(dir)];
	double tri_right = t1.X[t1.get_point(dir)];
	double z_val_left = t1.Z[t1.get_bot(dir)];
	double z_val_right = t1.Z[t1.get_point(dir)];
	double shading_left = t1.shading[t1.get_bot(dir)];
	double shading_right = t1.shading[t1.get_point(dir)];
	for (double j = std::max(ceil__441(tri_left), 0.0); j <= std::min(floor__441(tri_right), 999.0); j += 1.0f) {
		double tri_bot = t1.slope(j, t1.get_bot(dir), t1.get_point(dir));
		double tri_top = t1.slope(j, t1.get_top(dir), t1.get_point(dir));
		double z_val_bot = lin_interp(z_val_left, z_val_right, tri_left, tri_right, j);
		double z_val_top = lin_interp(t1.Z[t1.get_top(dir)], z_val_right, tri_left, tri_right, j);
		double shading_bot = lin_interp(shading_left, shading_right, tri_left, tri_right, j);
		double shading_top = lin_interp(t1.shading[t1.get_top(dir)], shading_right, tri_left, tri_right, j);
		double c_bot[3];
		for (int i = 0; i < 3; i++) {
			c_bot[i] = lin_interp(t1.colors[t1.get_bot(dir)][i], t1.colors[t1.get_point(dir)][i], tri_left, tri_right, j);
		}

		double c_top[3];
		for (int i = 0; i < 3; i++) {
			c_top[i] = lin_interp(t1.colors[t1.get_top(dir)][i], t1.colors[t1.get_point(dir)][i], tri_left, tri_right, j);
		}

		for (double k = std::max(ceil__441(tri_bot), 0.0); k <= std::min(floor__441(tri_top), 999.0); k += 1.0f) {
			double zval = lin_interp(z_val_bot, z_val_top, tri_bot, tri_top, k);
			double shading = lin_interp(shading_bot, shading_top, tri_bot, tri_top, k);
			double colors[3];
			for (int i = 0; i < 3; i++) {
				colors[i] = lin_interp(c_bot[i], c_top[i], tri_bot, tri_top, k);
			}
			int index = k * 1000 + j;
			if (zval > zbuffer[index]) {
				buffer[3 * index] = (unsigned char)ceil__441(std::max(0.0, std::min(1.0, shading*colors[0]))*255);
				buffer[3 * index + 1] = (unsigned char)ceil__441(std::max(0.0, std::min(1.0, shading*colors[1]))*255);
				buffer[3 * index + 2] = (unsigned char)ceil__441(std::max(0.0, std::min(1.0, shading*colors[2]))*255);
				zbuffer[index] = zval;
			}		
		}
	}
}

void fill_left(Triangle t1, unsigned char* buffer, double* zbuffer) {
	int dir = 0;
	double tri_left = t1.X[t1.get_point(dir)];
	double tri_right = t1.X[t1.get_top(dir)];
	double z_val_left = t1.Z[t1.get_point(dir)];
	double z_val_right = t1.Z[t1.get_bot(dir)];
	double shading_left = t1.shading[t1.get_point(dir)];
	double shading_right = t1.shading[t1.get_bot(dir)];
	for (double j = std::max(ceil__441(tri_left), 0.0); j <= std::min(floor__441(tri_right), 999.0); j += 1.0f) {
		double tri_bot = t1.slope(j, t1.get_bot(dir), t1.get_point(dir));
		double tri_top = t1.slope(j, t1.get_top(dir), t1.get_point(dir));
		double z_val_bot = lin_interp(z_val_left, z_val_right, tri_left, tri_right, j);
		double z_val_top = lin_interp(z_val_left, t1.Z[t1.get_top(dir)], tri_left, tri_right, j);
		double shading_bot = lin_interp(shading_left, shading_right, tri_left, tri_right, j);
		double shading_top = lin_interp(shading_left, t1.shading[t1.get_top(dir)], tri_left, tri_right, j);
		double c_bot[3];
		for (int i = 0; i < 3; i++) {
			c_bot[i] = lin_interp(t1.colors[t1.get_point(dir)][i], t1.colors[t1.get_bot(dir)][i], tri_left, tri_right, j);
		}

		double c_top[3];
		for (int i = 0; i < 3; i++) {
			c_top[i] = lin_interp(t1.colors[t1.get_point(dir)][i], t1.colors[t1.get_top(dir)][i], tri_left, tri_right, j);
		}

		for (double k = std::max(ceil__441(tri_bot), 0.0); k <= std::min(floor__441(tri_top), 999.0); k += 1.0f) {
			double zval = lin_interp(z_val_bot, z_val_top, tri_bot, tri_top, k);
			double shading = lin_interp(shading_bot, shading_top, tri_bot, tri_top, k);
			double colors[3];
			for (int i = 0; i < 3; i++) {
				colors[i] = lin_interp(c_bot[i], c_top[i], tri_bot, tri_top, k);
			}
			
			int index = k * 1000 + j;
			
			

			if (zval > zbuffer[index]) {
				buffer[3 * index] = (unsigned char)ceil__441(std::max(0.0, std::min(1.0, shading*colors[0]))*255);
				buffer[3 * index + 1] = (unsigned char)ceil__441(std::max(0.0, std::min(1.0, shading*colors[1]))*255);
				buffer[3 * index + 2] = (unsigned char)ceil__441(std::max(0.0, std::min(1.0, shading*colors[2]))*255);
				zbuffer[index] = zval;
			}
		}
	}
}


void render_triangles(std::vector<Triangle> triangles, unsigned char* buffer, double* zbuffer) {
	for (int i = 0; i < triangles.size(); i++) {
		if (triangles[i].isArb()) {

			int m = triangles[i].find_mid();
			int l = triangles[i].get_point(0);
			int r = triangles[i].get_point(1);
			double bot_y = triangles[i].slope(triangles[i].X[m], l, r);
			Triangle t1;
			Triangle t2;
			t1.X[0] = triangles[i].X[l];
			t1.Y[0] = triangles[i].Y[l];
			t1.Z[0] = triangles[i].Z[l];
			t1.shading[0] = triangles[i].shading[l];
			
			for (int q = 0; q < 3; q++) {
				t1.colors[0][q] = triangles[i].colors[l][q];
			}
			t1.X[1] = triangles[i].X[m];
			t1.Y[1] = triangles[i].Y[m];
			t1.Z[1] = triangles[i].Z[m];
			t1.shading[1] = triangles[i].shading[m];
			
			for (int w = 0; w < 3; w++) {
				t1.colors[1][w] = triangles[i].colors[m][w];
			}
			t1.X[2] = triangles[i].X[m];
			t1.Y[2] = bot_y;
			t1.Z[2] = lin_interp(triangles[i].Z[l], triangles[i].Z[r], triangles[i].X[l], triangles[i].X[r], triangles[i].X[m]);
			t1.shading[2] = lin_interp(triangles[i].shading[l], triangles[i].shading[r], triangles[i].X[l], triangles[i].X[r], triangles[i].X[m]);
			
			for (int p = 0; p < 3; p++) {
				t1.colors[2][p] = lin_interp(triangles[i].colors[l][p], triangles[i].colors[r][p], triangles[i].X[l], triangles[i].X[r], triangles[i].X[m]);
			}

			t2.X[0] = triangles[i].X[r];
			t2.Y[0] = triangles[i].Y[r];
			t2.Z[0] = triangles[i].Z[r];
			t2.shading[0] = triangles[i].shading[r];
			
			for (int q = 0; q < 3; q++) {
				t2.colors[0][q] = triangles[i].colors[r][q];
			}
			t2.X[1] = triangles[i].X[m];
			t2.Y[1] = triangles[i].Y[m];
			t2.Z[1] = triangles[i].Z[m];
			t2.shading[1] = triangles[i].shading[m];
			for (int w = 0; w < 3; w++) {
				t2.colors[1][w] = triangles[i].colors[m][w];
			}
			t2.X[2] = triangles[i].X[m];
			t2.Y[2] = bot_y;
			t2.Z[2] = lin_interp(triangles[i].Z[l], triangles[i].Z[r], triangles[i].X[l], triangles[i].X[r], triangles[i].X[m]);
			t2.shading[2] = lin_interp(triangles[i].shading[l], triangles[i].shading[r], triangles[i].X[l], triangles[i].X[r], triangles[i].X[m]);
			for (int p = 0; p < 3; p++) {
				t2.colors[2][p] = lin_interp(triangles[i].colors[l][p], triangles[i].colors[r][p], triangles[i].X[l], triangles[i].X[r], triangles[i].X[m]);
			}

			fill_left(t1, buffer, zbuffer);
			fill_right(t2, buffer, zbuffer);
		}
		else {
			if (triangles[i].isRight()) {
				fill_right(triangles[i], buffer, zbuffer);

			}
			else {
				fill_left(triangles[i], buffer, zbuffer);
			}
		}
	}
}

int main()
{
	
	std::vector<Triangle> triangles = GetTriangles();
    for (int i = 0; i < 1000; i++) {

		vtkImageData* image = NewImage(1000, 1000);
		unsigned char* buffer =
			(unsigned char*)image->GetScalarPointer(0, 0, 0);
		double* zbuffer;
		zbuffer = new double[1000 * 1000];
		int npixels = 1000 * 1000;
	    for (int i = 0; i < npixels * 3; i++) {
		    if (i < npixels) {
			    zbuffer[i] = -1.0;
		    }
		    buffer[i] = 0;
	    }

	    Screen screen;
	    screen.buffer = buffer;
	    screen.width = 1000;
	    screen.height = 1000;
		
		Camera c1 = GetCamera(i, 1000);
		LightingParameters lp = GetLighting(c1);

		std::vector<double> O(3);
		std::vector<double> W(3);
		std::vector<double> Up(3);
		std::vector<double> d(3);
		
		for (int j = 0; j < 3; j++) {
			O[j] = c1.position[j];
			d[j] = -1.0 * c1.position[j];
			Up[j] = c1.up[j];
		}
		for (int y = 0; y < 3; y++) {
			W[y] = O[y] - c1.focus[y];
		}
		
		W = normalize(W);
		std::vector<double> U = cross_product(Up, W);
		U = normalize(U);
		std::vector<double> V = cross_product(W, U);
		V = normalize(V);
		
		Matrix V_MATRIX = c1.ViewTransform(c1.angle, c1.near, c1.far);
		Matrix C_MATRIX = c1.CameraTransform(U, V, W, d);
		Matrix D_MATRIX = c1.DeviceTransform(1000, 1000);
		Matrix M = Matrix::ComposeMatrices(Matrix::ComposeMatrices(C_MATRIX, V_MATRIX), D_MATRIX);
		std::vector<Triangle> tris(triangles);

		for (int p = 0; p < tris.size(); p++) {
			for (int t = 0; t < 3; t++) {
				double ptin[4] = { tris[p].X[t], tris[p].Y[t],tris[p].Z[t], 1};
				double ptout[4];

				std::vector<double> viewDir(3);
				viewDir[0] = c1.position[0] - tris[p].X[t];
				viewDir[1] = c1.position[1] - tris[p].Y[t];
				viewDir[2] = c1.position[2] - tris[p].Z[t];
				viewDir = normalize(viewDir);

				M.TransformPoint(ptin, ptout);
				tris[p].X[t] = ptout[0]/ptout[3];
				tris[p].Y[t] = ptout[1]/ptout[3];
				tris[p].Z[t] = ptout[2]/ptout[3];
				tris[p].shading[t] = std::max(0.0, calc_shading(lp, viewDir, triangles[p].normals[t]));
			}
		}
		render_triangles(tris, buffer, zbuffer);
		char image_name[10];
		sprintf(image_name, "frame%03d", i);
		WriteImage(image, image_name);
   }
}
