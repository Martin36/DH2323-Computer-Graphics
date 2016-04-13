#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <algorithm>

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::vec2;
// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

SDL_Surface* screen;
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
int t;
vector<Triangle> triangles;
vec3 cameraPos(0, 0, -2.001);
float focalLength = SCREEN_WIDTH / 2;
mat3 R(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
float yaw = 0;  // Yaw angle controlling camera rotation around y-axis
float dyaw = M_PI / 40;
vec3 rightDir(R[0][0], R[0][1], R[0][2]);
vec3 downDir(R[1][0], R[1][1], R[1][2]);
vec3 forwardDir(R[2][0], R[2][1], R[2][2]);
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

vec3 lightPos(0, -0.5, -0.9);
vec3 lightPower = 10.1f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
vec2 reflectance(0.5, 0.5);

vec3 currentNormal;
vec3 currentReflectance;
// ----------------------------------------------------------------------------

struct Pixel
{
	int x;
	int y;
	float zinv;
	vec3 illumination;
	vec3 pos3d;
};

struct Vertex
{
	vec3 position;
	vec3 normal;
	vec3 reflectance;
};
// FUNCTIONS

void Update();
void Draw();

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b - a) / float(glm::max(N - 1, 1));
	vec2 current(a);
	for (int i = 0; i < N; ++i)
	{
		result[i] = current;
		current += step;
	}
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
	int N = result.size();
	vec3 step((b.x - a.x), (b.y - a.y), (b.zinv - a.zinv));
	step /= float(glm::max(N - 1, 1));
	vec3 current(a.x, a.y, a.zinv);
//	vec3 illuminationStep = (b.illumination - a.illumination) / float(glm::max(N - 1, 1));
//	vec3 currentIllumination = a.illumination;
	vec3 pos3DStep = (b.pos3d - a.pos3d) / float(glm::max(N - 1, 1));
	vec3 current3DPos = a.pos3d;

	for (int i = 0; i < N; i++){
		result[i].x = current.x;
		result[i].y = current.y;
		result[i].zinv = current.z;
		result[i].pos3d = current3DPos;
//		result[i].illumination = currentIllumination;
		current += step;
		current3DPos += pos3DStep;
//		currentIllumination += illuminationStep;
	}
}

void VertexShader(const vec3& v, ivec2& p){
	vec3 P = (v - cameraPos)*R;
	p.x = (focalLength * P.x / P.z) + SCREEN_WIDTH / 2;
	p.y = (focalLength * P.y / P.z) + SCREEN_HEIGHT / 2;
}

void VertexShader(const vec3& v, Pixel& p){
	vec3 P = (v - cameraPos)*R;
	p.zinv = 1.0f / abs(P.z);
	p.x = (focalLength * P.x / P.z) + SCREEN_WIDTH / 2;
	p.y = (focalLength * P.y / P.z) + SCREEN_HEIGHT / 2;
}
//Transforms the vertex to a pixel
void VertexShader(const Vertex& v, Pixel& p){
	/*
	//Task 6.1
	//Compute illumination
	vec3 r(lightPos - v.position);
	float distanceToLight = glm::distance(lightPos, v.position);
	float r2 = glm::dot(r, r);
	vec3 dirToLight = glm::normalize(r);
	vec3 n = v.normal;
	float maxrn = glm::max(glm::dot(dirToLight, n), 0.0f);
	vec3 D = (lightPower * maxrn) / float(4 * M_PI*r2);
	vec3 illumination = v.reflectance * (D + indirectLightPowerPerArea);
	p.illumination = illumination;
	*/
	vec3 P = (v.position - cameraPos)*R;
	p.zinv = 1.0f / abs(P.z);
	p.x = (focalLength * P.x / P.z) + SCREEN_WIDTH / 2;
	p.y = (focalLength * P.y / P.z) + SCREEN_HEIGHT / 2;
	p.pos3d = v.position;
}
void PixelShader(const Pixel& p)
{
	int x = p.x;
	int y = p.y;
	
	if (p.zinv > depthBuffer[y][x])
	{
		//Compute illumination
		vec3 r(lightPos - p.pos3d);
		float distanceToLight = glm::distance(lightPos, p.pos3d);
		float r2 = glm::dot(r, r);
		vec3 dirToLight = glm::normalize(r);
		vec3 n = currentNormal;
		float maxrn = glm::max(glm::dot(dirToLight, n), 0.0f);
		vec3 D = (lightPower * maxrn) / float(4 * M_PI*r2);
		vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea);
//		p.illumination = illumination;

		depthBuffer[y][x] = p.zinv;
		PutPixelSDL(screen, x, y, illumination);
	}
	
}
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color){
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	vector<ivec2> line(pixels);
	Interpolate(a, b, line);
	for (int i = 0; i < pixels; i++){
		ivec2 pos = line[i];
		PutPixelSDL(surface, pos.x, pos.y, color);
	}
}

void DrawPolygonEdges(const vector<vec3>& vertices)
{
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<ivec2> projectedVertices(V);
	for (int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}
	// Loop over all vertices and draw the edge from it to the next vertex:
	for (int i = 0; i < V; ++i)
	{
		int j = (i + 1) % V; // The next vertex
		vec3 color(1, 1, 1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j],
			color);
	}
}

void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels,
	vector<ivec2>& rightPixels)
{
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int N = vertexPixels.size();
	vector<int> xValues(N);
	vector<int> yValues(N);
	for (int i = 0; i < N; i++){
		xValues[i] = vertexPixels[i].x;
		yValues[i] = vertexPixels[i].y;
	}
	auto yMinMax = minmax_element(yValues.begin(), yValues.end());		//returns min and max value as a pair<ForwardIterator, ForwardIterator>
	int yMin = *yMinMax.first;
	int yMax = *yMinMax.second;
	int nrOfRows = yMax - yMin + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels.resize(nrOfRows);
	rightPixels.resize(nrOfRows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < nrOfRows; i++)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
		leftPixels[i].y = i + yMin;
		rightPixels[i].y = i + yMin;
	}
	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (int i = 0; i < N; i++){
		ivec2 start = vertexPixels[i];
		ivec2 end = vertexPixels[(i + 1) % N];
		int yDist = abs(end.y - start.y);
		vector<ivec2> edge(yDist + 1);
		Interpolate(start, end, edge);
		for (int j = 0; j <= yDist; j++){
			//Loop through the points on the egde
			int x = edge[j].x;
			int y = edge[j].y;
			int ind = y - yMin;			//y - yMin is the position corresponding to the coordinate y
			leftPixels[ind].x = (leftPixels[ind].x > x) ? x : leftPixels[ind].x;
			rightPixels[ind].x = (rightPixels[ind].x < x) ? x : rightPixels[ind].x;
		}

	}

}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels,
	vector<Pixel>& rightPixels)
{
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int N = vertexPixels.size();
	vector<int> xValues(N);
	vector<int> yValues(N);
	for (int i = 0; i < N; i++){
		xValues[i] = vertexPixels[i].x;
		yValues[i] = vertexPixels[i].y;
	}
	auto yMinMax = minmax_element(yValues.begin(), yValues.end());		//returns min and max value as a pair<ForwardIterator, ForwardIterator>
	int yMin = *yMinMax.first;
	int yMax = *yMinMax.second;
	int nrOfRows = yMax - yMin + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels.resize(nrOfRows);
	rightPixels.resize(nrOfRows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < nrOfRows; i++)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
		leftPixels[i].y = i + yMin;
		rightPixels[i].y = i + yMin;
	}
	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (int i = 0; i < N; i++){
		Pixel start = vertexPixels[i];
		Pixel end = vertexPixels[(i + 1) % N];
		int yDist = abs(end.y - start.y);
		vector<Pixel> edge(yDist + 1);

		//Interpolating p/z instead of just p
		start.pos3d *= start.zinv;
		end.pos3d *= end.zinv;
		Interpolate(start, end, edge);
		start.pos3d /= start.zinv;
		end.pos3d /= end.zinv;
		for (int j = 0; j < edge.size(); j++){
			edge[j].pos3d /= edge[j].zinv;
		}

		for (int j = 0; j <= yDist; j++){
			//Loop through the points on the egde
			int ind = edge[j].y - yMin;			//y - yMin is the position corresponding to the coordinate y
			if (leftPixels[ind].x > edge[j].x){
				leftPixels[ind] = edge[j];
			}
			if (rightPixels[ind].x < edge[j].x){
				rightPixels[ind] = edge[j];
			}

		}

	}

}

void DrawRows(const vector<ivec2>& leftPixels,
	const vector<ivec2>& rightPixels)
{
	for (int i = 0; i < leftPixels.size(); i++){
		vector<ivec2> row(rightPixels[i].x - leftPixels[i].x + 1);
		Interpolate(leftPixels[i], rightPixels[i], row);
		for (int j = 0; j < row.size(); j++){
			PutPixelSDL(screen, row[j].x, row[j].y, currentColor);
			
		}
	}
}

void DrawRows(const vector<Pixel>& leftPixels,
	const vector<Pixel>& rightPixels)
{
	for (int i = 0; i < leftPixels.size(); i++){
		vector<Pixel> row(rightPixels[i].x - leftPixels[i].x + 1);
		Interpolate(leftPixels[i], rightPixels[i], row);
		for (int j = 0; j < row.size(); j++){
			PixelShader(row[j]);
			/*
			if (row[j].zinv >= depthBuffer[row[j].y][row[j].x]){		//Its in front of current pixel at that point
				depthBuffer[row[j].y][row[j].x] = row[j].zinv;
				PutPixelSDL(screen, row[j].x, row[j].y, currentColor);
			}
			*/
		}
	}
}

//Using ivec2 as coordinate	
void DrawPolygonVec2(const vector<vec3>& vertices)
{
	int V = vertices.size();
	vector<ivec2> vertexPixels(V);
	for (int i = 0; i<V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(leftPixels, rightPixels);
}


// Using Pixel as coordinate

void DrawPolygon(const vector<vec3>& vertices)
{
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for (int i = 0; i<V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(leftPixels, rightPixels);
}

void DrawPolygon(const vector<Vertex>& vertices)
{
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for (int i = 0; i<V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(leftPixels, rightPixels);
}

int main(int argc, char* argv[])
{
	LoadTestModel(triangles);
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.
	/*		Test code for Task 4
	vector<ivec2> vertexPixels(3);
	vertexPixels[0] = ivec2(10, 5);
	vertexPixels[1] = ivec2(5, 10);
	vertexPixels[2] = ivec2(15, 15);
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for (int row = 0; row<leftPixels.size(); ++row)
	{
		cout << "Start: ("
			<< leftPixels[row].x << ","
			<< leftPixels[row].y << "). "
			<< "End: ("
			<< rightPixels[row].x << ","
			<< rightPixels[row].y << "). " << endl;
	}
	*/
	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState(0);

	if (keystate[SDLK_UP])
		cameraPos += forwardDir / 10.f;


	if (keystate[SDLK_DOWN])
		cameraPos -= forwardDir / 10.f;

	if (keystate[SDLK_RIGHT])
		cameraPos += rightDir / 10.f;

	if (keystate[SDLK_LEFT])
		cameraPos -= rightDir / 10.f;

	if (keystate[SDLK_RSHIFT])
		;

	if (keystate[SDLK_RCTRL])
		;

	if (keystate[SDLK_w])
		lightPos += forwardDir;

	if (keystate[SDLK_s])
		lightPos -= forwardDir;

	if (keystate[SDLK_d])
		lightPos += rightDir;

	if (keystate[SDLK_a])
		lightPos -= rightDir;


	if (keystate[SDLK_e]){
		// Rotate camera clockwise around y-axis
		yaw += dyaw;
		vec3 col1(cosf(yaw), 0, -sinf(yaw));
		vec3 col3(sinf(yaw), 0, cosf(yaw));
		R = mat3(col1, vec3(0, 1, 0), col3);
		rightDir = vec3(R[0][0], R[0][1], R[0][2]);
		downDir = vec3(R[1][0], R[1][1], R[1][2]);
		forwardDir = vec3(R[2][0], R[2][1], R[2][2]);
	}

	if (keystate[SDLK_q]){
		// Rotate camera anti-clockwise around y-axis
		yaw -= dyaw;
		//		cout << yaw;
		vec3 col1(cosf(yaw), 0, -sinf(yaw));
		vec3 col3(sinf(yaw), 0, cosf(yaw));
		R = mat3(col1, vec3(0, 1, 0), col3);
		rightDir = vec3(R[0][0], R[0][1], R[0][2]);
		downDir = vec3(R[1][0], R[1][1], R[1][2]);
		forwardDir = vec3(R[2][0], R[2][1], R[2][2]);
	}
		
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);
	
	//Reset depth buffer
	for (int y = 0; y < SCREEN_HEIGHT; y++){
		for (int x = 0; x < SCREEN_WIDTH; x++){
			depthBuffer[y][x] = 0.0f;
		}
	}

	for (int i = 0; i < triangles.size(); ++i)
	{
		currentColor = triangles[i].color;
		/*
		//Using vec3
		vector<vec3> vertices(3);
		//Extract the vertices of current triangle
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		*/
		

		//Using Vertex struct		
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		/*
		//Task 6.1
		//Set the vertex normals
		vertices[0].normal = triangles[i].normal;
		vertices[1].normal = triangles[i].normal;
		vertices[2].normal = triangles[i].normal;
		//Set the reflectance
		vertices[0].reflectance = currentColor;
		vertices[1].reflectance = currentColor;
		vertices[2].reflectance = currentColor;
		*/
		//Task 6.2
		currentNormal = triangles[i].normal;
		currentReflectance = currentColor;

		DrawPolygon(vertices);
		//DrawPolygonEdges(vertices);		//Task 3
		//		DrawLineSDL(screen, ivec2(0, 0), ivec2(500, 500), vec3(1, 1, 1));
		// Add drawing
		/*
		for (int v = 0; v<3; ++v)
		{
		ivec2 projPos;
		VertexShader(vertices[v], projPos);
		vec3 color(1, 1, 1);
		PutPixelSDL(screen, projPos.x, projPos.y, color);

		}
		*/
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}