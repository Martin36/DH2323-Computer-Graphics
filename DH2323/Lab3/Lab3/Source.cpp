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

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
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

// ----------------------------------------------------------------------------
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

void VertexShader(const vec3& v, ivec2& p){
	vec3 P = (v - cameraPos)*R;
	p.x = (focalLength * P.x / P.z) + SCREEN_WIDTH / 2;
	p.y = (focalLength * P.y / P.z) + SCREEN_HEIGHT / 2;
}

void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color){
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	vector<ivec2> line(pixels);
	Interpolate(a, b, line);
	for (int i = 0; i < pixels; i++){
		ivec2 pos = line[i];
		//		pos.x = pos.x + SCREEN_WIDTH / 2;
		//		pos.y = pos.y + SCREEN_HEIGHT / 2;
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
	vector<int> xValues;
	vector<int> yValues;
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
	for (int i = 0; i < nrOfRows; ++i)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
	}
	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.

}

int main(int argc, char* argv[])
{
	LoadTestModel(triangles);
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

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
		;

	if (keystate[SDLK_s])
		;

	if (keystate[SDLK_d]){
		// Rotate camera clockwise around y-axis
		yaw += dyaw;
		vec3 col1(cosf(yaw), 0, -sinf(yaw));
		vec3 col3(sinf(yaw), 0, cosf(yaw));
		R = mat3(col1, vec3(0, 1, 0), col3);
		rightDir = vec3(R[0][0], R[0][1], R[0][2]);
		downDir = vec3(R[1][0], R[1][1], R[1][2]);
		forwardDir = vec3(R[2][0], R[2][1], R[2][2]);
	}


	if (keystate[SDLK_a]){
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


	if (keystate[SDLK_e])
		;

	if (keystate[SDLK_q])
		;
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);
	for (int i = 0; i < triangles.size(); ++i)
	{
		vector<vec3> vertices(3);
		//Extract the vertices of current triangle
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		DrawPolygonEdges(vertices);
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