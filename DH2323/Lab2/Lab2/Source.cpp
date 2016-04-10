#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
float focalLength = SCREEN_WIDTH/2;
vec3 cameraPos(0, 0, -2);
mat3 R(vec3(1,0,0),vec3(0,1,0),vec3(0,0,1));
float yaw = 0;
float dyaw = M_PI / 10;
vec3 rightDir(R[0][0], R[0][1], R[0][2]);
vec3 downDir(R[1][0], R[1][1], R[1][2]);
vec3 forwardDir(R[2][0], R[2][1], R[2][2]);
vec3 lightPos(0, -0.5, -0.7);
vec3 lightColor = 14.f * vec3(1, 1, 1);
vec3 indirectLight = 0.5f*vec3(1, 1, 1);

struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
};
// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();

bool ClosestIntersection(vec3 start, vec3 dir,
	const vector<Triangle>& triangles,
	Intersection& closestIntersection){
	//Intersection intersection;
	vec3 s = start;
	vec3 d = dir;
	bool state = false;
	float m = std::numeric_limits<float>::max();
	
	for (int i = 0; i < triangles.size(); i++){
		Triangle triangle = triangles[i];
		vec3 v0 = triangle.v0;
		vec3 v1 = triangle.v1;
		vec3 v2 = triangle.v2;
		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = s - v0;
		mat3 A(-d, e1, e2);
		vec3 x = glm::inverse(A) * b;
		if (x.x >= 0 && x.y >= 0 && x.z >= 0 && x.y + x.z <= 1){	//Then we have an intersection
			vec3 r = s + x.x * d;
			float dist = glm::distance(s,r);
			if (dist < m){
				state = true;
				m = dist;
				closestIntersection.distance = dist;
				closestIntersection.position = r;
				closestIntersection.triangleIndex = i;
			}
		}

	}
	return state;
	

}
vec3 DirectLight(const Intersection& i){
	vec3 P = lightColor;
	vec3 rhat = glm::normalize(lightPos - i.position);
	vec3 nhat = triangles[i.triangleIndex].normal;
	float r = glm::distance(lightPos,i.position);
//	vec3 dir = (lightPos - i.position);
	//The shadow effect
	Intersection newI;
	if (ClosestIntersection(i.position + rhat*0.01f, rhat, triangles, newI)){
		if (newI.distance < r){
			vec3 D(0, 0, 0);
			return D;
		}
	}
	
	float y = glm::dot(rhat, nhat);
	float max = glm::max(y, 0.f);
	vec3 D = P * float(max / (4 * M_PI*r*r));
	return D;
	
}

int main(int argc, char* argv[])
{
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.
	LoadTestModel(triangles);

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
	if (keystate[SDLK_w]){
		lightPos += forwardDir;
	}
	if (keystate[SDLK_s]){
		lightPos -= forwardDir;
	}
	if (keystate[SDLK_a]){
		lightPos -= rightDir;
	}
	if (keystate[SDLK_d]){
		lightPos += rightDir;
	}
	if (keystate[SDLK_q]){
		lightPos -= downDir;
	}
	if (keystate[SDLK_e]){
		lightPos += downDir;
	}

	if (keystate[SDLK_UP])
	{
		// Move camera forward
		cameraPos += forwardDir;
		
	}
	if (keystate[SDLK_DOWN])
	{
		// Move camera backward
		cameraPos -= forwardDir;
	}
	if (keystate[SDLK_LEFT])
	{
/*		// Move camera to the left
		cameraPos.x--;
*/
		// Rotate camera anti-clockwise around y-axis
		yaw += dyaw;
		cout << yaw;
		vec3 col1(cosf(yaw), 0, -sinf(yaw));
		vec3 col3(sinf(yaw), 0, cosf(yaw));
		R = mat3(col1, vec3(0, 1, 0), col3);
		vec3 rightDir(R[0][0], R[0][1], R[0][2]);
		vec3 downDir(R[1][0], R[1][1], R[1][2]);
		vec3 forwardDir(R[2][0], R[2][1], R[2][2]);
	}
	if (keystate[SDLK_RIGHT])
	{
/*
		// Move camera to the right
		cameraPos.x++;
*/
		// Rotate camera clockwise around y-axis
		yaw -= dyaw;
		vec3 col1(cosf(yaw), 0, -sinf(yaw));
		vec3 col3(sinf(yaw), 0, cosf(yaw));
		R = mat3(col1, vec3(0, 1, 0), col3);
		vec3 rightDir(R[0][0], R[0][1], R[0][2]);
		vec3 downDir(R[1][0], R[1][1], R[1][2]);
		vec3 forwardDir(R[2][0], R[2][1], R[2][2]);
	}

}

void Draw()
{
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (int y = 0; y<SCREEN_HEIGHT; ++y)
	{
		for (int x = 0; x<SCREEN_WIDTH; ++x)
		{
			vec3 d = vec3(x - SCREEN_WIDTH / 2, y - SCREEN_HEIGHT / 2, focalLength) * R;
//		vec3 d(x - SCREEN_WIDTH / 2, y - SCREEN_HEIGHT / 2, focalLength);
			vec3 s = cameraPos;
			Intersection closestIntersection;
			if (ClosestIntersection(s, d, triangles, closestIntersection)){
				vec3 start = closestIntersection.position;
				vec3 dirToLight = lightPos - start;
				vec3 D = DirectLight(closestIntersection);
				vec3 rho = triangles[closestIntersection.triangleIndex].color;
//			vec3 color = rho;		//Plain color
//			vec3 color = D;			//Illumination without color
//			vec3 color = rho * D;		//Illumination from direct light with color
			vec3 color = rho * (D + indirectLight);		//Illumination from direct and indirect light with color


				PutPixelSDL(screen, x, y, color);

//				vec3 color = triangles[closestIntersection.triangleIndex].color;
			}
			else{
				vec3 color(0, 0, 0);
				PutPixelSDL(screen, x, y, color);

			}
		}
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}