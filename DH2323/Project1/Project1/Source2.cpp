// Introduction lab that covers:
// * C++
// * SDL
// * 2D graphics
// * Plotting pixels
// * Video memory
// * Color representation
// * Linear interpolation
// * glm::vec3 and std::vector

#include "SDL.h"
#include <iostream>
#include <C:\Users\martin\Documents\CgLab1\glm\glm\glm.hpp>
#include <vector>
#include "C:\Users\martin\Documents\CgLab1\SDLauxiliary.h"

using namespace std;
using glm::vec3;

// --------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
SDL_Surface* screen;
const vector<vec3> stars(1000);


// --------------------------------------------------------
// FUNCTION DECLARATIONS

void Draw();

void Interpolate(float a, float b, vector<float>& result){
	int N = result.capacity();
	if (N == 1){
		result[0] = a;
	}
	else{
		float d = (b - a) / (N - 1);
		for (int i = 0; i < N; i++){
			result[i] = a + i*d;
		}
	}
}

void Interpolate(vec3 a, vec3 b, vector<vec3>& result){
	int N = result.capacity();
	if (N == 1){
		result[0] = a;
	}
	else{
		float dx = (b.x - a.x) / (N - 1);
		float dy = (b.y - a.y) / (N - 1);
		float dz = (b.z - a.z) / (N - 1);
		for (int i = 0; i < N; i++){
			vec3 temp(a.x + i*dx, a.y + i*dy, a.z + i*dz);
			result[i] = temp;
		}
	}
}

// --------------------------------------------------------
// FUNCTION DEFINITIONS

int main(int argc, char* argv[])
{
	/*
	vector<float> result(1); // Create a vector width 10 floats
	Interpolate(5, 4, result); // Fill it with interpolated values
	for (int i = 0; i<result.size(); ++i)
		cout << result[i] << " "; // Print the result to the terminal
		vector<vec3> result(4);
		vec3 a(1, 4, 9.2);
		vec3 b(4, 1, 9.8);
		Interpolate(a, b, result);
		for (int i = 0; i<result.size(); ++i)
		{
		cout << "( "
		<< result[i].x << ", "
		<< result[i].y << ", "
		<< result[i].z << " ) ";
		}

	*/
	for (i = 0; i < stars.capacity(); i++){
		float x = -1 + float(rand()) / float(2);
		float y = -1 + float(rand()) / float(2);
		float z = float(rand()) / float(RAND_MAX);
		vec3 temp = (x, y, z);
		stars[i] = temp;
	}
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	while (NoQuitMessageSDL())
	{
		Draw();
	}
	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}
/*void Draw()
{
	vec3 topLeft(1, 0, 0); // red
	vec3 topRight(0, 0, 1); // blue
	vec3 bottomLeft(0, 1, 0); // green
	vec3 bottomRight(1, 1, 0); // yellow
	vector<vec3> leftSide(SCREEN_HEIGHT);
	vector<vec3> rightSide(SCREEN_HEIGHT);
	Interpolate(topLeft, bottomLeft, leftSide);
	Interpolate(topRight, bottomRight, rightSide);
	for (int y = 0; y<SCREEN_HEIGHT; ++y)
	{
		vector<vec3> temp(SCREEN_WIDTH);
		Interpolate(leftSide[y], rightSide[y], temp);

		for (int x = 0; x<SCREEN_WIDTH; ++x)
		{
			vec3 color = temp[x];
			PutPixelSDL(screen, x, y, color);
		}
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}*/
void Draw()
{
	SDL_FillRect(screen, 0, 0);
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);
	for (size_t s = 0; s<stars.size(); ++s)
	{
		vec3 coord = stars[s];
		vec3 color = (1, 1, 1);
		PutPixelSDL(screen, coord.x, coord.y, color);

		// Add code for projecting and drawing each star
	}
	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);
	SDL_UpdateRect(screen, 0, 0, 0, 0);
}
