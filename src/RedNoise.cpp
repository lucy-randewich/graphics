#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <glm/glm.hpp>
#include "glm/ext.hpp"
#include <Utils.h>
#include <fstream>
#include <vector>

#define WIDTH 320
#define HEIGHT 240

using namespace std;


vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    vector<float> vec;
    float increment = (to-from)/(numberOfValues-1);
    for(size_t i = 0; i < numberOfValues; i++){
        vec.push_back(from);
        from += increment;
    }
    return vec;
}

vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
    vector<glm::vec3> vec;

    vector<float> val1 = interpolateSingleFloats(from[0], to[0], numberOfValues);
    vector<float> val2 = interpolateSingleFloats(from[1], to[1], numberOfValues);
    vector<float> val3 = interpolateSingleFloats(from[2], to[2], numberOfValues);

    for(size_t i=0; i<numberOfValues; i++) {
        glm::vec3 tmp(val1[i], val2[i], val3[i]);
        vec.push_back(tmp);
    }
    return vec;
}

void drawGreyGradient(DrawingWindow &window) {
	window.clearPixels();
    vector<float> widthPixelValues = interpolateSingleFloats(255, 0, window.width);
    for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = widthPixelValues[x];
			float green = widthPixelValues[x];
			float blue = widthPixelValues[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void drawColourGradient(DrawingWindow &window){
    window.clearPixels();

    glm::vec3 topLeft(255, 0, 0);        // red
    glm::vec3 topRight(0, 0, 255);       // blue
    glm::vec3 bottomRight(0, 255, 0);    // green
    glm::vec3 bottomLeft(255, 255, 0);   // yellow

    vector<glm::vec3> redToYellow = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
    vector<glm::vec3> blueToGreen = interpolateThreeElementValues(topRight, bottomRight, window.height);

    for (size_t y = 0; y < window.height; y++) {
        vector<glm::vec3> rowVals = interpolateThreeElementValues(redToYellow[y], blueToGreen[y], window.width);
        for (size_t x = 0; x < window.width; x++) {
            glm::vec3 thisRowVals = rowVals[x];

            float red = thisRowVals[0];
            float green = thisRowVals[1];
            float blue = thisRowVals[2];
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);

        }
    }
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;
    for (float i=0.0; i<numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        window.setPixelColour(round(x), round(y), (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue));
    }
}

void strokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    drawLine(window, triangle.v0(), triangle.v1(), colour);
    drawLine(window, triangle.v1(), triangle.v2(), colour);
    drawLine(window, triangle.v0(), triangle.v2(), colour);
}


void filledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    // Sort vertices into top, middle, bottom
    if(triangle.v0().y > triangle.v1().y) swap(triangle.v0(), triangle.v1());
    if(triangle.v0().y > triangle.v2().y) swap(triangle.v0(), triangle.v2());
    if(triangle.v1().y > triangle.v2().y) swap(triangle.v1(), triangle.v2());
    cout << triangle << endl;

    // Divide triangle in half horizontally - find "extra" point
    int extraY = triangle.v1().y;
    int a_x = triangle.v0().x;
    int a_y = triangle.v0().y;
    int c_x = triangle.v2().x;
    int c_y = triangle.v2().y;
    int b_x = triangle.v1().x;
    int b_y = triangle.v1().y;

    float ratio = (float) (extraY-a_y)/(float) (c_y-a_y);
    float extraX = (c_x - a_x) * ratio + a_x;

    // Top triangle
    for(int y = a_y; y < extraY; y++){
        float ratioi = (float) (y-a_y)/(float) (c_y-a_y);
        float i = (c_x - a_x) * ratioi + a_x;

        float ratioj = (float) (y-a_y)/(float) (b_y-a_y);
        float j = (b_x - a_x) * ratioj + a_x;

        drawLine(window, CanvasPoint(i, y), CanvasPoint(j, y), colour);
    }

    // Bottom triangle
    for(int y = extraY; y < c_y; y++){
        float ratioi = (float) (y-a_y)/(float) (c_y-a_y);
        float i = (c_x - a_x) * ratioi + a_x;

        float ratioj = (float) (y-b_y)/(float) (c_y-b_y);
        float j = (c_x - b_x) * ratioj + b_x;

        drawLine(window, CanvasPoint(i, y), CanvasPoint(j, y), colour);
    }



    // Draw white border
    strokedTriangle(window, triangle, Colour(255, 255, 255));
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u) {
            CanvasTriangle triangle(CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height));
            Colour colour(rand()%256, rand()%256, rand()%256);
            strokedTriangle(window, triangle, colour);
        }else if (event.key.keysym.sym == SDLK_f) {
            CanvasTriangle triangle(CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height));
            Colour colour(rand()%256, rand()%256, rand()%256);
            filledTriangle(window, triangle, colour);
        }
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}


int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);

        CanvasPoint from(window.width/2, 0);
        CanvasPoint to(window.width/2, window.height);
        Colour colour(0, 255, 0);
        //drawLine(window, from, to, colour);
        //drawColourGradient(window);
        CanvasTriangle triangle(CanvasPoint(10, 10), CanvasPoint(10, 200), CanvasPoint(30, 100));
        //strokedTriangle(window, triangle, colour);
        filledTriangle(window, triangle, colour);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}

}
