#include <CanvasTriangle.h>
#include <DrawingWindow.h>
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

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}


int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

    glm::vec3 from(1, 4, 9.2);
    glm::vec3 to(4, 1, 9.8);
    vector<glm::vec3> result = interpolateThreeElementValues(from, to, 4);
    for(size_t i=0; i<result.size(); i++) {
        std::cout << glm::to_string(result[i]);
    }

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		drawColourGradient(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}

}
