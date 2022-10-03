#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
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
    float numberOfSteps = max(abs(to.x - from.x), abs(to.y - from.y));
    float xStepSize = (to.x - from.x)/numberOfSteps;
    float yStepSize = (to.y - from.y)/numberOfSteps;
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

vector<CanvasPoint> interpolatePoints(CanvasPoint from, CanvasPoint to, float numberOfSteps){
    vector<CanvasPoint> linePoints;
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;
    for (float i=0.0; i<numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        linePoints.push_back(CanvasPoint(round(x), round(y)));
    }
    return linePoints;
}

vector<CanvasPoint> interpolateTexturePoints(TexturePoint from, TexturePoint to, float numberOfSteps){
    vector<CanvasPoint> linePoints;
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;
    for (float i=0.0; i<numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        linePoints.push_back(CanvasPoint(round(x), round(y)));
    }
    return linePoints;
}

void filledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    // Sort vertices into top, middle, bottom
    if(triangle.v0().y > triangle.v1().y) swap(triangle.v0(), triangle.v1());
    if(triangle.v0().y > triangle.v2().y) swap(triangle.v0(), triangle.v2());
    if(triangle.v1().y > triangle.v2().y) swap(triangle.v1(), triangle.v2());

    // Divide triangle in half horizontally - find "extra" point
    int extraY = triangle.v1().y;
    int a_x = triangle.v0().x;
    int a_y = triangle.v0().y;
    int c_x = triangle.v2().x;
    int c_y = triangle.v2().y;

    float ratio = (float) (extraY-a_y)/(float) (c_y-a_y);
    float extraX = (c_x - a_x) * ratio + a_x;

    // Top triangle
    float maxLeftLine = max(abs(triangle.v0().x - extraX), abs(triangle.v0().y - extraY));
    float maxRightLine = max(abs(triangle.v0().x - triangle.v1().x), abs(triangle.v0().y - triangle.v1().y));
    float numberOfSteps = max(maxLeftLine, maxRightLine);
    vector<CanvasPoint> left_line = interpolatePoints(triangle.v0(), CanvasPoint(extraX, extraY), numberOfSteps);
    vector<CanvasPoint> right_line = interpolatePoints(triangle.v0(), triangle.v1(), numberOfSteps);
    for(int i = 0; i < left_line.size(); i++){
        drawLine(window, left_line[i], right_line[i], colour);
    }

    // Bottom triangle
    float maxLeftLine2 = max(abs(extraX - triangle.v2().x), abs(extraY - triangle.v2().y));
    float maxRightLine2 = max(abs(triangle.v1().x - triangle.v2().x), abs(triangle.v1().y - triangle.v2().y));
    float numberOfSteps2 = max(maxLeftLine2, maxRightLine2);
    vector<CanvasPoint> left_line_2 = interpolatePoints(CanvasPoint(extraX, extraY), triangle.v2(), numberOfSteps2);
    vector<CanvasPoint> right_line_2 = interpolatePoints(triangle.v1(), triangle.v2(), numberOfSteps2);
    for(int i = 0; i < left_line_2.size(); i++){
        drawLine(window, left_line_2[i], right_line_2[i], colour);
    }

    // White border
    strokedTriangle(window, triangle, Colour(255, 255, 255));
}

void drawSubTextureTriangle(DrawingWindow &window, CanvasPoint a, CanvasPoint b, CanvasPoint c, TextureMap textureFile) {
    float maxLeftLine = max(abs(a.x - b.x), abs(a.y - b.y));
    float maxRightLine = max(abs(a.x - c.x), abs(a.y - c.y));
    float numberOfSteps = max(maxLeftLine, maxRightLine);
    vector<CanvasPoint> left_line = interpolatePoints(a, b, numberOfSteps);
    vector<CanvasPoint> right_line = interpolatePoints(a, c, numberOfSteps);

    float maxLeftLine_texture = max(abs(a.texturePoint.x - b.texturePoint.x), abs(a.texturePoint.y - b.texturePoint.y));
    float maxRightLine_texture = max(abs(a.texturePoint.x - c.texturePoint.x), abs(a.texturePoint.y - c.texturePoint.y));
    float numberOfSteps_texture = max(maxLeftLine_texture, maxRightLine_texture);
    vector<CanvasPoint> left_line_texture = interpolateTexturePoints(a.texturePoint, b.texturePoint, numberOfSteps_texture);
    vector<CanvasPoint> right_line_texture = interpolateTexturePoints(a.texturePoint, c.texturePoint, numberOfSteps_texture);

    for(int i = 0; i < left_line.size(); i++){
        float fraction_down_line = (float)i/left_line.size();
        int position_on_left_texture_line = round(fraction_down_line * left_line_texture.size());
        int position_on_right_texture_line = round(fraction_down_line * right_line_texture.size());

        // Draw line
        CanvasPoint to = left_line[i];
        CanvasPoint from = right_line[i];
        float numberOfSteps = max(abs(to.x - from.x), abs(to.y - from.y));
        float xStepSize = (to.x - from.x)/numberOfSteps;
        float yStepSize = (to.y - from.y)/numberOfSteps;
        CanvasPoint to_texture = left_line_texture[position_on_left_texture_line];
        CanvasPoint from_texture = right_line_texture[position_on_right_texture_line];

        float numberOfSteps_texture = max(abs(to_texture.x - from_texture.x), abs(to_texture.y - from_texture.y));
        float xStepSize_texture = (to_texture.x - from_texture.x)/numberOfSteps_texture;
        float yStepSize_texture = (to_texture.y - from_texture.y)/numberOfSteps_texture;
        for (float i=0.0; i<numberOfSteps; i++) {
            float x = from.x + (xStepSize * i);
            float y = from.y + (yStepSize * i);
            float x_texture = from_texture.x + (xStepSize_texture * i);
            float y_texture = from_texture.y + (yStepSize_texture * i);
            window.setPixelColour(round(x), round(y), textureFile.pixels[y_texture * textureFile.width + x_texture]);
        }
    }
}

void textureMapTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap textureFile) {
    // Sort vertices into top, middle, bottom
    if(triangle.v0().y > triangle.v1().y) swap(triangle.v0(), triangle.v1());
    if(triangle.v0().y > triangle.v2().y) swap(triangle.v0(), triangle.v2());
    if(triangle.v1().y > triangle.v2().y) swap(triangle.v1(), triangle.v2());

    // Find "extra" point
    float ratio = (float) (triangle.v1().y-triangle.v0().y)/(float) (triangle.v2().y-triangle.v0().y);
    float extraX = (triangle.v2().x - triangle.v0().x) * ratio + triangle.v0().x;
    CanvasPoint extra = CanvasPoint(extraX, triangle.v1().y);
    float extraX_texture = (triangle.v2().texturePoint.x - triangle.v0().texturePoint.x) * ratio + triangle.v0().texturePoint.x;
    extra.texturePoint = TexturePoint(extraX_texture, triangle.v1().texturePoint.y);

    // Top triangle
    drawSubTextureTriangle(window, triangle.v0(), extra, triangle.v1(), textureFile);

    // Bottom triangle
    drawSubTextureTriangle(window, triangle.v2(), extra, triangle.v1(), textureFile);

    // White border
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

        CanvasPoint cp1 = CanvasPoint(160, 10);
        cp1.texturePoint = TexturePoint(195, 5);
        CanvasPoint cp2 = CanvasPoint(300, 230);
        cp2.texturePoint = TexturePoint(395, 380);
        CanvasPoint cp3 = CanvasPoint(10, 150);
        cp3.texturePoint = TexturePoint(65, 330);
        CanvasTriangle triangle(cp1, cp2, cp3);
        TextureMap textureFile("texture.ppm");
        //CanvasTriangle textureTriangle(CanvasPoint(195, 5), CanvasPoint(395, 380), CanvasPoint(65, 330));

        //drawLine(window, from, to, colour);
        //drawColourGradient(window);
        //strokedTriangle(window, textureTriangle, colour);
        //filledTriangle(window, triangle1, colour);
        textureMapTriangle(window, triangle, textureFile);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}

}
