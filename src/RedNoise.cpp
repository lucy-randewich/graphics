#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <glm/glm.hpp>
#include "glm/ext.hpp"
#include <Utils.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define WIDTH 320.0f*3
#define HEIGHT 240.0f*3
#define FOCAL_LENGTH 2.0f
#define SCALER WIDTH

using namespace std;

string renderer = "ray_traced";

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
    float numberOfSteps = max(abs(to.x - from.x), abs(to.y - from.y));
    float xStepSize = (to.x - from.x)/numberOfSteps;
    float yStepSize = (to.y - from.y)/numberOfSteps;
    for (float i=0.0; i<=numberOfSteps; i++) {
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

vector<CanvasPoint> interpolatePoints(CanvasPoint from, CanvasPoint to, float numberOfSteps) {
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
    //strokedTriangle(window, triangle, Colour(255, 255, 255));
}

vector<ModelTriangle> readOBJFile(string objfile, float scale_factor, vector<Colour> colour_library) {
    ifstream file(objfile);
    char character;
    float x,y,z;
    string tmpv1, tmpv2, tmpv3;
    int v1, v2, v3;
    string line;
    string colour_name;
    Colour colour;
    glm::vec3 p1, p2, p3;

    vector<glm::vec3> vertices;
    vector<glm::vec3> faces;
    vector<ModelTriangle> triangles;

    // Add all vertices to vector and triangles to vector
    while(getline(file, line)) {
        istringstream stream(line);
        stream >> character;
        switch (character) {
            case 'v':
                stream >> x >> y >> z;
                vertices.push_back({x * scale_factor, y * scale_factor, z * scale_factor});
                break;
            case 'f':
                stream >> tmpv1 >> tmpv2 >> tmpv3;
                v1 = stoi(tmpv1);
                v2 = stoi(tmpv2);
                v3 = stoi(tmpv3);

                p1 = vertices[v1-1];
                p2 = vertices[v2-1];
                p3 = vertices[v3-1];
                triangles.push_back(ModelTriangle(p1, p2, p3, colour));
                break;
            case 'u':
                stream >> tmpv1 >> colour_name;
                for (Colour curr_colour: colour_library){
                    if (curr_colour.name == colour_name) {
                        colour = curr_colour;
                    }
                }
        }
    }
    return triangles;
}

vector<Colour> readMTLFile(string mtlfile) {
    ifstream file(mtlfile);
    vector<Colour> colours;
    string keyword, name, line, line2;
    float r, g, b;

    while(getline(file, line)) {
        // Get colour name from first line
        istringstream stream(line);
        stream >> keyword;
        stream >> name;

        // Get colour values from second line
        getline(file, line2);
        istringstream stream2(line2);
        stream2 >> keyword;
        stream2 >> r >> g >> b;
        colours.push_back(Colour(name, r*255, g*255, b*255));

        // Iterate over blank line and ignore
        getline(file, line2);
    }
    return colours;
}

CanvasPoint getCanvasIntersectionPointWithOrientation(DrawingWindow &window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, glm::vec3 vertexPosition, float focalLength) {
    glm::vec3 projectedVertex = (vertexPosition - cameraPosition) * cameraOrientation;
    //glm::vec3 projectedVertex = cameraOrientation * (cameraPosition - vertexPosition);

    float u = focalLength * -projectedVertex.x/projectedVertex.z;
    float v = focalLength * projectedVertex.y/projectedVertex.z;
    float z = focalLength * projectedVertex.z;

    u *= window.width;
    v *= window.width;
    //z *= window.width;

    u += window.width/2.0f;
    v += window.height/2.0f;

    return CanvasPoint(round(u), round(v), z);
}

vector<CanvasPoint> interpolatePointsWithDepth(CanvasPoint from, CanvasPoint to, float numberOfSteps) {
    vector<CanvasPoint> linePoints;
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float zDiff = to.depth - from.depth;
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;
    float zStepSize = zDiff/numberOfSteps;
    for (float i=0.0; i<numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        float z = from.depth + (zStepSize * i);
        linePoints.push_back(CanvasPoint(round(x), round(y), z));
    }
    return linePoints;
}

void drawLineWithDepth(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour, float **depth_buffer, glm::vec3 cameraPosition) {
    float numberOfSteps = max(abs(to.x - from.x), abs(to.y - from.y));
    float xStepSize = (to.x - from.x)/numberOfSteps;
    float yStepSize = (to.y - from.y)/numberOfSteps;
    float zStepSize = (to.depth - from.depth)/numberOfSteps;
    for (float i=0.0; i<=numberOfSteps; i++) {
        int x = int(round(from.x + (xStepSize * i)));
        int y = int(round(from.y + (yStepSize * i)));
        float point_depth = 1/((1+exp(-(from.depth + (zStepSize * i)))));
        //float point_depth = 1/((from.depth + (zStepSize * i)));

        if(x < window.width && x > 0 && y < window.height && y > 0){
                if (depth_buffer[x][y] <= point_depth) {
                    window.setPixelColour(x, y, (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue));
                    depth_buffer[x][y] = point_depth;
                }
        }
    }
}

void filledTriangleWithDepth(DrawingWindow &window, CanvasTriangle triangle, Colour colour, float **depth_buffer, glm::vec3 cameraPosition) {
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
    float c_z = triangle.v2().depth;
    float a_z = triangle.v0().depth;

    float ratio = (float) (extraY-a_y)/(float) (c_y-a_y);
    float extraX = (c_x - a_x) * ratio + a_x;
    float extraZ = (c_z - a_z) * ratio + a_z;

    // Top triangle
    float maxLeftLine = max(abs(triangle.v0().x - extraX), abs(triangle.v0().y - extraY));
    float maxRightLine = max(abs(triangle.v0().x - triangle.v1().x), abs(triangle.v0().y - triangle.v1().y));
    float numberOfSteps = max(maxLeftLine, maxRightLine);
    vector<CanvasPoint> left_line = interpolatePointsWithDepth(triangle.v0(), CanvasPoint(extraX, extraY, extraZ), numberOfSteps);
    vector<CanvasPoint> right_line = interpolatePointsWithDepth(triangle.v0(), triangle.v1(), numberOfSteps);
    for(int i = 0; i < left_line.size(); i++){
        drawLineWithDepth(window, left_line[i], right_line[i], colour, depth_buffer, cameraPosition);
    }

    // Bottom triangle
    float maxLeftLine2 = max(abs(extraX - triangle.v2().x), abs(extraY - triangle.v2().y));
    float maxRightLine2 = max(abs(triangle.v1().x - triangle.v2().x), abs(triangle.v1().y - triangle.v2().y));
    float numberOfSteps2 = max(maxLeftLine2, maxRightLine2);
    vector<CanvasPoint> left_line_2 = interpolatePointsWithDepth(CanvasPoint(extraX, extraY, extraZ), triangle.v2(), numberOfSteps2);
    vector<CanvasPoint> right_line_2 = interpolatePointsWithDepth(triangle.v1(), triangle.v2(), numberOfSteps2);
    for(int i = 0; i < left_line_2.size(); i++){
        drawLineWithDepth(window, left_line_2[i], right_line_2[i], colour, depth_buffer, cameraPosition);
    }

}

void rasteriseObj(DrawingWindow &window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, vector<Colour> colour_library, vector<ModelTriangle> triangles, bool wireframe){
    // Set up matrix of floats for 1/z depth buffer
    float **depth_buffer;
    depth_buffer = new float *[window.width];
    for(int i = 0; i <window.width; i++)
        depth_buffer[i] = new float[window.height];

    for(int i = 0; i < window.width; i++){
        for(int j = 0; j < window.height; j++){
            depth_buffer[i][j] = 0;
        }
    }

    for (ModelTriangle triangle : triangles) {
        Colour pixel_colour = triangle.colour;

        // Get canvas intersection points of triangle vertices
        CanvasPoint p1 = getCanvasIntersectionPointWithOrientation(window, cameraPosition, cameraOrientation, triangle.vertices[0], 2);
        CanvasPoint p2 = getCanvasIntersectionPointWithOrientation(window, cameraPosition, cameraOrientation, triangle.vertices[1], 2);
        CanvasPoint p3 = getCanvasIntersectionPointWithOrientation(window, cameraPosition, cameraOrientation, triangle.vertices[2], 2);
        CanvasTriangle ctriangle(p1, p2, p3);

        if(wireframe) strokedTriangle(window, ctriangle, pixel_colour);
        else filledTriangleWithDepth(window, ctriangle, pixel_colour, depth_buffer, cameraPosition);
    }
    //drawDepth(window, depth_buffer);
}

void lookAt(glm::vec3 pointToLookAt, glm::mat3 &cameraOrientation, glm::vec3 &cameraPosition) {
    glm::vec3 forward = glm::normalize(cameraPosition - pointToLookAt);
    glm::vec3 vertical = glm::vec3(0, 1, 0);
    glm::vec3 right = glm::normalize(glm::cross(vertical, forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
    cameraOrientation = glm::mat3(right, up, forward);
}

RayTriangleIntersection getClosestIntersection(glm::vec3 &cameraPosition, glm::vec3 rayDirection, vector<ModelTriangle> &triangles) {
    RayTriangleIntersection closestIntersection = RayTriangleIntersection();
    closestIntersection.intersectionFound = false;
    float smallest_t = 999999;
    int index = 0;
    for (ModelTriangle triangle : triangles) {
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;     // In form (t,u,v)
        float t = possibleSolution[0];
        float u = possibleSolution[1];
        float v = possibleSolution[2];

        if((u>=0.0) && (u<=1.0) && (v>=0.0) && (v<=1.0) && ((u+v)<=1.0)){
            if((t < smallest_t) && t>0.0005){
                glm::vec3 intersection = (triangle.vertices[0] + (u * e0) + (v * e1));
                //glm::vec3 intersection = (t * rayDirection) + cameraPosition;
                smallest_t = t;
                closestIntersection = RayTriangleIntersection(intersection, t, triangle, index, true);
            }
        }
        index++;
    }
    return closestIntersection;
}

void rayTraceObj(DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, vector<Colour> &colour_library, vector<ModelTriangle> &triangles, glm::vec3 &lightsource) {
    for (float x = 0; x < window.width; x++){
        for (float y = 0; y < window.height; y++){
            glm::vec3 rayDirection = glm::vec3(x, y, -1.0f);
            rayDirection[0] = ((rayDirection[0] - WIDTH/2.0f)/(SCALER*FOCAL_LENGTH));
            rayDirection[1] = -((rayDirection[1] - HEIGHT/2.0f)/(SCALER*FOCAL_LENGTH));

            rayDirection = rayDirection * glm::inverse(cameraOrientation);

            RayTriangleIntersection intersectionTriangle = getClosestIntersection(cameraPosition, glm::normalize(rayDirection), triangles);

            if(intersectionTriangle.intersectionFound) {
                glm::vec3 shadowRay = lightsource - intersectionTriangle.intersectionPoint;
                float distance = glm::distance(lightsource,intersectionTriangle.intersectionPoint);

                RayTriangleIntersection shadow_intersection = getClosestIntersection(intersectionTriangle.intersectionPoint, glm::normalize(shadowRay), triangles);

                if ((shadow_intersection.intersectionFound) && (shadow_intersection.distanceFromCamera < distance) && (intersectionTriangle.triangleIndex != shadow_intersection.triangleIndex)){
                    window.setPixelColour(x, y, 0);
                }else {     // No shadow
                    Colour colour = intersectionTriangle.intersectedTriangle.colour;
                    window.setPixelColour(x, y, (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue));
                }

            }
        }
    }
}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, glm::vec3 &lightsource) {
	if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_1) renderer = "wireframe";
        else if (event.key.keysym.sym == SDLK_2) renderer = "rasterised";
        else if (event.key.keysym.sym == SDLK_3) renderer = "ray_traced";
        else if (event.key.keysym.sym == SDLK_q) {
            lightsource[1] = lightsource[1] + 0.1;
            cout << "lightsource is " << lightsource[0] << " " << lightsource[1] << " " << lightsource[2] << endl;
        }
        else if (event.key.keysym.sym == SDLK_a) {
            lightsource[1] = lightsource[1] - 0.1;
            cout << "lightsource is " << lightsource[0] << " " << lightsource[1] << " " << lightsource[2] << endl;
        }
		else if (event.key.keysym.sym == SDLK_LEFT) cameraPosition[0] = cameraPosition[0] - 0.2;
		else if (event.key.keysym.sym == SDLK_RIGHT) cameraPosition[0] = cameraPosition[0] + 0.2;
		else if (event.key.keysym.sym == SDLK_UP) cameraPosition[1] = cameraPosition[1] + 0.2;
		else if (event.key.keysym.sym == SDLK_DOWN) cameraPosition[1] = cameraPosition[1] - 0.2;
        else if (event.key.keysym.sym == SDLK_w) cameraPosition[2] = cameraPosition[2] + 0.2;
        else if (event.key.keysym.sym == SDLK_s) cameraPosition[2] = cameraPosition[2] - 0.2;
        else if (event.key.keysym.sym == SDLK_j) {             // Pan camera (y axis)
            float theta = glm::radians(2.5);
            glm::mat3 rotate_matrix = glm::mat3(cos(theta), 0.0, -sin(theta),
                                                0.0, 1.0, 0.0,
                                                sin(theta), 0.0, cos(theta));
            cameraPosition = cameraPosition * rotate_matrix;
            lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
        }else if (event.key.keysym.sym == SDLK_l) {            // Tilt camera (x axis)
            float theta = glm::radians(2.5);
            glm::mat3 rotate_matrix = glm::mat3(1, 0, 0,
                                                0, cos(theta), sin(theta),
                                                0, -sin(theta), cos(theta));
            cameraPosition = cameraPosition * rotate_matrix;
            lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
        }else if (event.key.keysym.sym == SDLK_u) {
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

    glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 3.5);
    glm::vec3 lightsource = glm::vec3(0, 0.4, 0);
    glm::mat3 cameraOrientation = glm::mat3(1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1);

    // Read obj data from files
    vector<Colour> colour_library = readMTLFile("cornell-box.mtl");
    vector<ModelTriangle> triangles = readOBJFile("cornell-box.obj", 0.15, colour_library);

    while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) {
            handleEvent(event, window, cameraPosition, cameraOrientation, lightsource);
        }

        window.clearPixels();
        if(renderer == "wireframe") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, true);
        else if (renderer == "rasterised") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, false);
        else if (renderer == "ray_traced") rayTraceObj(window, cameraPosition, cameraOrientation, colour_library, triangles, lightsource);

        // Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
