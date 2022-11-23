#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <glm/glm.hpp>
#include "glm/ext.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <thread>

#define WIDTH (320.0f)*3.0f
#define HEIGHT (240.0f)*3.0f
#define FOCAL_LENGTH 2.0f
#define SCALER WIDTH

using namespace std;

string renderer = "wireframe";
bool soft_shadows = false;
bool sphere = false;
bool mirror = false;
bool phong = true;
bool glass = true;
vector<glm::vec3> lightsources;

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

vector<ModelTriangle> readOBJFile(string objfile, float scale_factor, vector<Colour> colour_library, vector<glm::vec3> &vertices, float shift) {
    ifstream file(objfile);
    char character;
    float x,y,z;
    string tmpv1, tmpv2, tmpv3;
    int v1, v2, v3;
    string line;
    string colour_name;
    Colour colour;
    colour.name = " ";
    glm::vec3 p0, p1, p2;
    glm::vec3 face_normal;

    vector<ModelTriangle> triangles;

    // Add all vertices to vector and triangles to vector
    while(getline(file, line)) {
        istringstream stream(line);
        stream >> character;
        switch (character) {
            case 'v':
                stream >> x >> y >> z;
                vertices.push_back({(x * scale_factor)+shift, (y * scale_factor)-shift, z * scale_factor});
                break;
            case 'f':
                stream >> tmpv1 >> tmpv2 >> tmpv3;
                v1 = stoi(tmpv1);
                v2 = stoi(tmpv2);
                v3 = stoi(tmpv3);

                p0 = vertices[v1-1];
                p1 = vertices[v2-1];
                p2 = vertices[v3-1];

                if (colour.name == " "){
                    colour = Colour("Red", 0xFF, 0x00, 0x00);
                }
                face_normal = glm::cross(p0-p2, p1-p2);
                triangles.push_back(ModelTriangle(p0, p1, p2, colour, face_normal));
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

vector<glm::vec3> getVertexNormals(vector<glm::vec3> vertices, vector<ModelTriangle> triangles){
    vector<glm::vec3> vertex_normals;

    for (glm::vec3 vertex : vertices){
        // Find normals of all triangles which use this vertex
        vector<glm::vec3> normalsToAverage;
        for (ModelTriangle triangle : triangles){
            for(glm::vec3 triangleVertex : triangle.vertices){
                if (vertex == triangleVertex){
                    normalsToAverage.push_back(triangle.normal);
                }
            }
        }
        // Average normals of these triangles
        float count = 0.0f;
        glm::vec3 sum = {0, 0, 0};
        for (glm::vec3 normal : normalsToAverage){
            sum += normal;
            count ++;
        }
        vertex_normals.push_back(glm::normalize(sum/count));
    }
    return vertex_normals;
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
            if((t < smallest_t) && t>0.000001){
                glm::vec3 intersection = (triangle.vertices[0] + (u * e0) + (v * e1));
                smallest_t = t;
                closestIntersection = RayTriangleIntersection(intersection, t, u, v, triangle, index, true, triangle.material);
            }
        }
        index++;
    }
    return closestIntersection;
}

glm::vec3 getNormal(RayTriangleIntersection intersectionTriangle, vector<glm::vec3> vertices, vector<glm::vec3> vertex_normals){
    // Get three normals of vertices of intersected triangle
    glm::vec3 v0 = intersectionTriangle.intersectedTriangle.vertices[0];
    glm::vec3 v1 = intersectionTriangle.intersectedTriangle.vertices[1];
    glm::vec3 v2 = intersectionTriangle.intersectedTriangle.vertices[2];
    glm::vec3 n0, n1, n2;
    for(size_t i = 0; i < vertices.size(); i++){
        if (vertices[i].x == v0.x && vertices[i].y == v0.y && vertices[i].z == v0.z){
            n0 = vertex_normals[i];
        }else if(vertices[i].x == v1.x && vertices[i].y == v1.y && vertices[i].z == v1.z){
            n1 = vertex_normals[i];
        }else if(vertices[i].x == v2.x && vertices[i].y == v2.y && vertices[i].z == v2.z){
            n2 = vertex_normals[i];
        }
    }

    // Interpolate normal from vertex normals
    float proportion_n0 = 1 - (intersectionTriangle.u + intersectionTriangle.v);
    float proportion_n1 = intersectionTriangle.u;
    float proportion_n2 = intersectionTriangle.v;
    glm::vec3 normal = glm::normalize((proportion_n0 * n0) + (proportion_n1 * n1) + (proportion_n2 * n2));
    return normal;
}

float getBrightness(float &distance, RayTriangleIntersection &intersectionTriangle, RayTriangleIntersection &shadow_intersection, glm::vec3 &lightsource, glm::vec3 &lightRay, glm::vec3 &rayDirection, glm::vec3 &normal){
    // PROXIMITY LIGHTING
    float proximity_brightness = glm::min(float((5.0f/(7.0f * pow(distance, 2.0f) * M_PI))), 1.0f);

    // INCIDENCE LIGHTING
    float incidence_angle = glm::dot(glm::normalize(normal), glm::normalize(lightsource - intersectionTriangle.intersectionPoint));
    float incident_brightness = glm::max(incidence_angle, 0.0f);

    // REFLECTION
    glm::vec3 reflectedRay = lightRay - (2 * normal) * (glm::dot(glm::normalize(lightRay), glm::normalize(normal)));
    float reflection_angle = glm::dot(glm::normalize(rayDirection), glm::normalize(reflectedRay));
    float specular_spread = pow(reflection_angle, 256);

    // Combine lighting types
    float brightness = glm::clamp((proximity_brightness * incident_brightness), 0.0f, 1.0f);

    if (specular_spread > brightness){
        brightness = specular_spread;
    }

    // SHADOW
    if ((shadow_intersection.intersectionFound) && (shadow_intersection.distanceFromCamera < distance) && (intersectionTriangle.triangleIndex != shadow_intersection.triangleIndex) && shadow_intersection.intersectedTriangle.material != "glass"){
        brightness = 0;
    }

    // AMBIENT LIGHTING
    brightness = glm::clamp(brightness, 0.15f, 1.0f);

    return brightness;
}

glm::vec3 refractt(const glm::vec3 &I, const glm::vec3 &N, const float &ior){
    float cosi = glm::clamp(glm::dot(I, N), -1.0f, 1.0f);
    float etai = 1, etat = ior;
    glm::vec3 n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    glm::vec3 direction =  k < 0 ? glm::vec3(0,0,0) : eta * I + (eta * cosi - sqrtf(k)) * n;
    return direction;
}

void fresnel(const glm::vec3 &I, const glm::vec3 &N, const float &ior, float &kr){
    float cosi = glm::clamp(glm::dot(I, N), -1.0f, 1.0f);
    float etai = 1, etat = ior;
    if (cosi > 0) { std::swap(etai, etat); }
    float sint = etai / etat * sqrtf(max(0.f, 1 - cosi * cosi));
    if (sint >= 1) {
        kr = 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
}

glm::vec3 reflectt(const glm::vec3 &I, const glm::vec3 &N)
{
    return I - 2 * glm::dot(I, N)* N;
}

Colour castRay(glm::vec3 startPoint, vector<glm::vec3> &vertices, vector<glm::vec3> &vertex_normals, vector<ModelTriangle> &triangles, glm::vec3 rayDirection, float &brightness){
    RayTriangleIntersection intersectionTriangle = getClosestIntersection(startPoint, rayDirection, triangles);
    Colour colour = intersectionTriangle.intersectedTriangle.colour;
    glm::vec3 normal = getNormal(intersectionTriangle, vertices, vertex_normals);

    // Find shadow intersection
    brightness = 0;
    for (glm::vec3 light : lightsources){
        glm::vec3 lightRay = light - intersectionTriangle.intersectionPoint;
        float distance = glm::distance(light, intersectionTriangle.intersectionPoint);
        RayTriangleIntersection shadow_intersection = getClosestIntersection(intersectionTriangle.intersectionPoint, glm::normalize(lightRay), triangles);
        if(phong){
            brightness += getBrightness(distance, intersectionTriangle, shadow_intersection, light, lightRay, rayDirection, normal);
        }else{
            brightness += getBrightness(distance, intersectionTriangle, shadow_intersection, light, lightRay, rayDirection, intersectionTriangle.intersectedTriangle.normal);
        }
    }
    brightness = glm::clamp(brightness/lightsources.size(), 0.2f, 1.0f);

    if(intersectionTriangle.intersectedMaterial == "mirror"){
        glm::vec3 reflected_ray = glm::normalize(rayDirection) - (2 * glm::normalize(intersectionTriangle.intersectedTriangle.normal)) * (glm::dot(glm::normalize(rayDirection), glm::normalize(intersectionTriangle.intersectedTriangle.normal)));
        colour = castRay(intersectionTriangle.intersectionPoint, vertices, vertex_normals, triangles, reflected_ray, brightness);
    }

    if(intersectionTriangle.intersectedMaterial == "glass"){

        float refractive_index = 2.0f;

        bool outside = glm::dot(intersectionTriangle.intersectedTriangle.normal, glm::normalize(rayDirection)) < 0.0f;
        float kr;
        fresnel(glm::normalize(rayDirection), glm::normalize(intersectionTriangle.intersectedTriangle.normal), refractive_index, kr);
        glm::vec3 bias = 0.0001f * intersectionTriangle.intersectedTriangle.normal;
        glm::vec3 refractionRayOrig = outside ? intersectionTriangle.intersectionPoint - bias : intersectionTriangle.intersectionPoint + bias;

        if(kr < 1){ // refract
            glm::vec3 refractionDirection = glm::normalize(refractt(glm::normalize(rayDirection), glm::normalize(intersectionTriangle.intersectedTriangle.normal), refractive_index));
            colour = castRay(refractionRayOrig, vertices, vertex_normals, triangles, refractionDirection, brightness);
        }else{  // total internal reflection
            glm::vec3 reflectionRayOrig = outside ? intersectionTriangle.intersectionPoint + bias : intersectionTriangle.intersectionPoint - bias;
            glm::vec3 reflectionDirection = glm::normalize(reflectt(glm::normalize(rayDirection), intersectionTriangle.intersectedTriangle.normal));
            colour = castRay(reflectionRayOrig, vertices, vertex_normals, triangles, reflectionDirection, brightness);
        }

    }

    if(!intersectionTriangle.intersectionFound){
        colour.red = 0;
        colour.blue = 0;
        colour.green = 0;
    }

    return colour;
}

void rayTraceObj(DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, vector<Colour> &colour_library, vector<ModelTriangle> &triangles, vector<glm::vec3> &lightsources, vector<glm::vec3> &vertices, vector<glm::vec3> &vertex_normals) {
    for (size_t x = 0; x < window.width; x++){
        for (size_t y = 0; y < window.height; y++){
            // Calculate ray from camera to pixel
            glm::vec3 rayDirection = glm::vec3(x, y, -1.0f);
            rayDirection[0] = ((rayDirection[0] - WIDTH/2.0f)/(SCALER*FOCAL_LENGTH));
            rayDirection[1] = -((rayDirection[1] - HEIGHT/2.0f)/(SCALER*FOCAL_LENGTH));
            rayDirection = rayDirection * glm::inverse(cameraOrientation);

            float brightness = 0;
            Colour colour = castRay(cameraPosition, vertices, vertex_normals, triangles, rayDirection, brightness);
            window.setPixelColour(x, y, (255 << 24) + (int(colour.red * brightness) << 16) + (int(colour.green * brightness) << 8) + int(colour.blue * brightness));
        }
    }
}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, glm::vec3 &lightsource) {
	if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_1) renderer = "wireframe";
        else if (event.key.keysym.sym == SDLK_2) renderer = "rasterised";
        else if (event.key.keysym.sym == SDLK_3) renderer = "ray_traced";
        else if (event.key.keysym.sym == SDLK_SPACE) phong = !phong;
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
    glm::vec3 lightsource = glm::vec3(0.0, 0.3, 0.0);
    lightsources.push_back(lightsource);

    if(soft_shadows) {
        for (int i = 0; i < 20; i++){
            lightsource.x += 0.005;
            lightsource.y += 0.01;
            lightsources.push_back(lightsource);
        }
    }

    glm::mat3 cameraOrientation = glm::mat3(1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1);
    float scaleFactor = 0.15;
    // USE 0.001 scale for hackspace logo model

    // Read obj data from files
    vector<glm::vec3> vertices;
    vector<glm::vec3> sphere_vertices;
    vector<Colour> colour_library = readMTLFile("cornell-box.mtl");
    vector<ModelTriangle> triangles = readOBJFile("cornell-box.obj", scaleFactor, colour_library, vertices, 0.0f);
    //vector<ModelTriangle> triangles = readOBJFile("fox.obj", 0.01, colour_library, vertices, 0.0f);

    // Assign glass to red block and mirror to front face of blue box
    for(size_t i = 0; i < triangles.size(); i++){
        if(glass && (triangles[i].colour.name == "Red")) {
            triangles[i].material = "glass";
        }else if(mirror && (i == 33 || i == 38)){
            triangles[i].material = "mirror";
        }else{
            triangles[i].material = "plastic";
        }
    }

    if(sphere){
        vector<ModelTriangle> sphere_triangles = readOBJFile("sphere.obj", scaleFactor-0.02f, colour_library, sphere_vertices, 0.2f);
        for(ModelTriangle triangle : sphere_triangles){triangles.push_back(triangle);}
        for(glm::vec3 vertex : sphere_vertices){vertices.push_back(vertex);}
    }

    vector<glm::vec3> vertexNormals = getVertexNormals(vertices, triangles);

    while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) {
            handleEvent(event, window, cameraPosition, cameraOrientation, lightsource);
        }

        window.clearPixels();
        if(renderer == "wireframe") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, true);
        else if (renderer == "rasterised") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, false);
        else if (renderer == "ray_traced") rayTraceObj(window, cameraPosition, cameraOrientation, colour_library, triangles, lightsources, vertices, vertexNormals);

        // Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
