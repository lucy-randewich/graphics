#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <RayTriangleIntersection.h>
#include <glm/glm.hpp>
#include "glm/ext.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <math.h>

#define WIDTH (320.0f)*2.0f
#define HEIGHT (240.0f)*2.0f
#define FOCAL_LENGTH 2.0f
#define SCALER WIDTH

using namespace std;

string renderer = "wireframe";
bool soft_shadows = true;
bool phong = true;
bool bump = false;
vector<glm::vec3> lightsources;

TextureMap cobblesTexture("brickwall.ppm");
TextureMap marbleTexture("marble.ppm");
TextureMap bumpFile("brickwall-normal.ppm");
TextureMap environmentMap("spacebox3.ppm");
vector<vector<uint32_t>> environmentPixelMatrix(environmentMap.width);
int frameIndex = 0;
float sphereLocation = -0.1f;
float sphere2Location = -0.2f;
float sphere3Location = 0.0f;
float sphereDirection = -2.5;
float sphere2Direction = -2.5;
float sphere3Direction = -2.5;
int numdirchanges = 0;
int numdirchanges2 = 0;
int numdirchanges3 = 0;
float changeLocation = -0.35;
float changeLocation2 = -0.4;
float changeLocation3 = -0.35;

int frameIndex2 = 1;

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
    float numberOfSteps = max(abs(to.x - from.x), abs(to.y - from.y));
    float xStepSize = (to.x - from.x)/numberOfSteps;
    float yStepSize = (to.y - from.y)/numberOfSteps;
    for (float i=0.0; i<=numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        if(round(x)<window.width && round(y)<window.height && round(x)>0 && round(y)>0){
            window.setPixelColour(round(x), round(y), (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue));
        }
    }
}

void strokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    drawLine(window, triangle.v0(), triangle.v1(), colour);
    drawLine(window, triangle.v1(), triangle.v2(), colour);
    drawLine(window, triangle.v0(), triangle.v2(), colour);
}

vector<ModelTriangle> readTextureOBJFile(string objfile, float scale_factor, vector<Colour> colour_library, vector<glm::vec3> &vertices, glm::vec3 shift) {
    ifstream file(objfile);
    string character;
    float x,y,z;
    string tmpv1, tmpv2, tmpv3;
    int v1, v2, v3;
    string line;
    string colour_name;
    Colour colour;
    colour.name = " ";
    glm::vec3 p0, p1, p2;
    glm::vec3 face_normal;
    vector<tuple<float, float>> vertexTextures;
    vector<ModelTriangle> triangles;

    // Add all vertices to vector and triangles to vector
    while(getline(file, line)) {
        istringstream stream(line);
        stream >> character;

        if(character == "v"){
            stream >> x >> y >> z;
            vertices.push_back({(x * scale_factor)+shift.x, (y * scale_factor)+shift.y, (z * scale_factor)+shift.z});
        }else if(character == "f"){
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
            face_normal = glm::normalize(glm::cross(p0-p2, p1-p2));
            ModelTriangle triangle = ModelTriangle(p0, p1, p2, colour, face_normal, "plastic", "box");
            if(colour.name == "Grass" || colour.name == "Leopard" || colour.name == "Cobbles" || colour.name == "Marble"){
                TextureMap textureFile;
                if(colour.name == "Cobbles"){
                    textureFile = cobblesTexture;
                }else if(colour.name == "Marble"){
                    textureFile = marbleTexture;
                }
                tuple<float, float> texturepoints1 = vertexTextures[stoi(tmpv1.substr(tmpv1.find("/") + 1))-1];
                tuple<float, float> texturepoints2 = vertexTextures[stoi(tmpv2.substr(tmpv2.find("/") + 1))-1];
                tuple<float, float> texturepoints3 = vertexTextures[stoi(tmpv3.substr(tmpv3.find("/") + 1))-1];
                TexturePoint tp1 = TexturePoint(get<0> (texturepoints1) * textureFile.width, get<1> (texturepoints1) * textureFile.height);
                TexturePoint tp2 = TexturePoint(get<0> (texturepoints2) * textureFile.width, get<1> (texturepoints2) * textureFile.height);
                TexturePoint tp3 = TexturePoint(get<0> (texturepoints3) * textureFile.width, get<1> (texturepoints3) * textureFile.height);
                triangle.texturePoints = {tp1, tp2, tp3};
                triangle.material = colour.name;
            }
            triangles.push_back(triangle);
        }else if(character == "usemtl"){
            stream >> colour_name;
            for (Colour curr_colour: colour_library){
                if (curr_colour.name == colour_name) {
                    colour = curr_colour;
                }
            }
        }else if(character == "vt"){
            float vt1, vt2;
            stream >> vt1 >> vt2;
            vertexTextures.push_back(make_tuple(vt1, vt2));
        }
    }
    return triangles;
}

vector<ModelTriangle> readTextureOBJFileWallAnimation(string objfile, float scale_factor, vector<Colour> colour_library, vector<glm::vec3> &vertices, glm::vec3 shift) {
    ifstream file(objfile);
    string character;
    float x,y,z;
    string tmpv1, tmpv2, tmpv3;
    int v1, v2, v3;
    string line;
    string colour_name;
    Colour colour;
    colour.name = " ";
    glm::vec3 p0, p1, p2;
    glm::vec3 face_normal;
    vector<tuple<float, float>> vertexTextures;
    vector<ModelTriangle> triangles;
    int vertexcount = 0;

    // Add all vertices to vector and triangles to vector
    while(getline(file, line)) {
        istringstream stream(line);
        stream >> character;

        if(character == "v"){
            stream >> x >> y >> z;
            if(vertexcount == 2 || vertexcount == 3){
                y -= (frameIndex2 * 0.05);
                z -= (frameIndex2 * 0.05);
                vertices.push_back({(x * scale_factor)+shift.x, (y * scale_factor)+shift.y, (z * scale_factor)+shift.z});
            }else {
                vertices.push_back({(x * scale_factor) + shift.x, (y * scale_factor) + shift.y, (z * scale_factor) + shift.z});
            }
            vertexcount++;
        }else if(character == "f"){
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
            face_normal = glm::normalize(glm::cross(p0-p2, p1-p2));
            ModelTriangle triangle = ModelTriangle(p0, p1, p2, colour, face_normal, "plastic", "box");
            if(colour.name == "Grass" || colour.name == "Leopard" || colour.name == "Cobbles" || colour.name == "Marble"){
                TextureMap textureFile;
                if(colour.name == "Cobbles"){
                    textureFile = cobblesTexture;
                }else if(colour.name == "Marble"){
                    textureFile = marbleTexture;
                }
                tuple<float, float> texturepoints1 = vertexTextures[stoi(tmpv1.substr(tmpv1.find("/") + 1))-1];
                tuple<float, float> texturepoints2 = vertexTextures[stoi(tmpv2.substr(tmpv2.find("/") + 1))-1];
                tuple<float, float> texturepoints3 = vertexTextures[stoi(tmpv3.substr(tmpv3.find("/") + 1))-1];
                TexturePoint tp1 = TexturePoint(get<0> (texturepoints1) * textureFile.width, get<1> (texturepoints1) * textureFile.height);
                TexturePoint tp2 = TexturePoint(get<0> (texturepoints2) * textureFile.width, get<1> (texturepoints2) * textureFile.height);
                TexturePoint tp3 = TexturePoint(get<0> (texturepoints3) * textureFile.width, get<1> (texturepoints3) * textureFile.height);
                triangle.texturePoints = {tp1, tp2, tp3};
                triangle.material = colour.name;
            }
            triangles.push_back(triangle);
        }else if(character == "usemtl"){
            stream >> colour_name;
            for (Colour curr_colour: colour_library){
                if (curr_colour.name == colour_name) {
                    colour = curr_colour;
                }
            }
        }else if(character == "vt"){
            float vt1, vt2;
            stream >> vt1 >> vt2;
            vertexTextures.push_back(make_tuple(vt1, vt2));
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
    float u = focalLength * -projectedVertex.x/projectedVertex.z;
    float v = focalLength * projectedVertex.y/projectedVertex.z;
    float z = focalLength * projectedVertex.z;

    u *= window.width;
    v *= window.width;

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

        if(x < window.width && x > 0 && y < window.height && y > 0){
                if (depth_buffer[x][y] <= point_depth) {
                    if(x<window.width && y<window.height && x>0 && y>0){
                        window.setPixelColour(x, y, (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue));
                    }
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
        glm::vec3 normalsToAverage = glm::vec3(0);
        for (ModelTriangle triangle : triangles){
            if (vertex == triangle.vertices[0] || vertex == triangle.vertices[1] || vertex == triangle.vertices[2]){
                normalsToAverage += triangle.normal;
            }
        }
        vertex_normals.push_back(glm::normalize(normalsToAverage));
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
                closestIntersection = RayTriangleIntersection(intersection, t, u, v, triangle, index, true, triangle.material, triangle.object);
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

glm::vec3 getPerterbedNormal(float x, float y){
    float F_x = bumpFile.pixels[floor(x+1) + bumpFile.width * floor(y)] - bumpFile.pixels[floor(x) + bumpFile.width * floor(y)];
    float F_y = bumpFile.pixels[floor(x) + bumpFile.width * floor(y+1)] - bumpFile.pixels[floor(x) + bumpFile.width * floor(y)];
    float N_x = -F_x/(sqrt(F_x * F_x + F_y * F_y + 1));
    float N_y = -F_y/(sqrt(F_x * F_x + F_y * F_y + 1));
    float N_z = 1/(sqrt(F_x * F_x + F_y * F_y + 1));
    return glm::vec3{N_x, N_y, N_z};
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
    brightness = glm::clamp(brightness, 0.20f, 1.0f);

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

glm::vec3 reflectt(const glm::vec3 &I, const glm::vec3 &N){
    return I - 2 * glm::dot(I, N)* N;
}

glm::vec3 getBumpedNormal(float x, float y){
    int pixel = bumpFile.pixels[round(x) + bumpFile.width * round(y)];
    int red = (pixel & 0x00FF0000) >> 16;
    int green = (pixel & 0x0000FF00) >> 8;
    int blue = (pixel & 0x000000FF);
    return glm::vec3{red, green, blue};
}

Colour castRay(glm::vec3 startPoint, vector<glm::vec3> vertices, vector<glm::vec3> vertex_normals, vector<ModelTriangle> triangles, glm::vec3 rayDirection, float &brightness, int depth){
    if (depth > 5) return Colour(0, 0, 0);
    RayTriangleIntersection intersectionTriangle = getClosestIntersection(startPoint, rayDirection, triangles);
    if(frameIndex >= 210 && frameIndex <= 240 && intersectionTriangle.intersectedTriangle.colour.name == "Grey"){
        intersectionTriangle = getClosestIntersection(intersectionTriangle.intersectionPoint, rayDirection, triangles);
    }
    Colour colour = intersectionTriangle.intersectedTriangle.colour;
    if(intersectionTriangle.intersectionFound){
        // GET NORMAL
        glm::vec3 normal = intersectionTriangle.intersectedTriangle.normal;
        if(intersectionTriangle.intersectedObject == "sphere"){
            normal = getNormal(intersectionTriangle, vertices, vertex_normals);
        }
        if(intersectionTriangle.intersectedMaterial == "Cobbles"){
            normal = getNormal(intersectionTriangle, vertices, vertex_normals);
            // Get barycentric coordinates
            float proportion0 = 1 - (intersectionTriangle.u + intersectionTriangle.v);
            float proportion1 = intersectionTriangle.u;
            float proportion2 = intersectionTriangle.v;
            TexturePoint tp0 = intersectionTriangle.intersectedTriangle.texturePoints[0];
            TexturePoint tp1 = intersectionTriangle.intersectedTriangle.texturePoints[1];
            TexturePoint tp2 = intersectionTriangle.intersectedTriangle.texturePoints[2];

            float x = (proportion0 * tp0.x) + (proportion1 * tp1.x) + (proportion2 * tp2.x);
            float y = (proportion0 * tp0.y) + (proportion1 * tp1.y) + (proportion2 * tp2.y);

            glm::vec3 perterbed_normal = getBumpedNormal(x, y);
            normal = perterbed_normal;
            normal = glm::normalize(normal);
        }
        if(bump && (intersectionTriangle.triangleIndex == 3 || intersectionTriangle.triangleIndex == 4)){
            normal = getNormal(intersectionTriangle, vertices, vertex_normals);
            // Get barycentric coordinates
            float proportion0 = 1 - (intersectionTriangle.u + intersectionTriangle.v);
            float proportion1 = intersectionTriangle.u;
            float proportion2 = intersectionTriangle.v;
            TexturePoint tp0 = intersectionTriangle.intersectedTriangle.texturePoints[0];
            TexturePoint tp1 = intersectionTriangle.intersectedTriangle.texturePoints[1];
            TexturePoint tp2 = intersectionTriangle.intersectedTriangle.texturePoints[2];

            float x = (proportion0 * tp0.x) + (proportion1 * tp1.x) + (proportion2 * tp2.x);
            float y = (proportion0 * tp0.y) + (proportion1 * tp1.y) + (proportion2 * tp2.y);

            glm::vec3 perterbed_normal = getPerterbedNormal(x, y);
            normal = normal + perterbed_normal;
            normal = glm::normalize(normal);
        }

        // SHADOWS
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

        // REFLECTION
        if(intersectionTriangle.intersectedMaterial == "mirror"){
            glm::vec3 reflected_ray = glm::normalize(rayDirection)-(2 * glm::normalize(normal))*(glm::dot(glm::normalize(rayDirection), glm::normalize(normal)));
            colour = castRay(intersectionTriangle.intersectionPoint, vertices, vertex_normals, triangles, reflected_ray, brightness, depth+1);
        }

        // METAL
        if(intersectionTriangle.intersectedMaterial == "gold"){
            glm::vec3 reflected_ray = glm::normalize(rayDirection)-(2 * glm::normalize(normal))*(glm::dot(glm::normalize(rayDirection), glm::normalize(normal)));
            colour = castRay(intersectionTriangle.intersectionPoint, vertices, vertex_normals, triangles, reflected_ray, brightness, depth+1);
            colour.green = glm::min(colour.green + 80, 255);
            colour.red = glm::min(colour.red + 100, 255);
            if(colour.blue > 200){
                colour.blue -= 150;
            }
        }

        // REFRACTION
        if(intersectionTriangle.intersectedMaterial == "glass"){
            float refractive_index = 2.5f;
            bool outside = glm::dot(normal, glm::normalize(rayDirection)) < 0.0f;
            float kr;
            fresnel(glm::normalize(rayDirection), glm::normalize(normal), refractive_index, kr);
            glm::vec3 bias = 0.0001f * normal;
            glm::vec3 refractionRayOrig = outside ? intersectionTriangle.intersectionPoint - bias : intersectionTriangle.intersectionPoint + bias;
            if(kr < 1){ // refract
                glm::vec3 refractionDirection = glm::normalize(refractt(glm::normalize(rayDirection), glm::normalize(normal), refractive_index));
                colour = castRay(refractionRayOrig, vertices, vertex_normals, triangles, refractionDirection, brightness, depth+1);
            }else{  // total internal reflection
                glm::vec3 reflectionRayOrig = outside ? intersectionTriangle.intersectionPoint + bias : intersectionTriangle.intersectionPoint - bias;
                glm::vec3 reflectionDirection = glm::normalize(reflectt(glm::normalize(rayDirection), normal));
                colour = castRay(reflectionRayOrig, vertices, vertex_normals, triangles, reflectionDirection, brightness, depth+1);
            }
            if(intersectionTriangle.intersectedObject == "Bear"){
                // Add red tint to glass
                colour.red = glm::min(colour.red + 80, 200);
            }
        }

        // TEXTURE
        TextureMap textureFile;
        bool applyTexture = false;
        if(intersectionTriangle.intersectedMaterial == "Marble"){
            textureFile = marbleTexture;
            applyTexture = true;
        }else if(intersectionTriangle.intersectedMaterial == "Cobbles"){
            textureFile = cobblesTexture;
            applyTexture = true;
        }
        if(applyTexture){
            float proportion0 = 1 - (intersectionTriangle.u + intersectionTriangle.v);
            float proportion1 = intersectionTriangle.u;
            float proportion2 = intersectionTriangle.v;
            TexturePoint tp0 = intersectionTriangle.intersectedTriangle.texturePoints[0];
            TexturePoint tp1 = intersectionTriangle.intersectedTriangle.texturePoints[1];
            TexturePoint tp2 = intersectionTriangle.intersectedTriangle.texturePoints[2];

            TexturePoint texturePoint;
            texturePoint.x = (proportion0 * tp0.x) + (proportion1 * tp1.x) + (proportion2 * tp2.x);
            texturePoint.y = (proportion0 * tp0.y) + (proportion1 * tp1.y) + (proportion2 * tp2.y);

            uint32_t colour_val = textureFile.pixels[floor(texturePoint.x) + (textureFile.width * floor(texturePoint.y))];
            int red = (colour_val & 0x00FF0000) >> 16;
            int green = (colour_val & 0x0000FF00) >> 8;
            int blue = (colour_val & 0x000000FF);
            colour.red = red;
            colour.green = green;
            colour.blue = blue;
        }

        // ENVIRONMENT MAP
        if(intersectionTriangle.intersectedMaterial == "envMap"){
            brightness = 1;
            glm::vec3 reflected_ray = glm::normalize(rayDirection)-(2 * glm::normalize(normal))*(glm::dot(glm::normalize(rayDirection), glm::normalize(normal)));
            float m = 2.0f * sqrt( pow(reflected_ray.x, 2) + pow(reflected_ray.y, 2) + pow(reflected_ray.z +1, 2));
            float u = ((reflected_ray.x/m) + 0.5f)*environmentMap.width;
            float v = ((reflected_ray.y/m) + 0.5f)*environmentMap.height;

            // Ensure that u and v lie within environmentMap
            if(u>environmentMap.width-1)u = environmentMap.width-1;
            if(v>environmentMap.height-1)v = environmentMap.height-1;
            if(v<0)v=0;
            if(u<0)u=0;
            uint32_t colour_val = environmentPixelMatrix[u][v];
            int red = (colour_val & 0x00FF0000) >> 16;
            int green = (colour_val & 0x0000FF00) >> 8;
            int blue = (colour_val & 0x000000FF);
            colour.red = red;
            colour.green = green;
            colour.blue = blue;
        }
    }

    return colour;
}

void rayTraceObj(DrawingWindow &window, glm::vec3 cameraPosition, glm::mat3 cameraOrientation, vector<Colour> colour_library, vector<ModelTriangle> triangles, vector<glm::vec3> lightsources, vector<glm::vec3> vertices, vector<glm::vec3> vertex_normals) {
    for (size_t x = 0; x < window.width; x++){
        for (size_t y = 0; y < window.height; y++) {
            // Calculate ray from camera to pixel
            glm::vec3 rayDirection = glm::vec3(x, y, -1.0f);
            rayDirection[0] = ((rayDirection[0] - WIDTH/2.0f)/(SCALER*FOCAL_LENGTH));
            rayDirection[1] = -((rayDirection[1] - HEIGHT/2.0f)/(SCALER*FOCAL_LENGTH));
            rayDirection = rayDirection * glm::inverse(cameraOrientation);

            float brightness = 0;
            Colour colour = castRay(cameraPosition, vertices, vertex_normals, triangles, rayDirection, brightness, 0);
            uint32_t colour_value = (255 << 24) + (int(colour.red * brightness) << 16) + (int(colour.green * brightness) << 8) + int(colour.blue * brightness);
            if(x<window.width && y<window.height && x>0 && y>0){
                window.setPixelColour(x, y, colour_value);
            }
        }
    }
}

bool handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, glm::vec3 &lightsource) {
	bool event_happened = false;
    if (event.type == SDL_KEYDOWN) {
        event_happened = true;
        if (event.key.keysym.sym == SDLK_1) renderer = "wireframe";
        else if (event.key.keysym.sym == SDLK_2) renderer = "rasterised";
        else if (event.key.keysym.sym == SDLK_3) renderer = "ray_traced";
		else if (event.key.keysym.sym == SDLK_LEFT) cameraPosition[0] = cameraPosition[0] - 0.1;
		else if (event.key.keysym.sym == SDLK_RIGHT) cameraPosition[0] = cameraPosition[0] + 0.1;
		else if (event.key.keysym.sym == SDLK_UP) cameraPosition[1] = cameraPosition[1] + 0.1;
		else if (event.key.keysym.sym == SDLK_DOWN) cameraPosition[1] = cameraPosition[1] - 0.1;
        else if (event.key.keysym.sym == SDLK_w) {
            cameraPosition[2] = cameraPosition[2] + 0.05;
            lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
        }
        else if (event.key.keysym.sym == SDLK_s) {
            cameraPosition[2] = cameraPosition[2] - 0.05;
            lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
        }
        else if (event.key.keysym.sym == SDLK_j) {             // Pan camera (y axis)
            float theta = glm::radians(3.0);
            glm::mat3 rotate_matrix = glm::mat3(cos(theta), 0.0, -sin(theta),
                                                0.0, 1.0, 0.0,
                                                sin(theta), 0.0, cos(theta));
            cameraPosition = cameraPosition * rotate_matrix;
            lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
        }else if (event.key.keysym.sym == SDLK_k) {             // Pan camera (y axis)
            float theta = glm::radians(-1.0);
            glm::mat3 rotate_matrix = glm::mat3(cos(theta), 0.0, -sin(theta),
                                                0.0, 1.0, 0.0,
                                                sin(theta), 0.0, cos(theta));
            cameraPosition = cameraPosition * rotate_matrix;
            lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
        }else if (event.key.keysym.sym == SDLK_l) {            // Tilt camera (x axis)
            float theta = glm::radians(1.0);
            glm::mat3 rotate_matrix = glm::mat3(1, 0, 0,
                                                0, cos(theta), sin(theta),
                                                0, -sin(theta), cos(theta));
            cameraPosition = cameraPosition * rotate_matrix;
            lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
        }
	}
    return event_happened;
}

vector<ModelTriangle> getTriangles(vector<glm::vec3> &vertices, float scaleFactor, vector<Colour> colour_library){
    vector<ModelTriangle> triangles = readTextureOBJFile("submission-box.obj", scaleFactor, colour_library, vertices, glm::vec3(0,0,0));

    // Add red ball
    vector<glm::vec3> sphere_vertices;
    vector<ModelTriangle> sphere_triangles;
    sphere_triangles = readTextureOBJFile("sphere.obj", scaleFactor-0.08f, colour_library, sphere_vertices, glm::vec3(0.2, sphereLocation, 0.2));
    for(ModelTriangle triangle : sphere_triangles){
        triangle.object = "sphere";
        triangle.material = "mirror";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : sphere_vertices){vertices.push_back(vertex);}

    // Add glass ball
    vector<glm::vec3> sphere2_vertices;
    vector<ModelTriangle> sphere2_triangles;
    sphere2_triangles = readTextureOBJFile("sphere.obj", scaleFactor-0.08f, colour_library, sphere2_vertices, glm::vec3(-0.3, sphere2Location, -0.2));
    for(ModelTriangle triangle : sphere2_triangles){
        triangle.object = "sphere";
        triangle.material = "glass";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : sphere2_vertices){vertices.push_back(vertex);}

    // Add env mapped ball
    vector<glm::vec3> sphere3_vertices;
    vector<ModelTriangle> sphere3_triangles;
    sphere3_triangles = readTextureOBJFile("sphere.obj", scaleFactor-0.08f, colour_library, sphere3_vertices, glm::vec3(0.2, sphere3Location, -0.4));
    for(ModelTriangle triangle : sphere3_triangles){
        triangle.object = "sphere";
        triangle.material = "envMap";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : sphere3_vertices){vertices.push_back(vertex);}

    // Add hackspace logo
    vector<glm::vec3> hackspace_vertices;
    vector<ModelTriangle> hackspace_triangles = readTextureOBJFile("hackspace-logo/logo.obj", 0.0005f, colour_library, hackspace_vertices, glm::vec3(-0.15, -0.23, -0.1));
    for(ModelTriangle triangle : hackspace_triangles){
        triangle.material = "gold";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : hackspace_vertices){vertices.push_back(vertex);}

    // Add plinth for hackspace logo to sit on
    vector<glm::vec3> plinth_vertices;
    vector<ModelTriangle> plinth_triangles = readTextureOBJFile("plinth_non_textured.obj", scaleFactor+0.05f, colour_library, plinth_vertices, glm::vec3(-0.1, -0.41, -0.2));
    for(ModelTriangle triangle : plinth_triangles){
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : plinth_vertices){vertices.push_back(vertex);}

    return triangles;
}

vector<ModelTriangle> getScene2Triangles2(vector<glm::vec3> &vertices, float scaleFactor, vector<Colour> colour_library, string scene){
    vector<ModelTriangle> triangles = readTextureOBJFile(scene, scaleFactor, colour_library, vertices, glm::vec3(0,0,0));

    // Add plinth for glass ball
    vector<glm::vec3> plinth_vertices1;
    vector<ModelTriangle> plinth_triangles1 = readTextureOBJFile("plinth.obj", scaleFactor+0.05f, colour_library, plinth_vertices1, glm::vec3(-0.5, -0.41, -0.2));
    for(ModelTriangle triangle : plinth_triangles1){
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : plinth_vertices1){vertices.push_back(vertex);}

    // Add plinth for mirror ball
    vector<glm::vec3> plinth_vertices2;
    vector<ModelTriangle> plinth_triangles2 = readTextureOBJFile("plinth.obj", scaleFactor+0.05f, colour_library, plinth_vertices2, glm::vec3(0.3, -0.41, -0.2));
    for(ModelTriangle triangle : plinth_triangles2){
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : plinth_vertices2){vertices.push_back(vertex);}

    // Add bear
    vector<glm::vec3> bear_vertices;
    vector<ModelTriangle> bear_triangles;
    bear_triangles = readTextureOBJFile("diamond.obj", 0.002, colour_library, bear_vertices, glm::vec3(-0.4, -0.21, -0.1));
    for(ModelTriangle triangle : bear_triangles){
        triangle.material = "glass";
        triangle.object = "Bear";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : bear_vertices){vertices.push_back(vertex);}

    // Add mirror ball
    vector<glm::vec3> sphere_vertices2;
    vector<ModelTriangle> sphere_triangles2;
    sphere_triangles2 = readTextureOBJFile("sphere.obj", scaleFactor-0.05f, colour_library, sphere_vertices2, glm::vec3(0.4, -0.26, -0.2));
    for(ModelTriangle triangle : sphere_triangles2){
        triangle.object = "sphere";
        triangle.material = "mirror";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : sphere_vertices2){vertices.push_back(vertex);}

    // Add hackspace logo
    vector<glm::vec3> hackspace_vertices;
    vector<ModelTriangle> hackspace_triangles = readTextureOBJFile("hackspace-logo/logo.obj", 0.0005f, colour_library, hackspace_vertices, glm::vec3(-0.15, -0.23, -0.1));
    for(ModelTriangle triangle : hackspace_triangles){
        triangle.material = "gold";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : hackspace_vertices){vertices.push_back(vertex);}

    // Add plinth for hackspace logo to sit on
    vector<glm::vec3> plinth_vertices;
    vector<ModelTriangle> plinth_triangles = readTextureOBJFile("plinth.obj", scaleFactor+0.05f, colour_library, plinth_vertices, glm::vec3(-0.1, -0.41, -0.2));
    for(ModelTriangle triangle : plinth_triangles){
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : plinth_vertices){vertices.push_back(vertex);}

    return triangles;
}

vector<ModelTriangle> getScene2Triangles3(vector<glm::vec3> &vertices, float scaleFactor, vector<Colour> colour_library, string scene){
    vector<ModelTriangle> triangles = readTextureOBJFileWallAnimation("submission-box-scene3.obj", scaleFactor, colour_library, vertices, glm::vec3(0,0,0));

    // Add plinth for glass ball
    vector<glm::vec3> plinth_vertices1;
    vector<ModelTriangle> plinth_triangles1 = readTextureOBJFile("plinth.obj", scaleFactor+0.05f, colour_library, plinth_vertices1, glm::vec3(-0.5, -0.41- (frameIndex2*0.002), -0.2));
    for(ModelTriangle triangle : plinth_triangles1){
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : plinth_vertices1){vertices.push_back(vertex);}

    // Add plinth for mirror ball
    vector<glm::vec3> plinth_vertices2;
    vector<ModelTriangle> plinth_triangles2 = readTextureOBJFile("plinth.obj", scaleFactor+0.05f, colour_library, plinth_vertices2, glm::vec3(0.3, -0.41- (frameIndex2*0.002), -0.2));
    for(ModelTriangle triangle : plinth_triangles2){
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : plinth_vertices2){vertices.push_back(vertex);}

    // Add bear
    vector<glm::vec3> bear_vertices;
    vector<ModelTriangle> bear_triangles;
    bear_triangles = readTextureOBJFile("diamond.obj", 0.002, colour_library, bear_vertices, glm::vec3(-0.4, -0.21- (frameIndex2*0.002), -0.1));
    for(ModelTriangle triangle : bear_triangles){
        triangle.material = "glass";
        triangle.object = "Bear";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : bear_vertices){vertices.push_back(vertex);}

    // Add mirror ball
    vector<glm::vec3> sphere_vertices2;
    vector<ModelTriangle> sphere_triangles2;
    sphere_triangles2 = readTextureOBJFile("sphere.obj", scaleFactor-0.05f, colour_library, sphere_vertices2, glm::vec3(0.4, -0.26- (frameIndex2*0.002), -0.2));
    for(ModelTriangle triangle : sphere_triangles2){
        triangle.object = "sphere";
        triangle.material = "mirror";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : sphere_vertices2){vertices.push_back(vertex);}

    // Add hackspace logo
    vector<glm::vec3> hackspace_vertices;
    vector<ModelTriangle> hackspace_triangles = readTextureOBJFile("hackspace-logo/logo.obj", 0.0005f, colour_library, hackspace_vertices, glm::vec3(-0.15, -0.23 - (frameIndex2*0.002), -0.1));
    for(ModelTriangle triangle : hackspace_triangles){
        triangle.material = "gold";
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : hackspace_vertices){vertices.push_back(vertex);}

    // Add plinth for hackspace logo to sit on
    vector<glm::vec3> plinth_vertices;
    vector<ModelTriangle> plinth_triangles = readTextureOBJFile("plinth.obj", scaleFactor+0.05f, colour_library, plinth_vertices, glm::vec3(-0.1, -0.41 - (frameIndex2*0.002), -0.2));
    for(ModelTriangle triangle : plinth_triangles){
        triangles.push_back(triangle);
    }
    for(glm::vec3 vertex : plinth_vertices){vertices.push_back(vertex);}

    return triangles;
}


int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

    glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 3.5);
    glm::vec3 lightsource = glm::vec3(0.0, 0.0 + 0.15, 0);
    lightsources.push_back(lightsource);

    if(soft_shadows) {
        glm::vec3 light1 = glm::vec3(0.0, 0.01 + 0.15, 0);
        lightsources.push_back(light1);
        glm::vec3 light2 = glm::vec3(0.01, 0.0 + 0.15, 0);
        lightsources.push_back(light2);
        glm::vec3 light3 = glm::vec3(0.01, 0.01 + 0.15, 0);
        lightsources.push_back(light3);
        glm::vec3 light4 = glm::vec3(0.0, -0.01 + 0.15, 0);
        lightsources.push_back(light4);
        glm::vec3 light5 = glm::vec3(-0.01, 0.0 + 0.15, 0);
        lightsources.push_back(light5);
        glm::vec3 light6 = glm::vec3(-0.01, -0.01 + 0.15, 0);
        lightsources.push_back(light6);
        glm::vec3 light7 = glm::vec3(0.0, 0.02 + 0.15, 0);
        lightsources.push_back(light7);
        glm::vec3 light8 = glm::vec3(0.02, 0.0 + 0.15, 0);
        lightsources.push_back(light8);
        glm::vec3 light9 = glm::vec3(0.02, 0.02 + 0.15, 0);
        lightsources.push_back(light9);
        glm::vec3 light10 = glm::vec3(0.0, -0.02 + 0.15, 0);
        lightsources.push_back(light10);
        glm::vec3 light11 = glm::vec3(-0.02, 0.0 + 0.15, 0);
        lightsources.push_back(light11);
        glm::vec3 light12 = glm::vec3(-0.02, -0.02 + 0.15, 0);
        lightsources.push_back(light12);
    }

    for(int i = 0; i < environmentMap.width; i++){environmentPixelMatrix[i] = vector<uint32_t>(environmentMap.height);}
    for(int y = 0; y < environmentMap.height; y++){
        for(int x = 0; x < environmentMap.width; x++){
            environmentPixelMatrix[x][y] = environmentMap.pixels[y * environmentMap.height + x];
        }
    }

    glm::mat3 cameraOrientation = glm::mat3(1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1);
    float scaleFactor = 0.15;

    vector<Colour> colour_library = readMTLFile("textured-cornell-box.mtl");

    vector<glm::vec3> vertices;
    vector<ModelTriangle> triangles = getTriangles(vertices, scaleFactor, colour_library);
    vector<glm::vec3> vertexNormals = getVertexNormals(vertices, triangles);

    if(renderer == "wireframe") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, true);
    else if (renderer == "rasterised") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, false);
    else if (renderer == "ray_traced") rayTraceObj(window, cameraPosition, cameraOrientation, colour_library, triangles, lightsources, vertices, vertexNormals);
    window.renderFrame();

    frameIndex = 

    while(true){
    //while (frameIndex <= 554) {
        if (window.pollForInputEvents(event)) {
            bool event_happened = handleEvent(event, window, cameraPosition, cameraOrientation, lightsource);
        }else{
            window.clearPixels();

            //if(frameIndex == 10){renderer = "rasterised";}
            //if(frameIndex == 17){renderer = "ray_traced";}
            //if(frameIndex == 356){renderer = "rasterised";}

            if(frameIndex <= 60){
                cameraPosition[2] = cameraPosition[2] - 0.05;
                if(sphereLocation <= -0.45 || (numdirchanges %2 != 0 && sphereLocation >= changeLocation)){numdirchanges++;changeLocation-=0.015;sphereDirection *= -1;}
                if(sphere2Location <= -0.45 || (numdirchanges2 %2 != 0 && sphere2Location >= changeLocation2)){numdirchanges2++;changeLocation2-=0.015;sphere2Direction *= -1;}
                if(sphere3Location <= -0.45 || (numdirchanges3 %2 != 0 && sphere3Location >= changeLocation3)){numdirchanges3++;changeLocation3-=0.015;sphere3Direction *= -1;}
                if(sphereLocation >= -0.45 || numdirchanges<=5){sphereLocation += sphereDirection * 0.01;}
                if(sphere2Location >= -0.45 || numdirchanges2<=5){sphere2Location += sphere2Direction * 0.01;}
                if(sphere3Location >= -0.45 || numdirchanges3<=5){sphere3Location += sphere3Direction * 0.01;}
                vertices.clear();
                triangles.clear();
                vertexNormals.clear();
                triangles = getTriangles(vertices, scaleFactor, colour_library);
                vertexNormals = getVertexNormals(vertices, triangles);
            }else if(frameIndex >= 61 && frameIndex <= 100){
                cameraPosition[2] = cameraPosition[2] + 0.05;
                vertices.clear();
                triangles.clear();
                vertexNormals.clear();
                triangles = getScene2Triangles2(vertices, scaleFactor, colour_library, "submission-box-scene2.obj");
                vertexNormals = getVertexNormals(vertices, triangles);
            }else if(frameIndex >= 101 && frameIndex <= 356){
                float theta = glm::radians(-1.4);
                glm::mat3 rotate_matrix = glm::mat3(cos(theta), 0.0, -sin(theta),
                                                    0.0, 1.0, 0.0,
                                                    sin(theta), 0.0, cos(theta));
                cameraPosition = cameraPosition * rotate_matrix;
                lookAt(glm::vec3(0, 0, 0), cameraOrientation, cameraPosition);
                cameraPosition[0] = cameraPosition[0] - 0.025;
                vertices.clear();
                triangles.clear();
                vertexNormals.clear();
                if(frameIndex >= 240){
                    triangles = getScene2Triangles2(vertices, scaleFactor, colour_library, "submission-box-scene3.obj");
                }else{
                    triangles = getScene2Triangles2(vertices, scaleFactor, colour_library, "submission-box-scene2.obj");
                }
                vertexNormals = getVertexNormals(vertices, triangles);
            }else if(frameIndex > 356 && frameIndex <= 485){
                vertices.clear();
                triangles.clear();
                vertexNormals.clear();
                triangles = getScene2Triangles3(vertices, scaleFactor, colour_library, "submission-box-scene3.obj");
                frameIndex2++;
            }
            if(frameIndex == 500){
                vertices.clear();
                triangles.clear();
                vertexNormals.clear();
            }
            if(frameIndex == 356){
                cameraOrientation = glm::mat3(1, 0, 0,
                                              0, 1, 0,
                                              0, 0, 1);
                cameraPosition.x = 0;
                cameraPosition.y = 0;
            }
            if(frameIndex > 355 && frameIndex <= 554){
                cameraPosition[2] = cameraPosition[2] - 0.01;
                sphereLocation = -0.1f;
                sphere2Location = -0.2f;
                sphere3Location = 0.0f;
                sphereDirection = -2.5;
                sphere2Direction = -2.5;
                sphere3Direction = -2.5;
                numdirchanges = 0;
                numdirchanges2 = 0;
                numdirchanges3 = 0;
                changeLocation = -0.35;
                changeLocation2 = -0.4;
                changeLocation3 = -0.35;
                vector<glm::vec3> vertices2;
                vector<ModelTriangle> triangles2 = getTriangles(vertices2, scaleFactor, colour_library);
                rasteriseObj(window, cameraPosition + glm::vec3(0, 0, 3), cameraOrientation, colour_library, triangles2, true);
            }

            if(renderer == "wireframe") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, true);
            else if (renderer == "rasterised") rasteriseObj(window, cameraPosition, cameraOrientation, colour_library, triangles, false);
            else if(renderer == "ray_traced") {rayTraceObj(window, cameraPosition, cameraOrientation, colour_library, triangles,lightsources, vertices, vertexNormals);}
            window.renderFrame();

            ostringstream filename;
            filename << "images/" << std::setw(5) << std::setfill('0') << frameIndex << ".ppm";
            window.savePPM(filename.str().c_str());

            frameIndex++;
        }
	}
}
