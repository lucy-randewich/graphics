#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include "ModelTriangle.h"

struct RayTriangleIntersection {
	glm::vec3 intersectionPoint;
	float distanceFromCamera;
    float u;
    float v;
	ModelTriangle intersectedTriangle;
	size_t triangleIndex;
    bool intersectionFound;

	RayTriangleIntersection();
	RayTriangleIntersection(const glm::vec3 &point, float distance, float u, float v, const ModelTriangle &triangle, size_t index, bool intersectionFound);
	friend std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection);
};
