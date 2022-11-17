#include "RayTriangleIntersection.h"

RayTriangleIntersection::RayTriangleIntersection() = default;
RayTriangleIntersection::RayTriangleIntersection(const glm::vec3 &point, float distance, float u, float v, const ModelTriangle &triangle, size_t index, bool intersectionFound) :
		intersectionPoint(point),
		distanceFromCamera(distance),
        u(u),
        v(v),
		intersectedTriangle(triangle),
		triangleIndex(index),
        intersectionFound(intersectionFound){}

std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection) {
	os << "Intersection is at [" << intersection.intersectionPoint[0] << "," << intersection.intersectionPoint[1] << "," <<
	   intersection.intersectionPoint[2] << "] on triangle " << intersection.intersectedTriangle <<
	   " at a distance of " << intersection.distanceFromCamera;
	return os;
}
