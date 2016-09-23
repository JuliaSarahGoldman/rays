#pragma once
#include <G3D/G3DAll.h>

/** */
class RayTracer {
    public:
        Radiance3 measureLight(const shared_ptr<Scene>& scene, const Ray ray, int numScatters, const TriTree& triangles, const CPUVertexArray& vertices);
        bool findIntersection(const shared_ptr<Surfel>& surfel, const Ray ray, const TriTree& triangles, const CPUVertexArray& vertices);
        bool findSphereIntersection(const shared_ptr<Surfel>& surfel, const Ray ray, const Point3 center, const float radius);
        bool findTriangleIntersection(const shared_ptr<Surfel>& surfel,const Ray ray, const Tri triangle, const CPUVertexArray& vertices);
        void rayTrace(const shared_ptr<Scene>& scene, const shared_ptr<Camera>& cam, const shared_ptr<Image>& image);
        Radiance3 shade(const Ray ray, const shared_ptr<Surfel>& surfel, const Array<shared_ptr<Light>>& lights);
};