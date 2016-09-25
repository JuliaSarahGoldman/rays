#pragma once
#include <G3D/G3DAll.h>

/** */
class RayTracer {
    protected:
         shared_ptr<TriTree> m_triangles;
    public:
        RayTracer();
        ~RayTracer();
        Radiance3 measureLight(const shared_ptr<Scene>& scene, const Ray ray, int numScatters);
        bool findIntersection(shared_ptr<Surfel>& surfel, const Ray ray);
        bool findSphereIntersection(const shared_ptr<Surfel>& surfel, const Ray ray, const Point3 center, const float radius);
        bool findTriangleIntersection(const Ray ray, const Tri triangle, const CPUVertexArray& vertices, float& t, float b[3], TriTree::Hit hit);
        void rayTrace(const shared_ptr<Scene>& scene, const shared_ptr<Camera>& cam, const shared_ptr<Image>& image);
        Radiance3 shade(const Ray ray, const shared_ptr<Surfel>& surfel, const Array<shared_ptr<Light>>& lights);
        bool isVisible(Point3 X, Point3 Y);
};