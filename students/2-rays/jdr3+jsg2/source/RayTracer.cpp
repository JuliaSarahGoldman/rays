/** \file RayTracer.cpp */
#include "RayTracer.h"

// Will calculate the radiance for a single ray of light by finding the intersection and recursively 
// scattering light
Radiance3 RayTracer::measureLight(const shared_ptr<Scene>& scene, const Ray ray, int numScatters, const TriTree& triangles, const CPUVertexArray& vertices){
  
    shared_ptr<Surfel> surfel = triangles.intersectRay(ray);
    Array<shared_ptr<Light>> lightArray = scene->lightingEnvironment().lightArray;

    if (surfel) return 1000.0f*shade(ray, surfel, lightArray);
    return Radiance3(1,0,0);
    // if (findIntersection(surfel, ray, triangles, vertices)) return Radiance3(1,1,1);
    //else return Radiance3(0,0,0);
}   

/*Will iterate through a TriTree, calling findTriangleIntersection for each Tri. If there are primitives, it will call findSphereIntersection for them. 
It should use the surfel it hits first (the one with the shortest distance from the camera.)*/
bool RayTracer::findIntersection(const shared_ptr<Surfel>& surfel, const Ray ray, const TriTree& triangles, const CPUVertexArray& vertices){
    Point3 P = ray.origin();
    Vector3 w = ray.direction();
    
    for (int i(0); i < triangles.size(); ++i) {
        if (findTriangleIntersection(surfel, ray, triangles[i], vertices)){
            return true;
        }
    }
    return false;
}

// Iterates through each pixel, creating a ray for each one from the pixel to the camera aperature, and using measureLight to find the radiance for each pixel
void RayTracer::rayTrace(const shared_ptr<Scene>& scene, const shared_ptr<Camera>& cam, const shared_ptr<Image>& image){
    
    // Get all the surfaces in a scene
    Array<shared_ptr<Surface>> sceneSurfaces;
    scene->onPose(sceneSurfaces);
    
    // Get all the triangles in a scene from all the surfaces
    TriTree sceneTris;
    CPUVertexArray vertices;
    sceneTris.setContents(sceneSurfaces);

    //Now iterate through all of the pixels
    for(int x(0); x < image->width();++x){
        for(int y(0); y<image->height();++y){
            float imWidth = image->width();
            float imHeight = image->height();
            Rect2D rect2D(Vector2(imWidth-1, imHeight-1));
            Ray ray = cam->worldRay(x+.5f,y+.5f,rect2D); //Maybe add .5 to x and y?

            Radiance3 radiance = measureLight(scene, ray, 2, sceneTris, vertices);
            image->set(Point2int32(x,y), radiance);
        }
    }
}

// Find the interesection of a ray to a sphere.
bool RayTracer::findSphereIntersection(const shared_ptr<Surfel>& surfel, const Ray ray, const Point3 center, const float radius){
    Point3 P = ray.origin();
    Vector3 w = ray.direction();

    return true;
}

// find the intersection of a ray to a triangle
bool RayTracer::findTriangleIntersection(const shared_ptr<Surfel>& surfel, const Ray ray, const Tri triangle, const CPUVertexArray& vertices){
    /*
    Point3 P = ray.origin();
    Vector3 w = ray.direction();

    Vector3 e1 = triangle.e1(vertices);
    Vector3 e2 = triangle.e2(vertices);
    Vector3 vertexNormal(triangle.normal(vertices));
    
    const Vector3 q = w.cross(e2);
    */

    Triangle realTriangle = triangle.toTriangle(vertices);
    float bar = 0.0f;
    float* baro = &bar;
    float dis = 100000.0f;
    return realTriangle.intersect(ray, dis, baro);
}

Radiance3 RayTracer::shade(const Ray ray, const shared_ptr<Surfel>& surfel, const Array<shared_ptr<Light>>& lights){
    Point3 P = ray.origin()+ray.direction()*abs(surfel->position-ray.origin());
    Vector3 w = ray.direction();
    Radiance3 surfrad = surfel->emittedRadiance(P); 
    Radiance3 L(surfel->emittedRadiance(P));
    float rVal(L.r);
    float gVal(L.g);
    float bVal(L.b);
    for(int i(0); i < lights.size(); ++i){
        Light  light = *lights[i];
        Vector3 lightVector = light.position().xzy() - P;
        float disToLight = lightVector.length();
        Vector3 w_i = lightVector/disToLight;
        if(light.inFieldOfView(surfel->position)){ 
            const Radiance3 L_i(light.emittedPower()/(4*pif()*disToLight*disToLight));
            float dotProd = abs(w_i.dot(surfel->shadingNormal));
            Radiance3 bir = light.biradiance(-1*surfel->position);
            Radiance3 surbi = surfel->finiteScatteringDensity(w_i, -1 * w);
            Radiance3 product = L_i * dotProd;
            product = product * bir;
            product = product * surbi;
            rVal += product.r;
            gVal += product.g;
            bVal += product.b;
        }
    }
    //return Radiance3(rVal, gVal, bVal);
    return Radiance3(rVal, gVal, bVal);
}