/** \file RayTracer.cpp */
#include "RayTracer.h"

// Will calculate the radiance for a single ray of light by finding the intersection and recursively 
// scattering light
Radiance3 RayTracer::measureLight(const shared_ptr<Scene>& scene, const Ray ray, int numScatters, const TriTree& triangles, const CPUVertexArray& vertices){
  
    //shared_ptr<Surfel> surfel = triangles.intersectRay(ray);
    shared_ptr<Surfel> surfel;
    bool doesIntersect = findIntersection(surfel, ray, triangles, vertices);
    Array<shared_ptr<Light>> lightArray = scene->lightingEnvironment().lightArray;

    if (doesIntersect) return shade(ray, surfel, lightArray);
    return Radiance3(1,0,0);
    // if (findIntersection(surfel, ray, triangles, vertices)) return Radiance3(1,1,1);
    //else return Radiance3(0,0,0);
}   

/*Will iterate through a TriTree, calling findTriangleIntersection for each Tri. If there are primitives, it will call findSphereIntersection for them. 
It should use the surfel it hits first (the one with the shortest distance from the camera.)*/
bool RayTracer::findIntersection(shared_ptr<Surfel>& surfel, const Ray ray, const TriTree& triangles, const CPUVertexArray& vertices){
    Point3 P = ray.origin();
    Vector3 w = ray.direction();
    
    float bar[3];
    float b[3];
    float t;
    float min = INFINITY;
    TriTree::Hit tempHit = TriTree::Hit();
    TriTree::Hit hit = TriTree::Hit();

    for (int i(0); i < triangles.size(); ++i) {
        bool intersects = findTriangleIntersection(ray, triangles[i], vertices, t, b, tempHit);
        if (intersects){
            if (t < min){
                min = t;
                bar[0] = b[0];
                bar[1] = b[1];
                bar[2] = b[2];
                hit = tempHit;
                hit.triIndex = i;
                hit.distance = t;
            }
        }
    }
    if (min == INFINITY){
        return false;
    }
    
    surfel = triangles.sample(hit);

    //return triangles[hit.triIndex].intersectionAlphaTest(vertices, u, v, 1.0f);
    
    return true;
    
}

// Iterates through each pixel, creating a ray for each one from the pixel to the camera aperature, and using measureLight to find the radiance for each pixel
void RayTracer::rayTrace(const shared_ptr<Scene>& scene, const shared_ptr<Camera>& cam, const shared_ptr<Image>& image){
    
    // Get all the surfaces in a scene
    Array<shared_ptr<Surface>> sceneSurfaces;
    scene->onPose(sceneSurfaces);
    
    // Get all the triangles in a scene from all the surfaces
    TriTree sceneTris;
    sceneTris.setContents(sceneSurfaces);
    CPUVertexArray vertices(sceneTris.vertexArray());
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
bool RayTracer::findTriangleIntersection(const Ray ray, const Tri triangle, const CPUVertexArray& vertices, float& t, float b[3], TriTree::Hit hit){
    /*
    Point3 P = ray.origin();
    Vector3 w = ray.direction();

    Vector3 e1 = triangle.e1(vertices);
    Vector3 e2 = triangle.e2(vertices);
    Vector3 vertexNormal(triangle.normal(vertices));
    
    const Vector3 q = w.cross(e2);
    */
    const Point3& P = ray.origin();
    const Vector3& w(ray.direction());
    const Vector3& e1(triangle.e1(vertices));
    const Vector3& e2(triangle.e2(vertices));
    const Vector3& n(triangle.normal(vertices));

    const Vector3& q(w.cross(e2));
    float a = e1.dot(q);

    if (n.dot(w) >= 0 || abs(a) <= 0.0001f) return false;

    const Vector3& s((P-triangle.vertex(vertices,0).position)/a);
    const Vector3& r = s.cross(e1);

    b[0] = s.dot(q);
    b[1] = r.dot(w);
    b[2] = 1.0f - b[0] - b[1];

    if(w.dot(n) > 0){
        hit.backface = true;
    }

    if(b[0] < 0.0f || b[1] < 0.0f || b[2] < 0.0f) return false;

    t = e2.dot(r);
    return (t >= 0.0f);

}

/*Radiance3 RayTracer::shade(const Ray ray, const shared_ptr<Surfel>& surfel, const Array<shared_ptr<Light>>& lights){
    Point3 P = ray.origin()+(.01*ray.direction()) + ray.direction()*abs(surfel->position-ray.origin()); 
    Vector3 w = ray.direction();
    Radiance3 L(surfel->emittedRadiance(P));
    float rVal(L.r);
    float gVal(L.g);
    float bVal(L.b);
    for(int i(0); i < lights.size(); ++i){
        Light  light = *lights[i];
        Vector3 lightVector = light.position().xzy() - P;
        float disToLight = lightVector.length();
        Vector3 w_i = lightVector.direction(); //lightVector/disToLight;
 //       if(light.inFieldOfView(surfel->position)){ 
            //const Radiance3 L_i(light.emittedPower()/(4*pif()*disToLight*disToLight));
            float dotProd = abs(w_i.dot(surfel->shadingNormal));
            Radiance3 bir = light.biradiance(surfel->position);
            Radiance3 surbi = surfel->finiteScatteringDensity(w_i, -1*w);
            //Radiance3 product = L_i * dotProd;
            //product = product * bir;
            Radiance3 product = dotProd * bir * surbi;
            rVal += product.r;
            gVal += product.g;
            bVal += product.b;
//        }
    }
    //return Radiance3(rVal, gVal, bVal);
    return Radiance3(rVal, gVal, bVal);
}*/

Radiance3 RayTracer::shade(const Ray ray, const shared_ptr<Surfel>& surfel, const Array<shared_ptr<Light>>& lights){
    Radiance3 L = surfel->emittedRadiance(-1*ray.direction());
    const Point3& X = surfel->position;
    const Vector3& n = surfel->shadingNormal;

    for(int i = 0; i < lights.size(); ++i) {
        const shared_ptr<Light> light(lights[i]);
        const Point3& Y = light->position().xyz();

        const Vector3& w_i = (Y-X).direction();
        Biradiance3& Bi = light->biradiance(X);

        const Color3& f = surfel->finiteScatteringDensity(w_i,-1*ray.direction());
        L+=Bi * f * abs(w_i.dot(n));
    }
    return L;
}