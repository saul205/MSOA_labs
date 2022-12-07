#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>
NORI_NAMESPACE_BEGIN
class PathTracing : public Integrator
{
public:
    PathTracing(const PropertyList &props)
    {
        /* No parameters this time */
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f &ray) const
    {
        Color3f Lo(0.);
        Intersection its;
        Ray3f iteRay(ray);
        Color3f Le(0.);
        Color3f bsdf(1.);
        while(true){

            if (!scene->rayIntersect(iteRay, its)) //Return environment
                return scene->getBackground(ray) * bsdf;
            else if(its.mesh->isEmitter()){
                const Emitter* em = its.mesh->getEmitter();
                EmitterQueryRecord emRecord(em, iteRay.o, its.p, its.shFrame.n, its.uv);
                Le = em->eval(emRecord);

                return Le * bsdf;
            }

            //Sample BSDF
            BSDFQueryRecord bsdfRecord(its.toLocal(-iteRay.d), its.uv);
            Color3f bsdf_aux = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
            bsdf *= bsdf_aux;

            float prob = bsdf_aux.maxCoeff(); 
            if(prob >= 1)
                prob = 0.9;
            if(sampler->next1D() < 1 - prob){
                return Color3f(0.);
            }

            iteRay = Ray3f(its.p, its.toWorld(bsdfRecord.wo));
            bsdf /= prob;
        }
        
        return Lo;
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(PathTracing, "path");
NORI_NAMESPACE_END