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
        // Find the surface that is visible in the requested direction
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return scene->getBackground(ray);
        else if (its.mesh->isEmitter()) {
            const Emitter* em = its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, ray.o, its.p, its.shFrame.n, its.uv);
            return em->eval(emRecord);
        }

        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.uv);
        Color3f bsdf = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
        Ray3f wo(its.p, its.toWorld(bsdfRecord.wo));

        float prob = 1- its.mesh->getBSDF()->eval(bsdfRecord).getLuminance();
        if(sampler->next1D() < 1 - prob){
            return Color3f(0.);
        }
        
        bsdf /= prob;

        return bsdf * Li(scene, sampler, wo);
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(PathTracing, "path");
NORI_NAMESPACE_END