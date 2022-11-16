#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>
NORI_NAMESPACE_BEGIN
class DirectMaterialSampling : public Integrator
{
public:
    DirectMaterialSampling(const PropertyList &props)
    {
        /* No parameters this time */
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f &ray) const
    {
        Color3f Lo(0.);
        // Find the surface that is visible in the requested direction
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return scene->getBackground(ray);

        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d));
        Color3f bsdf = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
        Color3f Le{0, 0, 0};
        Intersection next_its;
        Ray3f wo(its.p, bsdfRecord.wo);
        if (!scene->rayIntersect(wo, next_its)){
            // Return BSDF with Background
            Le = scene->getBackground(wo);
        }
        else if(next_its.mesh->isEmitter()){
            const Emitter* em = next_its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, wo.o, next_its.p, next_its.shFrame.n, next_its.uv);
            Le = em->eval(emRecord);
        }

        Lo = Le * bsdf;

        if (its.mesh->isEmitter()) {
            const Emitter* em = its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, ray.o, its.p, its.shFrame.n, its.uv);
            Lo += em->eval(emRecord);
        }

        return Lo;
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(DirectMaterialSampling, "direct_mats");
NORI_NAMESPACE_END