#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>
NORI_NAMESPACE_BEGIN
class PathTracingNEE : public Integrator
{
public:
    PathTracingNEE(const PropertyList &props)
    {
        /* No parameters this time */
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f &ray) const
    {
        return Li_rec(scene, sampler, ray, 0, false);
    }

    Color3f Li_rec(const Scene* scene, Sampler* sampler, const Ray3f &ray, int bounce, bool isSpecular) const
    {
        // Find the surface that is visible in the requested direction
        Intersection its;
        if (!scene->rayIntersect(ray, its)){

            if(isSpecular || bounce == 0)
                return scene->getBackground(ray);

            return Color3f(0.);
        }
        else if (its.mesh->isEmitter()) {

            if(isSpecular || bounce == 0){
                const Emitter* em = its.mesh->getEmitter();
                EmitterQueryRecord emRecord(em, ray.o, its.p, its.shFrame.n, its.uv);
                return em->eval(emRecord);
            }

            return Color3f(0.);
        }

        Color3f Le(0.);
        float pdflight;
        EmitterQueryRecord emitterRecord(its.p);
        const Emitter* emit = scene->sampleDirect(sampler->next1D(), pdflight);
        emitterRecord.emitter = emit;
        Color3f Le_emit = emit->sample(emitterRecord, sampler->next2D(), 0.);
        BSDFQueryRecord bsdfRecord_emit(its.toLocal(-ray.d),
            its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);

        Ray3f sray(its.p, emitterRecord.wi);
        Intersection it_shadow;
        if (scene->rayIntersect(sray, it_shadow))
            if (it_shadow.t >= (emitterRecord.dist - 1.e-5))
                Le = Le_emit * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord_emit) / (pdflight * emitterRecord.pdf);

        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.uv);
        Color3f bsdf = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
        Ray3f wo(its.p, its.toWorld(bsdfRecord.wo));

        float prob = 1 - its.mesh->getBSDF()->eval(bsdfRecord).getLuminance();
        if(sampler->next1D() < 1 - prob){
            return Color3f(0.);
        }
        
        bsdf;

        return (Le + bsdf * Li_rec(scene, sampler, wo, ++bounce, bsdfRecord.measure == EDiscrete)) / prob;
    }

    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(PathTracingNEE, "path_nee");
NORI_NAMESPACE_END