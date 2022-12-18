#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>
NORI_NAMESPACE_BEGIN
class DirectEmitterSampling : public Integrator
{
public:
    DirectEmitterSampling(const PropertyList &props)
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
        float pdflight;
        EmitterQueryRecord emitterRecord(its.p);
        
        float rnd = sampler->next1D();

        const Emitter* emit = scene->sampleEmitter(rnd, pdflight);
        emitterRecord.emitter = emit;

        Color3f Le = emit->sample(emitterRecord, sampler->next2D(), 0.);

        if (its.mesh->isEmitter()) {
            const Emitter* em = its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, ray.o, its.p, its.shFrame.n, its.uv);
            Lo += em->eval(emRecord);
        }
        
        Ray3f sray(its.p, emitterRecord.wi);
        Intersection it_shadow;
        if (scene->rayIntersect(sray, it_shadow))
            if (it_shadow.t < (emitterRecord.dist - 1.e-5))
                return Lo;

        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d),
            its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);

        Lo += Le * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord) / (pdflight * emitterRecord.pdf);

        return Lo;
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(DirectEmitterSampling, "direct_ems");
NORI_NAMESPACE_END
