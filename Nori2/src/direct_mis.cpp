#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>
NORI_NAMESPACE_BEGIN
class DirectMIS : public Integrator
{
public:
    DirectMIS(const PropertyList &props)
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

        Color3f emLe = SampleLights(scene, sampler, ray, its);
        Color3f matLe = SampleBSDF(scene, sampler, ray, its);

        Lo += matLe;
        Lo += emLe;

        if (its.mesh->isEmitter()) {
            const Emitter* em = its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, ray.o, its.p, its.shFrame.n, its.uv);
            Lo += em->eval(emRecord);
        }

        return Lo;
    }

    Color3f SampleLights(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& its) const {

        Color3f Loem(0.);
        float pdflight;

        //Sample a ligh source
        EmitterQueryRecord emitterRecord(its.p);
        const std::vector<Emitter*> lights = scene->getLights();
        float rnd = sampler->next1D();
        const Emitter* emit = scene->sampleEmitter(rnd, pdflight);

        Color3f Le = emit->sample(emitterRecord, sampler->next2D(), 0.);
        
        //Check visibility
        Ray3f sray(its.p, emitterRecord.wi);
        Intersection it_shadow;
        if (scene->rayIntersect(sray, it_shadow))
            if (it_shadow.t < (emitterRecord.dist - 1.e-5))
                return Color3f(0.);

        //Compute BSDF with light from sampled light source
        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d),
            its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);

        float emPdf = pdflight * emitterRecord.pdf;
        float matPdf = its.mesh->getBSDF()->pdf(bsdfRecord);
        float base = emPdf + matPdf;
        float w = 0;
        if(base != 0 && !isinf(emPdf))
            w = emPdf / base;

        Loem = Le * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord) / (pdflight * emitterRecord.pdf);

        return Loem * w;
    }

    Color3f SampleBSDF(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& its) const {

        Color3f Lomat(0.);

        //Sample BSDF, result is bsdf evaluation 
        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d));
        Color3f bsdf = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());

        float matPdf  = its.mesh->getBSDF()->pdf(bsdfRecord);

        //Cast ray towards sampled direction
        Color3f Le(0.);
        Intersection next_its;
        float emPdf = 0;
        Ray3f wo(its.p, its.toWorld(bsdfRecord.wo));
        if (!scene->rayIntersect(wo, next_its)){ // Not intersected
            // Background Light
            Le = scene->getBackground(wo);
            if(isnan(Le.x()))
                Le = Color3f(0.);
            emPdf = 1.f;   
        }
        else if(next_its.mesh->isEmitter()){ // Intersected Emitter
            // Emitter light
            const Emitter* em = next_its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, its.p, next_its.p, next_its.shFrame.n, next_its.uv);
            emPdf = em->pdf(emRecord);
            Le = em->eval(emRecord);
        }

        float base = matPdf + emPdf;
        float w = 0;
        if(base != 0)
            w = matPdf / base;


        Lomat = Le * bsdf * w;

        return Lomat;
    }

    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(DirectMIS, "direct_mis");
NORI_NAMESPACE_END