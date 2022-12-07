#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>
NORI_NAMESPACE_BEGIN
class PathTracingMIS : public Integrator
{
public:
    PathTracingMIS(const PropertyList &props)
    {
        /* No parameters this time */
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f &ray) const
    {
        return Li_rec(scene, sampler, ray, 0, 1);
    }

    float weight(float mainPdf, float auxPdf) const {
        float base = mainPdf + auxPdf;
        float w = 0;
        if(base != 0 && !isinf(mainPdf))
            w = mainPdf / base;

        return w;
    }

    Color3f Li_rec(const Scene* scene, Sampler* sampler, const Ray3f &ray, int bounce, float matPdf) const
    {
        // Find the surface that is visible in the requested direction
        Intersection its;
        if (!scene->rayIntersect(ray, its)){

            float emPdf = 0;
            const Emitter* emEnv = scene->getEnvironmentalEmitter();
            if(emEnv != nullptr){
                EmitterQueryRecord emitter_intersection(
					emEnv, ray.o, its.p, its.shFrame.n, its.uv);
                emPdf = scene->pdfEmitter(emEnv) * emEnv->pdf(emitter_intersection);  
            }
            return scene->getBackground(ray) * (bounce == 0 ? 1 : weight(matPdf, emPdf));
        }
        else if (its.mesh->isEmitter()) {

            const Emitter* em = its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, ray.o, its.p, its.shFrame.n, its.uv);
            float emPdf = em->pdf(emRecord) * scene->pdfEmitter(em);
            return em->eval(emRecord) * (bounce == 0 ? 1 : weight(matPdf, emPdf));
        }

        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.uv);
        Color3f bsdf = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());

        Color3f Le(0.);
        if(bsdfRecord.measure != EDiscrete){
            float pdflight;
            EmitterQueryRecord emitterRecord(its.p);
            const Emitter* emit = scene->sampleEmitter(sampler->next1D(), pdflight);
            emitterRecord.emitter = emit;
            Color3f Le_emit = emit->sample(emitterRecord, sampler->next2D(), 0.);
            BSDFQueryRecord bsdfRecord_emit(its.toLocal(-ray.d),
                its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);

            Ray3f sray(its.p, emitterRecord.wi);
            Intersection it_shadow;
            if (scene->rayIntersect(sray, it_shadow))
                if (it_shadow.t >= (emitterRecord.dist - 1.e-5)){
                    Le = Le_emit * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord_emit) / (pdflight * emitterRecord.pdf);
                    float emPdf = pdflight * emitterRecord.pdf;
                    float matPdf_em = its.mesh->getBSDF()->pdf(bsdfRecord_emit);
                    Le *= weight(emPdf, matPdf_em);
                }
        }

        Ray3f wo(its.p, its.toWorld(bsdfRecord.wo));

        float prob = bsdf.maxCoeff();
        if(prob >= 1){
            prob = 0.9;
        }
        if(sampler->next1D() < 1 - prob){
            return Color3f(0.);
        }

        return (Le + bsdf * Li_rec(scene, sampler, wo, ++bounce, its.mesh->getBSDF()->pdf(bsdfRecord))) / prob;
    }

    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(PathTracingMIS, "path_mis");
NORI_NAMESPACE_END