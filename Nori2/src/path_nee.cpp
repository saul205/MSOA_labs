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
        Color3f Lo(0.);
        Intersection its;
        Ray3f iteRay(ray);
        Color3f Le(0.);
        Color3f bsdf(1.);
        float probs = 1;

        if (!scene->rayIntersect(iteRay, its)) //Return environment
            return (Le + scene->getBackground(ray)) * bsdf;
        else if(its.mesh->isEmitter()){
            const Emitter* em = its.mesh->getEmitter();
            EmitterQueryRecord emRecord(em, iteRay.o, its.p, its.shFrame.n, its.uv);
            Le += em->eval(emRecord);

            return Le * bsdf;
        }

        while(true){
            
            //Sample BSDF
            BSDFQueryRecord bsdfRecord(its.toLocal(-iteRay.d), its.uv);
            Color3f new_bsdf = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());

            if(bsdfRecord.measure != EDiscrete){
                // Emitter sampling for NEE
                float pdflight;
                EmitterQueryRecord emitterRecord(its.p);
                float rnd = sampler->next1D();
                const Emitter* emit = scene->sampleEmitter(rnd, pdflight);
                emitterRecord.emitter = emit;
                Color3f Le_sample = emit->sample(emitterRecord, sampler->next2D(), 0.);

                Ray3f sray(its.p, emitterRecord.wi);
                Intersection it_shadow;
                if (scene->rayIntersect(sray, it_shadow) && it_shadow.t >= (emitterRecord.dist - 1.e-5)){
                    BSDFQueryRecord bsdfRecord_emit(its.toLocal(-iteRay.d),
                    its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                    Le += bsdf * Le_sample * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord_emit) / (pdflight * emitterRecord.pdf);
                }
            }

            float prob = 1 - its.mesh->getBSDF()->eval(bsdfRecord).getLuminance();
            // update ray
            iteRay = Ray3f(its.p, its.toWorld(bsdfRecord.wo));
            if(sampler->next1D() < 1 - prob){
                return Color3f(0.);
            }

            bsdf /= prob;
            bsdf *= new_bsdf;

            if (!scene->rayIntersect(iteRay, its)){ //Return environment

                if(bsdfRecord.measure == EDiscrete)
                    return (Le + scene->getBackground(ray)) * bsdf;

                return Le * bsdf;
            }
            else if(its.mesh->isEmitter()){
                if(bsdfRecord.measure == EDiscrete){
                    const Emitter* em = its.mesh->getEmitter();
                    EmitterQueryRecord emRecord(em, iteRay.o, its.p, its.shFrame.n, its.uv);
                    Le += em->eval(emRecord);

                    return Le * bsdf;
                }
                return Le * bsdf;
            }
        }
        
        return Lo;
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(PathTracingNEE, "path_nee");
NORI_NAMESPACE_END