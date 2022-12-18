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
        Intersection its;
        Ray3f iteRay(ray);
        Color3f Le(0.);
        Color3f bsdf(1.);

        bool isSpecular = false;

        for(int bounce = 0;;++bounce){

            if (!scene->rayIntersect(iteRay, its)){

                if(isSpecular || bounce == 0){

                    return Le + scene->getBackground(iteRay) * bsdf;
                }

                return Le;
            }
            else if (its.mesh->isEmitter()) {

                if(isSpecular || bounce == 0){
                    const Emitter* em = its.mesh->getEmitter();
                    EmitterQueryRecord emRecord(em, iteRay.o, its.p, its.shFrame.n, its.uv);
                    return Le + em->eval(emRecord) * bsdf;
                }

                return Le;
            }

            //Sample BSDF
            BSDFQueryRecord bsdfRecord(its.toLocal(-iteRay.d), its.uv);
            Color3f bsdf_aux = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
            isSpecular = bsdfRecord.measure == EDiscrete;
            
            Color3f Le_emiter(0.);
            if(!isSpecular){ // Emitter sampling for NEE
                
                float pdflight;
                EmitterQueryRecord emitterRecord(its.p);
                const Emitter* emit = scene->sampleEmitter(sampler->next1D(), pdflight);
                emitterRecord.emitter = emit;
                Color3f Le_em = emit->sample(emitterRecord, sampler->next2D(), 0.);

                Ray3f sray(its.p, emitterRecord.wi);
                Intersection it_shadow;
                if (scene->rayIntersect(sray, it_shadow))
                    if(it_shadow.t >= (emitterRecord.dist - 1.e-5)){
                    BSDFQueryRecord bsdfRecord_emit(its.toLocal(-iteRay.d),
                        its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                    Le_emiter = Le_em * bsdf * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord_emit) / (pdflight * emitterRecord.pdf);
                }
            }

            bsdf *= bsdf_aux;
            if(bounce > 3){

                float prob = std::min(bsdf_aux.maxCoeff(), 0.95f);
                if(sampler->next1D() < 1 - prob){
                    return Le;
                }

                Le_emiter /= prob;
                bsdf = bsdf / prob;
            }
            
            Le += Le_emiter;

            iteRay = Ray3f(its.p, its.toWorld(bsdfRecord.wo));
        }
        
        return Le;
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(PathTracingNEE, "path_nee");
NORI_NAMESPACE_END