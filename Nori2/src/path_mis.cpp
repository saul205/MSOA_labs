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

    float weight(float matPdf, float emPdf) const {
        float base = matPdf + emPdf;
        float w = 0;
        if(base != 0)
            w = matPdf / base;

        return w;
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f &ray) const
    {
        Intersection its;
        Ray3f iteRay(ray);
        Color3f Le(0.);
        Color3f bsdf(1.);

        float emPdf = 0, matPdf = 1;
        bool isSpecular = false;
        for(int bounce = 0;;++bounce){
            emPdf = 0;
            if (!scene->rayIntersect(iteRay, its)){

                if(isSpecular || bounce == 0){

                    const Emitter* emEnv = scene->getEnvironmentalEmitter();
                    if(emEnv != nullptr){
                        EmitterQueryRecord emitter_intersection(
                            emEnv, iteRay.o, its.p, its.shFrame.n, its.uv);
                        emPdf = scene->pdfEmitter(emEnv) * emEnv->pdf(emitter_intersection);  
                    }

                    return Le + scene->getBackground(iteRay) * bsdf * weight(matPdf, emPdf);
                }
                    
                return Le;
            }
            else if (its.mesh->isEmitter()) {

                if(isSpecular || bounce == 0){
                    const Emitter* em = its.mesh->getEmitter();
                    EmitterQueryRecord emRecord(em, iteRay.o, its.p, its.shFrame.n, its.uv);
                    emPdf = em->pdf(emRecord) * scene->pdfEmitter(em);
                    return Le + em->eval(emRecord) * bsdf * weight(matPdf, emPdf);
                }

                return Le;
            }
            
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

                    emPdf = pdflight * emitterRecord.pdf;
                    float matPdf_emit = its.mesh->getBSDF()->pdf(bsdfRecord_emit);
                    Le_emiter = Le_em * bsdf * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord_emit) * weight(matPdf_emit, emPdf)
                        / (pdflight * emitterRecord.pdf);
                }
            }

            //Sample BSDF
            BSDFQueryRecord bsdfRecord(its.toLocal(-iteRay.d), its.uv);
            Color3f bsdf_aux = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
            bsdf *= bsdf_aux;

            matPdf = its.mesh->getBSDF()->pdf(bsdfRecord);

            float prob = bsdf_aux.maxCoeff();
            if(sampler->next1D() < 1 - prob){
                return Le;
            }

            Le += Le_emiter / prob;
            bsdf = bsdf / prob;

            iteRay = Ray3f(its.p, its.toWorld(bsdfRecord.wo));
            isSpecular = bsdfRecord.measure == EDiscrete;
        }
        
        return Le;
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }
};
NORI_REGISTER_CLASS(PathTracingMIS, "path_mis");
NORI_NAMESPACE_END