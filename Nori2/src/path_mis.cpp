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

    float weight(float mainPdf, float auxPdf) const {

        if (isinf(mainPdf) || (mainPdf == 0 && auxPdf == 0))
            return 1;

        float base = mainPdf + auxPdf;
        return mainPdf / base;
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f &ray) const
    {
        Intersection its;
        Ray3f iteRay(ray);
        Color3f Le(0.);
        Color3f bsdf(1.);

        float emPdf = 0, matPdf = 0;
        bool isSpecular = false;
        for(int bounce = 0;;++bounce){

            // Hit a lightsource
            // If bounce == 0, first object intersected is lightsource so MIS weight is not taken into account
            if (!scene->rayIntersect(iteRay, its)){

                const Emitter* emEnv = scene->getEnvironmentalEmitter();
                if(emEnv != nullptr){
                    EmitterQueryRecord emitter_intersection(
                        emEnv, iteRay.o, its.p, its.shFrame.n, its.uv);
                    emPdf = scene->pdfEmitter(emEnv) * emEnv->pdf(emitter_intersection);  
                }

                // Return accumulated light + 
                //      Background * bsdf accumulated * weight
                //  The matPdf comes from bsdf evaluation on previous iteration
                return Le + scene->getBackground(iteRay) * bsdf 
                    * (bounce == 0 || isSpecular ? 1 : weight(matPdf, emPdf));
            }
            else if (its.mesh->isEmitter()) {

                const Emitter* em = its.mesh->getEmitter();
                EmitterQueryRecord emRecord(em, iteRay.o, its.p, its.shFrame.n, its.uv);
                emPdf = em->pdf(emRecord) * scene->pdfEmitter(em);

                // Return accumulated light + 
                //      light evaluation * bsdf accumulated * weight
                //  The matPdf comes from bsdf evaluation on previous iteration
                return Le + em->eval(emRecord) * bsdf 
                    * (bounce == 0 || isSpecular ? 1 : weight(matPdf, emPdf));
            }

            //Sample BSDF
            BSDFQueryRecord bsdfRecord(its.toLocal(-iteRay.d), its.uv);
            Color3f bsdf_aux = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
            isSpecular = bsdfRecord.measure == EDiscrete;
            
            // If its not specular, sample a light source
            // The reason is that direct sampling an emitter doesn't make sense with specular materials
            // as the probability of sampling that direction is approximated to 0
            Color3f Le_emiter(0.);
            if(!isSpecular){ // Emitter sampling
                
                float pdflight;
                EmitterQueryRecord emitterRecord(its.p);
                const Emitter* emit = scene->sampleEmitter(sampler->next1D(), pdflight);
                emitterRecord.emitter = emit;
                Color3f Le_em = emit->sample(emitterRecord, sampler->next2D(), 0.);

                // Visibility check
                Ray3f sray(its.p, emitterRecord.wi);
                Intersection it_shadow;
                if (scene->rayIntersect(sray, it_shadow))
                    if(it_shadow.t >= (emitterRecord.dist - 1.e-5)){
                    BSDFQueryRecord bsdfRecord_emit(its.toLocal(-iteRay.d),
                        its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);

                    // Calculate pdfs for the MIS
                    emPdf = pdflight * emitterRecord.pdf;
                    float matPdf_emit = its.mesh->getBSDF()->pdf(bsdfRecord_emit);
                    Le_emiter = Le_em * bsdf * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord_emit) * weight(emPdf, matPdf_emit)
                        / (pdflight * emitterRecord.pdf);
                }
            }
            
            // Accumulate bsdf for the different bounces / iterations
            bsdf *= bsdf_aux;
            // Compute material pdf for the next iteration (in case it intersects light emitter prepare the matPdf in advance)
            matPdf = its.mesh->getBSDF()->pdf(bsdfRecord);

            // Russian roulette
            // The probability of the ray dying is 1 - maxCoeff(bsdf)
            // We return light accumulated by direct emitter sampling as this is iterative
            float prob = bsdf_aux.maxCoeff();
            if(sampler->next1D() < 1 - prob){
                return Le;
            }

            // Account for the russian roulette probs
            Le += Le_emiter / prob;
            bsdf = bsdf / prob;

            // Update ray
            iteRay = Ray3f(its.p, its.toWorld(bsdfRecord.wo));
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