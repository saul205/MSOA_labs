#include <nori/warp.h>
#include <nori/KDTree.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <random>
#include <tbb/task_scheduler_init.h>
#include <thread>


NORI_NAMESPACE_BEGIN

/*
Process:

- Ray trace from lights to store photons
    - Caustics Map
    - Global Map(does not store caustics not to double sample caustics)

- Render:
    - Direct ilumination from path tracing
    - Indirect ilumination
        - Can render indirect ilumination with path tracing modified not to follow caustics
        - Use the global map (this allows to add caustics trivially) 
    - Add caustics visible on the first pass

Doubts: 
    - Can we use more than 1 bounce when rendering
        - If that is the case, can we average path tracing and pm to obtain smoother results
    - Can we not use the global map?

*/


class PhotonMapper : public Integrator
{
public:
    PhotonMapper(const PropertyList &props)
    {
        /* No parameters this time */
        n_rays = props.getInteger("rays", 100000);
        max_caustics = props.getInteger("caustics", 100000);
        max_global = props.getInteger("photons", 100000);
        m_bounces = props.getInteger("max_bounces", 25);
    }

    void preprocess_worker(const Scene *scene, 
            int n_r, int ind, int max_c, int max_g) {
        //sample a direction
        std::unique_ptr<Sampler> sampler(scene->getSampler()->clone());
        cout << max_c << " " << max_g << endl;
        //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        Intersection its;
        bool isCaustic = false;
        for(int i = 0; i < n_r; ++i){
            float pdflight = 0;
            int index;
            //select light source
            
            const Emitter* emit = scene->sampleDirect(sampler->next1D(), pdflight, index);
            lightCount[ind][index] ++;
            Color3f energy;
            Ray3f ray = emit->trace_ray(sampler.get(), energy);
            energy /= pdflight;
            for(int bounce = 0; bounce <= 50; ++bounce){
                if(!scene->rayIntersect(ray, its) || its.mesh->isEmitter()){
                    break;
                }

                BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.uv);
                Color3f bsdf_aux = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
                

                if(bsdfRecord.measure == EDiscrete){
                    isCaustic = true;
                }
                else if(bounce > 0){
                    Photon photon(its.p, ray.d, energy, its.shFrame.n, index);
                    storePhoton(photon, isCaustic, caustics[ind], photons[ind], max_c, max_g);
                    isCaustic = false;
                }

                energy *= bsdf_aux;
                float prob = std::min(bsdf_aux.sum() / 3, 0.95f);
                if(sampler->next1D() > prob){
                    break;
                }

                energy /= prob;
                
                ray = Ray3f(its.p, its.toWorld(bsdfRecord.wo));
            }

            if(caustics[ind].size() >= max_c && photons[ind].size() >= max_g)
                n_r = i;
        }
    }

    void preprocess(const Scene *scene) override {

        int n_threads = 16;

        lightCount = std::vector<std::vector<long long int>>(n_threads, std::vector<long long int>(scene->getLights().size()));
        caustics = std::vector<std::vector<Photon>>(n_threads, std::vector<Photon>());
        photons = std::vector<std::vector<Photon>>(n_threads, std::vector<Photon>());

        cout << lightCount.size() << endl;
        //sample a direction
        
        std::vector<thread> threads;

        for(int i = 0; i < n_threads; i++){

            threads.push_back(thread( [this, i, scene, n_threads] 
                {preprocess_worker(scene, n_rays / n_threads, i, max_caustics / n_threads, max_global / n_threads); }));
        }

        for(int i = 0; i < n_threads; i++){
            threads[i].join();
        }

        // Merge lists
        std::vector<long long int> light_count(scene->getLights().size());
        for(int i = 0; i < n_threads; i++){

            for(int j = 0; j < lightCount[i].size(); j++){
                light_count[j] += lightCount[i][j];
            }

            for(auto ph : caustics[i]){
                m_caustics.push_back(ph);
            }

            for(auto ph : photons[i]){
                m_photons.push_back(ph);
            }
        }

        caustics.clear();
        photons.clear();
        lightCount.clear();

        cout << m_caustics.size() << " " << m_photons.size() << " " << n_rays << endl;
        cout << light_count[0] << " " << light_count[1] << endl;
        for (Photon ph : m_photons) {
            ph.i /= light_count[ph.light_index];
            tree_global.store(std::vector<float>{ph.p[0], ph.p[1], ph.p[2]}, ph);
        }
        for (Photon ph : m_caustics) {
            ph.i /= light_count[ph.light_index];
            tree_caustics.store(std::vector<float>{ph.p[0], ph.p[1], ph.p[2]}, ph);
        }

        int a = m_caustics.size();
        int b = m_photons.size();

        m_caustics.clear();
        m_photons.clear();

        if (a != 0) {
            tree_caustics.balance();
        }
        if (b != 0) {
            tree_global.balance();
        }

        cout << "AA" << endl;
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f &ray) const
    {
        Intersection its;
        Ray3f iteRay(ray);
        Color3f Le(0.);
        Color3f bsdf(1.);

        bool isSpecular = false;
        for(int bounce = 0; bounce < 50 ;++bounce){

            // Hit a lightsource
            // If bounce == 0, first object intersected is lightsource so MIS weight is not taken into account
            if (!scene->rayIntersect(iteRay, its)){

                auto background = scene->getBackground(iteRay);
                if(!background.isValid()){
                    return Le;
                }

                return Le + scene->getBackground(iteRay) * bsdf;
            }
            else if (its.mesh->isEmitter()) {

                const Emitter* em = its.mesh->getEmitter();
                EmitterQueryRecord emRecord(em, iteRay.o, its.p, its.shFrame.n, its.uv);

                return Le + em->eval(emRecord) * bsdf;
            }else{

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
                        if(it_shadow.t < (emitterRecord.dist - 1.e-5))
                            return estimate_radiance(its) * bsdf_aux * bsdf;

                    BSDFQueryRecord bsdfRecord_emit(its.toLocal(-iteRay.d),
                        its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);

                    // Calculate pdfs for the MIS
                    //emPdf = pdflight * emitterRecord.pdf;
                    //float matPdf_emit = its.mesh->getBSDF()->pdf(bsdfRecord_emit);
                    if(emitterRecord.pdf == 0)
                        Le_emiter = Color3f(0.f);
                    else 
                    Le_emiter = Le_em * bsdf * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord_emit) //* weight(emPdf, matPdf_emit)
                        / (pdflight * emitterRecord.pdf);

                    return estimate_radiance(its) * bsdf_aux * bsdf + Le_emiter;
                }
                
                // Accumulate bsdf for the different bounces / iterations
                bsdf *= bsdf_aux;
                
                // Russian roulette
                // The probability of the ray dying is 1 - maxCoeff(bsdf)
                // We return light accumulated by direct emitter sampling as this is iterative
                // Update ray
                iteRay = Ray3f(its.p, its.toWorld(bsdfRecord.wo));
            }
        }
        
        return Le;
    }

    Color3f kernel(Photon ph, Intersection its, float dist) const {
        return ph.i / (dist*dist*M_PIf) * (ph.p - its.p).norm() / dist;
    }

    Color3f kernel_epanechnikov(Photon ph, Intersection its, float dist) const {
        auto l = ph.i / (dist*dist*M_PIf) * abs(1 - ((ph.p - its.p).norm() / dist) * ((ph.p - its.p).norm() / dist));
        if(isinf(l.x())){
            return Color3f(0.);
        }
        return l;
    }

    Color3f kernel_gaussian(Photon ph, Intersection its, float dist) const {
        float alfa = 0.918;
		float beta = 1.953;
		float coef = 1 / (dist * dist * 3.1416);
        float dp = ((its.p - ph.p).norm());
        float aux = (1. - exp(-beta * (dp * dp) / (2. * dist * dist))) / (1. - exp(-beta));
        return coef * ph.i * alfa * abs(1. - aux);
    }

    Color3f estimate_radiance(Intersection its) const {

        Color3f L{0.};
        float dist = 0.001, dist_caustics = 0.001;
        int progressive_pass = 1;

        std::vector<const KDTree<Photon, 3>::Node*> nodes;
        tree_caustics.find(std::vector<float>{its.p[0], its.p[1], its.p[2]}, 100, nodes, dist_caustics);
        if(dist_caustics != INFINITY)
        for(auto it = nodes.begin(); it != nodes.end(); ++it){
            if((*it)->data().n.dot(its.shFrame.n) >= 0.7)
                L += kernel_gaussian((*it)->data(), its, dist_caustics);
        }
        int n_ph_c = nodes.size();

        
        tree_global.find(std::vector<float>{its.p[0], its.p[1], its.p[2]}, 100, nodes, dist);
        if(dist != INFINITY)
        for(auto it = nodes.begin(); it != nodes.end(); ++it){
            if((*it)->data().n.dot(its.shFrame.n) >= 0.7)
                L += kernel_gaussian((*it)->data(), its, dist);
        }

        float delta = 0.6;
        int n_ph = nodes.size();

        for(int i = 1; i < progressive_pass; i++){

            //Compute the radiance given a kernel
            std::list<const KDTree<Photon, 3>::Node*> nodes;
            tree_caustics.find(std::vector<float>{its.p[0], its.p[1], its.p[2]}, dist_caustics, &nodes);
            for(auto it = nodes.begin(); it != nodes.end(); ++it){
                if((*it)->data().n.dot(its.shFrame.n) >= 0.7)
                    L += kernel_gaussian((*it)->data(), its, dist_caustics);
            }

            int m_ph_c = nodes.size();

            tree_global.find(std::vector<float>{its.p[0], its.p[1], its.p[2]}, dist, &nodes);
            for(auto it = nodes.begin(); it != nodes.end(); ++it){
                if((*it)->data().n.dot(its.shFrame.n) >= 0.7)
                    L += kernel_gaussian((*it)->data(), its, dist);
            }

            int m_ph = nodes.size();

            dist *= sqrt( (n_ph + delta * m_ph) / (n_ph + m_ph));
            dist_caustics *= sqrt((n_ph_c + delta * m_ph_c) / (n_ph_c + m_ph_c));
        }

        return L / progressive_pass;
    }
    std::string toString() const
    {
        return "Direct Emitter Sampler []";
    }

private:
    void storePhoton(Photon ph, bool isCaustic, std::vector<Photon> &caustics, std::vector<Photon> &photons, int max_c, int max_g){
        if(isCaustic && caustics.size() < max_c)
            caustics.push_back(ph);
        else if(!isCaustic && photons.size() < max_g)
            photons.push_back(ph);
    }

        
    std::vector<Photon> m_photons;
    std::vector<Photon> m_caustics;

    KDTree<Photon, 3U> tree_global;
    KDTree<Photon, 3U> tree_caustics;

    std::vector<std::vector<Photon>> caustics;
    std::vector<std::vector<Photon>> photons;
    std::vector<std::vector<long long int>> lightCount;

    long long int n_rays;
    long long int max_caustics; 
    long long int max_global;
    int m_bounces;

};
NORI_REGISTER_CLASS(PhotonMapper, "photon_mapper");
NORI_NAMESPACE_END