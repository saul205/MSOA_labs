#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN
class PointEmitter : public Emitter
{
public:
    PointEmitter(const PropertyList &props)
    {
        m_type = EmitterType::EMITTER_POINT;
        m_position = props.getPoint("position", Point3f(0., 100., 0.));
        m_radiance = props.getColor("radiance", Color3f(1.f));
    }
    virtual std::string toString() const
    {
        return tfm::format(
            "PointEmitter[\n"
            " position = %s,\n"
            " radiance = %s,\n"
            "]",
            m_position.toString(),
            m_radiance.toString());
    }
    virtual Color3f eval(const EmitterQueryRecord &lRec) const
    {
        // This function assumes that a ray have been traced towards
        // the light source. However, since the probability of randomly
        // sampling a point in space is 0, its evaluation returns 0.
        return 0.;
    }
    virtual Color3f sample(EmitterQueryRecord &lRec,
                           const Point2f &sample,
                           float optional_u) const
    {
        lRec.p = m_position;
        lRec.dist = (lRec.p - lRec.ref).norm();
        lRec.wi = (lRec.p - lRec.ref) / lRec.dist;
        // Note that the pdf should be infinite, but for numerical
        // reasons it is more convenient to just leave as 1
        lRec.pdf = 1.;
        // Note that here it is assumed perfect visibility; this means
        // that visibility should be taken care of in the integrator.
        return m_radiance / (lRec.dist*lRec.dist);
    } // Note that the pdf should be infinite, but for numerical reasons
    // it is more convenient to just leave as 1
    virtual float pdf(const EmitterQueryRecord &lRec) const
    {
        return 1.;
    }

    Ray3f trace_ray(Sampler* sampler, Color3f &energy) const override {
		energy = m_radiance;

		Vector3f d = Warp::squareToUniformSphere(sampler->next2D());

        energy /= Warp::squareToUniformSpherePdf(d);

		Ray3f ray(m_position, d);
		return ray;
	}

    Color3f getRadiance() override {
		return m_radiance;
	}

protected:
    Point3f m_position;
    Color3f m_radiance;
};
NORI_REGISTER_CLASS(PointEmitter, "pointlight")
NORI_NAMESPACE_END