#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN
class SpotlightEmitter : public Emitter
{
public:
    SpotlightEmitter(const PropertyList &props)
    {
        m_type = EmitterType::EMITTER_POINT;
        m_position = props.getPoint("position", Point3f(0., 100., 0.));
        m_radiance = props.getColor("radiance", Color3f(1.f));
        m_cameraToWorld = props.getTransform("toWorld", Transform());
        m_distance = props.getFloat("distance", 10.f);
    }
    virtual std::string toString() const
    {
        Point3f p(0, 0, m_distance);
        Vector3f dir = (p - Point3f(0, 0, 0)).normalized();
        Vector3f world = m_cameraToWorld * dir;

        Point3f target(-2.4, -0.268277, -1.4);
        Vector3f wi = (m_position - target).normalized();

        return tfm::format(
            "PointEmitter[\n"
            " position = %s,\n"
            " radiance = %s,\n"
            " matrix = %s, \n"
            " z = %s, \n"
            " dir = %s, \n"
            " dot = %s, \n"
            "]",
            m_position.toString(),
            m_radiance.toString(),
            m_cameraToWorld.toString(),
            world.toString(),
            wi.toString(),
            abs(wi.dot(world)));
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
        lRec.pdf = pdf(lRec);
        // Note that here it is assumed perfect visibility; this means
        // that visibility should be taken care of in the integrator.
        return m_radiance / (lRec.dist);
    } // Note that the pdf should be infinite, but for numerical reasons
    // it is more convenient to just leave as 1
    virtual float pdf(const EmitterQueryRecord &lRec) const
    {
        Point3f p(0, 0, m_distance);
        Vector3f dir = (p - Point3f(0.f, 0.f, 0.f)).normalized();
        Vector3f world = m_cameraToWorld * dir;
        world = world.normalized();

        float refAngle = m_distance / sqrt(m_distance*m_distance + 1);
        float angle = lRec.wi.dot(-world);

        return refAngle <= angle ? 1.f / M_PI : 0.;
    }

    Ray3f trace_ray(Sampler* sampler, Color3f &energy) const override {
		energy = m_radiance;

		Point2f d = Warp::squareToUniformDisk(sampler->next2D());
        Point3f p(d.x(), d.y(), m_distance);

        Vector3f dir = (m_cameraToWorld * p - m_position).normalized();

        energy /= Warp::squareToUniformDiskPdf(d);

		Ray3f ray(m_position, dir);
		return ray;
	}

protected:
    Point3f m_position;
    Color3f m_radiance;
    Transform m_cameraToWorld;
    float m_distance;
};
NORI_REGISTER_CLASS(SpotlightEmitter, "spotlight")
NORI_NAMESPACE_END