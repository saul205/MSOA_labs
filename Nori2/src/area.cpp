/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
	
	v1 - Dec 01 2020
	Copyright (c) 2020 by Adrian Jarabo

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/mesh.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter {
public:
	AreaEmitter(const PropertyList &props) {
		m_type = EmitterType::EMITTER_AREA;
		m_radiance = new ConstantSpectrumTexture(props.getColor("radiance", Color3f(1.f)));
		m_scale = props.getFloat("scale", 1.);
	}

	virtual std::string toString() const {
		return tfm::format(
			"AreaLight[\n"
			"  radiance = %s,\n"
			"  scale = %f,\n"
			"]",
			m_radiance->toString(), m_scale);
	}

	// We don't assume anything about the visibility of points specified in 'ref' and 'p' in the EmitterQueryRecord.
	virtual Color3f eval(const EmitterQueryRecord & lRec) const {
		if (!m_mesh)
			throw NoriException("There is no shape attached to this Area light!");

		if(lRec.n.dot(-lRec.wi) < 0)
		 	return Color3f{0, 0,0};

		return m_radiance->eval(lRec.uv);
	}

	virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample, float optional_u) const {
		if (!m_mesh)
			throw NoriException("There is no shape attached to this Area light!");

		m_mesh->samplePosition(sample, lRec.p, lRec.n, lRec.uv);
        lRec.dist = (lRec.p - lRec.ref).norm();
        lRec.wi = (lRec.p - lRec.ref) / lRec.dist;
		lRec.pdf = pdf(lRec);

		return eval(lRec);
	}

	// Returns probability with respect to solid angle given by all the information inside the emitterqueryrecord.
	// Assumes all information about the intersection point is already provided inside.
	// WARNING: Use with care. Malformed EmitterQueryRecords can result in undefined behavior. 
	//			Plus no visibility is considered.
	virtual float pdf(const EmitterQueryRecord &lRec) const {
		if (!m_mesh)
			throw NoriException("There is no shape attached to this Area light!");

		return m_mesh->pdf(lRec.p) * (lRec.dist * lRec.dist) / abs(lRec.n.dot(lRec.wi));
	}


	// Get the parent mesh
	void setParent(NoriObject *parent)
	{
		auto type = parent->getClassType();
		if (type == EMesh)
			m_mesh = static_cast<Mesh*>(parent);
	}

	// Set children
	void addChild(NoriObject* obj, const std::string& name = "none") {
		switch (obj->getClassType()) {
		case ETexture:
			if (name == "radiance")
			{
				delete m_radiance;
				m_radiance = static_cast<Texture*>(obj);
			}
			else
				throw NoriException("AreaEmitter::addChild(<%s>,%s) is not supported!",
					classTypeName(obj->getClassType()), name);
			break;

		default:
			throw NoriException("AreaEmitter::addChild(<%s>) is not supported!",
				classTypeName(obj->getClassType()));
		}
	}

	Color3f getRadiance() override {
		return m_radiance->eval(Point2f(0.f));
	}

	Ray3f trace_ray(Sampler* sampler, Color3f &energy) const override {
		
		Point3f p;
		Normal3f n;
		Point2f uv;
		float posPdf, dirPdf;
		m_mesh->samplePosition(sampler->next2D(), p, n, uv);
		n.normalize();
		posPdf = m_mesh->pdf(p);
		
		Frame shFrame(n);
		Vector3f d = Warp::squareToCosineHemisphere(sampler->next2D());
		d.normalize();
		dirPdf = Warp::squareToCosineHemispherePdf(d);

		energy = m_radiance->eval(Point2f()) / (posPdf * dirPdf) * abs(n.dot(shFrame.toWorld(d)));

		Ray3f ray(p, shFrame.toWorld(d));
		return ray;
	}
protected:
	Texture* m_radiance;
	float m_scale;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END
