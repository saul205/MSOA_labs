/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    throw NoriException("Warp::squareToTent() is not yet implemented!");
}

float Warp::squareToTentPdf(const Point2f &p) {
    throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float rho = sqrt(sample.x());
    float theta = 2 * M_PI * sample.y();
    return Point2f(rho * cos(theta), rho * sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return (p.x() * p.x() + p.y() * p.y() <= 1) ? 1.f / M_PI : 0.0f;
}

Point2f Warp::squareToUniformTriangle(const Point2f& sample) {

    Point2f p;
    if(sample.x() + sample.y() > 1){
        p = Point2f(1-sample.x(), 1-sample.y()); 
    }
    else
        p = sample;

    return p;
}

float Warp::squareToUniformTrianglePdf(const Point2f& p) {
    return ((p.array() >= 0).all() && (p.array() <= 1).all()) && p.x() + p.y() <= 1 ? 1.f / 0.5f : 0.0f;
}


Vector3f Warp::squareToUniformSphere(const Point2f &sample) {

    float theta = 2*M_PIf*sample.x();
    float phi = acos(2*sample.y()-1);
    return Vector3f(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));

}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {

    return ((v.array() >= -1).all() && (v.array() <= 1).all()) ? 1.f / (4.f * M_PIf) : 0.0f;

}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {

    float theta = 2*M_PIf*sample.x();
    float phi = acos(1 - sample.y());
    return Vector3f(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return (v.x() >= -1 && v.x() <= 1 && v.y() <= 1 && v.y() >= -1 && v.z() >= 0 && v.z() <= 1) ? 1.f / (2.f * M_PIf) : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {

    float theta = 2*M_PIf*sample.x();
    float phi = acos(sqrt(sample.y()));
    return Vector3f(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {

    float phi = atan(sqrt(v.x()*v.x() + v.y()*v.y()) / v.z());
    return (v.array() >= -1 && v.array() <= 1).all() && v.z() >= 0 ? cos(phi) / M_PIf : 0;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float alpha_squared = pow(alpha, 2);
    float theta = 2*M_PIf*sample.x();
    float phi = atan(sqrt(-alpha_squared*log(1 - sample.y())));
    return Vector3f(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {

    float alpha_squared = pow(alpha, 2);
    float phi = atan(sqrt(m.x()*m.x() + m.y()*m.y()) / m.z());
    float pdf = exp(-pow(tan(phi), 2) / alpha_squared) / ( M_PIf * alpha_squared * pow(cos(phi), 3));
    return (m.array() >= -1 && m.array() <= 1).all() && m.z() >= 0 ? pdf: 0;
}

NORI_NAMESPACE_END
