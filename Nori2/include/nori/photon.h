#pragma once

NORI_NAMESPACE_BEGIN

class Photon{
    public:
        Photon(){};
        Photon(Point3f p, Vector3f d, Color3f energy, Normal3f n, int light_index) : p(p), d(d), i(energy), n(n),
            light_index(light_index) {
        }

        int light_index;
        Normal3f n;
        Point3f p;
        Vector3f d;
        Color3f i;
};

NORI_NAMESPACE_END