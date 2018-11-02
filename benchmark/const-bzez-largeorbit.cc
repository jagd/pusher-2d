#include <pusher/pusher.h>
#include <iostream>

int main()
{
    const double Ez = -100e3 / 10e-2;
    const double Bz = 1.0; // f ~ 28 GHz
    const double u = 1e8; // f ~ 26 GHz considering gamma
    const double r = -M0/Q0/Bz * u;
    const double gamma = u2gamma(u);
    const double v = u/gamma;
    const double omega = v/r;

    const auto ef = std::make_shared<ConstEzField>(Ez);
    const auto mf = std::make_shared<ConstBzField>(Bz);
    auto aphi = APhiPusher(ef, mf);
    auto boris = BorisPusher(ef, mf);

    const double dz = 1e-5;

    for (double dt = 1e-13; dt > 1e-15; dt *= 0.5) {
        std::clog << "dt " << dt << '\n';
        const double borisStartAngle = std::atan(dt * omega/2);
        boris.setElectronInfo(r*std::cos(borisStartAngle), r*std::sin(borisStartAngle),0, 0, u, 0);
        aphi.setElectronInfo(0, r, 0, 0, aphi.pTheta(0, r, u), gamma);
        double nextZ = 0;
        while (true) {
            aphi.step(dt);
            boris.step(dt);
            const auto pa = aphi.pos();
            if (pa.z > nextZ) { // every ~10 degree
                const auto pb = boris.pos();
                nextZ = pa.z + dz;
                std::cout << pb.z << ' '
                          << std::sqrt(pb.x*pb.x + pb.y*pb.y) << ' '
                          << pa.z << ' '
                          << pa.r  << ' '
                          << aphi.gammaCurrent() << '\n';
            }
            if (pa.z > 100e-2) {
                break;
            }
        }
        std::cout << '\n';
    }
    return 0;
}
