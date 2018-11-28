#include <pusher/pusher.h>
#include <iostream>
#include <iomanip>


#ifdef DEMO_BORIS
static void demoBoris(BorisPusher &boris, double zInit)
{
    double every = 1;
    for (double dt = 1e-12; dt > 1e-16; dt *= 0.5) {
        boris.setElectronInfo(1, 0, zInit, 0, 0, 0);
        double counter = 0;
        int iter = 0;
        int64_t i = 0;
        while(true) {
            ++i;
            boris.step(dt);
            const auto p = boris.pos();
            if (++counter == every) {
                std::cout << i*dt << ' ' << p.z << '\n';
                counter = 0;
            }
            if (((iter % 2 == 0) && (p.z > 0)) || ((iter % 2 != 0) && (p.z < 0))) {
                if (++iter > 2) {
                    break;
                }
            }
        }
        std::cout << '\n';
        every += every;

    }
}
#endif


#ifdef DEMO_APHI
static void demoAPhi(LeapFrogPusher2D &aphi, double zInit)
{
    double every = 1;
    for (double dt = 1e-12; dt > 1e-16; dt *= 0.5) {
        aphi.setElectronInfo(zInit, 1, 0, 0, aphi.pTheta(zInit, 1, 0), 1.0);
        double counter = 0;
        int iter = 0;
        int64_t i = 0;
        while(true) {
            ++i;
            aphi.step(dt);
            const auto p = aphi.pos();
            if (++counter == every) {
                std::cout << i*dt << ' ' << p.z << '\n';
                counter = 0;
            }
            if (((iter % 2 == 0) && (p.z > 0)) || ((iter % 2 != 0) && (p.z < 0))) {
                if (++iter > 2) {
                    break;
                }
            }
        }
        std::cout << '\n';
        every += every;

    }
}
#endif



int main()
{
#ifdef VERY_RELATIVISTIC
    const auto ez = 1000e3/10e-2;
#else
    const auto ez = 100e3/10e-2;
#endif
    const auto ef = std::make_shared<MirroredEzField>(ez);
    const auto mf = std::make_shared<ConstBzField>(0);
    auto boris = BorisPusher(ef, mf);
    auto aphi = LeapFrogPusher2D(ef, mf);

    const double zInit = -1e-2;

#ifdef DEMO_BORIS
    demoBoris(boris, zInit);
#elif DEMO_APHI
    demoAPhi(aphi, zInit);
#else

    const int pseudoOrder = 30;
    int64_t maxSteps = static_cast<int64_t>(1) << pseudoOrder;
    const double minDt = 1e-18;
    const double r = 1e-2;
    /*
        // Calculate the pseudo solution
        boris.setElectronInfo(r, 0, zInit, 0, 0, 0);
        std::clog << "Calculating the pseudo-exact solution.\n";
        for (int i = 0; i < maxSteps; ++i) {
            boris.step(minDt);
        }
        const double refZ = boris.pos().z;
        std::clog << std::setprecision(20) <<  refZ << std::endl;
    */
#ifdef VERY_RELATIVISTIC
    const double refZ =  0.0083568139228484463604; // for 1 MV / 10 cm
#else
    const double refZ = -0.0030687973670488798497; // for 100 kV / 10 cm
#endif

    const double totalEnergy = Q0*ef->pot(zInit, 0);

    const int orderOffset = 8;
    int64_t steps = maxSteps >> orderOffset;
    double dt = minDt*(1 << orderOffset);
    for (int c = pseudoOrder-orderOffset; c >=0 ; --c) {
        boris.setElectronInfo(r, 0, zInit, 0, 0, 0);
        aphi.setElectronInfo(zInit, r, 0, 0, aphi.pTheta(zInit, r, 0), 1.0);
        std::clog << "dt = " << dt << '\n';
        for (int i = 0; i < steps; ++i) {
            aphi.step(dt);
            boris.step(dt);
        }
        const auto pBoris = boris.pos();
        const auto pAPhi = aphi.pos();
        const double potEnergyBoris = Q0*ef->pot(pBoris.z, 0);
        const double kinEnergyBoris = (boris.gammaCurrent()-1.0)*(M0*C0*C0);
        std::cout << std::setprecision(18)
                  << dt << ' '
                  << (pBoris.z - refZ)/refZ << ' '
                  << (pAPhi.z - refZ)/refZ  << ' '
                  << std::abs(potEnergyBoris + kinEnergyBoris - totalEnergy) / std::abs(totalEnergy) << '\n';
        dt *= 2;
        steps /= 2;
    }
#endif
    return 0;
}
