#include <pusher/pusher.h>
#include <iostream>
#include <iomanip>


#ifdef DEMO
static void demoBoris(BorisPusher &boris)
{
    double every = 1;
    for (double dt = 1e-13; dt > 1e-16; dt *= 0.5) {
        boris.setElectronInfo(1, 0, -10e-2, 0, 0, 0);
        double counter = 0;
        int64_t i = 0;
        while(true) {
            ++i;
            boris.step(dt);
            const auto p = boris.pos();
            if (++counter == every) {
                std::cout << i*dt << ' ' << p.z << '\n';
                counter = 0;
            }
            if (p.z > 0) {
                break;
            }
        }
        std::cout << '\n';
        every += every;

    }
}
#endif


int main()
{
    const auto ef = std::make_shared<MirOscZEField>(10e3/10e-2);
    const auto mf = std::make_shared<HomogeneousMagField>(0);
    auto boris = BorisPusher(ef, mf);

#ifdef DEMO
    demoBoris(boris);
#else
    auto aphi = APhiPusher(ef, mf);

    const int pseudoOrder = 30;
    int64_t maxSteps = static_cast<int64_t>(1) << pseudoOrder;
    const double minDt = 1e-18;
    const double r = 1e-2;
    /*
    boris.setElectronInfo(r, 0, -10e-2, 0, 0, 0);
    std::clog << "Calculating the pseudo-exact solution.\n";
    for (int i = 0; i < maxSteps; ++i) {
        boris.step(minDt);
    }
    const double refZ = boris.pos().z;
    std::clog << std::setprecision(20) <<  refZ << std::endl;
    */
    const double refZ = -0.089871131212472837868; // for 10 kV / 10 cm
//    const double refZ = 0.087713137568685900503; // for 1MV / 10cm

    const double totalEnergy = Q0*ef->pot(10e-2, 0);

    const int orderOffset = 8;
    int64_t steps = maxSteps >> orderOffset;
    double dt = minDt*(1 << orderOffset);
    for (int c = pseudoOrder-orderOffset; c >=0 ; --c) {
        boris.setElectronInfo(r, 0, -10e-2, 0, 0, 0);
        aphi.setElectronInfo(-10e-2, r, 0, 0, aphi.pTheta(-10e-2, r, 0));
        std::clog << "dt = " << dt << '\n';
        for (int i = 0; i < steps; ++i) {
            aphi.step(dt);
            boris.step(dt);
        }
        const auto pBoris = boris.pos();
        const auto pAPhi = aphi.pos();
        const double potEnergyBoris = Q0*ef->pot(pBoris.z, 0);
        const double kinEnergyBoris = (boris.gamma()-1.0)*(M0*C0*C0);
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
