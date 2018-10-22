#include <pusher/pusher.h>
#include <iostream>

int main()
{
	const double Ez = -10e3 / 10e-2;
    const double k = 1.0;
	const double r0 = 11e-3;
	const double z0 = -8e-3;

    const auto ef = std::make_shared<ConstEzField>(Ez);
    const auto mf = std::make_shared<LinearBzField>(0, k);
    auto aphi = APhiPusher(ef, mf);
	auto boris = BorisPusher(ef, mf);
	auto leapfrog = LeapFrog(ef, mf);

	const double flux0 = r0 * r0*mf->bz(z0, r0);

	double dtScale = 1;
	const double dt = 1e-15 * dtScale;
	const int64_t maxStep = (1 << 22)/dtScale;
	std::clog << "dt " << dt << '\n';
	leapfrog.setElectronInfo(r0, 0, z0, 0, 0, 0);
	boris.setElectronInfo(r0, 0, z0, 0, 0, 0);
	aphi.setElectronInfo(z0, r0, 0, 0, aphi.pTheta(z0, r0, 0), 1.0);
	for (int64_t i = 0; i < maxStep; ++i) {
		if (i % 4096 == 0) {
			const auto pa = aphi.pos();
			const auto pb = boris.pos();
			const auto pl = leapfrog.pos();
			std::clog << i << " / " << maxStep << '\n';
			std::cout << i * dt << ' '
				<< pb.z << ' '
				<< std::sqrt(pb.x*pb.x + pb.y*pb.y) << ' '
				<< pa.z << ' '
				<< pa.r << ' '
				<< pl.z << ' '
				<< std::sqrt(pl.x*pl.x + pl.y*pl.y) << ' '
				// field line:
				<< pa.z << ' '
				<< std::sqrt(std::abs(flux0 / mf->bz(pa.z, 0))) << '\n';
		}
		aphi.step(dt);
		boris.step(dt);
		leapfrog.step(dt);
	}

	return 0;
}
