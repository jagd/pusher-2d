#include <pusher/pusher.h>
#include <iostream>

int main()
{
	const double Ez = 0 / 10e-2;
    const double k = 100;
	const double r0 = 1e-3;
    const double u0 = 1e8;

    const auto ef = std::make_shared<ConstEzField>(Ez);
    const auto mf = std::make_shared<LinearBzField>(10e-3, k);
    auto aphi = APhiPusher(ef, mf);
	auto boris = BorisPusher(ef, mf);
	auto rk = RK4Pusher(ef, mf);
	auto leapfrog = LeapFrog(ef, mf);

	// const double flux0 = r0 * r0*mf->bz(0, r0);

	double dtScale = 16;
	const double dt = 1e-15 * dtScale;
	const int64_t maxStep = (1 << 20)/dtScale;
	std::clog << "dt " << dt << '\n';

	leapfrog.setElectronInfo(r0, 0, 0, 0, 0, u0);
	boris.setElectronInfo(r0, 0, 0, 0, 0, u0);
	rk.setElectronInfo(r0, 0, 0, 0, 0, u0);
	aphi.setElectronInfo(0, r0, u0/u2gamma(u0), 0, aphi.pTheta(0, r0, 0), u2gamma(u0));

	for (int64_t i = 0; i < maxStep; ++i) {
		if (i % static_cast<int64_t>(std::ceil(256/dtScale)) == 0) {
			const auto pl = leapfrog.pos();
			const auto pb = boris.pos();
			const auto pr = rk.pos();
			const auto pa = aphi.pos();
			std::clog << i*100.0/maxStep << " %\n" << '\n';
			std::cout << i * dt << ' '

				<< pl.z << ' '
				<< std::sqrt(pl.x*pl.x + pl.y*pl.y) << ' '

				<< pb.z << ' '
				<< std::sqrt(pb.x*pb.x + pb.y*pb.y) << ' '

				<< pa.z << ' '
				<< pa.r << ' '

				<< pr.z << ' '
				<< std::sqrt(pr.x*pr.x + pr.y*pr.y) << ' '

				// cross section
				<< pb.x << ' '
				<< pb.y << '\n';
		}
		aphi.step(dt);
		boris.step(dt);
		rk.step(dt);
		leapfrog.step(dt);
	}

	return 0;
}
