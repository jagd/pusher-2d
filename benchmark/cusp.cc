#include <pusher/pusher.h>
#include <iostream>

int main()
{
    const double b1 = -0.01;
    const double b2 = -b1;
    const double z0field = 1e-4;
	const double r0 = 1e-2;
    const double u0 = 1e8;
    const double z0 = -0.1;

    const auto ef = std::make_shared<ConstEzField>(0);
    const auto mf = std::make_shared<BiUniformMagField>(z0field, b1, b2);
    auto aphi = APhiPusher(ef, mf);
	auto boris = BorisPusher(ef, mf);
	auto rk = RK4Pusher(ef, mf);
	auto leapfrog = LeapFrogPusher(ef, mf);

	// const double flux0 = r0 * r0*mf->bz(0, r0);

    double dtScale = 1;
	const double dt = 1e-12 * dtScale;
	const int64_t maxStep = static_cast<int64_t>((1 << 13)/dtScale);
	std::clog << "dt " << dt << '\n';

	leapfrog.setElectronInfo(r0, 0, z0, 0, 0, u0);
	boris.setElectronInfo(r0, 0, z0, 0, 0, u0);
	rk.setElectronInfo(r0, 0, z0, 0, 0, u0);
	aphi.setElectronInfo(z0, r0, u0, 0, aphi.pTheta(z0, r0, 0), u2gamma(u0));

	for (int64_t i = 0; i < maxStep; ++i) {
		if (i % static_cast<int64_t>(std::ceil(64/dtScale)) == 0) {
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
