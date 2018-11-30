#ifndef PUSHER_H
#define PUSHER_H

#include "magfield.h"
#include "efield.h"
#include "pv.h"
#include <memory>

#define C0 2.99792458e8
#define Q0 -1.60217662e-19
#define M0 9.10938356e-31

class IPusher2D
{
public:
    IPusher2D(
        std::shared_ptr<IStaticEField>,
        std::shared_ptr<IStaticMagField>
    );
    virtual ~IPusher2D();
    virtual void step(double dt) = 0;
protected:
    std::shared_ptr<IStaticEField> efield_;
    std::shared_ptr<IStaticMagField> magfield_;
};

class RK4Pusher: public IPusher2D
{
public:
    RK4Pusher(
        std::shared_ptr<IStaticEField> e,
        std::shared_ptr<IStaticMagField> m
    ): IPusher2D(e, m) {}
    void setElectronInfo(
        double x, double y, double z,
        double ux, double uy, double uz
    );
    virtual void step(double dt) override;
    PV3D pos() const;
    PV3D uCurrent() const;
    double gammaCurrent() const; //! current gamma

private:
	PV3D b3D(const PV3D &pos) const;
	PV3D e3D(const PV3D &pos) const;
	PV3D a3D(const PV3D &pos, const PV3D &v) const;
    PV3D pos_;
    PV3D u_; //! u := gamma * v
};

class LeapFrogPusher: public IPusher2D
{
public:
    LeapFrogPusher(
        std::shared_ptr<IStaticEField> e,
        std::shared_ptr<IStaticMagField> m
    ): IPusher2D(e, m) {}
    void setElectronInfo(
        double x, double y, double z,
        double ux, double uy, double uz
    );
    virtual void step(double dt) override;
    PV3D pos() const;
    PV3D uLastHalf() const;
    double gammaLastHalf() const;

private:
    PV3D pos_;
    PV3D uLastHalf_; //! u := gamma * v
};

class BorisPusher: public IPusher2D
{
public:
    BorisPusher(
        std::shared_ptr<IStaticEField> e,
        std::shared_ptr<IStaticMagField> m
    ): IPusher2D(e, m) {}
    void setElectronInfo(
        double x, double y, double z,
        double ux, double uy, double uz
    );
    virtual void step(double dt) override;
    PV3D pos() const;
    PV3D uLastHalf() const;
    double gammaCurrent() const;

private:
    PV3D pos_;
    PV3D uLastHalf_; //! u := gamma * v
    double gammaAtPos_;
};

class LeapFrogPusher2D: public IPusher2D
{
public:
    LeapFrogPusher2D(
        std::shared_ptr<IStaticEField> e,
        std::shared_ptr<IStaticMagField> m
    ): IPusher2D(e, m) {}
    double pTheta(double z, double r, double uTheta) const;
    double pTheta() const {return pTheta_;}
    void step(double dt) override;
    void setElectronInfo(
        double z, double r,
        double uzLastHalf, double urLastHalf,
        double pTheta,
        double gamma
    );
    PV2D pos() const;
    PV2D vLastHalf() const;
    double gammaCurrent() const;
private:
    PV2D pos_;
    PV2D uLastHalf_;
    double pTheta_;
    double totalEnergy_; //! in Volt*Q, where Q is signed
};

class LeapFrogPusher2DSync: public IPusher2D
{
public:
    LeapFrogPusher2DSync(
        std::shared_ptr<IStaticEField> e,
        std::shared_ptr<IStaticMagField> m
    ): IPusher2D(e, m) {}
    double pTheta() const {return pTheta_;}
    void step(double dt) override;
    void setElectronInfo(
        double z, double r,
        double uz, double ur,
        double uTheta
    );
    PV2D pos() const { return pos_; }
    double gammaCurrent() const { return gammaAt(pos_); }
private:
    double pTheta(double z, double r, double uTheta) const;
    double gammaAt(const PV2D &pos) const;
    PV2D pos_;
    PV2D u_;
    double pTheta_;
    double totalEnergy_; //! in Volt*Q, where Q is signed
};

double v2gamma(double v);
double u2gamma(double u);
double usqr2gamma(double usqr);
double v2u(double v);
double gamma2v(double gamma);

#endif // PUSHER_H
