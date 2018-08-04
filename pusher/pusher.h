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
    virtual ~IPusher2D() = default;
    virtual void step(double dt) = 0;
protected:
    std::shared_ptr<IStaticEField> efield_;
    std::shared_ptr<IStaticMagField> magfield_;
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
    PV3D u() const;
    PV3D gamma() const;

private:
    PV3D pos_;
    PV3D uLastHalf_; //! u := gamma * v
    PV3D gammaAtPos_;
};

#endif // PUSHER_H
