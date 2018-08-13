#ifndef EFIELD_H
#define EFIELD_H

class IStaticEField
{
public:
    IStaticEField() = default;
    virtual double pot(double z, double r) const;
    virtual double ez(double z, double r) const = 0;
    virtual double er(double z, double r) const = 0;
    virtual ~IStaticEField();
};


class ZeroEField: public IStaticEField
{
public:
    virtual ~ZeroEField() override = default;
    virtual double pot(double, double) const override { return 0; }
    virtual double ez(double, double) const override { return 0; }
    virtual double er(double, double) const override { return 0; }
};


/// a mirrored symmetric field, which has inifnite length.
/// z=0 is the symmetric plane. Field strength Ez at z>=0,
/// otherwise -Ez.
class MirOscZEField: public IStaticEField
{
public:
    MirOscZEField(double ez);
    virtual double pot(double z, double r) const override;
    virtual double ez(double z, double r) const override;
    virtual double er(double z, double r) const override;
private:
    const double ez_;
};

#endif // EFIELD_H
