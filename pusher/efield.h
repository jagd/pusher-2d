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


class ConstEzField: public IStaticEField
{
public:
    ConstEzField(double ez);
    virtual ~ConstEzField() override = default;
    virtual double pot(double, double) const override;
    virtual double ez(double, double) const override;
    virtual double er(double, double) const override;
private:
    const double ez_;
};


class ConstErField: public IStaticEField
{
public:
    ConstErField(double er);
    virtual ~ConstErField() override = default;
    virtual double pot(double, double) const override;
    virtual double ez(double, double) const override;
    virtual double er(double, double) const override;
private:
    const double er_;
};

/// @brief  E(z) = Ez0 + k*z
class LinearEzField: public IStaticEField
{
public:
    LinearEzField(double ez0, double k);
    virtual double pot(double z, double r) const override;
    virtual double ez(double z, double r) const override;
    virtual double er(double z, double r) const override;
private:
    const double ez0_;
    const double k_;
};

/// a mirrored symmetric field, which has inifnite length.
/// z=0 is the symmetric plane. Field strength Ez at z>=0,
/// otherwise -Ez.
class MirroredEzField: public IStaticEField
{
public:
    MirroredEzField(double ez);
    virtual double pot(double z, double r) const override;
    virtual double ez(double z, double r) const override;
    virtual double er(double z, double r) const override;
private:
    const double ez_;
};

#endif // EFIELD_H
