#ifndef MAGFIELD_H
#define MAGFIELD_H

class IStaticMagField
{
public:
    IStaticMagField();
    virtual ~IStaticMagField();
    virtual double aTheta(double z, double r) const = 0;
    virtual double aTheta2z(double z, double r) const = 0;
    virtual double aTheta2r(double z, double r) const = 0;
    virtual double br(double z, double r) const = 0;
    virtual double bz(double z, double r) const = 0;
};


/// A zero field can be created if ctor(0)
class ConstBzField : public IStaticMagField
{
public:
    ConstBzField(double bz): bz_(bz) {}
    virtual double aTheta(double z, double r) const override;
    virtual double aTheta2z(double z, double r) const override;
    virtual double aTheta2r(double z, double r) const override;
    virtual double br(double z, double r) const override;
    virtual double bz(double z, double r) const override;
private:
    const double bz_;
};

#endif // MAGFIELD_H
