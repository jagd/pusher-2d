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


/// @brief a linear Bz, fulfilling Div(B)=0
///        Bz = b0 + k*z
///		   Br = -r*k/2;
class LinearBzField : public IStaticMagField
{
public:
    LinearBzField(double b0, double k): b0_(b0), k_(k) {}
    virtual double aTheta(double z, double r) const override;
    virtual double aTheta2z(double z, double r) const override;
    virtual double aTheta2r(double z, double r) const override;
    virtual double br(double z, double r) const override;
    virtual double bz(double z, double r) const override;
private:
    const double b0_;
    const double k_;
};

#endif // MAGFIELD_H
