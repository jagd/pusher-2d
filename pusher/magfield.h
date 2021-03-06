#ifndef MAGFIELD_H
#define MAGFIELD_H

class IStaticMagField
{
public:
    IStaticMagField();
    virtual ~IStaticMagField();
    virtual double aTheta(double z, double r) const = 0;
    virtual double br(double z, double r) const = 0;
    virtual double bz(double z, double r) const = 0;
};


/// A zero field can be created if ctor(0)
class ConstBzField : public IStaticMagField
{
public:
    ConstBzField(double bz): bz_(bz) {}
    virtual double aTheta(double z, double r) const override;
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
    virtual double br(double z, double r) const override;
    virtual double bz(double z, double r) const override;
private:
    const double b0_;
    const double k_;
};

/// A zero field can be created for ctor(0)
class ConstBrField : public IStaticMagField
{
public:
    ConstBrField(double br): br_(br) {}
    virtual double aTheta(double z, double r) const override;
    virtual double br(double z, double r) const override;
    virtual double bz(double z, double r) const override;
private:
    const double br_;
};


/// @brief Uniform B-fields for |z|>z0, inbetween there is a infinity
///        smooth transition
class BiUniformMagField : public IStaticMagField
{
public:
    BiUniformMagField(double z0, double b1, double b2);
    virtual double aTheta(double z, double r) const override;
    virtual double br(double z, double r) const override;
    virtual double bz(double z, double r) const override;

private:
    const double z0_;
    const double z02_;
    const double z03_;
    const double z04_;
    const double z05_;
    const double b1_;
    const double b2_;
    const double coef_;
};

#endif // MAGFIELD_H
