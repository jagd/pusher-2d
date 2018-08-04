#ifndef MAGFIELD_H
#define MAGFIELD_H

class IStaticMagField
{
public:
    IStaticMagField();
    virtual ~IStaticMagField();
    virtual double aphi(double z, double r) const = 0;
    virtual double aphi2z(double z, double r) const = 0;
    virtual double aphi2r(double z, double r) const = 0;
    virtual double br(double z, double r) const = 0;
    virtual double bz(double z, double r) const = 0;
};

class HomogeneousMagField : public IStaticMagField
{
public:
    HomogeneousMagField(double bz): bz_(bz) {}
    ~HomogeneousMagField() = default;
    virtual double aphi(double z, double r) const override;
    virtual double aphi2z(double z, double r) const override;
    virtual double aphi2r(double z, double r) const override;
    virtual double br(double z, double r) const override;
    virtual double bz(double z, double r) const override;
private:
    const double bz_;
};

#endif // MAGFIELD_H
