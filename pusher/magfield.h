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
    ~HomogeneousMagField() {}
    virtual double aphi(double z, double r) const;
    virtual double aphi2z(double z, double r) const;
    virtual double aphi2r(double z, double r) const;
    virtual double br(double z, double r) const;
    virtual double bz(double z, double r) const;
private:
    const double bz_;
};

#endif // MAGFIELD_H
