#ifndef EFIELD_H
#define EFIELD_H

class IStaticEField
{
public:
    IStaticEField();
    virtual double pot(double, double) const;
    virtual double ez(double, double) const = 0;
    virtual double er(double, double) const = 0;
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

#endif // EFIELD_H
