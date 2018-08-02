#ifndef MAGFIELD_H
#define MAGFIELD_H

class IMagField
{
public:
    IMagField();
    virtual double aphi(double z, double r) const;
    virtual double bz(double z, double r) const;
    virtual double br(double z, double r) const;
    virtual ~IMagField();
};

#endif // MAGFIELD_H
