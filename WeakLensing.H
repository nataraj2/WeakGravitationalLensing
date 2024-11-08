#include <iostream>     // For standard input/output
#include <cmath>    // For mathematical functions like sqrt, pow
#include <iomanip>      // For output formatting
#include <limits>       // For handling numeric limits
#include <vector>


inline double convert_ly_to_m(double& dist_ly) {
    return dist_ly*3e8*3600.0*24.0*365.0;
}

inline double convert_Gpc_to_ly(double& dist_pc) {
    return dist_pc*3.26e9;
}

inline double convert_Gpc_to_m(double& dist_pc) {
    double tmp = convert_Gpc_to_ly(dist_pc);
    tmp = convert_ly_to_m(tmp);
    return tmp;
}

class Triplet
{
    public:
        Triplet(double val0, double val1, double val2) :
                m_val0(val0), m_val1(val1), m_val2(val2) {}

    // Overload operator[] to provide array-like access to x, y, z
    double& operator[](size_t index) {
        if (index == 0) return m_val0;
        if (index == 1) return m_val1;
        if (index == 2) return m_val2;
        throw std::out_of_range("Index out of range for Triplet");
    }

    const double& operator[](size_t index) const {
        if (index == 0) return m_val0;
        if (index == 1) return m_val1;
        if (index == 2) return m_val2;
        throw std::out_of_range("Index out of range for Triplet");
    }

    private:
        double m_val0, m_val1, m_val2;
};

class PolyRay
{
    public:
        std::vector<Triplet> coord;
        std::vector<Triplet> dir;
        double theta_start, phi_start;
};


void write_polyray_bundle(const std::vector<PolyRay>& polyrayvec,
                          const std::string filename);
