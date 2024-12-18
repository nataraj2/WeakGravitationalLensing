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
