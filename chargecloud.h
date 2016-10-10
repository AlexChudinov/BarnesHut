#ifndef CHARGECLOUD_H
#define CHARGECLOUD_H

#include "binarytree.h"

#include "mathutility.h"
#include "matrixtemplate.h"

//Forward declaration
template <class Float>
class ChargeCloud;

/**
 * Charge distribution approximation using central charge
 * and quadrupole moments
 */
template <class Float> struct charge_approximation
{
    typedef Float charge_type;
    typedef math::Vector<Float, 3> position_type;
    typedef math::Matrix<Float, 3> matrix3D_type;

    charge_type total_charge;
    position_type central_point;
    position_type max_dispersion;

    charge_approximation() {}

    charge_approximation(const ChargeCloud<Float>& cloud)
        :
          total_charge(cloud.total_charge()),
          central_point(cloud.chargeCenter()),
          max_dispersion(cloud.maxDispersion())
    {}

};

/*!
 * Distribution of charges in space
 */
template <class Float>
class ChargeCloud
{
    typedef Float charge_type;
    typedef math::vector_c<Float, 3> position_type;
    typedef math::matrix_c<Float, 3> matrix3D_type;
    typedef std::vector<charge_type> vector_of_charges_type;
    typedef std::vector<position_type> vector_of_positions_type;
    typedef data_structs::triple<charge_approximation<Float>, ChargeCloud> triple;
    typedef std::pair<ChargeCloud,ChargeCloud> pair;

    vector_of_charges_type m_Charges;
    vector_of_positions_type m_Positions;

public:
    typedef data_structs::binary_tree<charge_approximation<Float>>
        barnes_hut_tree;

    ChargeCloud() : m_Charges(), m_Positions() {;}

    ChargeCloud
    (
            const vector_of_charges_type& charges,
            const vector_of_positions_type& positions)
        :
          m_Charges(charges),
          m_Positions(positions)
    { assert(m_Charges.size() == m_Positions.size()); }

    ChargeCloud& addParticle
    (
            const charge_type& charge,
            const position_type& position)
    {
        m_Charges.push_back(charge);
        m_Positions.push_back(position);
        return *this;
    }

    /*!
     * Computes charge center of a cloud
     */
    position_type chargeCenter() const
    {
        position_type result(0.0);
        math::mean(m_Positions, m_Charges, result);
        return result;
    }

    /*!
     * \brief maxDispersion looks for the direction in which
     * a charge dispersion is maximal
     * \param center A center of dispersion
     * \return direction of maximal charge dispersion
     */
    position_type maxDispersion(const position_type& center) const
    {
        //Shift the system into the central coordinate system
        vector_of_positions_type centralized(m_Positions);

        //Prepare unbiased system
        for (position_type& p: centralized) p -= center;

        return math::pc1(centralized,m_Charges);
    }

    /*!
     * \brief total_charge calculates total charge of a cloud
     * \return total charge of a cloud
     */
    charge_type total_charge() const
    {
        charge_type result(0.0);
        for(const charge_type& q : m_Charges) result += q;
        return result;
    }

    barnes_hut_tree create_barnes_hut_tree() const
    {
        struct create_tree_node_data
        {
            inline triple operator()
                    (const ChargeCloud& cloud) const
            {
                triple result;
                result.first = charge_approximation<Float>(cloud);
                result.second = cloud.split();
                return result;
            }
        };

        return barnes_hut_tree(create_tree_node_data(), *this);
    }

    /*!
     * \brief empty Checks if cloud is empty
     * \return true if there are no any charges
     */
    bool empty() const
    {
        return m_Charges.empty();
    }

    std::size_t size() const { return m_Charges.size(); }

    void reserve(std::size_t n)
    {
        m_Charges.reserve(n);
        m_Positions.reserve(n);
    }

    /*!
     * Splits charged cloud using plane that goes trough the
     * center of charge and orthogonal to a direction
     * of a max dispersion
     */
    pair split() const
    {
        ChargeCloud cloud1, cloud2;

        std::vector<bool> charge_alignment(this->size());
        std::size_t n1 = 0, n2 = 0;

        position_type charge_center = this->chargeCenter();
        position_type max_dispersion
                = this->maxDispersion(charge_center);

        if (math::abs(max_dispersion) == 0.0)
            return pair(cloud1, cloud2);

        for (std::size_t i = 0; i < this->size(); ++i)
        {
            if( (m_Positions[i] - charge_center)
                    * max_dispersion < 0.0)
            {
                charge_alignment[i] = true;
                ++n1;
            }
            else
            {
                charge_alignment[i] = false;
                ++n2;
            }
        }

        cloud1.reserve(n1); cloud2.reserve(n2);

        for (std::size_t i = 0; i < this->size(); ++i)
        {
            if (charge_alignment[i])
                cloud1.addParticle(m_Charges[i], m_Positions[i]);
            else
                cloud2.addParticle(m_Charges[i], m_Positions[i]);
        }

        return pair(cloud1, cloud2);
    }
};

#endif // CHARGECLOUD_H
