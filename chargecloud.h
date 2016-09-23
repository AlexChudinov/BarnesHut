#ifndef CHARGECLOUD_H
#define CHARGECLOUD_H

#include "mathutility.h"
#include <vector>
#include "vectortemplate.h"

/*!
 * Distribution of charges in space
 */
template <class Float>
class ChargeCloud
{
    typedef Float charge_type;
    typedef math::Vector<Float, 3> position_type;
    typedef std::vector<charge_type> vector_of_charges_type;
    typedef std::vector<position_type> vector_of_positions_type;

    vector_of_charges_type m_Charges;
    vector_of_positions_type m_Positions;

public:
    ChargeCloud() : m_Charges(), m_Positions() {;}

    ChargeCloud
    (
            const vector_of_charges_type& charges,
            const vector_of_positions_type& positions)
        :
          m_Charges(charges),
          m_Positions(positions)
    {}

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
        position_type result;
        math::mean(m_Positions, m_Charges, result);
        return result;
    }

    /*!
     * \brief maxDispersion finds the direction in which
     * a charge dispersion is maximal
     * \return direction of maximal charge dispersion
     */
    position_type maxDispersion() const
    {
        position_type result;

        //Shift the system into the central coordinate system
        vector_of_positions_type centralized(m_Positions);
        typename vector_of_positions_type::iterator it
                = centralized.begin();
        position_type center = this->chargeCenter();

        Float
                sum_xx = 0.0,
                sum_yy = 0.0,
                sum_zz = 0.0,
                sum_xy = 0.0,
                sum_xz = 0.0,
                sum_zy = 0.0;

        for(; it != centralized.end(); ++it)
        {
            *it -= center;
            sum_xx += (it->X()*it->X());
            sum_yy += (it->Y()*it->Y());
            sum_xy += (it->X()*it->Y());
        }

        //Find first two components of a projection vector
        if (sum_xx != 0.0)
        {
            result.X() = 1.0;
            result.Y() = sum_xy/sum_xx;
            result /= (result*result);
        }
        else result.Y() = 1.0;

        //Find last component of a projection vector
        sum_xx = 0.0; sum_xy = 0.0;
        for(; it != centralized.end(); ++it)
        {
            double x = result * *it;
            sum_xx += (x*x);
            sum_xy = (x*it->Z());
        }
        if(sum_xx != 0.0)
        {
            result.Z() = sum_xy/sum_xx;
            result /= (result*result);
        }
        else
        {
            result.Y() = 0.0;
            result.Z() = 1.0;
        }

        return result;
    }

};

#endif // CHARGECLOUD_H
