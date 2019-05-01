/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef LIBRARIES_RECONSTRUCTION_ASHIFT_ESTIMATOR_H_
#define LIBRARIES_RECONSTRUCTION_ASHIFT_ESTIMATOR_H_

#include "data/point3D.h"
#include <vector>
#include <cassert>
#include <limits>
#include <complex>
#include "data/hw.h"
#include "data/dimensions.h"
#include "core/xmipp_error.h"

namespace Alignment {

enum class AlignType { None, OneToN, NToM, Consecutive };

template<typename T>
class AShiftEstimator {
public:
    AShiftEstimator() {
        setDefault();
    }
    virtual ~AShiftEstimator() {
        release();
    }

    virtual void init2D(const HW &hw, AlignType type,
               const Dimensions &dims, size_t batch, Point3D<size_t> maxShift) = 0;

// FIXME DS add function to load reference (in spatial domain)

//    virtual void computeShift2DOneToN(
//        T *h_others) = 0; // FIXME DS uncomment

    inline std::vector<Point3D<T>> getShifts() {
        // FIXME DS add check that it's computed
        return m_shifts;
    }

    virtual void release();

protected:
    // various
    AlignType m_type;
    const Dimensions *m_dims;
    size_t m_batch;
    Point3D<size_t> m_maxShift;

    // computed shifts
    std::vector<Point3D<T>> m_shifts;

    // flags
    bool m_is_ref_loaded;
    bool m_is_shift_computed;
    bool m_isInit;

    virtual void setDefault();
    virtual void init2D(AlignType type, const Dimensions &dims,
               size_t batch, Point3D<size_t> maxShift);
    virtual void check();
};

} /* namespace Alignment */

#endif /* LIBRARIES_RECONSTRUCTION_ASHIFT_ESTIMATOR_H_ */
