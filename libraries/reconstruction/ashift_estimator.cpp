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
#include "ashift_estimator.h"

namespace Alignment {

template<typename T>
std::vector<Point2D<T>> AShiftEstimator<T>::computeShiftFromCorrelations2D(
        T *h_centers, MultidimArray<T> &helper, size_t nDim,
        size_t centerSize, size_t maxShift) {
    assert(centerSize == (2 * maxShift + 1));
    assert(helper.xdim == helper.ydim);
    assert(helper.xdim == centerSize);
    T x;
    T y;
    auto result = std::vector<Point2D<T>>();
    helper.setXmippOrigin(); // tell the array that the 'center' is in the center
    for (size_t n = 0; n < nDim; ++n) {
        helper.data = h_centers + n * centerSize * centerSize;
        bestShift(helper, x, y, nullptr, maxShift);
        result.emplace_back(x, y);
    }
    // avoid data corruption
    helper.data = nullptr;
    return result;
}

template class AShiftEstimator<float>;

} /* namespace Alignment */