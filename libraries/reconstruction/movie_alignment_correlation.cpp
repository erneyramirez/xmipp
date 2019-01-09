/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
 *             David Strelak (davidstrelak@gmail.com)
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

#include "reconstruction/movie_alignment_correlation.h"

template<typename T>
void ProgMovieAlignmentCorrelation<T>::defineParams() {
    AProgMovieAlignmentCorrelation<T>::defineParams();
    this->addExampleLine(
                "xmipp_movie_alignment_correlation -i movie.xmd --oaligned alignedMovie.stk --oavg alignedMicrograph.mrc");
    this->addSeeAlsoLine("xmipp_cuda_movie_alignment_correlation");
}

template<typename T>
AlignmentResult<T> ProgMovieAlignmentCorrelation<T>::computeGlobalAlignment(
        const MetaData &movie, const Image<T> &dark, const Image<T> &gain) {
    loadData(movie, dark, gain);
    size_t N = this->nlast - this->nfirst + 1; // no of images to process
    Matrix2D<T> A(N * (N - 1) / 2, N - 1);
    Matrix1D<T> bX(N * (N - 1) / 2), bY(N * (N - 1) / 2);

    if (this->verbose)
        std::cout << "Computing shifts between frames ..." << std::endl;
    // Now compute all shifts
    computeShifts(N, bX, bY, A);

    // Choose reference image as the minimax of shifts
    auto ref = core::optional<size_t>();
    return this->computeAlignment(bX, bY, A, ref, N);
}

template<typename T>
void ProgMovieAlignmentCorrelation<T>::applyShiftsComputeAverage(
            const MetaData& movie, const Image<T>& dark, const Image<T>& gain,
            Image<T>& initialMic, size_t& Ninitial, Image<T>& averageMicrograph,
            size_t& N, const LocalAlignmentResult<T> &alignment) {
    throw std::logic_error("Not implemented");
}

template<typename T>
LocalAlignmentResult<T> ProgMovieAlignmentCorrelation<T>::computeLocalAlignment(
        const MetaData &movie, const Image<T> &dark, const Image<T> &gain,
        const AlignmentResult<T> &globAlignment) {
    throw std::logic_error("Not implemented");
}

template<typename T>
void ProgMovieAlignmentCorrelation<T>::loadData(const MetaData& movie,
        const Image<T>& dark, const Image<T>& gain) {
    sizeFactor = this->computeSizeFactor();
    MultidimArray<T> filter;
    FourierTransformer transformer;
    bool firstImage = true;
    int n = 0;
    FileName fnFrame;
    Image<T> croppedFrame, reducedFrame;

    if (this->verbose) {
        std::cout << "Computing Fourier transform of frames ..." << std::endl;
        init_progress_bar(movie.size());
    }

    FOR_ALL_OBJECTS_IN_METADATA(movie)
    {
        if (n >= this->nfirst && n <= this->nlast) {
            this->loadFrame(movie, dark, gain, __iter.objId, croppedFrame);

            if (firstImage) {
                newXdim = croppedFrame().xdim * sizeFactor;
                newYdim = croppedFrame().ydim * sizeFactor;
            }

            // Reduce the size of the input frame
            scaleToSizeFourier(1, newYdim, newXdim, croppedFrame(),
                    reducedFrame());

            // Now do the Fourier transform and filter
            MultidimArray<std::complex<T> > *reducedFrameFourier =
                    new MultidimArray<std::complex<T> >;
            transformer.FourierTransform(reducedFrame(), *reducedFrameFourier,
                    true);
            if (firstImage) {
                firstImage = false;
                filter = this->createLPF(this->getTargetOccupancy(), newXdim,
                    newYdim);
            }
            for (size_t nn = 0; nn < filter.nzyxdim; ++nn) {
                T wlpf = DIRECT_MULTIDIM_ELEM(filter, nn);
                DIRECT_MULTIDIM_ELEM(*reducedFrameFourier,nn) *= wlpf;
            }
            frameFourier.push_back(reducedFrameFourier);
        }
        ++n;
        if (this->verbose)
            progress_bar(n);
    }
    if (this->verbose)
        progress_bar(movie.size());
}

template<typename T>
void ProgMovieAlignmentCorrelation<T>::computeShifts(size_t N,
        const Matrix1D<T>& bX, const Matrix1D<T>& bY, const Matrix2D<T>& A) {
    int idx = 0;
    assert(frameFourier.size() > 0);
    MultidimArray<T> Mcorr;
    Mcorr.resizeNoCopy(newYdim, newXdim);
    Mcorr.setXmippOrigin();
    CorrelationAux aux;
    for (size_t i = 0; i < N - 1; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            bestShift(*frameFourier[i], *frameFourier[j], Mcorr, bX(idx),
                    bY(idx), aux, NULL, this->maxShift * sizeFactor);
            bX(idx) /= sizeFactor; // scale to expected size
            bY(idx) /= sizeFactor;
            if (this->verbose)
                std::cerr << "Frame " << i + this->nfirst << " to Frame "
                        << j + this->nfirst << " -> ("
                        << bX(idx) << ","
                        << bY(idx) << ")\n";
            for (int ij = i; ij < j; ij++)
                A(idx, ij) = 1;

            idx++;
        }
        delete frameFourier[i];
    }
}

template<typename T>
void ProgMovieAlignmentCorrelation<T>::applyShiftsComputeAverage(
        const MetaData& movie, const Image<T>& dark, const Image<T>& gain,
        Image<T>& initialMic, size_t& Ninitial, Image<T>& averageMicrograph,
        size_t& N, const AlignmentResult<T> &globAlignment) {
    // Apply shifts and compute average
    Image<T> frame, croppedFrame, reducedFrame, shiftedFrame;
    Matrix1D<T> shift(2);
    FileName fnFrame;
    int j = 0;
    int n = 0;
    Ninitial = N = 0;
    FOR_ALL_OBJECTS_IN_METADATA(movie)
    {
        if (n >= this->nfirstSum && n <= this->nlastSum) {
            // get shifts
            XX(shift) = -globAlignment.shifts.at(n).x; // we want to compensate
            YY(shift) = -globAlignment.shifts.at(n).y; // the shift

            // load frame
            this->loadFrame(movie, dark, gain, __iter.objId, croppedFrame);
            if (this->bin > 0) {
                scaleToSizeFourier(1, floor(YSIZE(croppedFrame()) / this->bin),
                        floor(XSIZE(croppedFrame()) / this->bin),
                        croppedFrame(), reducedFrame());
                croppedFrame() = reducedFrame();
            }

            if ( ! this->fnInitialAvg.isEmpty()) {
                if (j == 0)
                    initialMic() = croppedFrame();
                else
                    initialMic() += croppedFrame();
                Ninitial++;
            }

            if (this->fnAligned != "" || this->fnAvg != "") {
                if (this->outsideMode == OUTSIDE_WRAP)
                    translate(this->BsplineOrder, shiftedFrame(),
                            croppedFrame(), shift, WRAP);
                else if (this->outsideMode == OUTSIDE_VALUE)
                    translate(this->BsplineOrder, shiftedFrame(),
                            croppedFrame(), shift, DONT_WRAP,
                            this->outsideValue);
                else
                    translate(this->BsplineOrder, shiftedFrame(),
                            croppedFrame(), shift, DONT_WRAP,
                            croppedFrame().computeAvg());
                if (this->fnAligned != "")
                    shiftedFrame.write(this->fnAligned, j + 1, true,
                            WRITE_REPLACE);
                if (this->fnAvg != "") {
                    if (j == 0)
                        averageMicrograph() = shiftedFrame();
                    else
                        averageMicrograph() += shiftedFrame();
                    N++;
                }
            }
            std::cout << "Frame " << std::to_string(j) << " processed." << std::endl;
            j++;
        }
        n++;
    }
}

template class ProgMovieAlignmentCorrelation<double> ;
