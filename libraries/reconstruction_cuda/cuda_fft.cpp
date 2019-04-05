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

#include "cuda_fft.h"
#include <cufft.h>
#include "cuda_utils.h"
#include "cuFFTAdvisor/advisor.h"

template<typename T>
void CudaFFT<T>::init(const FFTSettingsNew<T> &settings) {
    release();

    m_settings = settings;
    m_plan = createPlan(m_settings);

    // allocate input data storage
    gpuErrchk(cudaMalloc(&m_d_SD, m_settings.sBytesBatch()));
    if (m_settings.isInPlace()) {
        // input data holds also the output
        m_d_FD = (std::complex<T>*)m_d_SD;
    } else {
        // allocate also the output buffer
        gpuErrchk(cudaMalloc(&m_d_FD, m_settings.fBytesBatch()));
    }

    m_isInit = true;
}

template<typename T>
void CudaFFT<T>::setDefault() {
    m_settings = FFTSettingsNew<T>(0);
    m_isInit = false;
    m_d_SD = nullptr;
    m_d_FD = nullptr;
}

template<typename T>
void CudaFFT<T>::release() {
    gpuErrchk(cudaFree(m_d_SD));
    if ((void*)m_d_FD != (void*)m_d_SD) {
        gpuErrchk(cudaFree(m_d_FD));
    }
    if (m_isInit) { // avoid destroying empty plan
        gpuErrchkFFT(cufftDestroy(m_plan));
    }
    setDefault();
}

template<typename T>
std::complex<T>* CudaFFT<T>::fft(T *d_inOut) {
    return fft(d_inOut, (std::complex<T>*) d_inOut);
}

template<typename T>
std::complex<T>* CudaFFT<T>::fft(cufftHandle plan, T *d_inOut) {
    return fft(plan, d_inOut, (std::complex<T>*)d_inOut);
}

template<typename T>
T* CudaFFT<T>::ifft(std::complex<T> *d_inOut) {
    return ifft(d_inOut, (T*)d_inOut);
}

template<typename T>
T* CudaFFT<T>::ifft(cufftHandle plan, std::complex<T> *d_inOut) {
    return ifft(plan, d_inOut, (T*)d_inOut);
}

template<typename T>
std::complex<T>* CudaFFT<T>::fft(const T *h_in,
        std::complex<T> *h_out) {
    auto isReady = m_isInit && m_settings.isForward();
    if ( ! isReady) {
        REPORT_ERROR(ERR_LOGIC_ERROR, "Not ready to perform Fourier Transform. "
                "Call init function first");
    }

    // process signals in batches
    for (size_t offset = 0; offset < m_settings.sDim().n(); offset += m_settings.batch()) {
        // how many signals to process
        size_t toProcess = std::min(m_settings.batch(), m_settings.sDim().n() - offset);

        // copy memory
        gpuErrchk(cudaMemcpy(
                m_d_SD,
                h_in + offset * m_settings.sDim().xyzPadded(),
                toProcess * m_settings.sBytesSingle(),
                    cudaMemcpyHostToDevice));

        fft(m_plan, m_d_SD, m_d_FD);

        // copy data back
        gpuErrchk(cudaMemcpy(
                h_out + offset * m_settings.fDim().xyzPadded(),
                m_d_FD,
                toProcess * m_settings.fBytesSingle(),
                cudaMemcpyDeviceToHost));
    }
    return h_out;
}

template<typename T>
T* CudaFFT<T>::ifft(const std::complex<T> *h_in,
        T *h_out) {
    auto isReady = m_isInit && ( ! m_settings.isForward());
    if ( ! isReady) {
        REPORT_ERROR(ERR_LOGIC_ERROR, "Not ready to perform Inverse Fourier Transform. "
                "Call init function first");
    }

    // process signals in batches
    for (size_t offset = 0; offset < m_settings.fDim().n(); offset += m_settings.batch()) {
        // how many signals to process
        size_t toProcess = std::min(m_settings.batch(), m_settings.fDim().n() - offset);

        // copy memoryvim
        gpuErrchk(cudaMemcpy(
                m_d_FD,
                h_in + offset * m_settings.fDim().xyzPadded(),
                toProcess * m_settings.fBytesSingle(),
                    cudaMemcpyHostToDevice));

        ifft(m_plan, m_d_FD, m_d_SD);

        // copy data back
        gpuErrchk(cudaMemcpy(
                h_out + offset * m_settings.sDim().xyzPadded(),
                m_d_SD,
                toProcess * m_settings.sBytesSingle(),
                cudaMemcpyDeviceToHost));
    }
    return h_out;
}

template<typename T>
size_t CudaFFT<T>::estimatePlanBytes(const FFTSettingsNew<T> &settings) {
    size_t size = 0;
    auto f = [&] (int rank, int *n, int *inembed,
            int istride, int idist, int *onembed, int ostride,
            int odist, cufftType type, int batch) {
        gpuErrchkFFT(cufftEstimateMany(rank, n, inembed,
                istride, idist, onembed, ostride,
                odist, type, batch, &size));
    };
    manyHelper(settings, f);
    return size;
}

template<typename T>
size_t CudaFFT<T>::estimateTotalBytes(const FFTSettingsNew<T> &settings) {
    size_t planBytes = estimatePlanBytes(settings);
    size_t dataBytes = settings.sBytesBatch()
            + (settings.isInPlace() ? 0 : settings.fBytesBatch());
    return planBytes + dataBytes;
}

template<typename T>
std::complex<T>* CudaFFT<T>::fft(cufftHandle plan, const T *d_in,
        std::complex<T> *d_out) {
    if (std::is_same<T, float>::value) {
        gpuErrchkFFT(cufftExecR2C(plan, (cufftReal*)d_in, (cufftComplex*)d_out));
    } else if (std::is_same<T, double>::value){
        gpuErrchkFFT(cufftExecD2Z(plan, (cufftDoubleReal*)d_in, (cufftDoubleComplex*)d_out));
    } else {
        REPORT_ERROR(ERR_TYPE_INCORRECT, "Not implemented");
    }
    return d_out;
}

template<typename T>
T* CudaFFT<T>::ifft(cufftHandle plan, const std::complex<T> *d_in,
        T *d_out) {
    if (std::is_same<T, float>::value) {
        gpuErrchkFFT(cufftExecC2R(plan, (cufftComplex*)d_in, (cufftReal*)d_out));
    } else if (std::is_same<T, double>::value){
        gpuErrchkFFT(cufftExecZ2D(plan, (cufftDoubleComplex*)d_in, (cufftDoubleReal*)d_out));
    } else {
        REPORT_ERROR(ERR_TYPE_INCORRECT, "Not implemented");
    }
    return d_out;
}

template<typename T>
template<typename F>
void CudaFFT<T>::manyHelper(const FFTSettingsNew<T> &settings, F function) {
    auto n = std::array<int, 3>{(int)settings.sDim().z(), (int)settings.sDim().y(), (int)settings.sDim().x()};
    int idist;
    int odist;
    cufftType type;
    if (settings.isForward()) {
        idist = settings.sDim().xyzPadded();
        odist = settings.fDim().xyzPadded();
        type = std::is_same<T, float>::value ? CUFFT_R2C : CUFFT_D2Z;
    } else {
        idist = settings.fDim().xyzPadded();
        odist = settings.sDim().xyzPadded();
        type = std::is_same<T, float>::value ? CUFFT_C2R : CUFFT_Z2D;
    }
    int rank = 3;
    if (settings.sDim().z() == 1) rank--;
    if ((2 == rank) && (settings.sDim().y() == 1)) rank--;

    int offset = 3 - rank;
    function(rank, &n[offset], nullptr,
            1, idist, nullptr, 1, odist, type, settings.batch());
}

template<typename T>
cufftHandle CudaFFT<T>::createPlan(const FFTSettingsNew<T> &settings) {
    cufftHandle plan;
    auto f = [&] (int rank, int *n, int *inembed,
            int istride, int idist, int *onembed, int ostride,
            int odist, cufftType type, int batch) {
        gpuErrchkFFT(cufftPlanMany(&plan, rank, n, inembed,
                istride, idist, onembed, ostride,
                odist, type, batch));
    };
    manyHelper(settings, f);
    return plan;
}

template<typename T>
FFTSettingsNew<T> CudaFFT<T>::findMaxBatch(const GPU &gpu,
        const FFTSettingsNew<T> &settings,
        size_t reserveBytes) {
    size_t freeBytes = gpu.lastFreeBytes();
    if (reserveBytes > freeBytes) {
        REPORT_ERROR(ERR_GPU_MEMORY, "You want to reserve more memory than is currently available");
    }
    size_t availableBytes = freeBytes - reserveBytes;
    size_t singleBytes = settings.sBytesSingle() + (settings.isInPlace() ? 0 : settings.fBytesSingle());
    size_t batch = (availableBytes / singleBytes) + 1; // + 1 will be deducted in the while loop
    while (batch > 1) {
        batch--;
        auto tmp = FFTSettingsNew<T>(settings.sDim(), batch, settings.isInPlace(), settings.isForward());
        size_t totalBytes = estimateTotalBytes(tmp);
        if (totalBytes <= availableBytes) {
            return tmp;
        }
    }
    REPORT_ERROR(ERR_GPU_MEMORY, "Estimated batch size is 0(zero). "
            "This probably means you don't have enough GPU memory for even a single transformation.");
}

template<typename T>
core::optional<FFTSettingsNew<T>> CudaFFT<T>::findOptimal(GPU &gpu,
        const FFTSettingsNew<T> &settings,
        size_t reserveBytes, bool squareOnly, int sigPercChange,
        bool crop, bool verbose) {
    using cuFFTAdvisor::Tristate::TRUE;
    using cuFFTAdvisor::Tristate::FALSE;
    using core::optional;
    size_t freeBytes = gpu.lastFreeBytes();
    std::vector<cuFFTAdvisor::BenchmarkResult const *> *options =
            cuFFTAdvisor::Advisor::find(30, gpu.device(),
                    settings.sDim().x(), settings.sDim().y(), settings.sDim().z(), settings.sDim().n(),
                    TRUE, // use batch
                    std::is_same<T, float>::value ? TRUE : FALSE,
                    settings.isForward() ? TRUE : FALSE,
                    settings.isInPlace() ? TRUE : FALSE,
                    cuFFTAdvisor::Tristate::TRUE, // is real
                    sigPercChange, memoryUtils::MB(freeBytes - reserveBytes),
                    false, // allow transposition
                    squareOnly, crop);

    auto result = optional<FFTSettingsNew<T>>();
    if (0 != options->size()) {
        auto res = options->at(0);
        auto optSetting = FFTSettingsNew<T>(
                res->transform->X,
                res->transform->Y,
                res->transform->Z,
                settings.sDim().n(),
                res->transform->N / res->transform->repetitions,
                settings.isInPlace(),
                settings.isForward());
        result = optional<FFTSettingsNew<T>>(optSetting);
    }
    for (auto& it : *options) delete it;
    delete options;
    return result;
}

// explicit instantiation
template class CudaFFT<float>;
template class CudaFFT<double>;
