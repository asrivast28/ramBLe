/***
 *  File: bitvector.hpp
 *  Created: Oct 22, 2016
 *
 *    Authors: Matthew Eichhorn <maeichho@buffalo.edu>
 *             Blake Hurlburt <blakehur@buffalo.edu>
 *             Grant Iraci <grantira@buffalo.edu>
 *             Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2016-2017 SCoRe Group http://www.score-group.org/
 *  Distributed under the MIT License.
 */

#ifndef BITVECTOR_HPP
#define BITVECTOR_HPP

#include <algorithm>
#include <mm_malloc.h>
#if defined(HAVE_AVX2_INSTRUCTIONS) || defined(HAVE_AVX512_INSTRUCTIONS)
#include <immintrin.h>
#include "popcnt-lookup.hpp"
#endif // defined(HAVE_AVX2_INSTRUCTIONS) || defined(HAVE_AVX512_INSTRUCTIONS)
#if defined(HAVE_AVX512BW_INSTRUCTIONS)
#include "popcnt-avx512-harley-seal.hpp"
#elif defined(HAVE_AVX2_INSTRUCTIONS)
#include "popcnt-avx2-harley-seal.hpp"
#endif


namespace detail {

  inline int simd_intersect__(int c, const uint64_t* A, const uint64_t* B, uint64_t* out) {
      int w = 0;

      // currently the fastest code, SIMD challenging but possible :-)
      for (int i = 0; i < c; i += 1) {
          out[i] = A[i] & B[i];
#if !defined(HAVE_AVX2_INSTRUCTIONS) && !defined(HAVE_AVX512BW_INSTRUCTIONS)
          w += __builtin_popcountll(out[i]);
#endif
      }
#if defined(HAVE_AVX512BW_INSTRUCTIONS)
      w = popcnt_AVX512_harley_seal(reinterpret_cast<const uint8_t*>(out), c*8);
#elif defined(HAVE_AVX2_INSTRUCTIONS)
      w = popcnt_AVX2_harley_seal(reinterpret_cast<const uint8_t*>(out), c*8);
#endif
      return w;
  } // simd_intersect__

} // namespace detail


class bitvector {
public:
    static constexpr size_t unit() {
        return sizeof(uint64_t) * 8;
    }

    static constexpr size_t align() {
#if defined(HAVE_AVX512BW_INSTRUCTIONS)
        return 64;
#elif defined(HAVE_AVX2_INSTRUCTIONS)
        return 32;
#else
        return 1;
#endif
    }

    explicit bitvector(int n = 0, bool unsafe = false)
        : size_(std::max(0, n)), capacity_((size_ + unit() - 1) / unit()), weight_(0), state_(NORMAL) {
        if (capacity_ == 0) ptr_ = nullptr;
        else {
            ptr_ = (uint64_t*)_mm_malloc(sizeof(uint64_t)*capacity_, align());
            if (!unsafe) std::fill(ptr_, ptr_ + capacity_, 0);
        }
    } // bitvector

    bitvector(const bitvector& A) : size_(A.size_), capacity_(A.capacity_), weight_(A.weight_), state_(A.state_) {
        if (capacity_ == 0) ptr_ = nullptr;
        else {
            // we may have to use placement new with posix_memalign
            ptr_ = (uint64_t*)_mm_malloc(sizeof(uint64_t)*capacity_, align());
            // we may consider memcpy since uint64_t is POD
            // not sure if it is worth the effort
            std::copy(A.ptr_, A.ptr_ + capacity_, ptr_);
        }
    } // bitvector

    bitvector(bitvector&& A) : size_(A.size_), capacity_(A.capacity_), weight_(A.weight_), state_(A.state_) {
        if (capacity_ == 0) ptr_ = nullptr;
        else {
            ptr_ = A.ptr_;
            A.ptr_ = nullptr;
            A.size_ = 0;
            A.capacity_ = 0;
        }
    } // bitvector

    bitvector& operator=(const bitvector& A) {
        if (this == &A) return *this;
        _mm_free(ptr_);

        size_ = A.size_;
        capacity_ = A.capacity_;
        weight_ = A.weight_;
        state_ = A.state_;

        if (capacity_ == 0) ptr_ = nullptr;
        else {
            ptr_ = (uint64_t*)_mm_malloc(sizeof(uint64_t)*capacity_, align());
            std::copy(A.ptr_, A.ptr_ + capacity_, ptr_);
        }
        return *this;
    } // operator=

    bitvector& operator=(bitvector&& A) {
        if (this == &A) return *this;
        _mm_free(ptr_);

        size_ = A.size_;
        capacity_ = A.capacity_;
        weight_ = A.weight_;
        state_ = A.state_;
        ptr_ = A.ptr_;

        A.ptr_ = nullptr;
        A.size_ = 0;
        A.capacity_ = 0;

        return *this;
    } // operator=

    ~bitvector() { _mm_free(ptr_); }

    bool empty() const { return weight() == 0; }

    int size() const { return size_; }

    int weight() const {
        // if (state_ == IDENTITY) return size_;
        // if (state_ == EMPTY) return 0;
        return weight_;
    } // weight

    bool operator[](int n) const {
        if (state_ == IDENTITY) return true;
        if (state_ == EMPTY) return false;
        return (ptr_[n >> 6] >> (n & 63)) & 1;
    } // operator[]

    // caller MUST ensure that insert is never called on a bit that is already in the vector
    void insert(int n) {
        // if (operator[](n)) return; // bit already a 1
        ptr_[n >> 6] |= (1LL << (n & 63));
        ++weight_;
    } // insert


    friend void intersect(const bitvector& A, const bitvector& B, bitvector& out) {
        int isz = std::max(B.size_, A.size_);

        if (A.state_ == IDENTITY) {
            out = B;
            return;
        }

        if (B.state_ == IDENTITY) {
            out = A;
            return;
        }

        if (A.state_ == EMPTY || B.state_ == EMPTY) {
            out.state_ = EMPTY;
            return;
        }

        if (out.capacity_ * 64 < isz) out = std::move(bitvector(A.size_, true));

        int c = std::min(B.capacity_, A.capacity_);
        out.weight_ = detail::simd_intersect__(c, A.ptr_, B.ptr_, out.ptr_);
    } // intersect


    static bitvector identity(int size) { return bitvector(IDENTITY, size); }


private:
    enum bv_state { EMPTY, NORMAL, IDENTITY };

    explicit bitvector(bv_state s, int m) : ptr_(nullptr), size_(m), capacity_(0), weight_(s == IDENTITY ? m : 0), state_(s)  { }

    uint64_t* ptr_ = nullptr;

    int size_ = 0;
    int capacity_ = 0;
    int weight_ = 0;

    bv_state state_ = EMPTY;

}; // class bitvector

#endif // BITVECTOR_HPP
