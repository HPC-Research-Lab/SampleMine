#pragma once

#include <chrono>
#include <iostream>
#include <vector>
#include <set>
#include <array>

namespace euler::util {

  class Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> begin;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double duration;

  public:
    Timer() : duration(0.0) {}

    void start() { begin = std::chrono::high_resolution_clock::now(); duration = 0; }

    void stop() {
      end = std::chrono::high_resolution_clock::now();
      duration += (std::chrono::time_point_cast<std::chrono::microseconds>(end)
        .time_since_epoch()
        .count() -
        std::chrono::time_point_cast<std::chrono::microseconds>(begin)
        .time_since_epoch()
        .count()) *
        1e-6;
    }

    double get() { return duration; }
  };

  template <typename T>
  class span {
    const T* addr;
    size_t len;

  public:
    span() : addr(nullptr), len(0) {}
    span(const T* a, size_t l) : addr(a), len(l) {}
    const T* data() { return addr; }
    size_t size() { return len; }
    T& operator[](size_t idx) {return addr[idx];}
  };

  template <class T>
  void print_vec(const std::vector<T>& a) {
    for (T i : a) std::cout << i << " ";
    std::cout << std::endl;
  }

  template <class T, size_t S>
  void print_vec(const std::array<T, S>& a) {
    for (T i : a) std::cout << i << " ";
    std::cout << std::endl;
  }

  double random_number();



  /*// atomically adds 128-bit src to dst, returning the old dst.
  inline __uint128_t fetch_add(__uint128_t *dst, __uint128_t src)
  {
      __uint128_t dstval, olddst;

      dstval = *dst;

      do
      {
          olddst = dstval;
          dstval = __sync_val_compare_and_swap(dst, dstval, dstval + src);
      }
      while(dstval != olddst);

      return dstval;
  }

  inline __uint64_t fetch_add(__uint64_t *dst, __uint64_t src) {
    __sync_fetch_and_add(dst, src);
    return *dst;
  }*/

}  // namespace euler::util