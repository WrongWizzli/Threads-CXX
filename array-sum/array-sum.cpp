#include <thread>
#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <mutex>
#include <assert.h>
#include "omp.h"
#include <fstream>

int64_t CountSumSequential(const std::vector<int> &v) {
    int64_t sum = 0;
    for (int i = 0; i < v.size(); ++i) {
        sum += v[i];
    }
    return sum;
}

int64_t CountSumOMP(const std::vector<int> &v) {
    int64_t sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < v.size(); ++i) {
        sum += v[i];
    }
    return sum;
}

void CountSumOneThread(int64_t &thread_sum, int64_t thread_id, int64_t n_threads, const std::vector<int> &v) {
    std::size_t start_pos = thread_id * v.size() / n_threads;
    std::size_t end_pos = (thread_id + 1) * v.size() / n_threads;
    int64_t sum = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = start_pos; i < end_pos; ++i) {
        sum += v[i];
    }
    thread_sum = sum;
    auto t1 = std::chrono::steady_clock::now();
}

int64_t CountSumThreads(const std::vector<int> &v) {
    unsigned int n_threads = std::thread::hardware_concurrency();
    std::cout << "Max concurrency threads: " << n_threads << std::endl;

    if (v.size() < 1000 * n_threads) {
        return CountSumSequential(v);
    }
    int64_t sum = 0;
    std::vector<std::thread> threads;
    std::vector<int64_t> thread_sums(n_threads, 0);
    for (int i = 0; i < n_threads; ++i) {
        threads.push_back(
            std::thread(
                CountSumOneThread, 
                std::ref(thread_sums[i]), 
                i, n_threads, std::cref(v))
        );
    }
    for (int i = 0; i < n_threads; ++i) {
        threads[i].join();
    }
    for (int i = 0; i < n_threads; ++i) {
        sum += thread_sums[i];
    }
    return sum;
}

std::vector<int> GetRandomVector(std::size_t size) {
    std::vector<int> v(size);
    srand(std::chrono::steady_clock::now().time_since_epoch().count());
    for (int i = 0; i < v.size(); ++i) {
        v[i] = rand() % 100;
    }
    return v;
}

int main(void) {
    std::vector<int64_t> vec_sizes(
        {200, 2000, 20000, 200000, 2000000, 
         20000000, 2000000000, 2000000000L}
    );
    std::ofstream out("time_results.csv");
    out << "vector_size,seq_time,omp_time,thr_time\n";
    for (auto vec_size: vec_sizes) {
        auto v = GetRandomVector(vec_size);
        auto t0 = std::chrono::steady_clock::now();
        int64_t sum1 = CountSumThreads(v);
        auto t1 = std::chrono::steady_clock::now();
        int64_t sum2 = CountSumSequential(v);
        auto t2 = std::chrono::steady_clock::now();
        int64_t sum3 = CountSumOMP(v);
        auto t3 = std::chrono::steady_clock::now();

        assert(sum1 == sum2 && sum2 == sum3);
        auto thr_time = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        auto seq_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        auto omp_time = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
        out << vec_size << ",";
        out << seq_time / 1000. << ",";
        out << omp_time / 1000. << ",";
        out << thr_time / 1000. << ",\n";
    }
    return 0;
}