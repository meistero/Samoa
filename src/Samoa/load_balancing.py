#!/usr/bin/python

import random
import numpy
import argparse


def low_freq_gauss(mu, sigma, n, frequency = 13):
    noise_low = [max(1, int(random.gauss(mu, sigma))) for i in xrange(0, frequency)]
    return [noise_low[i * (frequency - 1) / (n - 1)] for i in xrange(0, n)]

def sine_noise(mu, sigma, n, amplitudes):
    octaves = int(math.log(n)) + 1
    phases = [random.random() for octave in xrange(0, octaves)]
    f = lambda x: sum([amplitudes(octave + 1) * math.sin(2.0 * math.pi * float(2 ** (octave - 1)) * float(x) / float(n) + phases[octave]) for octave in xrange(0, octaves)])
    return [max(1, int(mu + sigma * f(x))) for x in xrange(0, n)]

def dint(arr):
    arr_int = [0 for i in xrange(0, len(arr) + 1)]

    arr_int[0] = 0
    for i in xrange(0, len(arr)):
        arr_int[i + 1] = arr_int[i] + arr[i]

    return arr_int

def dder(arr_int):
    return [arr_int[i] - arr_int[i - 1] for i in xrange(1, len(arr_int))]

def exact(L, c):
    n = len(L) - 1

    sol_exact_int = [L[n] * i / c for i in xrange(0, c + 1)]
    return dder(sol_exact_int)

def start_cutoff(L, c):
    n = len(L) - 1
    j = 0
    l = dder(L)

    sol_approx = [0 for i in xrange(0, c)]

    for j in xrange(0, n):
        i = (c * L[j]) / L[n]

        sol_approx[i] += l[j]

    return sol_approx

def midpoint_cutoff(L, c):
    n = len(L) - 1
    j = 0
    l = dder(L)

    sol_approx = [0 for i in xrange(0, c)]

    for j in xrange(0, n):
        i = (c * (L[j] + L[j + 1])) / (2 * L[n])

        sol_approx[i] += l[j]

    return sol_approx

def iterative_binary(L, c):
    n = len(L) - 1

    sol_iterative_int = dint(midpoint_cutoff(L, c))
    current_max = max(dder(sol_iterative_int))
    current_min = max(L[n] / c, max(dder(L)))

    sol_test_int = [sol_iterative_int[i] for i in xrange(0, c + 1)]

    iters = 0
    test = (current_min + current_max) / 2

    while current_max > current_min:
        j = 0
        for i in xrange(1, c):
            #S_i := max(L_j; L_j \in L and L_j < S_{i-1} + test)
            #S_i < S_k + (i - k) * test

            while j < n and L[j + 1] < sol_test_int[i - 1] + test:
                j += 1

            sol_test_int[i] = L[j]

        if L[n] < sol_test_int[c - 1] + test:
            current_max = max(dder(sol_test_int))
            test = (current_min + current_max) / 2

            for i in xrange(1, c):
                sol_iterative_int[i] = sol_test_int[i]
        else:
            current_min = test
            test = current_max

        iters = iters + 1

    return dder(sol_iterative_int), iters

def iterative_linear(L, c):
    n = len(L) - 1

    sol_iterative_int = dint(midpoint_cutoff(L, c))
    current_max = max(dder(sol_iterative_int)) + 1

    sol_test_int = [sol_iterative_int[i] for i in xrange(0, c + 1)]

    iters = 0
    converged = False
    while(not converged):
        j = 0
        for i in xrange(1, c):
            while j < n and L[j + 1] < sol_test_int[i - 1] + current_max:
                j += 1

            sol_test_int[i] = L[j]

        if L[n] < sol_test_int[c - 1] + current_max:
            current_max = max(dder(sol_test_int))

            for i in xrange(1, c):
                sol_iterative_int[i] = sol_test_int[i]
        else:
            converged = True

        iters = iters + 1

    return dder(sol_iterative_int), iters

def brute_force(L, c):
    n = len(L) - 1

    def brute_force_recursive(d, i_p):
        if (d >= c):
            if max(dder(sol_test_int)) < max(dder(sol_brute_int)):
                for i in xrange(0, c + 1):
                    sol_brute_int[i] = sol_test_int[i]
        else:
            for i in xrange(i_p, n + 1):
                sol_test_int[d] = L[i]
                brute_force_recursive(d + 1, i + 1)

    sol_brute_int = [0 for i in xrange(0, c + 1)]
    sol_brute_int[c] = L[n]

    sol_test_int = [0 for i in xrange(0, c + 1)]
    sol_test_int[c] = L[n]

    brute_force_recursive(1, 0)
    return dder(sol_brute_int)

def main():
    parser = argparse.ArgumentParser(description='Tests load balancing strategies for distribution of n sections to c cores.')

    parser.add_argument('-c', metavar='<value>', type=int, default=8,
                       help='number of cores')

    parser.add_argument('-n', metavar='<value>', type=int, default=16,
                       help='number of sections')

    parser.add_argument('-mu', metavar='<value>', type=int, default=1000,
                       help='load average (load per section is gauss distributed around mu with a standard deviation of epsilon)')

    parser.add_argument('-sigma', metavar='<value>', type=int, default=100,
                       help='load deviation (load per section is gauss distributed around mu with a standard deviation of sigma)')

    args = parser.parse_args()

    c = args.c
    n = args.n
    mu = args.mu
    sigma = args.sigma

    i_tests = 0
    start_imbalance = 0
    midpoint_imbalance = 0
    optimal_imbalance = 0
    i_binary_iterations = 0
    i_linear_iterations = 0

    while(True):
        l = low_freq_gauss(mu, sigma, n)
        L = dint(l)

        print "Testing..", "avg(l_i):", sum(l) / n, "max(l_i):", max(l), "sum(l_i) / c:", L[n] / c

        sol_exact = exact(L, c)

        sol_approx = start_cutoff(L, c)
        assert(sum(sol_approx) == sum(sol_exact))
        print " Start cutoff         ", "maximum:", max(sol_approx), "imbalance:", 100 * (max(sol_approx) - 1) / max(sol_exact) - 99, "%"
        start_imbalance += 100 * (max(sol_approx) - 1) / max(sol_exact) - 99

        sol_approx = midpoint_cutoff(L, c)
        assert(sum(sol_approx) == sum(sol_exact))
        print " Midpoint cutoff      ", "maximum:", max(sol_approx), "imbalance:", 100 * (max(sol_approx) - 1) / max(sol_exact) - 99, "%"
        midpoint_imbalance += 100 * (max(sol_approx) - 1) / max(sol_exact) - 99

        sol_binary_search, iters = iterative_binary(L, c)
        assert(sum(sol_binary_search) == sum(sol_exact))
        print " Binary search        ", "maximum:", max(sol_binary_search), "imbalance:", 100 * (max(sol_binary_search) - 1) / max(sol_exact) - 99, "%"
        i_binary_iterations += iters
        optimal_imbalance += 100 * (max(sol_binary_search) - 1) / max(sol_exact) - 99

        if (n < 33000):
            sol_linear_search, iters = iterative_linear(L, c)
            assert(sum(sol_linear_search) == sum(sol_exact))
            print " Linear search        ", "maximum:", max(sol_linear_search), "imbalance:", 100 * (max(sol_linear_search) - 1) / max(sol_exact) - 99, "%"
            i_linear_iterations += iters

            #Check if the binary search maximum equals the linear search maximum
            if max(sol_binary_search) != max(sol_linear_search):
                print "Error: Binary search maximum is bigger than linear search maximum:"
                print sol_binary_search, "maxima: ", max(sol_binary_search), max(sol_linear_search)
                assert(False)

        if (n < 24):
            sol_brute = brute_force(L, c)
            assert(sum(sol_brute) == sum(sol_exact))
            print " Brute force solution ", "maximum:", max(sol_brute), "imbalance:", 100 * (max(sol_brute) - 1) / max(sol_exact) - 99, "%"

            #Check if the linear search maximum equals the brute force maximum
            if max(sol_linear_search) != max(sol_brute):
                print "Error: Linear search maximum is bigger than brute force solution maximum:"
                print sol_linear_search, "maxima: ", max(sol_linear_search), max(sol_brute)
                assert(False)


        i_tests += 1
        print "Test complete! (total tests:", i_tests, ", start imbalance:", start_imbalance / i_tests, "%, midpoint imbalance:", midpoint_imbalance / i_tests, "%, optimal imbalance:", optimal_imbalance / i_tests, "%, binary search iterations:", i_binary_iterations / i_tests, "on average", "linear search iterations:", i_linear_iterations / i_tests, "on average)"
main()

