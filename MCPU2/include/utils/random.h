#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <vector>
#include <stdexcept>

template <typename T>

T selectElementWithWeights(const std::vector<T>& elements, const std::vector<double>& weights) {
    if (elements.size() != weights.size()) {
        throw std::invalid_argument("Elements and weights must have the same size");
    }

    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a discrete distribution based on the weights
    std::discrete_distribution<> dist(weights.begin(), weights.end());

    // Select an element based on the distribution
    int index = dist(gen);
    return elements[index];
}

double generateGaussian(double mean = 0.0, double stddev = 1.0) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::normal_distribution<> d(mean, stddev);
    return d(gen);
}

#endif // RANDOM_H