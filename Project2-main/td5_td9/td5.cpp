#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <cstring>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Vector.h"

std::normal_distribution<double> normalDist(0.0, 1.0);
static std::mt19937_64 engine{ std::random_device{}() };

void normalizeVector(double &x, double &y, double &z) {
    double magnitude = std::sqrt(x * x + y * y + z * z);
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
}

void generateRandomUnitVector(double &x, double &y, double &z) {
    x = normalDist(engine);
    y = normalDist(engine);
    z = normalDist(engine);
    normalizeVector(x, y, z);
    return ;
}

void projectAndSortPixels(const double *sourceImage, const double *targetImage, double *outputImage, int numPixels, double x, double y, double z,
                            std::vector<std::pair<double, int>> &sortedSource, std::vector<double> &sortedTarget) {
    for (int i = 0; i < numPixels; ++i) {

        double sourceProjection = outputImage[i * 3] * x +  outputImage[i * 3 + 1] * y + outputImage[i * 3 + 2] * z;
        double targetProjection = targetImage[i * 3] * x + targetImage[i * 3 + 1] * y +  targetImage[i * 3 + 2] * z;
        sortedSource[i] = std::make_pair(sourceProjection , i);
        sortedTarget[i] = targetProjection;
    }

    std::sort(sortedSource.begin(), sortedSource.end());
    std::sort(sortedTarget.begin(), sortedTarget.end());
}

void adjustImagePixels(double *outputImage, const std::vector<std::pair<double, int>> &sortedSource, const std::vector<double> &sortedTarget, int numPixels, double x, double y, double z) {


    for (int i = 0; i < numPixels; ++i) {

        double difference = sortedTarget[i] - sortedSource[i].first;
        int pixelIndex = sortedSource[i].second;


        outputImage[pixelIndex * 3] += difference * x;
        outputImage[pixelIndex * 3 + 1] += difference * y;
        outputImage[pixelIndex * 3 + 2] += difference * z;
    }
}

void sliced_matching(const double *sourceImage, double *targetImage, int width, int height, double *outputImage) {
    int numPixels = width * height;

    std::memcpy(outputImage, sourceImage, numPixels * 3 * sizeof(double));

    for (int iteration = 0; iteration < 100; ++iteration) {
        double x, y, z;
        generateRandomUnitVector(x, y, z);


        std::vector<std::pair<double, int>> sortedSource(numPixels);
        std::vector<double> sortedTarget(numPixels);

        projectAndSortPixels(sourceImage, targetImage, outputImage, numPixels, x, y, z, sortedSource, sortedTarget);
        adjustImagePixels(outputImage, sortedSource, sortedTarget, numPixels, x, y, z);
    }
}

unsigned char clampToByte(double value, double minValue, double maxValue) {
    return static_cast<unsigned char>(std::max(minValue, std::min(value, maxValue)));

}