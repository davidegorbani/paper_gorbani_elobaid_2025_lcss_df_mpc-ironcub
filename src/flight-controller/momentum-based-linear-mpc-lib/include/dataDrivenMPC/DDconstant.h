#ifndef CONSTANT_DD_H
#define CONSTANT_DD_H

#include <array>

#define N_JOINTS 8
#define N_THRUSTS 4

constexpr std::array<int, 3> CoMPosIdx = {0, 1, 2};
constexpr std::array<int, 3> linMomIdx = {3, 4, 5};
constexpr std::array<int, 3> rpyIdx = {6, 7, 8};
constexpr std::array<int, 3> angMomIdx = {9, 10, 11};

constexpr std::array<int, N_THRUSTS> thrustHatIdx = {0, 1, 2, 3};
constexpr std::array<int, N_THRUSTS> thrustHatDotIdx = {4, 5, 6, 7};

constexpr std::array<int, N_JOINTS> createJointIdx()
{
    std::array<int, N_JOINTS> temp = {};
    for (int i = 0; i < N_JOINTS; ++i)
    {
        temp[i] = i;
    }
    return temp;
}

constexpr std::array<int, N_JOINTS> deltaJointIdx = createJointIdx();

constexpr std::array<int, N_THRUSTS> createThrustIdx()
{
    std::array<int, N_THRUSTS> temp = {};
    for (int i = 0; i < N_THRUSTS; ++i)
    {
        temp[i] = i;
    }
    return temp;
}

constexpr std::array<int, N_THRUSTS> createThrottleIdx()
{
    std::array<int, N_THRUSTS> temp = {};
    for (int i = 0; i < N_THRUSTS; ++i)
    {
        temp[i] = i;
    }
    return temp;
}

constexpr std::array<int, N_THRUSTS> thrustIdx = createThrustIdx();

constexpr std::array<int, N_THRUSTS> throttleIdx = createThrottleIdx();

#endif // CONSTANT_DF_H
