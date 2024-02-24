"""
Metodi di Analisi Applicata.
Andrea Di Antonio, 858798.
Assignment 3.
Version 2.

For code description, refer to the README.md file.
"""

print("Andrea Di Antonio, 858798.")


# LIBRARIES


from time import time as timer
import numpy as np

import matplotlib.pyplot as plt


# NUMERICS.


# Numerical method.
def finiteVolume(
    func: callable,
    startingCondition: callable,
    spaceRange: list,
    timeRange: list,
    spaceSteps: int = 2**10,
    timeSteps: int = 2**12,
    last: bool = False,
) -> np.array:
    global spaceStep
    global timeStep

    space, spaceStep = np.linspace(
        spaceRange[0], spaceRange[1], spaceSteps, retstep=True
    )
    time, timeStep = np.linspace(timeRange[0], timeRange[1], timeSteps, retstep=True)

    assert maxDerivative * timeStep / spaceStep < 1  # CFL.

    solution = np.zeros([timeSteps, startingCondition(space).size])
    solution[0, :] = startingCondition(space)

    for i in range(1, timeSteps):
        upSolution = np.concatenate([solution[i - 1, :], [solution[i - 1, -1]]])
        downSolution = np.concatenate([[solution[i - 1, 0]], solution[i - 1, :]])

        solution[i, :] = solution[i - 1, :] - timeStep / spaceStep * (
            func(solution[i - 1, :], upSolution[1:])
            - func(downSolution[:-1], solution[i - 1, :])
        )

    return (solution[-1, :] if last else solution), space, time


# Numerical fluxes.
def laxFriedrichs(x: float, y: float) -> float:  # Works with numpy.array.
    return 0.5 * (flux(x) + flux(y)) - 0.5 * spaceStep / timeStep * (y - x)


def godunov(x: float, y: float) -> float:  # Works with numpy.array.
    return (
        (x == y) * flux(x)
        + (x < y) * ((flux(x) > flux(y)) * flux(y) + (flux(x) <= flux(y)) * flux(x))
        + (x > y)
        * (
            (0 <= fluxPrime(x)) * flux(x)
            + (0 >= fluxPrime(y)) * flux(y)
            + (fluxPrime(x) < 0) * (0 < fluxPrime(y)) * flux(invFluxPrimeZero)
        )
    )


# Parabolic flux.
def parabolicFlux(x: float) -> float:  # Works with numpy.array.
    return x * (1 - x)


def parabolicFluxPrime(x: float) -> float:  # Works with numpy.array.
    return 1 - 2 * x


# L1 norm.
def normL1(x: np.array, dx: float) -> float:
    return np.sum(np.abs(x)) * dx


# Total variation.
def totalVariation(x: np.array) -> float:
    return np.sum(np.abs(x[1:] - x[:-1]))


# DATA AND FUNCTIONS.


# Variables.
bar_x = -1

# Flux data.
flux = parabolicFlux
fluxPrime = parabolicFluxPrime
invFluxPrimeZero = 0.5
maxDerivative = 1


# Analytical functions.
def analytical_1(t: float, x: float) -> float:
    gamma = lambda t: t - 2 * np.sqrt(-bar_x * t)
    return (1.0 * (x >= bar_x) * (x <= -t) + 0.5 * (1 - x / t) * (abs(x) <= t)) * (
        t <= -bar_x
    ) + (0.5 * (1 - x / t) * (x >= gamma(t)) * (x <= t)) * (t >= -bar_x)


def analytical_2(t: float, x: float) -> float:
    gamma = lambda t: -1 + np.sqrt(2 * t)
    return (
        1.0 * (x <= -1 - t)
        + 0.5 * (1 - (x + 1) / t) * (x > -1 - t) * (x < -1 + t) * (t < 2)
        + 0.5 * (1 - (x + 1) / t) * (x >= -1 - t) * (x <= gamma(t)) * (t >= 2)
        + 0.5 * (x >= -1 + t) * (t < 2)
        + 0.5 * (x >= gamma(t)) * (t >= 2)
    )


# Starting conditions.
numerical_1 = lambda x: x * 0.0 + 1.0 * (x >= bar_x) * (x <= 0)
numerical_2 = lambda x: x * 0.0 + 1.0 * (x <= -1) + 0.5 * (x >= 0)


# SECTION 1.


execution = timer()

# Space and time ranges.
timeRange = [0, 10]
spaceRange = [-25, 25]

# Tests variables.
tests = 10
spaceStepsNumbers = (2 ** np.linspace(7, 13, tests)).astype(int)
timeStepsNumbers = (2 ** np.linspace(9, 15, tests)).astype(int)

results_1_1 = np.zeros((tests, 2))
results_1_2 = np.zeros((tests, 2))

start = timer()

for j in range(tests):  # L1 norm comparison.
    LF_1, _, _ = finiteVolume(
        laxFriedrichs,
        numerical_1,
        spaceRange,
        timeRange,
        spaceStepsNumbers[j],
        timeStepsNumbers[j],
        last=True,
    )
    G_1, _, _ = finiteVolume(
        godunov,
        numerical_1,
        spaceRange,
        timeRange,
        spaceStepsNumbers[j],
        timeStepsNumbers[j],
        last=True,
    )

    LF_2, _, _ = finiteVolume(
        laxFriedrichs,
        numerical_2,
        spaceRange,
        timeRange,
        spaceStepsNumbers[j],
        timeStepsNumbers[j],
        last=True,
    )
    G_2, space, time = finiteVolume(
        godunov,
        numerical_2,
        spaceRange,
        timeRange,
        spaceStepsNumbers[j],
        timeStepsNumbers[j],
        last=True,
    )

    results_1_1[j, 0] = normL1(
        LF_1 - analytical_1(timeRange[1], space), space[1] - space[0]
    )
    results_1_1[j, 1] = normL1(
        G_1 - analytical_1(timeRange[1], space), space[1] - space[0]
    )

    results_1_2[j, 0] = normL1(
        LF_2 - analytical_2(timeRange[1], space), space[1] - space[0]
    )
    results_1_2[j, 1] = normL1(
        G_2 - analytical_2(timeRange[1], space), space[1] - space[0]
    )

stop = timer() - start

print(
    "\nTime taken for section 1's {} test(s): {} second(s).\n\nSection 1, first test results:\n\nSpace\tTime\tLax-Friedrichs\tGodunov\n".format(
        tests, round(stop, 2)
    )
)

for j in range(tests):
    print(
        "{}\t{}\t{:.4e}\t{:.4e}".format(
            spaceStepsNumbers[j],
            timeStepsNumbers[j],
            results_1_1[j, 0],
            results_1_1[j, 1],
        )
    )

print("\nSection 1, second test results:\n\nSpace\tTime\tLax-Friedrichs\tGodunov\n")

for j in range(tests):
    print(
        "{}\t{}\t{:.4e}\t{:.4e}".format(
            spaceStepsNumbers[j],
            timeStepsNumbers[j],
            results_1_2[j, 0],
            results_1_2[j, 1],
        )
    )


# SECTION 2.


# Tests variables.
spaceSteps = spaceStepsNumbers[-1]
timeSteps = timeStepsNumbers[-1]

start = timer()

LF_1, _, _ = finiteVolume(
    laxFriedrichs, numerical_1, spaceRange, timeRange, spaceSteps, timeSteps
)
G_1, _, _ = finiteVolume(
    godunov, numerical_1, spaceRange, timeRange, spaceSteps, timeSteps
)

LF_2, _, _ = finiteVolume(
    laxFriedrichs, numerical_2, spaceRange, timeRange, spaceSteps, timeSteps
)
G_2, space, time = finiteVolume(
    godunov, numerical_2, spaceRange, timeRange, spaceSteps, timeSteps
)

results_2_1 = np.zeros((time.size, 3))
results_2_2 = np.zeros((time.size, 3))

for j in range(1, time.size):  # Skips t = 0.
    results_2_1[j, 0] = totalVariation(analytical_1(time[j], space))
    results_2_1[j, 1] = totalVariation(LF_1[j, :])
    results_2_1[j, 2] = totalVariation(G_1[j, :])

    results_2_2[j, 0] = totalVariation(analytical_2(time[j], space))
    results_2_2[j, 1] = totalVariation(LF_2[j, :])
    results_2_2[j, 2] = totalVariation(G_2[j, :])

stop = timer() - start

print(
    "\nTime taken for section 2's test: {} second(s)\n\nSection 2, some first test results:\n\nTime\tAnalytical 1\tLax-Friedrichs\tGodunov\n".format(
        round(stop, 2)
    )
)

print("Space steps: {}\nTime steps: {}\n".format(spaceSteps, timeSteps))

for j in range(1, time.size, time.size // 9):
    print(
        "{:.4f}\t{:.4e}\t{:.4e}\t{:.4e}".format(
            time[j], results_2_1[j, 0], results_2_1[j, 1], results_2_1[j, 2]
        )
    )

print(
    "\nSection 2, some second test results:\n\nTime\tAnalytical 2\tLax-Friedrichs\tGodunov\n"
)

for j in range(1, time.size, time.size // 9):
    print(
        "{:.4f}\t{:.4e}\t{:.4e}\t{:.4e}".format(
            time[j], results_2_2[j, 0], results_2_2[j, 1], results_2_2[j, 2]
        )
    )


# PLOTTING.

fig, ax = plt.subplots(2, 2)

ax[0, 0].set_title("Section 1, first test.")
ax[0, 0].set_ylim(-0.25, 1.25)
ax[0, 0].plot(space, G_1[0, :], label="Time: {}".format(time[0]))
ax[0, 0].plot(
    space,
    analytical_1(time[-1], space),
    label="Analytical, time: {}.".format(time[-1]),
)
ax[0, 0].plot(space, LF_1[-1, :], label="Lax-Friedrichs, time: {}.".format(time[-1]))
ax[0, 0].plot(space, G_1[-1, :], label="Godunov: {}.".format(time[-1]))

ax[0, 1].set_title("Section 1, second test.")
ax[0, 1].set_ylim(-0.25, 1.25)
ax[0, 1].plot(space, G_2[0, :], label="Time: {}".format(time[0]))
ax[0, 1].plot(
    space,
    analytical_2(time[-1], space),
    label="Analytical, time: {}.".format(time[-1]),
)
ax[0, 1].plot(space, LF_2[-1, :], label="Lax-Friedrichs, time: {}.".format(time[-1]))
ax[0, 1].plot(space, G_2[-1, :], label="Godunov: {}.".format(time[-1]))

ax[1, 0].set_title("Section 2, first test.")
ax[1, 0].plot(time, results_2_1[:, 0], label="Analytical's total variation.")
ax[1, 0].plot(time, results_2_1[:, 1], label="Lax-Friedrichs' total variation.")
ax[1, 0].plot(time, results_2_1[:, 2], label="Godunov's total variation.")
ax[1, 0].legend(loc="best")

ax[1, 1].set_title("Section 2, second test.")
ax[1, 1].plot(time, results_2_2[:, 0], label="Analytical's total variation.")
ax[1, 1].plot(time, results_2_2[:, 1], label="Lax-Friedrichs' total variation.")
ax[1, 1].plot(time, results_2_2[:, 2], label="Godunov's total variation.")
ax[1, 1].legend(loc="best")

ax[0, 0].legend(loc="best")
ax[0, 1].legend(loc="best")


# END.


print(
    "\nFinal results have been plotted.\n\nTime taken: {} second(s).".format(
        round(timer() - execution, 2)
    )
)
plt.show()
