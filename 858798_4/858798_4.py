"""
Metodi di Analisi Applicata.
Andrea Di Antonio, 858798.
Consegna 4.
Versione 1.

Per la descrizione del codice si veda il file README.md

Chiamare il programma con 'python 858798_4.py save' per salvare l'output grafico.
"""

print("Andrea Di Antonio, 858798.\n")


# LIBRARIES


from typing import Iterable, Callable
from time import time as timer
import numpy as np
import sys

import matplotlib.pyplot as plt
import matplotlib.animation as animation


# CLASSES


class road:
    def __init__(self, name: str, span: Iterable[float]) -> None:
        assert isinstance(name, str)
        assert name

        assert isinstance(span, Iterable)
        assert len(span) == 2
        assert span[0] < span[1]

        assert isinstance(span[0], (int, float))
        assert isinstance(span[1], (int, float))

        assert isinstance(sSteps, int)
        assert sSteps > 0

        self.name = name

        self.span = span

        self.space, self.sStep = np.linspace(span[0], span[1], sSteps, retstep=True)
        self.solution = np.zeros_like(self.space)

        self.extension: list[list[float]] = [[], []]

    def __repr__(self) -> str:
        return self.name

    def startingCondition(self, func: Callable) -> None:
        self.solution = func(self.space)

    def extend(self) -> np.array:
        solution = self.solution.copy()

        if not self.extension[0]:
            downSolution = np.concatenate([[self.solution[0]], solution])

        else:
            downSolution = np.concatenate([self.extension[0], solution])

        if not self.extension[1]:
            upSolution = np.concatenate([solution, [self.solution[-1]]])

        else:
            upSolution = np.concatenate([solution, self.extension[1]])

        return downSolution, upSolution


class intersection:
    def __init__(
        self, name: str, inbound: list[road], outbound: list[road], coeff: float = 0.5
    ) -> None:
        assert isinstance(name, str)
        assert name

        assert isinstance(inbound, list)
        assert isinstance(outbound, list)

        assert len(inbound) in {1, 2}
        assert len(outbound) in {1, 2}

        # Only supports 1x2, 2x1 intersections.
        assert len(inbound) + len(outbound) == 3

        for r in inbound + outbound:
            assert isinstance(r, road)

        for rIn in inbound:
            for rOut in outbound:
                assert rIn.span[1] == rOut.span[0]

        assert isinstance(coeff, float)
        assert 0 < coeff < 1

        self.name = name

        self.inbound = inbound
        self.outbound = outbound

        # Distribution/priority coefficient.
        self.coeff = coeff

        # Intersection type.
        if len(self.inbound) == 2 and len(self.outbound) == 1:
            self.family = "2x1"

        elif len(self.inbound) == 1 and len(self.outbound) == 2:
            self.family = "1x2"

    def __repr__(self) -> str:
        return "{}: {} -> {} [{:.2f}]".format(
            self.name, self.inbound, self.outbound, self.coeff
        )

    def getBoundaries(self) -> tuple[float]:
        if self.family == "1x2":
            road_1 = self.inbound[0]
            sol_1 = road_1.solution

            road_2 = self.outbound[0]
            road_3 = self.outbound[1]
            sol_2 = road_2.solution
            sol_3 = road_3.solution

            gm_1 = flux(sigma) if sol_1[-1] >= sigma else flux(sol_1[-1])

            gm_2 = flux(sigma) if sol_2[0] <= sigma else flux(sol_2[0])
            gm_3 = flux(sigma) if sol_3[0] <= sigma else flux(sol_3[0])

            alpha = self.coeff
            gm_1 = min({gm_1, gm_2 / alpha, gm_3 / (1 - alpha)})

            gm_2 = alpha * gm_1
            gm_3 = (1 - alpha) * gm_1

            # Boundary densities.
            sol_1_b = invFlux(gm_1)[1]
            sol_2_b = invFlux(gm_2)[0]
            sol_3_b = invFlux(gm_3)[0]

        elif self.family == "2x1":
            road_1 = self.inbound[0]
            road_2 = self.inbound[1]
            sol_1 = road_1.solution
            sol_2 = road_2.solution

            road_3 = self.outbound[0]
            sol_3 = road_3.solution

            gm_1 = flux(sigma) if sol_1[-1] >= sigma else flux(sol_1[-1])
            gm_2 = flux(sigma) if sol_2[-1] >= sigma else flux(sol_2[-1])

            gm_3 = flux(sigma) if sol_3[0] <= sigma else flux(sol_3[0])

            if gm_1 + gm_2 > gm_3:  # Priority rules apply which redefine gm_1 and gm_2.
                q = self.coeff

                if q * gm_3 <= gm_1 and (1 - q) * gm_3 <= gm_2:
                    gm_1 = q * gm_3
                    gm_2 = (1 - q) * gm_3

                else:
                    if (1 - q) * gm_3 > gm_2:
                        gm_1 = q * gm_3

                    if q * gm_3 > gm_1:
                        gm_2 = (1 - q) * gm_3

            gm_3 = gm_1 + gm_2

            # Boundary densities.
            sol_1_b = invFlux(gm_1)[1]
            sol_2_b = invFlux(gm_2)[1]
            sol_3_b = invFlux(gm_3)[0]

        # Returns boundary densities.
        return sol_1_b, sol_2_b, sol_3_b


# GRAPHICS


def animate(index: int) -> None:
    plt.suptitle("Solution at time {:.2f}".format(time[index]))
    for j in range(len(roads)):
        plots[j].set_ydata(solution[index, :, j])


# NUMERICS.


# Numerical method.
def networkSolver(
    roads: list[road],
    intersections: list[intersection],
    timeRange: list,
) -> list:
    """
    Network Solver.

    Returns the numerical solution given multiple roads connected through intersections.
    """

    # Numerical solution.
    time, tStep = np.linspace(timeRange[0], timeRange[1], tSteps, retstep=True)
    solution = np.zeros((tSteps, sSteps, len(roads)))

    # Loads starting conditions.
    for j in range(len(roads)):
        solution[0, :, j] = roads[j].solution

    counter = timer()
    for j in range(1, tSteps):
        # Numerical solver "loader".
        if j in range(0, tSteps, tSteps // 25):
            print(
                "Numerical solution [{}/{} - {:.2f}s]".format(
                    j, tSteps, timer() - counter
                )
            )
            counter = timer()

        # Sets boundary densities.
        for inter in intersections:
            boundaries = inter.getBoundaries()

            if inter.family == "1x2":
                inter.inbound[0].extension[1] = [boundaries[0]]

                inter.outbound[0].extension[0] = [boundaries[1]]
                inter.outbound[1].extension[0] = [boundaries[2]]

            elif inter.family == "2x1":
                inter.inbound[0].extension[1] = [boundaries[0]]
                inter.inbound[1].extension[1] = [boundaries[1]]

                inter.outbound[0].extension[0] = [boundaries[2]]

        # Evolves roads.
        for k in range(len(roads)):
            # CFL on every road.
            assert maxDerivative * tStep / roads[k].sStep < 1

            downSolution, upSolution = roads[k].extend()

            solution[j, :, k] = solution[j - 1, :, k] - tStep / roads[k].sStep * (
                godunov(solution[j - 1, :, k], upSolution[1:])
                - godunov(downSolution[:-1], solution[j - 1, :, k])
            )

            # Updates roads for next extension.
            roads[k].solution = solution[j, :, k]
            roads[k].extension = [[], []]

    return solution, time


# Numerical flux.
def godunov(x: float, y: float) -> float:  # Works with numpy.array.
    return (
        (x == y) * flux(x)
        + (x < y) * ((flux(x) > flux(y)) * flux(y) + (flux(x) <= flux(y)) * flux(x))
        + (x > y)
        * (
            (0 <= fluxPrime(x)) * flux(x)
            + (0 >= fluxPrime(y)) * flux(y)
            + (fluxPrime(x) < 0) * (0 < fluxPrime(y)) * flux(sigma)
        )
    )


# Parabolic flux.
def parabolicFlux(x: float) -> float:  # Works with numpy.array.
    return x * (1 - x)


# Parabolic flux's first derivative.
def parabolicFluxPrime(x: float) -> float:  # Works with numpy.array.
    return 1 - 2 * x


# Inverse parabolic flux.
def invParabolicFlux(y: float) -> tuple[float]:
    return 0.5 * (1 - np.sqrt(1 - 4 * y)), 0.5 * (1 + np.sqrt(1 - 4 * y))


# DATA AND SIMULATION


# Simulation data.
# Slightly longer and more precise simulation on "long" as a parameter.
tSteps = 10**4 * (10 if "long" in sys.argv else 1)
sSteps = 10**3 * (10 if "long" in sys.argv else 1)

# Flux data.
flux = parabolicFlux
fluxPrime = parabolicFluxPrime
invFlux = invParabolicFlux

sigma = 0.5
maxDerivative = 1

# Roads.
I1 = road("I1", [0, 2])
I2 = road("I2", [2, 3])
I3 = road("I3", [2, 4])
I4 = road("I4", [3, 4])
I5 = road("I5", [3, 6])
I6 = road("I6", [4, 6])

# Intersection.
J1 = intersection("J1", [I1], [I2, I3], 0.5)
J2 = intersection("J2", [I2], [I4, I5], 0.4)
J3 = intersection("J3", [I3, I4], [I6], 0.6)

# Starting conditions.
I1.startingCondition(lambda x: 0.75 * (x >= 1.5))

# Roads and intersections lists.
roads = [I1, I2, I3, I4, I5, I6]
intersections = [J1, J2, J3]

# Numerical solution.
start = timer()
solution, time = networkSolver(roads, intersections, [0, 5])
stop = timer() - start


# PLOTTING


fig, ax = plt.subplots(2, 3)

plots = [None] * 6

h = 0
k = 0

for j in range(len(roads)):
    (plots[j],) = ax[k, h].plot(roads[j].space, solution[0, :, j])

    ax[k, h].set_title("Road {}".format(j + 1))

    ax[k, h].set_xticks(roads[j].span)

    ax[k, h].set_xlim(roads[j].span[0] - 0.1, roads[j].span[1] + 0.1)
    ax[k, h].set_ylim(-0.1, 1.1)

    k += 1

    if k == 2:
        k = 0
        h += 1

# Final output.
print("\nParameters.\n")

print("\tTime steps: {}".format(tSteps))
print("\tSpace steps: {}".format(sSteps))

print("\nRoads.\n")

for rd in roads:
    print("\t{}: {}".format(rd, rd.span))

print("\nIntersections.\n")

for inter in intersections:
    print("\t{}".format(inter))

print("\nTime taken: {:.2f}s".format((stop)))

ani = animation.FuncAnimation(
    fig, animate, interval=1000 // 60, frames=range(1, tSteps, tSteps // 300)
)

# Saves the animation.
try:
    if "save" in sys.argv:
        ani.save("858798_4_video.gif", fps=60)
except:
    pass

# Shows the animation.
plt.show()
