"""
Metodi di Analisi Applicata.
Andrea Di Antonio, 858798.
Assignment 1.

Simulation of the second "Follow The Leader" model.

The solution of the system of ODEs is done using the Euler method.
The system is solved as a whole rather than acting on each pair of equations individually.
The position-velocity vector is constructed as [x1, ..., xN, v1, ..., vN], in order to have a vector field equal to [v1, ..., vN, V'(x(j+1) - xj)(v(j+1) - vj) + (V(x(j+1) - xj) - vj)/t, ..., (Vinf - vN)/t], where "t" represents the reaction time, "Vinf" the limit speed, "V()" the optimal speed depending on the distance between cars and "V'()" its derivative.
At the end of the simulation, the trajectories of the car positions and velocities over time are displayed.

The default values set (carsNumber, carsLength, infSpeed, reactionTime, endingTime) represent only an example and can be modified.
"""

import numpy as np
import matplotlib.pyplot as plt


# Solver: Euler.
def euler(vField, startingCondition: np.array, tSpan: list, steps=4096) -> np.array:
    trajectory = np.array([np.zeros_like(startingCondition) for _ in range(steps)])
    time = np.linspace(tSpan[0], tSpan[1], steps)
    step = time[1] - time[0]

    trajectory[0] = startingCondition

    for t in range(steps - 1):
        trajectory[t + 1] = trajectory[t] + vField(time[t], trajectory[t]) * step

    return time, trajectory


print("Follow the leader, second model simulation.\nAndrea Di Antonio, 858798.")

# Variables.
carsNumber = 4
carsLength = 0.5
infSpeed = 10
reactionTime = 0.025
endingTime = 5

# Assertions
assert carsNumber > 1
assert carsLength >= 0
assert infSpeed > 0
assert reactionTime > 0
assert endingTime > 0

print(
    "\nCurrent configuration:\n\nNumber of cars: {}\nLength: {}\nReaction time: {}\nSpeed limit: {}\nTime span: {}".format(
        carsNumber, carsLength, reactionTime, infSpeed, endingTime
    )
)


# Optimal speed.
def optimal(distance: float) -> float:  # Works with numpy.array.
    # It returns 0 should the distance between cars be lower or equal than their lengths.
    return (
        infSpeed
        * np.arctan(distance - carsLength)
        * 2
        / np.pi
        * (distance - carsLength > 0)
    )


# Optimal speed's first derivative.
def optimalPrime(distance: float) -> float:  # Works with numpy.array.
    # It returns 0 should the distance between cars be lower than their lengths.
    return (
        infSpeed
        * (2 / np.pi)
        / (1 + (distance - carsLength) ** 2)
        * (distance - carsLength >= 0)
    )


# Velocity Field for the ODE system.
def velocityField(time: float, point: np.array) -> np.array:
    # point: [X_1...N, V_1...N]

    distances = np.array(
        [point[index + 1] - point[index] for index in range(int(point.size / 2) - 1)]
    )  # Evaluates the distances between the cars.

    vField = np.zeros_like(point)
    vField[0 : int(point.size / 2)] = point[
        int(point.size / 2) : point.size
    ]  # [x1', ..., xN']

    vField[int(point.size / 2) : point.size - 1] = (
        optimalPrime(distances)
        * (
            point[int(point.size / 2) + 1 : point.size]
            - point[int(point.size / 2) : point.size - 1]
        )
        + (optimal(distances) - point[int(point.size / 2) : point.size - 1])
        / reactionTime
    )  # [v1', ..., v(N-1)']
    vField[point.size - 1] = (infSpeed - point[point.size - 1]) / reactionTime  # [vN']

    return vField


# Starting conditions.
carsArray = np.array([index + carsLength for index in range(carsNumber)])
speedsArray = np.array([0.0] * carsNumber)  # Stopped cars.

startingCondition = np.array([carsArray, speedsArray])
startingCondition = startingCondition.reshape(
    startingCondition.size,
)  # This simulation works with [x1, ..., xN, v1, ..., vN] vectors.

# Solutions.
time, trajectories = euler(
    velocityField, startingCondition, [0, endingTime]
)  # ODE solution.
trajectories = trajectories.T

# Plotting.
positionsTrajectories = trajectories[0 : int(len(trajectories) / 2)]
speedsTrajectories = trajectories[int(len(trajectories) / 2) : int(len(trajectories))]

fig, ax = plt.subplots(2)

ax[0].set_title("Cars positions")
ax[1].set_title("Cars speeds")

carIndex = 1
for trajectory in positionsTrajectories:
    label = label = (
        "Car {}".format(carIndex)
        if carIndex != len(positionsTrajectories)
        else "Leader"
    )
    ax[0].plot(time, trajectory, label=label)
    carIndex += 1

carIndex = 1
for trajectory in speedsTrajectories:
    label = label = (
        "Car {}".format(carIndex) if carIndex != len(speedsTrajectories) else "Leader"
    )
    ax[1].plot(time, trajectory, label=label)
    carIndex += 1

if carsNumber <= 10:  # Avoids showing the legend on a high number of cars.
    ax[0].legend(loc="best")
    ax[1].legend(loc="best")

plt.show()
