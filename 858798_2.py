"""
Metodi di Analisi Applicata.
Andrea Di Antonio, 858798.
Assignment 2.
Version 2.

Simulation of the second "Follow The Leader" model in the presence of an arbitrary number of traffic lights.

The solution of the system of ODEs is done using the RK4 method.
The system is solved as a whole rather than acting on each pair of equations individually.

Before the simulation, the road segments delimited by the traffic lights are evaluated, for which the leaders will then be identified.

The dynamics of the cars are initialized as follows: the position-velocity vector is constructed as [x1, ..., xN, v1, ..., vN], in order to have a vector field equal to [v1, ..., vN, V'(x(j+1) - xj)(v(j+1) - vj) + (V(x(j+1) - xj) - vj)/t, ..., (Vinf - vN)/t], where "t" represents the reaction time, "Vinf" the limit speed, "V()" the optimal speed depending on the distance between cars and "V'()" its derivative, equivalently to assignment 1.
In the absence of traffic lights or in the presence of a green light, this dynamics is preserved.

In the presence of a yellow or red light, the dynamics of the leaders are updated based on the four cases depending on their distance and position*:
i**) x(t_g) + v(t_g)(t_y) >= l: The leader crosses the light before it turns red without changing its dynamics.
ii) x(t_g) + v(t_g)(t_y) <= l, x(t_g) + v(t_g)(t_y t_r) >= l: The leader, proceeding at constant speed, manages to cross the light after it turns green without slowing down; in this case, the speed is kept constant.
iii***) x(t_g) + v(t_g)(t_y) <= l, x(t_g) + v(t_g)(t_y t_r) <= l, x(t_g) + v(t_g)(t_y t_r)/2 <= l: The leader is forced to slow down and stop before the light turns green; in this case, the dynamics correspond to: v' = -v(t_g)^2 / (l - x(t_g)) assuming t_g <= t <= t_g + 2(l - x(t_g)) / v(t_g), otherwise v' = 0.
iv***) x(t_g) + v(t_g)(t_y) <= l, x(t_g) + v(t_g)(t_y t_r) <= l, x(t_g) + v(t_g)(t_y t_r)/2 >= l: The leader is forced to slow down but not to stop in order to pass the light with a generally non-zero velocity when it turns green; in this case, the dynamics correspond to: v' = -2(x(t_g) + v(t_g)(t_y + t_r) - l) / (t_y + t_r)^2.

*: Position and velocity of the cars at t_g are determined by interpolating the trajectories for greater precision.
**: In this case, the traffic light dynamics apply to the first vehicle in a position lower than that of the leader such that condition 2 to 4 is met.  
***: To avoid problems due to the discretization of the model by computer calculation, a correction is introduced so that the cars do not pass the traffic light. Through the "physics" flag, this correction is equal to the maximum between the standard correction, of order equal to the step of the algorithm used (RK4), and half the length of the cars, as in the real case.

At the end of the simulation, the trajectories of the car positions, together with the positions of the traffic lights, and the velocity over time are displayed.

The default values set (carsNumber, carsLength, infSpeed, reactionTime, endingTime, lightsNumber, lightsDistances, green, yellow, red) represent only an example and can be modified.
"""

from time import time as timer
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt


# Solver: RK4.
def rk4(vField, startingCondition: np.array, tSpan: list, steps=4096) -> np.array:
    trajectory = np.array([np.zeros_like(startingCondition) for _ in range(steps)])
    time = np.linspace(tSpan[0], tSpan[1], steps)

    global step
    step = time[1] - time[0]

    trajectory[0] = startingCondition

    for t in range(1, steps):
        k1 = vField(time[t - 1], trajectory[t - 1], time[0:t], trajectory.T[:, 0:t])
        k2 = vField(
            time[t - 1] + step / 2,
            trajectory[t - 1] + step * k1 / 2,
            time[0:t],
            trajectory.T[:, 0:t],
        )
        k3 = vField(
            time[t - 1] + step / 2,
            trajectory[t - 1] + step * k2 / 2,
            time[0:t],
            trajectory.T[:, 0:t],
        )
        k4 = vField(
            time[t - 1] + step,
            trajectory[t - 1] + step * k3,
            time[0:t],
            trajectory.T[:, 0:t],
        )

        trajectory[t] = trajectory[t - 1] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return time, trajectory


print(
    "Follow the leader, second model simulation with traffic lights.\nAndrea Di Antonio, 858798."
)

# Flags
physiscs = True
spaceTime = False  # Changes the way the final plot gets presented.

# Variables.
steps = 2**13

carsNumber = 6
carsLength = 1
infSpeed = 0.75
reactionTime = 0.05
endingTime = 300

# Traffic lights parameters. Lights share the same distance between them.
lightsNumber = 6
lightsDistances = 20
green = 25
yellow = 10
red = 20

# Assertions
assert steps > 0

assert lightsNumber > 0
assert lightsDistances > 0
assert green > 0
assert yellow > 0
assert red > 0

assert carsNumber > 1
assert carsLength >= 0
assert infSpeed > 0
assert reactionTime > 0
assert endingTime > green

print(
    "\nCurrent configuration:\n\nNumber of cars: {}\nLength: {}\nReaction time: {}\nReaction time: {}\nSpeed limit: {}\nTime span: {}".format(
        carsNumber, carsLength, lightsNumber, reactionTime, infSpeed, endingTime
    )
)
print(
    "\nNumber of lights: {}\nDistances: {}\nGreen duration: {}\nYellow duration: {}\nRed duration: {}".format(
        lightsNumber, lightsDistances, green, yellow, red
    )
)

# Traffic lights positions.
lights = [index * lightsDistances for index in range(1, lightsNumber + 1)]

# Segments of road ending with traffic lights.
segments = [[0, lights[0]]] + [
    [lights[index], lights[index + 1]] for index in range(lightsNumber - 1)
]


# Optimal speed.
def optimal(distance: float) -> float:  # Works with numpy.array.
    # It returns 0 should the distance between cars be lower than their lengths.
    return (
        infSpeed
        * (2 / np.pi)
        * np.arctan(distance - carsLength)
        * (distance - carsLength >= 0)
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


# Traffic lights status.
def lightStatus(time: float) -> int:
    # 0: Green.
    # 1: Yellow.
    # 2: Red.

    cycle = green + yellow + red

    if time % cycle < green:  # Green light.
        return 0

    elif time % cycle <= green + yellow:  # Yellow light.
        return 1

    else:
        return 2  # Red light.


# Last time the lights went yellow.
def lastYellow(time: float) -> float:
    cycle = green + yellow + red
    return int((time - green) / cycle) * cycle + green


# Leader index.
def leaderIndex(
    cars: np.array, span: list
) -> int:  # Returns the leader of a specific segment.
    try:
        return list(cars).index(max(car for car in cars if span[0] <= car <= span[1]))

    except:
        return -1


# Velocity Field for the ODE system.
def velocityField(
    time: float, point: np.array, times: np.array, trajectories: np.array
) -> np.array:
    # point: [X_1...N, V_1...N]
    leaders = [
        leaderIndex(point[0 : int(point.size / 2)], span) for span in segments
    ]  # Leaders for every segment.

    distances = np.array(
        [point[index + 1] - point[index] for index in range(int(point.size / 2) - 1)]
    )  # Evaluates the distances between the cars.

    vField = np.zeros_like(point)
    vField[0 : int(point.size / 2)] = point[
        int(point.size / 2) : point.size
    ]  # [x1', ..., xN']

    # Initializes cars as if there were no traffic lights.
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

    # Updates the vector field based on traffic lights.
    if lightStatus(time) in [1, 2]:
        vField = lightDynamics(time, point, times, trajectories, leaders, vField)

    return vField


# Traffic lights dynamics.
def lightDynamics(
    time: float,
    point: np.array,
    times: np.array,
    trajectories: np.array,
    leaders: list,
    vField: np.array,
) -> np.array:
    for leader in leaders:
        if leader == -1:  # Segment has no leader.
            continue

        # Leader's speed index.
        speed = leader + int(point.size / 2)

        # Last time the traffic lights went yellow.
        y = lastYellow(time)

        # Evaluates position and speed at the yellow time with scipy.interpolate.interp1d.
        yellowPosition = interp1d(
            times, trajectories[leader], fill_value="extrapolate"
        )(y)
        yellowSpeed = interp1d(times, trajectories[speed], fill_value="extrapolate")(y)

        light = min([light for light in lights if light >= point[leader]])

        # Fixes roundoff errors by letting cars stop sooner either by half of their lengths or by a small epsilon given by the step of the solver.
        epsilon = max([step, carsLength / 2]) if physiscs else step

        # First case.
        if yellowPosition + yellowSpeed * yellow > light:
            vField = lightDynamics(
                time, point, times, trajectories, [leader - 1], vField
            )  # Update the vField for the following car(s).

        # Second case.
        elif (
            yellowPosition + yellowSpeed * yellow <= light
            and yellowPosition + yellowSpeed * (yellow + red) < light
        ):
            vField[speed] = 0.0

        # Third case.
        elif (
            yellowPosition + yellowSpeed * yellow <= light
            and yellowPosition + yellowSpeed * (yellow + red) >= light
            and yellowPosition + yellowSpeed * (yellow + red) / 2 >= light
        ):
            vField[speed] = (
                -(yellowSpeed**2) / (light - epsilon - yellowPosition) / 2
                if time <= y + 2 * (light - epsilon - yellowPosition) / yellowSpeed
                else 0.0
            )

        # Fourth case.
        elif (
            yellowPosition + yellowSpeed * yellow <= light
            and yellowPosition + yellowSpeed * (yellow + red) >= light
            and yellowPosition + yellowSpeed * (yellow + red) / 2 < light
        ):
            vField[speed] = (
                -2
                * (yellowPosition + yellowSpeed * (yellow + red) - light + epsilon)
                / ((yellow + red) ** 2)
            )

    return vField


# Starting conditions.
carsArray = np.array([index + carsLength for index in range(carsNumber)])
speedsArray = np.array([0.0] * carsNumber)  # Stopped cars.

startingCondition = np.array([carsArray, speedsArray])
startingCondition = startingCondition.reshape(
    startingCondition.size,
)  # This simulation works with [x1, ..., xN, v1, ..., vN] vectors.

# Solutions.
start = timer()
time, trajectories = rk4(
    velocityField, startingCondition, [0, endingTime], steps
)  # ODE solution, RK4.
stop = timer() - start
trajectories = trajectories.T

# Plotting.
positionsTrajectories = trajectories[0 : int(len(trajectories) / 2)]
speedsTrajectories = trajectories[int(len(trajectories) / 2) : int(len(trajectories))]

fig, ax = plt.subplots(1, 2)

ax[0].set_title("Positions and traffic lights")
ax[1].set_title("Speeds")

# Traffic lights.
colours = ["green", "yellow", "red"]
statuses = np.array([lightStatus(t) for t in time])

splits = [
    index + 1
    for index in range(statuses.size - 1)
    if statuses[index] != statuses[index + 1]
]

statuses = np.split(statuses, splits)
times = np.split(time, splits)

if spaceTime:
    carIndex = 1
    for trajectory in positionsTrajectories:
        label = label = (
            "Car {}".format(carIndex)
            if carIndex != len(positionsTrajectories)
            else "Leader"
        )
        ax[0].plot(time, trajectory, label=label, linewidth=1)
        carIndex += 1

    carIndex = 1
    for trajectory in speedsTrajectories:
        label = label = (
            "Car {}".format(carIndex)
            if carIndex != len(speedsTrajectories)
            else "Leader"
        )
        ax[1].plot(time, trajectory, label=label, linewidth=0.25)
        carIndex += 1

    for light in lights:
        for index in range(len(times)):
            ax[0].plot(
                times[index],
                light * np.ones_like(times[index]),
                color=colours[statuses[index][0]],
                linewidth=1,
            )

else:
    carIndex = 1
    for trajectory in positionsTrajectories:
        label = label = (
            "Car {}".format(carIndex)
            if carIndex != len(positionsTrajectories)
            else "Leader"
        )
        ax[0].plot(trajectory, time, label=label, linewidth=1)
        carIndex += 1

    carIndex = 1
    for trajectory in speedsTrajectories:
        label = label = (
            "Car {}".format(carIndex)
            if carIndex != len(speedsTrajectories)
            else "Leader"
        )
        ax[1].plot(trajectory, time, label=label, linewidth=0.25)
        carIndex += 1

    for light in lights:
        for index in range(len(times)):
            ax[0].plot(
                light * np.ones_like(times[index]),
                times[index],
                color=colours[statuses[index][0]],
                linewidth=1,
            )

if carsNumber <= 8:  # Avoids showing the legend on a high number of cars.
    ax[0].legend(loc="best")
    ax[1].legend(loc="best")

# Timing.
print("\nDone in {} second(s) with {} steps".format(round(stop, 2), steps))

# Shows the final plot.
plt.show()
