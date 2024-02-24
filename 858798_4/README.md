# Metodi di Analisi Applicata, 858798_4

Simulation of a simple road network using the finite volume method.

Road network composed of the following roads:

- **R1**[^1]
- **R2**
- **R3**
- **R4**
- **R5**
- **R6**

[^1]: Renamed compared to the original assignment.

and the following intersections:

- **J1**: intersection *1x2*
    - Incoming road: **R1**
    - Outgoing roads: **R2**, **R3**
- **J2**: intersection *1x2*
    - Incoming road: **R2**
    - Outgoing roads: **R4**, **R5**
- **J3**: intersection *2x1*
    - Incoming roads: **R3**, **R4**
    - Outgoing road: **R6**

## Algorithm

The code uses the Godunov numerical flux for evaluating the numerical solution on various roads connected by intersections.
At each time step, the algorithm updates the boundary conditions for the various intersections present and extends the numerical solutions on each road to simulate the intersections.

The classes handling the intersections require that the incoming roads terminate at the same spatial coordinate, and the same goes for the outgoing roads.
Additionally, the terminal spatial coordinate of the incoming roads must coincide with the initial spatial coordinate of the outgoing roads.
This is not crucial for the functioning of the algorithm but represents a mere stylistic choice.

## Example Output

Test executed on Apple M2 CPU.

### Textual Output

[Textual Output](./858798_4_output.txt) (excerpt) with example settings.

Reproducible via `python 858798_4.py long`[^2]

[^2]: The `long` option significantly increases execution time.

```
Andrea Di Antonio, 858798.

Numerical solution [4000/100000 - 17.97s]
Numerical solution [8000/100000 - 17.98s]
...
Numerical solution [92000/100000 - 17.93s]
Numerical solution [96000/100000 - 17.92s]

Parameters.

    Time steps: 100000
    Space steps: 10000

Roads.

    R1: [0, 2]
    R2: [2, 3]
    R3: [2, 4]
    R4: [3, 4]
    R5: [3, 6]
    R6: [4, 6]

Intersections.

    J1: [R1] -> [R2, R3] [0.50]
    J2: [R2] -> [R4, R5] [0.40]
    J3: [R3, R4] -> [R6] [0.60]

Time taken: 448.79s
```

### Animated Output

[Animated Output](./858798_4_gif.gif) with example settings[^3].

[^3]: When saved, the gif does not display y-axis ticks.

![Test result](./858798_4_gif.gif)