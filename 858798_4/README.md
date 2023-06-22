# Metodi di Analisi Applicata, 858798_4

Simulazione di una rete stradale tramite metodo a volume finito.

Rete stradale composta dalle seguenti strade:

- **I1**
- **I2**
- **I3**
- **I4**
- **I5**
- **I6**

e dai seguenti incroci:

- **J1**: incrocio *1x2*
	- Strada entrante: **I1**
	- Strade uscenti: **I2**, **I3**
- **J2**: incrocio *1x2*
	- Strada entrante: **I2**
	- Strade uscenti: **I4**, **I5**
- **J3**: incrocio *2x1*
	- Strade entranti: **I3**, **I4**
	- Strada uscente: **I6**

## Algoritmo

Il codice usa il flusso numerico di Godunov per la valutazione della soluzione numerica su varie strade connesse da incroci.  
Le condizioni al bordo sono valutate a ogni ciclo temporale e aggiornate.

## Output di esempio

Test eseguito su CPU Apple M2.

### Output testuale

[Output testuale](./858798_4_output.txt) ridotto con impostazioni di esempio.

Riproducibile tramite `python 858798_4.py long`

```
Andrea Di Antonio, 858798.

Numerical solution [4000/100000 - 17.84s]
Numerical solution [8000/100000 - 17.81s]
...
Numerical solution [92000/100000 - 17.89s]
Numerical solution [96000/100000 - 17.87s]

Parameters.

	Time steps: 100000
	Space steps: 10000

Roads.

	I1: [0, 2]
	I2: [2, 3]
	I3: [2, 4]
	I4: [3, 4]
	I5: [3, 6]
	I6: [4, 6]

Intersections.

	J1: [I1] -> [I2, I3] [0.50]
	J2: [I2] -> [I4, I5] [0.40]
	J3: [I3, I4] -> [I6] [0.60]

Time taken: 448.06s

```

### Output grafico

[Output grafico](./858798_4_gif.gif) con impostazioni di esempio.

![Test result](./858798_4_gif.gif)