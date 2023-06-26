# Metodi di Analisi Applicata, 858798_4

Simulazione di una semplice rete stradale tramite metodo a volume finito.

Rete stradale composta dalle seguenti strade:

- **R1**[^1]
- **R2**
- **R3**
- **R4**
- **R5**
- **R6**

[^1]: Cambio di nomenclatura rispetto alla consegna originale.

e dai seguenti incroci:

- **J1**: incrocio *1x2*
	- Strada entrante: **R1**
	- Strade uscenti: **R2**, **R3**
- **J2**: incrocio *1x2*
	- Strada entrante: **R2**
	- Strade uscenti: **R4**, **R5**
- **J3**: incrocio *2x1*
	- Strade entranti: **R3**, **R4**
	- Strada uscente: **R6**

## Algoritmo

Il codice usa il flusso numerico di Godunov per la valutazione della soluzione numerica su varie strade connesse da incroci.  
A ogni ciclo temporale l'algoritmo aggiorna le condizioni al bordo per i vari incroci presenti e estende le soluzioni numeriche su ogni strada per simulare le intersezioni.

Le classi che gestiscono gli incroci richiedono che le strade entranti terminino alla stessa cordinata spaziale e lo stesso vale per le strade uscenti.
Inoltre, la coordinata spaziale terminale delle strade entranti deve coincidere con la coordinata spaziale iniziale delle strade uscenti.  
Questo non è cruciale per il funzionamento dell'algoritmo, ma rappresenta una mera scelta stilistica.

## Output di esempio

Test eseguito su CPU Apple M2.

### Output testuale

[Output testuale](./858798_4_output.txt) (ridotto) con impostazioni di esempio.

Riproducibile tramite `python 858798_4.py long`[^2]

[^2]: L'opzione `long` aumenta notevolmente il tempo di esecuzione.

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

### Output animato

[Output animato](./858798_4_gif.gif) con impostazioni di esempio[^3].

[^3]: Quando salvata, la gif non presenta i *ticks* sull'asse y.

![Test result](./858798_4_gif.gif)