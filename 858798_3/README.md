# Metodi di Analisi Applicata, 858798_3

Confronto della soluzione numerica per due problemi di conservazione tramite flussi
numerici di Lax-Friedrichs (LF) e Godunov (G).

I test eseguiti da questo codice sono divisi in due sezioni.

### Sezione 1.

Confronto numerico dei flussi numerici di Godunov e Lax-Friedrichs.

I test della sezione 1 consistono nella valutazione della soluzione numerica,
usando flussi numerici di Lax-Friedrichs e Godunov, per diversi valori di 
∆t e ∆x e la valutazione della distanza, in norma L1, con la soluzione analitica.

### Sezione 2.

Valutazione della variazione totale delle soluzioni.

I test della sezione 2 consistono nella valutazione, fissati ∆t e ∆x, della
variazione totale delle tre soluzioni (LF, G, Analitica) a tempi diversi con
annesso confronto grafico.

## Conclusioni

Dai test si evince, sia numericamente che qualitativamente, la superiorità del flusso numerico di Godunov.

## Test di esempio

Test eseguito su CPU Apple M2.

### Output testuale

[Output testuale](./858798_3_output.txt) del test di esempio.

```
Andrea Di Antonio, 858798.

Time taken for section 1's 10 test(s): 47.59 second(s).

Section 1, first test results:

Space	Time	Lax-Friedrichs	Godunov

128		512		1.6025e+00	3.6992e-01
203		812		1.4986e+00	2.8479e-01
322		1290	1.2299e+00	2.2247e-01
512		2048	1.1149e+00	1.3857e-01
812		3250	9.7566e-01	9.6350e-02
1290	5160	8.2158e-01	5.8368e-02
2048	8192	6.6999e-01	4.5120e-02
3250	13003	5.1653e-01	2.9056e-02
5160	20642	3.8527e-01	2.2176e-02
8192	32768	2.7427e-01	1.3126e-02

Section 1, second test results:

Space	Time	Lax-Friedrichs	Godunov

128		512		3.5469e+00	4.8367e-01
203		812		2.8775e+00	4.2263e-01
322		1290	2.2628e+00	2.8174e-01
512		2048	1.7877e+00	1.9203e-01
812		3250	1.4290e+00	1.3953e-01
1290	5160	1.1431e+00	1.0010e-01
2048	8192	8.6981e-01	6.7191e-02
3250	13003	6.3918e-01	4.8503e-02
5160	20642	4.4986e-01	3.2452e-02
8192	32768	3.1037e-01	2.2001e-02

Time taken for section 2's test: 31.54 second(s)

Section 2, some first test results:

Time	Analytical 1	Lax-Friedrichs	Godunov

Space steps: 8192
Time steps: 32768

0.0003	2.0000e+00	2.0000e+00	2.0000e+00
1.1112	1.8927e+00	1.4590e+00	1.8703e+00
2.2221	1.3393e+00	1.0690e+00	1.3437e+00
3.3329	1.0943e+00	8.7957e-01	1.0978e+00
4.4438	9.4849e-01	7.6409e-01	9.4034e-01
5.5547	8.4780e-01	6.8454e-01	8.4017e-01
6.6655	7.7426e-01	6.2551e-01	7.6638e-01
7.7764	7.1702e-01	5.7950e-01	7.0927e-01
8.8873	6.7065e-01	5.4233e-01	6.6322e-01
9.9982	6.3215e-01	5.1150e-01	6.2511e-01

Section 2, some second test results:

Time	Analytical 2	Lax-Friedrichs	Godunov

0.0003	5.0000e-01	1.5000e+00	1.5000e+00
1.1112	1.4961e+00	1.2706e+00	1.5000e+00
2.2221	1.4486e+00	1.0768e+00	1.4169e+00
3.3329	1.2735e+00	9.7307e-01	1.2625e+00
4.4438	1.1708e+00	9.0797e-01	1.1630e+00
5.5547	1.0993e+00	8.6264e-01	1.0944e+00
6.6655	1.0470e+00	8.2889e-01	1.0432e+00
7.7764	1.0066e+00	8.0255e-01	1.0033e+00
8.8873	9.7415e-01	7.8130e-01	9.7112e-01
9.9982	9.4711e-01	7.6369e-01	9.4440e-01

Final results have been plotted.

Time taken: 79.25 second(s).

```

### Output grafico

[Output grafico](./858798_3_image.png) del test di esempio.

![Test result](./858798_3_image.png)