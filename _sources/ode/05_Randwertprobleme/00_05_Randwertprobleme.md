(numerische-lösung-von-randwertproblemen)=
Numerische Lösung von Randwertproblemen
===

In diesem Kapitel der Vorlesung wollen wir uns mit der Lösung von
Randwertproblemen für **lineare gewöhnliche Differentialgleichungen
zweiter Ordnung** beschäftigen. Insbesondere beschäftigen wir uns mit
Lösungen $u \in C^2([0,1])$ für Differentialgleichungen der Gestalt
```{math}
:label: eq:dgl_zweite_ordnung
 -u''(x) - p(x) \cdot u'(x) + q(x) \cdot u(x) \ = \ g(x),
```
für $x\in (0,1)$ mit den zusätzlichen Randbedingungen
```{math}
:label: eq:randwert_rb1
\alpha_0 \cdot u'(0) + \beta_0 \cdot u(0) \ = \ g_0,

```
```{math}
:label: eq:randwert_rb2
\alpha_1 \cdot u'(1) + \beta_1 \cdot u(1) \ = \ g_1. 

```
Im Fall $\alpha_i = 0$ und $\beta_i = 1$ für $i=0,1$ spricht man von
**Dirichlet-Randbedingungen** für die gilt:
```{math}
u(0) \, = \, g_0, \qquad u(1) \, = \, g_1.
```
Im Fall $\alpha_i=1$ und $\beta_i=0$ für $i=0,1$ hingegen sprechen wir
von **Neumann-Randbedingungen** für die gilt:
```{math}
u'(0) \, = \, g_0, \qquad u'(1) \, = \, g_1.
```
In allen anderen Fällen erhalten wir gemischte Randwertbedingungen, die
auch *Robin-Randbedingungen* genannt werden.

Solche linearen Randwertprobleme treten beispielsweise auf, wenn man die
Auslenkung eines an zwei Seiten festgebundenen Seils beschreiben will.
Hierbei wird klar, dass die Auslenkung des Seils am linken und rechten
Rand verschwindet muss, d.h. wir haben Dirichlet-Randbedingungen
vorliegen mit $u(0) = u(1) = 0$. Ein ähnliches Problem untersucht man
bei der Berechnung der Temperaturverteilung eines Metallstabs, der an
beiden Enden erhitzt wird.

Nehmen wir an $P(x)$ sei eine Stammfunktion der Funktion $p(x)$, so dass
$P'(x) = p(x)$ gilt. Dann können wir beide Seiten der Gleichung
{eq}`eq:dgl_zweite_ordnung` mit dem Term $a(x) \coloneqq e^{P(x)}$
multiplizieren und erhalten somit:
```{math}
-a(x) \cdot u''(x) - e^{P(x)} \cdot p(x) \cdot u'(x) + e^{P(x)} \cdot q(x) \cdot u(x) \ = \ e^{P(x)} \cdot g(x).
```
Wenn wir nun die Produktregel für Differentiation anwenden sehen wir,
dass gilt
```{math}
-(a(x) \cdot u'(x))' \ = \ -a(x) \cdot u''(x) - a'(x) \cdot u'(x) \ = \ -a(x) \cdot u''(x) - e^{P(x)} \cdot p(x) \cdot u'(x).
```
Zusammen mit den Hilfsfunktionen $c(x) \coloneqq e^{P(x)} \cdot q(x)$
und $f(x) \coloneqq e^{P(x)} \cdot g(x)$ können wir das Problem
{eq}`eq:dgl_zweite_ordnung` schließlich in eine sogenannte
**Divergenzform** bringen mit:
```{math}
:label: eq:divergenzform
-(a(x) \cdot u'(x))' + c(x) \cdot u(x) \ = \ f(x).
```
Die Divergenzform {eq}`eq:divergenzform` der ursprünglichen
Differentialgleichung hat einige Vorteile für die Analyse der
Eigenschaften der Differentialgleichung und ebenfalls bei ihrer
numerischen Lösung, so dass wir diese im Folgenden weiter verwenden
werden.

