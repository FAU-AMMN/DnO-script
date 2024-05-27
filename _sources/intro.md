(einleitung)=
Einleitung
===

In dieser Vorlesung werden wir einige weiterführende Aspekte der
numerischen Mathematik diskutieren, nämlich numerische Verfahren zur
Lösung von Optimierungsproblemen und von (gewöhnlichen)
Differentialgleichungen.

Im ersten Teil der Vorlesung beschäftigen wir uns mit
**Optimierungsproblemen**, welche in vielen mathematischen
Anwendungsbereichen auftreten, von klassischen ökonomischen Problemen
über Materialoptimierung bis hin zu modernen Problemen in der
mathematischen Bildverarbeitung und im maschinellen Lernen. Hierbei
konzentrieren wir uns auf hauptsächlich auf *unbeschränkte
Minimierungsprobleme* der Form
```{math}
\min_{x \in \Omega} F(x),
```
wobei $F$ eine gegebene, zu minimierende Funktion ist und
$\Omega \subset \R^n$ eine geeignete Menge von möglichen
Eingabevektoren. Methodisch knüpfen wir im Teil zur Optimierung an die
*iterativen Methoden zur Lösung von Gleichungssystemen* an, allerdings
kommen hier noch einige Aspekte dazu: Mit einem Optimierungsproblem im
Hintergrund können wir die Iterationsverfahren geeignet anpassen um
tatsächlich mit jeder Iteration die Werte einer gegebenen Funktion zu
verkleinern. Darüber hinaus werden wir geeignete *Wahlen der
Schrittweite* der Iterationsverfahren kennenzulernen, um die Konvergenz
gegen Minimierer oder zumindest stationäre Punkte gewährleisten zu
können. Ein weiterer Aspekt ist die *Optimierung unter Nebenbedingung*
und die *Minimierung konvexer nicht-differenzierbarer Probleme*, die in
vielen modernen Anwendungen auftreten. Dazu werden wir exemplarisch ein
spezielles Verfahren kennenlernen.  \
 \
Der zweite Teil der Vorlesung beschäftigt sich mit **Lösungen von
gewöhnlichen Differentialgleichungen**. Wir beginnen mit *einfachen
Anfangswertproblemen* der Form
```{math}
u'(t) \ = \ F(u(t),t), \qquad u(0) \, = \, u_0,
```
für eine unbekannte Funktion $u \colon [0, T] \rightarrow \R^m$ und die
nur für spezielle Formen von
$F \colon \R^{m} \times [0,T] \rightarrow \R^m$ gelöst werden können. In
einer Vielzahl von Anwendungen treten jedoch *allgemeinere Funktionen*
$F$ auf, für die eine analytische Lösung nicht mehr möglich ist, wie zum
Beispiel bei den Newtonschen Gesetzen für die Dynamik von mehr als
$N > 2$ Teilchen (siehe das *$N$-Körper-Problem* {cite:p}`N_koerper_problem`.
Die entstehenden Gleichungssysteme können dann auch beliebig groß
werden, z.B. in der Molekulardynamik, wo $u(t)$ die räumlichen
Koordinaten $k\in\N$ verschiedener Teilchen im Zeitverlauf beschreibt
und man somit Gleichungssysteme der Größe $n=3k$ erhält. Andere
klassische Anwendungsgebiete gewöhnlicher Differentialgleichungen sind
die Modellierung von Populationsdynamiken oder auch von Aktienmärkten,
wo meist noch eine zufällige Komponente hinzugefügt wird und man somit
*stochastische Differentialgleichungen* erhält. Für solche Anwendungen
lassen sich im Allgemeinen nur Lösungen mittels numerischer Verfahren
approximieren, die wir in dieser Vorlesung herleiten wollen.

Die numerischen Verfahren zur Lösung von Anfangswertproblemen lassen
sich unterschiedlich gestalten. Sie sind einerseits ähnlich zu bereits
bekannten Iterationsverfahren, nämlich dann wenn die Ableitung durch
*Differenzenquotienten* auf einem Gitter approximiert wird. Andererseits
lässt sich ein Bezug zur numerischen Integration herstellen, wenn man
die äquivalente Formulierung als *Volterra-Integralgleichung*
```{math}
u(t) \ = \ u_0 + \int_0^t F(u(s),s)\, \mathrm{d}s
```
benutzt und anschließend numerische Quadraturformeln auf das Integral
anwendet. Ein wichtiger Aspekt ist in beiden Fällen die
*Diskretisierung*, d.h. wir approximieren das Problem in einem
endlich-dimensionalen Lösungsraum, z.B. durch Werte auf einem Gitter.
Mathematisch stellt sich dann natürlich die Frage ob und in welchem
Sinne das diskretisierte Problem gegen das ursprüngliche Problem
konvergiert.

Weiterhin werden wir auch *Randwertprobleme* betrachten, die in späteren
Vorlesungen zu partiellen Differentialgleichungen führen werden. Ein
einfaches Beispiel ist die numerische Lösung einer gewöhnlichen
Differentialgleichung der Form
```{math}
- (a(x) u'(x))' + c(x) u(x) \ = \ f(x),  \quad x \in (0,1),
```
mit vorgegebenen Randwerten $u(0)=u(1)=0$. Hier müssen wir die
Diskretisierung für das gesamte Intervall $(0,1)$ auf einmal durchführen
und nicht von einem Schritt zum Nächsten wie bei den oben beschriebenen
Anfangswertproblemen. Diese Diskretisierung liefert uns ein lineares
Gleichungssystem, das wir anschließend mit bekannten Methoden der
Numerik lösen müssen. Die Abschätzung des Diskretisierungsfehlers
erfordert weiterführende Methoden, welche wir im Laufe der Vorlesung
kennenlernen werden.
