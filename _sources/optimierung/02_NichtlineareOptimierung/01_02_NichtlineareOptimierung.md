(s:opt_grundlagen)=
# Mathematische Grundlagen

Im Folgenden wollen wir die mathematischen Grundlagen zur Untersuchung
von allgemeinen Optimierungsproblemen einführen. Wir folgen hierbei zu
großen Teilen der Notation von Nocedal und Wright in
{cite:p}`nocedal_1999`. Wir beginnen mit der Definition des allgemeinen
Optimierungsproblems, welches wir im Verlauf der Vorlesung noch weiter
konkretisieren werden.

````{prf:definition} Allgemeines Optimierungsproblem
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Funktion, welche wir *Zielfunktion* nennen. Unser Ziel ist es einen
unbekannten Vektor $x \in \Omega$, auch *Parametervektor* genannt, zu
finden, welcher das folgende **allgemeine Optimierungsproblem** löst:
```{math}
:label: eq:optimierungsproblem_allgemein
\min_{x \in \Omega} F(x) \qquad \text{ mit } \qquad
\begin{cases}
c_i(x) \ = \ 0, \quad &i \in \mathcal{E},\\
c_i(x) \ \geq \ 0, \quad &i \in \mathcal{I}.
\end{cases}
```
Die reellwertigen Funktionen
$c_i \colon \Omega \rightarrow \mathbb{R}, i=1,\ldots,m$ bilden einen
Vektor von *Nebenbedingungen*, welcher das Optimierungsproblem
restringiert. Die Indexmengen $\mathcal{E}$ und $\mathcal{I}$ legen
hierbei fest, ob es ich bei der jeweiligen Nebenbedingung um eine
Gleichung oder eine Ungleichung handelt.

````

Zur Veranschaulichung betrachten wir ein zweidimensionales Beispiel für
ein beschränktes, nichtlineares Optimierungsproblem.

````{prf:example} Beschränktes Optimierungsproblem
:label: ex:optimierungsproblem_restringiert
Wir betrachten das folgende mathematische Problem:
```{math}
\min_{x \in \mathbb{R}^2} (x_1 - 2)^2 + (x_2 - 1)^2
```
unter den Nebenbedingungen $x_1^2 - x_2 \leq 0$ und $x_1 + x_2 \leq 2$.
Wir können dieses Problem in die allgemeine Form des
Optimierungsproblems {eq}`eq:optimierungsproblem_allgemein` umschreiben
als:
```{math}
\min_{x \in \mathbb{R}^2} F(x) \ = \ \min_{x \in \mathbb{R}^2} (x_1 - 2)^2 + (x_2 - 1)^2, \quad \text{ mit } \quad
\begin{cases}
c_1(x) = -x_1^2 + x_2 \geq 0,\\
c_2(x) = -x_1 - x_2 + 2 \geq 0.
\end{cases}
```
Hierbei gilt für die Indexmengen $\mathcal{I} = \lbrace 1,2 \rbrace$ und
$\mathcal{E} = \emptyset$. Visualisiert man die Niveaulinien der
Zielfunktion $F$ zusammen mit den Nebenbedingungen, so erkennt man
direkt, dass der globale Minimierer der quadratischen Zielfunktion $F$,
nämlich $x = (x_1, x_2)^T = (2,1)^T$, nicht in der erlaubten Menge der
Parameter liegt, welche durch die Nebenbedingungen beschrieben ist.
Trotzdem existiert ein eindeutiges globales Minimum des beschränkten
Optimierungsproblems wie man in {numref}`fig:optimierung_restringiert`
sehen kann, nämlich $F(x^*) = 1$ für den Minimierer
$x^* = (x^*_1, x^*_2)^T = (1,1)^T$.

````

<figure id="fig:optimierung_restringiert">
<p>ToDo!</p>
<figcaption>Visualisierung der Niveaumengen der Zielfunktion <span
class="math inline"><em>F</em></span> und den Nebenbedingungen <span
class="math inline"><em>c</em><sub>1</sub></span> und <span
class="math inline"><em>c</em><sub>2</sub></span> des restringierten
Optimierungsproblems in .</figcaption>
</figure>

Im Rahmen dieser Vorlesung wollen wir uns zunächst auf eine bestimmte
Klasse von allgemeinen Optimierungsproblemen konzentrieren, den
*unbeschränkten* oder *unrestringierten* Optimierungsproblemen.

````{prf:definition} Unbeschränkte Optimierung
Liegt ein allgemeines Optimierungsproblem der Form
{eq}`eq:optimierungsproblem_allgemein` ohne Nebenbedingungen vor, d.h.,
für die Indexmengen gilt $\mathcal{E} = \mathcal{I} = \emptyset$, so
sprechen wir von einem **unbeschränkten** oder **unrestringierten
Optimierungsproblem**.

````

````{prf:remark} Relaxation
Häufig lassen sich restringierte Optimierungsprobleme in unrestringierte
Optimierungsprobleme überführen, indem man zusätzliche Strafterme zur
Zielfunktion hinzufügt, die eine Verletzung der ursprünglichen
Nebenbedingungen zwar mit Kosten belegt, diese jedoch grundsätzlich
erlaubt. Hierbei spricht man auch von *relaxierten
Optimierungsproblemen*.

Hat man beispielsweise das folgende restringierte Optimierungsproblem
vorliegen
```{math}
\min_{x \in \R} \left\lbrace F(x) \coloneqq e^x \right \rbrace \qquad \text{ mit } \qquad c(x) \coloneqq x \, \geq \, 0,
```
so lässt sich stattdessen auch folgendes relaxiertes Optimierungsproblem
ohne Nebenbedingungen betrachten:
```{math}
\min_{x \in \R} \left\lbrace G(x) \coloneqq e^x + \lambda \cdot (\min(0,x))^2 \right \rbrace,
```
wobei $\lambda \in \R^+$ den Strafterm zur Einhaltung der Nebenbedingung
gewichtet.

````

Daher gehen wir für den weiteren Verlauf der Vorlesung immer (wenn nicht
anders beschrieben) von einem unrestringierten Optimierungsproblem aus.

Neben der Unterscheidung von Optimierungsproblemen in beschränkte und
unbeschränkte Formulierungen, lassen sich noch weitere Kriterien zur
Charakterisierung eines Optimierungsproblems heran ziehen:

-   **Anzahl der unbekannten Parameter**, z.B. groß oder klein

-   **Eigenschaften der Zielfunktion**, z.B. Linearität, Beschränktheit,
    Konvexität, Stetigkeit, Differenzierbarkeit

-   **Charakteristik des Optimums**, z.B. Sattelpunkt, lokales oder
    globales Optimum

-   **Modelleigenschaften**, z.B. stochastisch oder deterministisch

Da wir uns intensiv mit der Bestimmung und numerischen Approximation von
Optima beschäftigen werden, macht es Sinn diese zuerst formal zu
beschreiben.

````{prf:definition} Lokales und globales Minimum
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Zielfunktion. Wir nennen einen Punkt $x^* \in \Omega$ einen **lokalen
Minimierer** der Zielfunktion $F$, falls es eine lokale Umgebung
$U \subset \Omega$ von $x^* \in U$ gibt, so dass für alle $x \in U$
gilt:
```{math}
:label: eq:lokales_minimum
F(x^*) \, \leq \, F(x), \quad \forall x \in U.
```
Den Funktionswert $F(x^*) \in \R$ nennen wir in diesem Fall ein
**lokales Minimum** der Zielfunktion $F$.

Wir nennen $x^* \in \Omega$ einen **globalen Minimierer** von $F$, falls
die Ungleichung {eq}`eq:lokales_minimum` für jede beliebige Umgebung
$U \subset \Omega$ gilt und somit insbesondere für $U = \Omega$. Den
Funktionswert $F(x^*) \in \R$ nennen wir in diesem Fall ein **globales
Minimum** der Zielfunktion $F$.

````

````{prf:definition} Striktes Minimum
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Zielfunktion. Sei $x^* \in \Omega$ ein lokaler Minimierer der
Zielfunktion $F$ in einer offenen Umgebung $U \subset \Omega$.

Wir nennen $F(x^*)$ ein **striktes Minimum** von $F$, falls gilt
```{math}
F(x^*) \, < \, F(x) \qquad \forall \, x\in U \setminus \lbrace x^* \rbrace.
```

````

````{prf:remark} Äquivalenz von Minimierung und Maximierung
In obiger Definition sprechen wir nur von Minima, jedoch ist klar, dass
sich jedes Maximierungsproblem durch einen Vorzeichenwechsel leicht in
ein Minimierungsproblem umformulieren lässt, d.h.,
```{math}
\max_{x\in\Omega} F(x) \quad \Leftrightarrow \quad \min_{x\in\Omega} -F(x) \ =: \ \min_{x_\in\Omega} G(x).
```
Formal lassen sich folgende Gleichungen zeigen:
```{math}
\begin{split}
\max_{x\in\Omega} F(x) \ &= \ - \min_{x\in\Omega} -F(x),\\
\operatorname{argmax}_{x\in\Omega} F(x) \ &= \ \operatorname{argmin}_{x\in\Omega} -F(x).
\end{split}
```

````

Da wir nun eine Charakterisierung von lokalen Minima haben, können wir
mit folgendem Satz die notwendigen Bedingungen für solch ein lokales
Minimum angeben.

````{prf:theorem} Notwendige Optimalitätsbedingungen 1. Ordnung
:label: thm:minimum_notwendig
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Zielfunktion. Sei $x^*\in \Omega$ ein lokaler Minimierer von $F$ in
$\Omega$ und die Zielfunktion $F$ sei stetig differenzierbar in einer
lokalen, offenen Umgebung $U \subset \Omega$ des Punkts $x^*$.

Dann gilt $\nabla F(x^*) = 0$.

````

````{prf:proof} 
Wir führen einen Beweis durch Widerspruch. Nehmen wir also an, dass
$x^* \in \mathbb{R}$ ein lokaler Minimierer von $F$ sei, jedoch aber
$\nabla F(x^*) \neq 0$ gelte. Wir wählen den Richtungsvektor
$\vec{p} \in \R^n$ mit $\vec{p} \coloneqq -\nabla F(x^*) \neq 0$. Es ist
somit klar, dass
```{math}
\langle \vec{p}, \nabla F(x^*) \rangle \ = \ - \langle \nabla F(x^*), \nabla F(x^*) \rangle \ = \ - ||\nabla F(x^*)||^2 \ < \ 0.
```
Da $\nabla F$ nach Vorraussetzung stetig in einer lokalen Umgebung
$U \subset \Omega$ von $x^*$ ist existiert ein $T > 0$, so dass auch
gilt:
```{math}
\langle \vec{p}, \nabla F(x^* + t\vec{p}) \rangle  \, < \, 0, \quad \text{ für alle } t\in [0,T].
```
Nach dem Satz von Taylor gilt aber auch für jedes $\tilde{t} \in (0,T]$:
```{math}
F(x^* + \tilde{t}\vec{p}) \ = \ F(x^*) +~\underbrace{\tilde{t}\langle \vec{p}, \nabla F(x^* + t\vec{p}) \rangle}_{<~0}, \quad \text{ für ein } t \in (0,\tilde{t}).
```
Somit gilt also $F(x^* + \tilde{t}\vec{p}) < F(x^*)$ für alle
$\tilde{t} \in (0, T]$ und wir haben offenbar eine Richtung
$\vec{p} \in \mathbb{R}^n / \lbrace 0\rbrace$ gefunden in der die
Funktionswerte von $F$ abnehmen. Also ist $x^* \in \Omega$ kein lokaler
Minimierer von $F$.

Das ist aber ein Widerspruch zur Annahme und somit ist die Behauptung
bewiesen. ◻

````

Die für die Optimierung interessanten Punkte $x^* \in \Omega$, die die
notwendigen Optimalitätsbedingungen aus {prf:ref}`thm:minimum_notwendig`
erfüllen, nennen wir stationäre Punkte.

````{prf:definition} Stationärer Punkt
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Zielfunktion. Wir nennen einen Punkt $x^* \in \Omega$ **stationären
Punkt** von $F$, falls $F$ in einer lokalen, offenen Umgebung
$U \subset \Omega$ von $x^*$ stetig differenzierbar ist und der Punkt
$x^*$ die Bedingung $\nabla F(x^*) = 0$ erfüllt.

````

Mit der Definition von stationären Punkten lässt sich folgendes Korollar
direkt ableiten.

````{prf:corollary} 
Jeder lokale Minimierer $x^* \in \Omega$ einer Zielfunktion
$F \colon \Omega \rightarrow \mathbb{R}$ ist ein stationärer Punkt.

````

````{prf:remark} Sattelpunkte
Die Umkehrung der Aussage in {prf:ref}`thm:minimum_notwendig` gilt im
Allgemeinen nicht. Man betrachte zum Beispiel die Zielfunktion
$F(x) \coloneqq -x^3$. Diese besitzt einen stationären Punkt in
$x^* = 0$, d.h., es gilt $\nabla F(0) = 0$. Dennoch handelt es sich
hierbei nicht um einen lokales Minimierer, sondern lediglich um einen
**Sattelpunkt**. Aus diesem Grund handelt es sich nur um notwendige und
nicht hinreichende Bedingungen.

````

Bei der Suche nach lokalen Minima einer Zielfunktion $F$ lässt sich ein
weiteres Kriterium anwenden, welches die zweite Ableitung der Funktion
verwendet.

````{prf:theorem} Notwendige Optimalitätsbedingungen 2. Ordnung
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Zielfunktion. Sei $x^*\in \Omega$ ein lokaler Minimierer von $F$ in
$\Omega$ und $F$ sei zweimal stetig differenzierbar in einer lokalen
Umgebung $U \subset \Omega$ von $x^*$, d.h., die Hessematrix
$\nabla^2 F$ von $F$ ist stetig in der offenen Umgebung
$U \subset \Omega$ von $x^*$.

Dann gilt $\nabla F(x^*) = 0$ und $\nabla^2 F(x^*)$ ist positiv
semidefinit, d.h., es gilt
```{math}
\langle \vec{p}, \nabla^2F(x^*) \vec{p} \rangle \, \geq \, 0 \qquad \forall \, \vec{p} \in \mathbb{R}^n.
```

````

````{prf:proof} 
In den Übungsaufgaben zu zeigen. ◻

````

Schlussendlich wollen wir auch eine hinreichende Bedingung für das
Vorliegen eines lokalen Minimums angeben.

````{prf:theorem} Hinreichende Optimalitätsbedingungen 2. Ordnung
:label: thm:minimum_hinreichend
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Zielfunktion. Sei $x^*\in \Omega$ ein lokaler Minimierer von $F$ in
$\Omega$ und $F$ sei zweimal stetig differenzierbar in einer lokalen
Umgebung $U \subset \Omega$ des Punkts $x^*$, d.h., die Hessematrix
$\nabla^2 F$ von $F$ sei stetig in einer offenen Umgebung
$U \subset \Omega$ von $x^*$. Außerdem gelte

1.  $\nabla F(x^*) \ = \ 0$,

2.  $\nabla^2 F(x^*)$ ist positiv definit.

Dann ist $F(x^*)$ ein striktes lokales Minimum von $F$.

````

````{prf:proof} 
Da die Hessematrix $\nabla^2 F$ von $F$ stetig und positiv definit in
$x^* \in \Omega$ ist nach Voraussetzung können wir einen Radius $r > 0$
finden, so dass $\nabla^2 F(x)$ positiv definit ist für alle
$x \in B_r(x^*)$. Für jeden Vektor
$\vec{p} \in \mathbb{R}^n / \lbrace 0\rbrace$ mit $||\vec{p}|| < r$ gilt
dann nach dem Satz von Taylor:
```{math}
F(x^* + \vec{p}) \ = \ F(x^*) +~\underbrace{\langle \vec{p}, \nabla F(x^*) \rangle}_{=~0} + \frac{1}{2} \langle \vec{p}, \nabla^2F(x^* + t\vec{p})\vec{p} \rangle, \quad \text{ für ein } t\in (0,1).
```
Da $||t\vec{p}|| < r$ ist nach Konstruktion wissen wir, dass
```{math}
\langle \vec{p}, \nabla^2F(x^* + t\vec{p})\vec{p} \rangle \ > \ 0
```
gilt und somit schon $F(x^* + \vec{p}) > F(x^*)$ gelten muss. Da
$\vec{p} \in \mathbb{R}^n / \lbrace 0 \rbrace$ mit $||\vec{p}|| < r$
beliebig gewählt war handelt es sich bei $F(x^*)$ um ein striktes
lokales Minimum der Zielfunktion $F$. ◻

````

````{prf:remark} Definitheit der Hessematrix
Die in {prf:ref}`thm:minimum_hinreichend` genannten Bedingungen sind
hinreichend, jedoch nicht notwendig für das Vorliegen eines strikten
lokalen Minimums. Dies sieht man ein, wenn man beispielsweise die
Zielfunktion $F(x) \coloneqq x^4$ betrachtet. $F$ besitzt ein striktes
lokales Minimum in $x^*=0$ und es gilt $\nabla F(0) = 0$, jedoch
verschwindet die zweite Ableitung $\nabla^2F(0) = 0$ und ist somit nicht
positiv definit.

````

Eine äußerst wertvolle Eigenschaft bei der Optimierung ist die
Konvexität einer Zielfunktion, da jedes lokale Optimum einer konvexen
Funktion bereits ein globales Optimum ist.

````{prf:definition} Konvexität
Sei $\Omega \subset \mathbb{R}^n$ ein offenes, zusammenhängendes Gebiet
und sei $F \colon \Omega \rightarrow \mathbb{R}$ eine reellwertige
Zielfunktion. Wir nennen $F$ **konvex** wenn für beliebige Vektoren
$x,y \in \Omega$ die folgende Ungleichung für alle
$0 \leq \alpha \leq 1$ gilt:
```{math}
F(\alpha x + (1-\alpha)y) \ \leq \ \alpha F(x) + (1-\alpha)F(y).
```
Anschaulich bedeutet Konvexität einer Funktion $F$, dass jede
Verbindungsgerade zwischen zwei Punkten $F(x)$ und $F(y)$ in diesem
Intervall oberhalb des Graphen der Funktion $F$ verläuft.

````

