(s:cg_verfahren)=
# Verfahren der konjugierten Gradienten

Im Folgenden wollen wir uns mit einem besonders eleganten Verfahren der
Optimierung beschäftigen: dem Verfahren der konjugierten Gradienten.
Ursprünglich wurde das Verfahren von Hestenes und Stiefel in
{cite:p}`hestenes_1952` im Jahr $1952$ vorgeschlagen. Obwohl das
Verfahren im Allgemeinen für die nichtlineare Optimierung eingesetzt
werden kann, wird es insbesondere zur Lösung von großen linearen
Gleichungssystemen $Ax = b$ mit symmetrischer, dünn besetzter, positiv
definiter Matrix $A \in \R^{n\times n}$ eingesetzt. Solche
Gleichungssysteme treten zum Beispiel bei der numerischen Modellierung
und Lösung partieller Differentialgleichungen auf. Das Verfahren lässt
sich in diesem Fall besonders anschaulich motivieren und herleiten.
Darum wollen wir uns im Folgenden zunächst auf das Lösen von großen
linearen Gleichungssystemen $Ax=b$ konzentrieren. Wir folgen bei der
Herleitung des Verfahren der konjugierten Gradienten der didaktisch sehr
gelungenen Arbeit von Jonathan Shewchuk in {cite:p}`shewchuk_1994`. Für
eine ansprechende, interaktive Visualisierung des Verfahren der
konjugierten Gradienten empfehlen wir den Mathematik Blog von Philipp
Wacker {cite:p}`wacker`.

(ss:cg_problemstellung)=
## Problemstellung

Sei im Folgenden also $A \in \mathbb{R}^{n \times n}$ eine sehr große
Matrix und $b \in \mathbb{R}^n$ ein reeller Vektor. Wir suchen einen
unbekannten Vektor $x \in \mathbb{R}^n$, der das lineare
Gleichungssystem
```{math}
:label: eq:LGS
Ax \ = \ b
```
löst. Wir suchen also nach denjenigen Koeffizienten, mit denen sich der
Vektor $b$ als Linearkombination aus Spaltenvektoren der Matrix $A$
darstellen lässt. Diese Koeffizienten entsprechen den Einträgen des
unbekannten Vektors $x$.

Aus der Vorlesung Einführung in die Numerik  in {cite:p}`numerik1` ist
bekannt, dass es genau dann eine eindeutige Lösung $x \in \mathbb{R}^n$
für die Gleichung {eq}`eq:LGS` gibt, falls die Determinante
$\operatorname{det}(A) \neq 0$ ist. Eine hinreichende Bedingung für die
Eindeutigkeit des Lösungsvektors $x \in \mathbb{R}^n$ ist es also zu
fordern, dass die Matrix $A$ symmetrisch und positiv definit ist.

Wir gehen aus diesem Grund im Folgenden immer davon aus, dass
$A \in \mathbb{R}^{n\times n}$ eine *symmetrische* und *positiv
definite* Matrix ist. In diesem Fall ist das Bestimmen einer Lösung von
{eq}`eq:LGS` ein gut-gestelltes Problem und die Lösung lässt sich direkt
angeben als:
```{math}
x \ = \ A^{-1}b.
```
Wie wir jedoch ebenfalls aus {cite:p}`numerik1` wissen ist die Inversion
einer Matrix $A \in \mathbb{R}^{n\times n}$ numerisch sehr aufwänding
und selbst unter Ausnutzung der Symmetrie lässt sich höchstens ein
Verfahren mit Rechenaufwand $\mathcal{O}(\frac{1}{6}n^3)$ angeben.
Sollte die Dimension des Problems jedoch sehr groß sein (d.h. wir nehmen
$n >> 1$ an), so ist eine direkte Lösung von {eq}`eq:LGS` mittels
Inversion nicht durchführbar. Glücklicherweise liefert uns das Verfahren
der konjugierten Gradienten (neben anderen iterativen
Lösungsalgorithmen) eine Möglichkeit das lineare Gleichungssystem
numerisch zu lösen. Man beachte, dass wir explizit darauf verzichten zu
fordern, dass die Matrix $A$ *dünnbesetzt* ist. In diesem Fall könnten
wir nämlich ebenfalls die numerischen Iterationsverfahren aus
{cite:p}`numerik1` anwenden.

Wir betrachten zunächst das folgende konvexe, quadratische
Optimierungsproblem der Form
```{math}
:label: eq:quadratic_problem
\min_{x\in \mathbb{R}^n} \, \frac{1}{2} \langle x, Ax \rangle - \langle b, x \rangle + c,
```
wobei $A$ und $b$ wie im Fall des linearen Gleichungssystems in
{eq}`eq:LGS` gewählt sind und $c \in \mathbb{R}$ eine beliebige, reelle
Konstante ist. Der folgende Satz liefert uns eine hilfreiche Aussage zur
Lösung des ursprünglichen Problems.

````{prf:theorem} Äquivalenzaussage für lineares Gleichungssystem
:label: thm:LGS_equivalent
Das konvexe, quadratische Minimierungsproblem in
{eq}`eq:quadratic_problem` ist äquivalent zum ursprünglichen linearen
Gleichungssystem in {eq}`eq:LGS`, d.h., jede Lösung von
{eq}`eq:quadratic_problem` ist schon Lösung von {eq}`eq:LGS` und anders
herum.

````

````{prf:proof} 
Für die erste Richtung des Beweises nehmen wir an, dass
$x^* \in \mathbb{R}^n$ eine Lösung des linearen Gleichungssystems
$Ax = b$ sei. Wir betrachten die hinreichenden Optimalitätsbedingungen
zweiter Ordnung aus {prf:ref}`thm:minimum_hinreichend` für das
Minimierungsproblem {eq}`eq:quadratic_problem`. Hierzu suchen wir
zunächst die stationären Punkte der Funktion
$F \colon \mathbb{R}^n \rightarrow \mathbb{R}$ mit
```{math}
F(x) \ \coloneqq \ \frac{1}{2} \langle x, Ax \rangle - \langle b, x \rangle + c.
```
Der Gradient von $F$ lässt sich wegen der Symmetrie von $A$ bestimmen
als
```{math}
\nabla F(x) \ = \ \frac{1}{2}(A + A^T)x - b \ = \ Ax - b \ \overset{!}{=} \ 0.
```
Alle stationären Punkte $x \in \mathbb{R}^n$ von $F$ mit
$\nabla F(x) = 0$ sind also gerade die Lösungen des linearen
Gleichungssystems $Ax = b$. Damit ist $x^*$ nach Vorraussetzung also
einziger stationärer Punkt von $F$. Um zu zeigen, dass
$x^* \in \mathbb{R}^n$ auch schon ein lokales Minimum von $F$ ist müssen
wir noch die Hessematrix von $F$ betrachten, welche gegeben ist durch:
```{math}
\nabla^2 F(x) \ = \ A.
```
Da $A$ nach Vorraussetzung positiv definit ist, ist auch die Hessematrix
$\nabla^2 F(x) = A$ positiv definit und somit sind die hinreichenden
Kriterien für das Vorliegen eines lokalen Minimums von $F$ im Punkt
$x^* \in \mathbb{R}^n$ erfüllt.

Für die Rückrichtung des Beweises nehmen wir, dass
$x^* \in \mathbb{R}^n$ ein lokales Minimum der Zielfunktion $F$ ist.
Damit folgt direkt, dass $x^*$ ein stationärer Punkt von $F$ ist und
somit muss gelten:
```{math}
\nabla F(x^*) \ = \ A x^* - b \ = \ 0.
```
Das bedeutet aber schon, dass $x^* \in \mathbb{R}^n$ Lösung des linearen
Gleichungssystems $Ax = b$ ist. ◻

````

Der {prf:ref}`thm:LGS_equivalent` erlaubt es uns also ein quadratisches
Optimierungsproblem der Form {eq}`eq:quadratic_problem` numerisch zu
lösen anstatt einen unbekannten Lösungsvektor für ursprüngliche lineare
Gleichungssystem {eq}`eq:LGS` zu finden.

Wir interessieren uns nun also für ein iteratives Verfahren, welches
eine Folge von Punkten $x_0,x_1,\ldots \in \mathbb{R}^n$ konstruiert,
die gegen ein Minimum von {eq}`eq:quadratic_problem` und somit gegen die
eindeutige Lösung des linearen Gleichungssystems {eq}`eq:LGS`
konvergiert. Hierfür benötigen wir noch zusätzliche Notation, um das
angestrebte Verfahren vernünftig zu beschreiben.

````{prf:definition} Fehler und Residuum
:label: def:residuum
Sei $x_{k+1} = G(x_k)$ ein Iterationsverfahren, dass gegen ein lokales
Minimum $x^* \in \mathbb{R}^n$ der quadratischen Funktion $F$ in
{eq}`eq:quadratic_problem` konvergiert, d.h., $x_k \rightarrow x^*$ für
$k \rightarrow \infty$. Dann können wir die beiden folgenden Begriffe
definieren:

1.  Wir bezeichnen den Vektor $e_k \in \mathbb{R}^n$ mit
    ```{math}
    e_k \ \coloneqq \ x_k - x^*
    ```
    als den aktuellen **Fehler**, den man durch den aktuellen Punkt
    $x_k \in \mathbb{R}^n$ macht.

2.  Wir bezeichnen den Vektor $r_k \in \mathbb{R}^n$ mit
    ```{math}
    r_k \ \coloneqq \ b - Ax_k
    ```
    als das aktuelle **Residuum**, das man durch den aktuellen Punkt
    $x_k \in \mathbb{R}^n$ erhält.

````

````{prf:remark} Fehler und Residuum
:label: rem:fehler_residuum
In Bezug auf {prf:ref}`def:residuum` lassen sich folgende Aussagen
festhalten:

1.  Der Fehler $e_k \in \mathbb{R}^n$ ist eher abstrakter Natur und
    dient zur besseren Analyse des Verfahrens der konjugierten
    Gradienten. Explizit werden wir diesen Vektor jedoch nie bestimmen
    können innerhalb des Iterationsverfahrens, da wir dann schon fertig
    wären mit einem einfach Update der Form $x^* = x_k - e_k$.

2.  Wie wir bereits im Beweis von {prf:ref}`thm:LGS_equivalent` gesehen
    haben, lässt sich das Residuum $r_k \in \mathbb{R}^n$ außerdem wie
    folgt umschreiben:
    ```{math}
    r_k \ = \ \underbrace{b - Ax_k}_{= -\nabla F(x_k)} \ = \ Ax^* - Ax_k \ = \ A (x^* - x_k) \ = \ -Ae_k.
    ```
    Daher lässt sich das Residuum $r_k$ auch als die Richtung des
    stärksten Abstiegs interpretieren und es ist klar, dass $r_k$ immer
    orthogonal zu den Niveaulinien der Funktion $F$ steht.

````

{numref}`fig:residuum` illustriert anschaulich die geometrische
Bedeutung der beiden in {prf:ref}`def:residuum` eingeführten Vektoren.
Wie man unschwer erkennt zeigen Fehler und Residuum im Allgemeinen nicht
in die selbe Richtung. Das erklärt auch warum das
Gradientenabstiegsverfahren in {ref}`ss:gradient_descent` selbst bei
optimaler Schrittweite $\alpha_k > 0$ nicht in einem Schritt die
gesuchte Lösung $x^* \in \mathbb{R}^n$ erreicht.

```{figure} ../atelier/img/residuum.png
---
name: "fig:residuum"
---
Visualisierung des Fehlers $e_0 \in \mathbb{R}^n$ und des Residuums
$r_0 \in \mathbb{R}^n$ für einen Startpunkt $x_0 \in \mathbb{R}^n$.
```
(ss:motivation)=
## Motivation

Um das Vorgehen beim Verfahren der konjugierten Gradienten zu motivieren
rufen wir uns noch einmal das Gradientenabstiegsverfahren aus
{ref}`ss:gradient_descent` in Erinnerung. Nehmen wir an wir befinden uns
im $k$-Schritt des Gradientenabstiegsverfahrens in Algorithmus
{prf:ref}`alg:gradient_descent_adaptive` in einem Punkt
$x_k \in \mathbb{R}^n$ und es sei eine Schrittweite $\alpha_k > 0$
gegeben. Dann erhalten wir den nächsten Punkt $x_{k+1} \in \mathbb{R}^n$
der Iterationsfolge durch folgendes Update:
```{math}
x_{k+1} \ = \ x_k - \alpha_k \nabla F(x_k) \ = \ x_k + \alpha_k r_k,
```
wobei $r_k \in \mathbb{R}^n$ das aktuelle Residuum des Punktes $x_k$
bezeichnet. Wir machen also in Richtung des steilsten Gradientenabstiegs
einen Schritt der Länge $\alpha_k > 0$. Da die Abstiegsrichtung in jedem
Schritt $x_k \rightarrow x_{k+1}$ orthogonal zu den Niveaulinien von $F$
steht, erhält man typischerweise einen Zickzack-Pfad durch das
Gradientenabstiegsverfahren (vgl.
{numref}`fig:gradient_descent_adaptive`).

Um dieses typische Verhalten besser zu verstehen können wir eine
Vorüberlegung zur Schrittweitenwahl für das quadratische
Optimierungsproblem in {eq}`eq:quadratic_problem` machen. Hierzu gehen
wir analog zur Bestimmung der optimalen Schrittrichtung in
{ref}`ss:gradient_descent` vor, nur dass wir diesmal die Schrittrichtung
$p_k \coloneqq -\nabla F(x_k)$ festhalten und bezüglich der unbekannten
Schrittweite optimieren.

Wir gehen davon aus, dass wir das lokale Minimum von $F$ noch nicht
erreicht haben, denn dann wäre $\alpha_k = 0$. Wir suchen also eine
Schrittweite $\alpha > 0$, so dass der Funktionswert $F(x_{k+1})$
entlang der Linie $x_k - \alpha \nabla F(x_k)$ minimal wird. Da $F$ eine
quadratische Funktion ist, wissen wir, dass ein eindeutiges Minimum
$\alpha_k$ entlang dieser Linie existieren muss. Wir nutzen also die
notwendigen Optimalitätsbedingungen aus {prf:ref}`thm:minimum_notwendig`
für das totale Differential, um folgenden Zusammenhang herzustellen:
```{math}
\begin{split}
\frac{\mathrm{d}}{\mathrm{d}\alpha} F(x_{k+1}) \ &= \ \Bigl\langle \nabla F(x_{k+1}), \frac{\mathrm{d}x_{k+1}}{\mathrm{d}\alpha} \Bigr\rangle \ = \ \Bigl\langle \nabla F(x_{k+1}), \frac{\mathrm{d} (x_k - \alpha \nabla F(x_k))}{\mathrm{d}\alpha} \Bigr\rangle \\ 
\ &= \ \langle \nabla F(x_{k+1}), -\nabla F(x_k) \rangle \ = \  \langle \nabla F(x_{k+1}), p_k \rangle \ \overset{!}{=} \ 0.
\end{split}
```
Das Ergebnis ist durchaus interessant. Die optimale Schrittweite
$\alpha > 0$ muss so gewählt werden, dass der nächste Punkt
$x_{k+1} \in \mathbb{R}^n$ der Iterationsfolge an der Stelle liegt an
der unsere Abstiegsrichtung orthognal auf den Gradienten der Funktion
$\nabla F(x_{k+1})$ trifft. Das bedeutet, dass die optimale Abfolge der
Abstiegsrichtungen im quadratischen Fall eine Menge von $90$ Grad
Zickzack-Linien ergibt, was zu unseren Beobachtungen in
{numref}`fig:gradient_descent_adaptive` passt. Da jedoch der Punkt
$x_{k+1} \in \mathbb{R}^n$ bislang noch unbekannt ist, können wir das
optimale $\alpha_k$ nicht in dieser Form angeben. Das folgende Lemma
bestimmt die optimale Schrittweite im Fall der quadratischen Optimierung
in {eq}`eq:quadratic_problem`.

````{prf:lemma} Optimale Schrittweite
:label: lem:optimale_schrittweite
Sei $F \colon \mathbb{R}^n \rightarrow \mathbb{R}$ die quadratische
Funktion aus {eq}`eq:quadratic_problem`. Wir betrachten das
Gradientenabstiegsverfahren im $k$-ten Iterationsschritt mit einer
unbekannten Schrittweite $\alpha_k > 0$, die jedoch so gewählt werden
muss, dass
```{math}
\langle \nabla F(x_{k+1}), \nabla F(x_k) \rangle \ \overset{!}{=} \ 0.
```
Sei außerdem $r_k = b - Ax_k$ das Residuum im aktuellen
Iterationsschritt. Dann lässt sich die optimale Schrittweite $\alpha_k$
berechnen als:
```{math}
:label: eq:optimal_step-size
\alpha_k \ = \ \frac{\langle r_k, r_k \rangle}{\langle r_k, Ar_k \rangle}.
```

````

````{prf:proof} 
Wir erinnern uns daran, dass $r_k = -\nabla F(x_k) = b - Ax_{k}$ ist und
somit können wir folgern:
```{math}
\begin{split}
0 \ &\overset{!}{=} \ \langle \nabla F(x_{k+1}), \nabla F(x_k) \rangle \ = \ \langle r_{k+1}, r_k \rangle \ = \ \langle b - Ax_{k+1}, r_k \rangle \\
 \ &= \ \langle b - A(x_k + \alpha_kr_k), r_k \rangle \ = \  \langle b - Ax_k, r_k \rangle - \alpha_k \langle Ar_k, r_k \rangle \\
 \ &= \langle r_k, r_k \rangle - \alpha_k \langle r_k, Ar_k \rangle
\end{split}
```
Da wir $A$ als positiv definit angenommen haben, können wir die
Gleichung umstellen und erhalten so die behauptete Berechnungsformel für
$\alpha_k$ in {eq}`eq:optimal_step-size`. ◻

````

Obwohl wir die optimale Schrittweite $\alpha_k$ in
{eq}`eq:optimal_direction` für das quadratische Optimierungsproblem
{eq}`eq:quadratic_problem` bestimmen konnten ist das
Gradientenabstiegsverfahren weit davon entfernt optimal zu sein. Trotz
optimaler Schrittweiten und optimaler Abstiegsrichtungen erhalten wir
eine Folge von Richtungsvektoren, die immer wieder in die gleiche
Richtung zeigen (siehe {numref}`fig:zig-zag`). Das ist numerisch gesehen
äußerst ineffizient. Man könnte sich also fragen, warum man nicht
einfach nur zwei orthogonale Schritte macht und die Schrittweiten als
Summe der optimalen Schrittweiten der geraden bzw. ungeraden
Iterationsschritte $k \in \mathbb{N}$ wählt. In der Tat würde man für
$N \in \mathbb{N}$ Schritte des Gradientenabstiegsverfahren im selben
Punkt $x_N \in \mathbb{R}^n$ mit nur zwei Iterationen landen, wie in
{numref}`fig:zig-zag` illustriert ist.

```{figure} ../atelier/img/residuum.png
---
name: "fig:zig-zag"
---
Vergleich des Gradientenabstiegsverfahrens mit optimaler Schrittweite
$\alpha_k > 0$ aus {prf:ref}`lem:optimale_schrittweite` (links) mit
einem idealen Abstiegsverfahren (rechts), bei dem alle orthogonalen
Teilschritte zusammengefasst sind.
```

Leider können wir nicht alle Schrittweiten aufaddieren, da wir zur
Berechnung der optimalen Schrittlänge $\alpha_k > 0$ bereits alle
vorangegangen Schritte $k=0,\ldots,k-1$ kennen müssten. Außerdem würde
ein großer, zusammengefasster Schritt in die erste der Richtungen
eventuell dazu führen, dass man keinen Abstieg der Funktionswerte von
$F$ mehr realisiert, sondern einen Aufstieg. Diese Beobachtung ist in
{numref}`fig:two-step` illustriert.

```{figure} ../atelier/img/two-step.png
---
name: "fig:two-step"
---
Illustration eines idealen Abstiegsverfahrens mit zwei orthogonalen
Richtungen. Man beachte, dass die Schrittweite $\alpha_0 > 0$ so gewählt
werden muss, dass man im ersten Schritt nicht in einem Punkt
$x_1 \in \mathbb{R}^2$ mit minimalen Funktionswert $F(x_1)$ entlang der
Richtung $x_0 - \alpha_0 \nabla F(x_0)$ endet.
```

Die Ideallösung wäre natürlich von einem Startpunkt
$x_0 \in \mathbb{R}^n$ in nur einem Schritt zum lokalen Optimum
$x^* \in \mathbb{R}^n$ zu gelangen. Da wir aber den Punkt $x^*$ a-priori
nicht kennen ist das eine unrealistische Forderung. Dennoch lässt sich
zeigen, dass das Gradientenabstiegsverfahren mit der optimalen
Schrittweite $\alpha_k$ in {eq}`eq:optimal_step-size` im Fall des
quadratischen Optimierungsproblems {eq}`eq:quadratic_problem` in genau
einem Schritt zum lokalen Minimum $x^* \in \mathbb{R}^n$ führt, wenn der
Fehler $e_0 = x_0 - x^*$ ein Eigenvektor von $A$ ist. Wir wissen nämlich
aus {prf:ref}`rem:fehler_residuum`, dass $r_k = -Ae_k$ gilt und somit
erhalten wir für das Gradientenabstiegsverfahren im ersten Schritt:
```{math}
x_1 \ = \ x_0 - \alpha_0 \nabla F(x_0) \ = \ x_0 + \alpha_0 r_0 \ = \ x_0 - \alpha_0 A e_0 \ = \ x_0 - \alpha_0 \lambda e_0.
```

Man müsste bei der Wahl des Startpunktes $x_0\in \R^n$ jedoch viel Glück
haben, um diese Forderung zu erfüllen. Darum wollen wir uns mit
alternativen Ideen beschäftigen.

(ss:orthogonal_descent)=
## Orthogonale Abstiegsrichtungen

Wir wünschen uns einen Algorithmus, der ähnlich dem
Gradientenabstiegsverfahren nur orthogonale Richtungen
$\lbrace d_0, \ldots, d_{n-1} \rbrace$ mit $d_i \in \mathbb{R}^n$ für
$0 \leq i \leq n-1$ verwendet, jedoch mit der Einschränkung, dass diese
nur ein einziges Mal genutzt werden können. Ziel dieses Verfahrens soll
es außerdem sein durch $n$ Schritte in die jeweils $n$ orthogonalen
Richtungen $\lbrace d_0, \ldots, d_{n-1} \rbrace$ das lokale Minimum der
Funktion zu erreichen. Damit hätten wir ein Itrerationsverfahren der
Form
```{math}
:label: eq:orthogonal_descent
x_{k+1} \ = \ x_k + \alpha_k d_k, \quad \alpha_k > 0, \ k = 0,\ldots,n-1
```
gewonnen. Wir könnten dies erzwingen indem wir im $k$-ten Schritt des
Iterationsverfahrens fordern, dass ein Schritt in Richtung
$d_k \in \mathbb{R}^n$ dazu führt, dass der Fehler
$e_{k+1} \in \mathbb{R}^n$ keinerlei Komponenten dieser Richtung mehr
enthält, d.h. wir fordern
```{math}
:label: eq:req_error
\langle e_{k+1}, d_k \rangle \ \overset{!}{=} \ 0.
```
Da wir den zu erwartenden Fehler $e_{k+1}$ in Bezug auf den aktuellen
Punkt $x_k \in \mathbb{R}^n$ folgendermaßen umschreiben können:
```{math}
e_{k+1} \ = \ x_{k+1} - x^* \ = \ x_k + \alpha_k d_k - x^* \ = \ e_k + \alpha_k d_k,
```
können wir die Forderung {eq}`eq:req_error` umformulieren zu:
```{math}
\langle e_k + \alpha_k d_k, d_k \rangle \ \overset{!}{=} \ 0.
```
Hieraus können wir die optimale Schrittweite $\alpha_k > 0$ in Richtung
$d_k \in \mathbb{R}^n$ ableiten als
```{math}
:label: eq:optimal_step-size_orthogonal
\alpha_k \ = \ - \frac{\langle e_k, d_k \rangle}{\langle d_k, d_k \rangle}.
```
Obwohl wir in {eq}`eq:optimal_step-size_orthogonal` eine optimale
Schrittweite $\alpha_k$ für das Verfahren mit orthogonalen
Abstiegsrichtungen in {eq}`eq:orthogonal_descent` bestimmen konnten,
hilft und diese nicht in der praktischen Anwendung des Verfahrens, da
sie von dem unbekannten Fehlervektor $e_k \in \mathbb{R}^n$. Dieser
hängt natürlich von der unbekannten Lösung $x^* \in \mathbb{R}^n$ ab und
wenn wir diese kennen würden, so müssten wir kein iteratives Verfahren
konstruieren. Selbst wenn man den Fehler $e_k$ weiter rekursiv
umschreibst, so würde man schlussendlich doch bei einer Abhängigkeit des
initialen Fehlers $e_0$ landen. Wir müssen uns also vorerst von dieser
Idee verabschieden und nach einer alternativen Möglichkeit suchen.

(ss:conjugated_descent)=
## Konjugierte Abstiegsrichtungen

Obwohl unsere Idee von orthogonalen Abstiegsrichtungen in
{ref}`ss:orthogonal_descent` nicht zum Ziel geführt hat, so war die Idee
gar nicht schlecht. Das Hauptproblem an der ursprünglichen Idee liegt in
der Forderung {eq}`eq:req_error`, nämlich dass der Fehlervektor
$e_{k+1}$ orthogonal zur aktuellen Richtung $d_k$ stehen soll. Diese
Forderung führt nämlich dazu, dass man orthogonale Vektoren erhält, die
nicht an die Geometrie des quadratischen Minimierungsproblems
{eq}`eq:quadratic_problem` angepasst sind.

Wenn man sich die Niveaulinien der Funktion $F$ genauer anschaut (siehe
zum Beispiel Abbildung {numref}`fig:two-step`), so erkennt man, dass es
Richtungen gibt entlang derer die Abstiegsrichtung zum lokalen Minimum
$x^* \in \mathbb{R}^n$ steiler verläuft als entlang der anderen
Richtungen. Die geometrischen Eigenschaften des Graphen von $F$ sind
maßgeblich durch die Gestalt der Matrix $A$, genauer gesagt durch deren
Eigenvektoren bestimmt. Daher wollen wir diese Eigenschaften bei der
Konstruktion eines iterativen Abstiegsverfahren berücksichtigen. Hierzu
führen wir folgendes hilfreiche Konzept ein.

````{prf:definition} Konjugierte Vektoren
Sei $A \in \mathbb{R}^{n\times n}$ eine symmetrische, positiv definite
Matrix und $u,v \in \mathbb{R}^n / \lbrace 0\rbrace$ zwei Vektoren. Wir
nennen $v$ und $w$ **konjugiert bezüglich $\mathbf{A}$** oder auch
**$\mathbf{A}$-orthogonal** falls gilt
```{math}
\langle v, Aw \rangle \ = \ \langle w, Av \rangle \ = \ 0.
```

````

Anstatt nun also die Orthogonalität unserer Richtungsvektoren
$\lbrace d_0, \ldots, d_{n-1}\rbrace$ zu erzwingen wie in
{ref}`ss:orthogonal_descent`, fordern wir nun, dass diese Vektoren
konjugiert bezüglich der Matrix $A$ und damit besser an das Problem
angepasst sind.

````{prf:remark} 
Es ist leicht einzusehen, dass $A$-orthogonal und orthogonal die selbe
Eigenschaft beschreiben, falls die Matrix $A$ ein Vielfaches der
Einheitsmatrix $I_n \in \mathbb{R}^{n \times n}$ ist. In diesem Fall ist
das quadratische Optimierungsproblem {eq}`eq:quadratic_problem`
symmetrisch in alle Richtungen.

````

Anschaulich lässt sich die Forderung nach $A$-Orthogonalität auch so
deuten, dass wir ein Paar von Vektoren
$v,w \in \mathbb{R}^n / \lbrace 0\rbrace$ suchen, welche in einem Winkel
so zueinander stehen, dass wenn man die Niveaulinien der Funktion $F$
symmetrisch reskaliert, diese Vektoren anschließend orthogonal
zueinander stehen. Diese Idee ist in Abbildung {numref}`fig:conjugacy`
dargestellt.

```{figure} ../atelier/img/conjugacy.png
---
name: "fig:conjugacy"
---
Illustration der Geometrie von konjugierten Vektoren im Referenzsystem
$\mathbb{R}^2$ (links) und der selben Vektoren in einem symmetrisierten
System bezüglich der Matrix $A$ (rechts).
```

Anstatt also ein Abstiegsverfahren der Form {eq}`eq:orthogonal_descent`
mit orthogonalen Vektoren zu verwenden, wollen wir ein Abstiegsverfahren
mit $A$-orthogonalen Vektoren $\lbrace d_0, \ldots, d_{n-1}\rbrace$
konstruieren, d.h., wir verwenden das Iterationsschema
```{math}
:label: eq:conjugate_descent
x_{k+1} \ = \ x_k + \alpha_k d_k, \quad \alpha_k > 0, \ k=0,\ldots,n-1
```
wobei für die Abstiegsrichtungen $d_k \in \mathbb{R}^n$ gelten soll:
```{math}
\langle d_i, Ad_j \rangle = 0 \qquad \text{ für alle } i\neq j.
```

Wir nehmen für den Moment an, dass wir einen numerischen Algorithmus
kennen mit dem wir eine Menge von $A$-orthogonalen Vektoren
$\lbrace d_0, \ldots, d_{n-1} \rbrace$ konstruieren können. Wie man
diese Menge konkret erhält werden wir uns im Anschluss erschließen. Sei
also nun im Folgenden $\lbrace d_0, \ldots, d_{n-1} \rbrace$ eine
gegebene Menge von $A$-orthogonalen Vektoren. Dann stellen wir uns die
Frage, wie die optimalen Schrittweiten $\alpha_k > 0$ in
{eq}`eq:conjugate_descent` gewählt werden müssen, um in $n$ Schritten
das lokale Minimum $x^* \in \mathbb{R}^n$ der Funktion $F$ zu erhalten.
Man beachte hierbei, dass wir nicht nur daran interessiert sind den
Punkt $x^*$ genügend gut zu approximieren, sondern wir fordern die
eindeutige Lösung des linearen Gleichungssystems $Ax = b$ in $n$
Schritten zu finden, d.h., wir nehmen explizit $x_n = x^*$ an.

Um das lokale Minimum wirklich in $n$ Schritten zu erreichen müssen wir
fordern, dass wir in jede Richtung $d_k$ nur einmal gehen und der
entstehende Fehler $e_{k+1}$ $A$-orthogonal hierzu ist. Das entspricht
der Forderung, dass man im entzerrten Problem auf der rechten Seite von
{numref}`fig:conjugacy` nur orthogonale Richtungen verwendet. Wir wollen
also folgende Eigenschaft erzwingen:
```{math}
:label: eq:req_error_conjugated
\langle Ae_{k+1}, d_k \rangle \ = \ 0.
```

Analog zur Idee der orthogonalen Richtungen in
{ref}`ss:orthogonal_descent` können wir den Fehler $e_{k+1}$ in
{eq}`eq:req_error_conjugated` wieder entwickeln, um die optimale
Schrittweitenlänge $\alpha_k > 0$ zu bestimmen
```{math}
\begin{split}
0 \ &\overset{!}{=} \ \langle d_k, A e_{k+1} \rangle \ = \ \langle d_k, A(x_{k+1} - x^*) \rangle \ = \ \langle d_k, A(x_k + \alpha_k d_k - x^*) \rangle \\
\ &= \ \langle d_k, A(e_k + \alpha_k d_k) \rangle \ = \ \langle d_k, -r_k + \alpha_k Ad_k \rangle \ = \ \alpha_k \langle d_k, Ad_k \rangle - \langle d_k, r_k\rangle.
\end{split}
```
Da wir $A$ als positiv definit vorausgesetzt haben, können wir die
folgende Gleichung umstellen zu
```{math}
:label: eq:optimal_step-size_conjugated
\alpha_k \ = \ \frac{\langle d_k, r_k \rangle}{\langle d_k, A d_k \rangle}.
```
Im Gegensatz zur Idee der orthogonalen Richtungen in
{eq}`eq:optimal_step-size_orthogonal` lässt sich der Ausdruck in
{eq}`eq:optimal_step-size_conjugated` explizit berechnen und hängt nicht
von dem unbekannten lokalen Minimum $x^* \in \mathbb{R}^n$ ab.
Zusammenfassend heißt das, dass wir aus der Bedingung, dass die
Abstiegsrichtung $d_k \in \mathbb{R}^n$ $A$-orthogonal zum Fehlervektor
$e_{k+1} \in \mathbb{R}^n$ sein soll, eine Schrittlänge $\alpha_k > 0$
finden konnten, welche diese Bedingung erfüllt.

Andersherum könnte man fragen, welche Bedingung man aus der Optimalität
einer unbekannten Schrittlänge $\alpha > 0$ folgern könnte. Dazu
betrachen wir wieder die notwendigen Optimalitätsbedingungen im totalen
Differential
```{math}
\begin{split}
0 \ &\overset{!}{=} \ \frac{\mathrm{d}}{\mathrm{d}\alpha}F(x_{k+1}) \ = \ \langle  \nabla F(x_{k+1}), \frac{\mathrm{d}}{\mathrm{d}\alpha}x_{k+1} \rangle \ = \ \langle -r_{k+1}, \frac{\mathrm{d}}{\mathrm{d}\alpha} (x_k + \alpha d_k) \rangle \\
\ &= \ \langle -r_{k+1}, d_k \rangle \ = \ \langle Ae_{k+1}, d_k \rangle.
\end{split}
```
Wir erhalten also für die Optimalität der unbekannten Schrittlänge
$\alpha > 0$, dass die Abstiegsrichtung $d_k \in \mathbb{R}^n$ und der
Fehlervektor $e_{k+1}$ konjugiert bezüglich der Matrix $A$ sein müssen.
Das ist aber genau die Eigenschaft, die wir bereits in
{eq}`eq:req_error_conjugated` gefordert hatten. Wegen der positiven
Definitheit der Matrix $A$ können wir ebenfalls folgern, dass die
Forderung, dass $e_{k+1}$ und $d_k$ konjugiert bezüglich $A$ sind, bei
optimaler Schrittweite $\alpha_k$ aus
{eq}`eq:optimal_step-size_conjugated` in jedem Schritt zu einem Abstieg
in Richtung $d_k$ führt. Wir erhalten also für eine gegebene Menge von
$A$-orthogonalen Vektoren $\lbrace d_0, \ldots, d_{n-1} \rbrace$ ein
Iterationsverfahren mit optimalen Schrittlangen $\alpha_k > 0$, die wir
in {eq}`eq:optimal_step-size_conjugated` angeben können und die uns
einen Abstieg garantieren.

Zunächst benötigen wir die Einsicht aus folgendem Lemma, die es uns
ermöglicht eine Basis aus $A$-orthogonalen Vektoren des $\R^n$ zu
betrachten.

````{prf:lemma} Basis von A-orthogonalen Vektoren
:label: lem:A-orthogonal_basis
Sei $A \in \mathbb{R}^{n\times n}$ eine symmetrische, positiv definite
Matrix und sei $\lbrace d_0, \ldots, d_{n-1}\rbrace$ eine Menge von
$A$-orthogonalen Vektoren mit $d_k \in \R^n \setminus \lbrace 0\rbrace$,
d.h., es gilt
$\langle d_i, A d_j \rangle = \langle A d_i, d_j \rangle = 0$ für alle
$i \neq j$.

Dann bilden die Vektoren
$d_0, \ldots, d_{n-1} \in \R^n \setminus \lbrace 0\rbrace$ eine Basis
des $\R^n$.

````

````{prf:proof} 
In den Übungsaufgaben zu zeigen. ◻

````

Folgender Satz zeigt uns, dass das Verfahren für eine gegebene Menge von
$A$-konjugierten Vektoren in der Tat in $n$ Schritten das lokalen
Minimum $x^* \in \mathbb{R}^n$ von $F$ erreicht.

````{prf:theorem} Konvergenz des CG-Verfahrens
:label: thm:cg_convergence
Sei eine Menge von $A$-konjugierten Vektoren
$\lbrace d_0, \ldots, d_{n-1} \rbrace$ mit
$d_k \in \mathbb{R}^n / \lbrace 0 \rbrace$ gegeben. Dann konvergiert das
Abstiegsverfahren in konjugierten Richtungen
```{math}
:label: eq:conjugate_descent_optimal
x_{k+1} \ = \ x_k + \alpha_k d_k, \qquad \alpha_k \ = \ \frac{\langle r_k, d_k \rangle}{\langle d_k, Ad_k\rangle}
```
in genau $n$ Schritten gegen die Lösung $x^* \in \mathbb{R}^n$ des
quadratischen Optimierungsproblems {eq}`eq:quadratic_problem`.

````

````{prf:proof} 
Für den Beweis der Konvergenz des Iterationsverfahrens
{eq}`eq:conjugate_descent_optimal` betrachten wir zunächst den initialen
Fehler $e_0\in \R^n$ durch die Wahl eines beliebigen Startpunktes
$x_0 \in \mathbb{R}^n$. Da wir $A$ als positiv definit angenommen haben
folgt mit {prf:ref}`lem:A-orthogonal_basis`, dass die Menge
$\lbrace d_k \rbrace_{k=0,\ldots,n-1}$ eine Basis des $\mathbb{R}^n$
bildet. Daher können wir den initialen Fehler $e_0 \in \mathbb{R}^n$ als
Linearkombination in dieser Basis darstellen als:
```{math}
:label: eq:error_basis
e_0 \ = \ \sum_{k=0}^{n-1} \delta_k d_k.
```
Um die unbekannten Koeffizienten $\delta_k \in \mathbb{R}^n$ zu
bestimmen können wir obige Gleichung nun jeweils von links mit einem
Vektor $d_i^T A, i=0,\ldots,n-1$ multiplizieren und erhalten so für
jeden Index eine Gleichung
```{math}
\langle d_i^T A, e_0 \rangle \ = \ \langle d_i^T A, \sum_{k=0}^{n-1} \delta_k d_k \rangle \ = \ \sum_{k=0}^{n-1} \delta_k \langle d_i^T A, d_k \rangle \ = \ \delta_i \langle d_i^T A, d_i \rangle.
```
Hierbei haben wir die Linearität des Skalarproduktes in $\mathbb{R}^n$
ausgenutzt und verwendet, dass die Vektoren
$\lbrace d_k \rbrace_{k=0,\ldots,n-1}$ konjugiert bezüglich der Matrix
$A$ sind. Damit können wir nach den unbekannten Koeffizienten
$\delta_i \in \R$ in jeder Gleichung auflösen und erhalten so einen
Ausdruck für die unbekannten Koeffizienten:
```{math}
\delta_i \ = \ \frac{\langle d_i^T A, e_0 \rangle}{\langle d_i^T A, d_i \rangle}.
```
Man beachte, dass dieser Ausdruck wohldefiniert ist, da wir angenommen
haben, dass die Matrix $A$ positiv definit ist. Wir addieren eine Null
hinzu, indem wir Terme hinzufügen, die $A$-konjugiert zur Richtung
$d_i \in \R^n$ sind:
```{math}
:label: eq:delta_coefficients
\delta_i \ = \ \frac{\langle d_i^T A, e_0 \rangle}{\langle d_i^T A, d_i \rangle} + \underbrace{\frac{\langle d_i^T A, \sum_{k=0}^{i-1} \alpha_k d_k \rangle}{\langle d_i^T A, d_i \rangle}}_{=\: 0} \ = \ \frac{\langle d_i^T A, e_0 +\sum_{k=0}^{i-1} \alpha_k d_k \rangle}{\langle d_i^T A, d_i \rangle}.
```
Wir verwenden wieder den Trick, dass sich der Fehlervektor
$e_{i+1} \in \mathbb{R}^n$ entwickeln lässt zu
$e_{i+1} = e_i + \alpha_i d_i$ und somit können wir rekursiv herleiten,
dass
```{math}
:label: eq:error_recursive
e_i \ = \ e_0 + \sum_{k=0}^{i-1} \alpha_k d_k.
```
Nun können wir die Gleichung {eq}`eq:error_recursive` in die Darstellung
der Koeffizienten $\delta_i$ in {eq}`eq:delta_coefficients` einsetzen
und erhalten:
```{math}
\delta_i \ = \ \frac{\langle d_i^T A, e_0 +\sum_{k=0}^{i-1} \alpha_k d_k \rangle}{\langle d_i^T A, d_i \rangle} \ = \ \frac{\langle d_i^T A, e_i \rangle}{\langle d_i^T A, d_i \rangle} \ = \ \frac{\langle d_i, Ae_i \rangle}{\langle d_i^T A, d_i \rangle} \ = \ - \frac{\langle d_i, r_i \rangle}{\langle d_i^T A, d_i \rangle} \ = \ -\alpha_i.
```
Das bedeutet, dass die Koeffizienten $\delta_i$ in
{eq}`eq:delta_coefficients` gerade den negativen optimalen Schrittweiten
$\alpha_i$ in {eq}`eq:optimal_step-size_conjugated` entsprechen, d.h.,
$\delta_i$ = $- \alpha_i$. Aus der Basisdarstellung des initialen
Fehlers $e_0 = x_0 - x^*$ in {eq}`eq:error_basis` können wir somit die
Behauptung des Satzes folgern:
```{math}
x^* \ = \ x_0 - e_0 \ = \ x_0 - \sum_{k=0}^{n-1} \delta_k d_k \ = \ x_0 + \sum_{k=0}^{n-1} \alpha_k d_k \ = \ x_n.
```
 ◻

````

````{prf:remark} Veränderung des Fehlers im Iterationsverfahren
:label: rem:error_cg
 \
Anstatt im Beweis von {prf:ref}`thm:cg_convergence` zu zeigen, dass sich
das lokale Minimum $x^* \in \mathbb{R}^n$ durch das Iterationsverfahren
zerlegen lässt, hätte man auch zeigen können, dass der Fehlervektor
$e_i \in \mathbb{R}^n$ in jedem Schritt des Iterationsverfahren kleiner
wird. Es gilt nämlich nach {eq}`eq:error_recursive`:
```{math}
e_i \ = \ e_0 + \sum_{k=0}^{i-1} \alpha_k d_k \ = \ \sum_{k=0}^{n-1}\delta_k d_k + \sum_{k=0}^{i-1}-\delta_k d_k \ = \ \sum_{k=i}^{n-1} \delta_k d_k.
```
Man sieht also, dass für eine wachsende Anzahl an Iterationen
$i=0,\ldots,n-1$ der Fehlerterm $e_i \in \mathbb{R}^n$ immer weniger
Terme hat, bis er schlussendlich ganz verschwindet.

Außerdem sagt es uns, dass der Abstieg mit konjugierten Richtungen in
dem Sinne optimal ist, als dass der Fehlerterm
$e_i = \sum_{k=i}^{n-1} \delta_k d_k$ keine Anteile der Richtungen
$\lbrace d_j \rbrace_{j=0,\ldots,k-1}$ mehr besitzt. Wir müssen also
nicht mehr entlang dieser Richtungen gehen, um zum lokalen Minimum
$x^* \in \mathbb{R}^n$ von $F$ zu gelangen. Aus Sicht der Numerik ist
das eine sehr schöne Eigenschaft, da wir nicht gezwungenermaßen $n$
Iterationen des Abstiegsverfahrens {eq}`eq:conjugate_descent_optimal`
durchführen müssen, sondern bereits nach $k < n$ abbrechen können, um
eventuell eine gute Approximation des lokalen Minimums
$x_k \approx x^* \in \mathbb{R}^n$ zu erhalten. Dies spielt insbesondere
bei sehr großen Dimensionen $n >\!\!> 1$ eine wichtige Rolle.

````

````{prf:example} 
Wir wollen im Folgenden ein Beispiel zur Durchführung eines
Abstiegsverfahrens mit gegebenen konjugierten Richtungen angeben. Seien
folgende Werte für das lineare Gleichungssystem $Ax = b$ gegeben:
```{math}
A \ = \
\begin{pmatrix}
3 & 2\\
2 & 6
\end{pmatrix},
\quad b \ = \ 
\begin{pmatrix}
2\\
-8
\end{pmatrix}.
```
Wir nehmen eine Menge von zwei $A$-orthogonalen Vektoren
$d_0, d_1 \in \mathbb{R}^2 / \lbrace 0 \rbrace$ als gegeben an mit:
```{math}
d_0 \ = \
\begin{pmatrix}
0\\
1
\end{pmatrix},
\quad d_1 \ = \ 
\begin{pmatrix}
3\\
-1
\end{pmatrix}.
```
Wir sehen ein, dass die Vektoren $d_0$ und $d_1$ konjugiert bezüglich
der Matrix $A$ sind, denn es gilt:
```{math}
\langle d_0, Ad_1 \rangle \ = \ (0,1)\cdot
\begin{pmatrix}
3 & 2\\
2 & 6
\end{pmatrix}
\begin{pmatrix}
3\\
-1
\end{pmatrix} \ = \ (0,1) \cdot
\begin{pmatrix}
7\\
0
\end{pmatrix} \ = \ 0.
```

Als Startwert für unser Iterationsverfahren wählen wir
$x_0 = (-2, 2)^T$. Für den ersten Schritt des Iterationsverfahren
berechnen wir zuerst das aktuelle Residuum
```{math}
r_0 \ = \ b - Ax_0 \ = \ 
\begin{pmatrix}
2 \\
-8
\end{pmatrix} - 
\begin{pmatrix}
3 & 2\\
2 & 6
\end{pmatrix}
\begin{pmatrix}
-2\\
2
\end{pmatrix}
\ = \
\begin{pmatrix}
2 \\
-8
\end{pmatrix} - 
\begin{pmatrix}
-2 \\
8
\end{pmatrix}
\ = \ 
\begin{pmatrix}
4 \\
-16
\end{pmatrix}.
```
Nun können wir die optimale Schrittweite $\alpha_0 > 0$ für den ersten
Schritt durch den Ausdruck {eq}`eq:optimal_step-size_conjugated`
bestimmen mit:
```{math}
\alpha_0 \ = \ \frac{\langle d_0, r_0\rangle}{\langle d_0, A d_0 \rangle} \ = \
\frac{4}{3}.
```
Hiermit können wir den ersten Abstieg durchführen und erhalten so den
nächsten Iterationspunkt
```{math}
x_1 \ = \ x_0 + \alpha_0 d_0 \ = \ 
\begin{pmatrix}
-2 \\
- 2/3
\end{pmatrix}.
```
Wir wollen nur den zweiten Schritt des Verfahrens angehen und benötigen
wiederum das aktuelle Residuum
```{math}
r_1 \ = \ b - Ax_1 \ = \ 
\begin{pmatrix}
28/3 \\
0
\end{pmatrix}.
```
Wir berechnen wieder die neue optimale Schrittweite mittels
{eq}`eq:optimal_step-size_conjugated`:
```{math}
\alpha_1 \ = \ \frac{\langle d_1, r_1\rangle}{\langle d_1, A d_1 \rangle} \ = \
\frac{28}{21}.
```
Mit dieser können wir den letzten Abstiegsschritt für $n=2$ berechnen
und erhalten somit:
```{math}
x_2 \ = \ x_1 + \alpha_1 d_1 \ = \
\begin{pmatrix}
2\\
-2
\end{pmatrix}
\ = \ x^*.
```

````

Der folgende Satz hilft uns zu verstehen, warum ein Abstiegsverfahren
mit konjugierten Richtungen besser funktioniert als das
Gradientenabstiegsverfahren in {ref}`ss:gradient_descent`.

````{prf:theorem} Orthogonalität des Residuums
:label: thm:residual_orthogonal
Das Residuum $r_{i+1} = b - Ax_{i+1}$ des Abstiegsverfahren mit
konjugierten Richtungen in {eq}`eq:conjugate_descent_optimal` ist
orthogonal zu allen bisherigen Abstiegsrichtungen $d_j, j=0,\ldots,i$,
d.h.
```{math}
\langle r_{i+1}, d_j \rangle \ = \ 0, \quad \text{ für alle } j=0,\ldots,i.
```

````

````{prf:proof} 
Aus {prf:ref}`rem:error_cg` wissen wir, dass wir den Fehler $e_{i+1}$
nach $i$ Iterationen des Abstiegsverfahrens angeben können als
```{math}
e_{i+1} \ = \ \sum_{k=i+1}^{n-1} \delta_k d_k.
```
Wir können beide Seiten der Gleichung mit einem Zeilenvektor
$-d_j^TA \in \R^n$ für einen Index $0 \leq j \leq i$ multiplizieren und
erhalten damit:
```{math}
-\langle d_j, Ae_{i+1} \rangle \ = \ - \sum_{k=i+1}^{n-1} \delta_k \underbrace{d_j Ad_k}_{=~0} \quad 
\Rightarrow \ \langle d_j, r_{i+1} \rangle \ = \ 0 \quad \text{ für alle } \ 0 \leq j \leq i.
```
 ◻

````

Man beachte, dass die Eigenschaft optimal bezüglich **aller vorherigen
Abstiegsrichtungen** nur für den Fall von konjugierten Richtungen
funktioniert und nicht im Fall des Gradientenabstiegsverfahren, wie wir
in {ref}`ss:motivation` gesehen haben. Hier war man nur optimal
bezüglich der **letzten Abstiegsrichtung** und nicht bezüglich aller
vorherigen Richtungen. Das resultiert in dem typischen Zickzack-Pfad
beim Abstieg, wie wir ihn in {numref}`fig:zig-zag` gesehen haben.

(ss:conjugated_gradients)=
## Konjugierte Gradienten

Wir haben in {ref}`ss:conjugated_descent` gesehen, dass wir ein
iteratives Abstiegsverfahren mit konjugierten Abstiegsrichtungen
$\lbrace d_j \rbrace_{j=0,\ldots,n-1}$ verwenden können, um in $n$
Iterationen die eindeutige Lösung des quadratischen Minimierungsproblems
{eq}`eq:quadratic_problem` und somit die Lösung des linearen
Gleichungssystems $Ax = b$ zu erhalten. Bisher sind wir jedoch davon
ausgegangen, dass wir die Menge der konjugierten Vektoren
$\lbrace d_0, \ldots, d_{n-1} \rbrace$ bereits kennen. Um einen
Algorithmus angeben zu können müssen wir also noch ergründen, wie sich
diese Menge mit möglichst geringen numerischen Aufwand finden lässt.

Eine naheliegende Idee wäre es das **Gram-Schmidtsche
Orthogonalisierungsverfahren** so umzugestalten, dass wir eine Menge von
linear unabhängigen Vektoren $\lbrace u_0,\ldots, u_{n-1} \rbrace$ mit
$u_k \in \R^n, k=0,\ldots,n-1$ konjugieren bezüglich der Matrix $A$.
Hierzu würde man die erste Abstiegsrichtung $d_0 \in \R^n$ des
Abstiegsverfahrens mit konjugierten Richtungen als den ersten Vektor der
Menge wählen, d.h., wir setzen $d_0 = u_0$. Anschließend konstruieren
wir die nächste Abstiegsrichtung $d_1$ indem wir alle Komponenten von
$u_1$ entfernen, die nicht $A$-orthogonal zu $d_0$ sind. Für die nächste
Abstiegsrichtung $d_2$ gehen wir analog vor, nur müssen wir darauf
achten alle Komponenten von $u_2$ zu entfernen, die nicht $A$-orthogonal
zu $d_0$ und $d_1$ sind. Dieses Vorgehen lässt sich iterativ bis zum
Vektor $d_{n-1}$ fortführen und man erhält eine Menge von konjugierten
Vektoren $\lbrace d_0, \ldots, d_{n-1} \rbrace$. Diese lassen sich in
geschlossener Form angeben als:
```{math}
:label: eq:directions_gram
d_i \ = \ u_i + \sum_{k=0}^{i-1} \beta_{i,k} d_k, \quad i=1,\ldots,n-1.
```
Wir müssen jedoch die Koeffizienten $\beta_{i,k}$ so bestimmen, dass die
Vektoren $d_i$ konjugiert zu allen vorherigen Richtungsvektoren
$d_j \in \R^n, 0 \leq j < i$ sind. Um diese Koeffizienten zu bestimmen
multiplizieren wir {eq}`eq:directions_gram` wieder von links mit einem
Zeilenvektor $d_j^TA \in \R^n$ für ein
$j \in \lbrace 0,\ldots,i-1\rbrace$ und erhalten
```{math}
\begin{split}
&\langle d_j, A d_i \rangle \ = \ \langle d_j, A u_i \rangle + \sum_{k=0}^{i-1} \beta_{i,k} \langle d_j, A d_k \rangle \\
\Rightarrow \quad & 0 \ = \ \langle d_j, Au_i \rangle + \beta_{i,j}\langle d_j, A d_j \rangle \\\
\Rightarrow \quad & \beta_{i,j} \ = \ - \frac{\langle d_j, Au_i \rangle}{\langle d_j, Ad_j \rangle}.
\end{split}
```
Der Ausdruck für die Koeffizienten $\beta_{i,j}$ ist wohldefiniert, da
wir angenommen haben, dass $A$ eine symmetrische, positiv definite
Matrix ist. Eigentlich könnten wir jetzt zufrieden sein, da wir ein
Verfahren angeben können mit dem sich ein Abstiegsverfahren mit
konjugierten Richtungen konstruieren lässt.

Leider haben wir durch die Verwendung des Gram-Schmidtschen
Orthogonalisierungsverfahrens nichts gewonnen, da der numerische Aufwand
zur Berechnung der unbekannten Koeffizienten im besten Fall in
$\mathcal{O}(n^3)$ liegt, was genau so teuer ist wie eine Invertierung
der Matrix $A$, zum Beispiel mit dem Eliminationsverfahren von Gauss
{cite:p}`numerik1`.

Glücklicherweise gibt es eine Möglichkeit eine Menge von konjugierten
Abstiegsrichtungen im Laufe des Iterationsverfahren
{eq}`eq:conjugate_descent_optimal` zu konstruieren ohne den numerischen
Rechenaufwand des modifizierten Gram-Schmidtschen
Orthogonalisierungsverfahren zu benötigen. Hierzu werden wir ähnlich mit
einer Menge von initialen Richtungsvektoren
$\lbrace u_0, \ldots, u_{n-1}\rbrace$ beginnen und diese geeignet
anpassen. Es stellt sich nämlich heraus, dass in jeder Iteration ein
Vektor existiert, der bereits $A$-orthogonal zu allen vorherigen
Abstiegsrichtungen ist mit Ausnahme der letzten. Diese Aussage wird
durch folgendes Lemma präzisiert.

````{prf:lemma} Eigenschaften des Residuums
:label: lem:residual_direction
Für $i \in \lbrace 1,\ldots,n-1 \rbrace$ befinden wir uns im $(i+1)$-ten
Schritt des Abstiegsverfahrens mit konjugierten Richtungen in
{eq}`eq:conjugate_descent_optimal` und
$\lbrace d_0, \ldots, d_{i-1} \rbrace$ ist eine Menge von
$A$-orthogonalen Vektoren und $u_i = r_i$ sei ein initialer
Richtungsvektor für den aktuellen Iterationsschritt. Dann gilt für
$r_{i+1} = b - A x_{i+1} = b - A (x_i + \alpha_i r_i)$:
```{math}
\langle r_{i+1}, r_j \rangle \ = \ 0, \quad \text{ für alle } \ j=0,\ldots,i,
```
und außerdem auch
```{math}
\langle r_{i+1}, Ad_j \rangle \ = \ 0, \quad \text{ für alle } \ j=0,\ldots,i-1.
```

````

````{prf:proof} 
Wir definieren zuerst die lineare Hülle, die durch die $A$-orthogonalen
Vektoren aufgespannt wird durch:
```{math}
\mathcal{D}_i \ \coloneqq \ \operatorname{span}\lbrace d_0, \ldots, d_{i-1} \rbrace \subset \R^n.
```
Wir gehen in diesem Beweis konstruktiv vor und wollen in jedem Schritt
den unbekannten $A$-orthogonalen Vektor $d_i \in \R^n$ aus dem aktuellen
Residuum $r_i \in \R^n$ und den vorigen Richtungsvektoren
$\lbrace d_0, \dots, d_{i-1}\rbrace$ konstruieren.

Wir wählen als erste Abstiegsrichtung $d_0 = r_0$ und erhalten damit:
```{math}
\operatorname{span}\lbrace r_0 \rbrace \ = \ \operatorname{span}\lbrace d_0 \rbrace \ = \ \mathcal{D}_1.
```
Aus {prf:ref}`thm:residual_orthogonal` wissen wir, dass $r_1$ orthogonal
zu $d_0$ steht und somit folgt auch schon:
```{math}
\langle r_0, r_1 \rangle \ = \ \langle d_0, r_1 \rangle \ = \ 0.
```
Da wir den nächsten $A$-orthogonalen Richtungsvektor $d_1 \in \R^n$ aus
dem aktuellen Residuum $r_1$ und dem Unterraum $\mathcal{D}_1$
konstruieren wollen, können wir damit folgern:
```{math}
\mathcal{D}_2 \ = \ \operatorname{span}\lbrace \mathcal{D}_1, d_1 \rbrace \ = \ \operatorname{span}\lbrace \mathcal{D}_1, r_1 \rbrace \ = \ \operatorname{span}\lbrace r_0, r_1 \rbrace.
```
Analog können wir nun für beliebiges $i \in \lbrace 1,\ldots,n-1\rbrace$
folgern, dass
```{math}
\mathcal{D}_i \ = \ \operatorname{span}\lbrace r_0,\ldots, r_{i-1} \rbrace.
```
Aus {prf:ref}`thm:residual_orthogonal` wissen wir, dass
$r_{i+1} \perp \mathcal{D}_{i+1}$ und damit wissen wir schon, dass die
erste Aussage des Satzes gilt:
```{math}
\langle r_{i+1}, r_j \rangle \ = \ 0, \quad \text{ für alle } j=0,\ldots,i.
```

Für die zweite Aussage des Lemmas drücken wir nun das Residuum $r_i$
durch den Fehler $e_i$ aus und erhalten:
```{math}
r_i \ = \ -Ae_i \ = \ -A(e_{i-1} + \alpha_{i-1}d_{i-1}) \ = \ \underbrace{r_{i-1}}_{\in \mathcal{D}_i} - \underbrace{\alpha_{i-1}Ad_{i-1}}_{\in A\mathcal{D}_i}.
```
Wir sehen also, dass
$r_i \in \operatorname{span}\lbrace r_{i-1}, A d_i \rbrace$. Außerdem
wissen wir durch unsere Folgerungen oben, dass
$r_{i-1} \in \mathcal{D}_i$ und $A d_{i-1} \in A\mathcal{D}_i$ gilt.
Damit gilt aber schon
```{math}
\mathcal{D}_{i+1} \ = \ \operatorname{span}\lbrace \mathcal{D}_i, r_i \rbrace \ = \ \operatorname{span}\lbrace \mathcal{D}_i, A \mathcal{D}_i \rbrace.
```
Wenn wir dies rekursiv entwickeln sehen wir interessanterweise ein, dass
```{math}
\mathcal{D}_i \ = \ \operatorname{span}\lbrace d_0, Ad_0, \ldots, A^{i-1}d_0\rbrace
```
Nach {prf:ref}`thm:residual_orthogonal` wissen wir jedoch auch, dass
$r_{i+1} \perp \mathcal{D}_{i+1}$ und somit muss gelten auch für den
Unterraum gelten $r_{i+1} \perp A\mathcal{D}_i$. Und damit haben wir die
zweite Aussage des Satzes gezeigt, nämlich dass
$r_{i+1} \perp_A \mathcal{D}_i$ oder
```{math}
\langle r_{i+1}, A d_j \rangle, \quad \text{ für alle } j=0,\ldots,i-1.
```
 ◻

````

{prf:ref}`lem:residual_direction` sagt aus, dass das Residuum $r_i$ ein
guter Ausgangspunkt für einen weiteren $A$-orthogonalen Richtungsvektor
$d_i \in \R^n$ im Punkt $x_i$ ist, da es bereits zu allen bisherigen
$A$-orthogonalen Richtungen $d_0,\ldots,d_{i-2}$ konjugiert bezüglich
der Matrix $A$ ist. Wir müssen also nur noch dafür sorgen, dass der
initiale Richtungsvektor $r_i \in \R^n$ $A$-orthogonal zur letzten
Suchrichtung $d_{i-1}$ wird. Dies ist numerisch wesentlich günstiger als
einen beliebig gewählten Richtungsvektor
$u_i \in \mathbb{R}^n / \lbrace 0 \rbrace$ $A$-orthogonal zu machen
bezüglich aller vorigen Richtungsvektoren $d_0,\ldots,d_{i-1}$.

Wir wollen also im Folgenden das vollständige Abstiegsverfahren mit
konjugierten Richtungen angeben. Da die initialen Richtungen nun als
$r_i = -\nabla F(x_i)$ für alle Iterationen $i=0,\ldots, n-1$ gewählt
werden, nennt man dieses Verfahren auch das **Abstiegsverfahren der
konjugierten Gradienten**.

````{prf:theorem} Konjugation der Residuen bezüglich A
Wir befinden uns im $(i+1)$-ten Schritt des Verfahrens der konjugierten
Gradienten für $i=0,\ldots, n-1$ in einem Punkt
$x_{i+1} \in \mathbb{R}^n$ und wir wählen als initialen Richtungsvektor
das aktuelle Residuum $r_{i+1} \in \R^n$, welches nach
{prf:ref}`lem:residual_direction` bereits $A$-orthogonal zu fast allen
vorherigen Richtungsvektoren $\lbrace d_0, \ldots, d_{i-1}\rbrace$ ist.
Indem wir den neuen Richtungsvektor $d_{i+1} \in \R^n$ definieren als
```{math}
:label: eq:conjugated_gradient
d_{i+1} \ \coloneqq \ r_{i+1} + \beta_{i+1} d_i, \qquad \beta_{i+1} \ \coloneqq \ \frac{\langle r_{i+1}, r_{i+1} \rangle}{\langle r_i, r_i \rangle},
```
erhalten wir die Eigenschaft, dass dieser Richtungsvektor nun
$A$-orthogonal zu allen bisherigen Richtungsvektoren des Verfahrens ist,
d.h.,
```{math}
d_{i+1} \: \perp_A \: d_j, \quad j = 0,\ldots,i.
```

````

````{prf:proof} 
Wir müssen den initialen Richtungsvektor
$u_{i+1} \in \mathbb{R}^n / \lbrace 0 \rbrace$ so modifizieren, dass der
resultierende Vektor $d_{i+1}$ konjugiert zu $d_i$ bezüglich der Matrix
$A$ ist. Mit dem modifizierten Gram-Schmidtschen
Orthogonalisierungsverfahren erhalten wir die Form
```{math}
d_{i+1} \ = \ r_{i+1} + \beta_{i+1} d_i, \qquad \beta_{i+1} \ = \ - \frac{\langle r_{i+1}, Ad_i\rangle}{\langle d_i, Ad_i \rangle}.
```
Die neue Abstiegsrichtung $d_{i+1} \in \mathbb{R}^n / \lbrace 0 \rbrace$
ist nach Konstruktion $A$-orthogonal zu allen vorherigen
Abstiegsrichtungen $\lbrace d_0, \ldots, d_i \rbrace$, jedoch wollen wir
den Koeffizienten $\beta_{i+1}$ noch näher charakterisieren im
Folgenden.

Aus dem Beweis von {prf:ref}`lem:residual_direction` wissen wir, dass
wir das aktuelle Residuum ausdrücken können als
```{math}
r_{i+1} \ = \ r_i - \alpha_iAd_i.
```
Wir multiplizieren diese Gleichung von links mit dem Zeilenvektor
$r_{i+1}^T \in \R^n$ und erhalten
```{math}
\langle r_{i+1}, r_{i+1} \rangle \ = \ \langle r_{i+
1}, r_i \rangle - \alpha_i\langle r_{i+1}, Ad_i \rangle.
```
Wir wissen aus {prf:ref}`lem:residual_direction` jedoch auch, dass
$\langle r_{i+1}, r_i \rangle = 0$ gilt und damit erhalten wir den
Ausdruck
```{math}
-\frac{1}{\alpha_i}\langle r_{i+1}, r_{i+1} \rangle \ = \ \langle r_{i+1}, Ad_i \rangle.
```
Wenn wir nun die optimale Schrittweite
$\alpha_i = \frac{\langle r_i, d_i \rangle}{\langle d_i, Ad_i \rangle}$
aus {eq}`eq:conjugate_descent_optimal` einsetzen erhalten wir für den
Koeffizienten $\beta_{i+1}$:
```{math}
\begin{split}
& \ \langle r_{i+1}, Ad_i \rangle  \ = \ -\frac{1}{\alpha_i}\langle r_{i+1}, r_{i+1} \rangle \ = \ -\frac{\langle d_i, Ad_i \rangle}{\langle r_i, d_i \rangle} \langle r_{i+1}, r_{i+1} \rangle \\
\Rightarrow & \ \frac{\langle r_{i+1}, r_{i+1}\rangle }{\langle r_i, d_i \rangle} \ = \ -\frac{\langle r_{i+1}, Ad_i \rangle}{\langle d_i, Ad_i \rangle} \ = \ \beta_{i+1}.
\end{split}
```

Schlussendlich können wir den Nenner in diesem Ausdruck noch weiter
umschreiben, da der letzte Richtungsvektor $d_i \in \R^n$ mit dem
modifizierten Gram-Schmidt Orthogonalisierungsverfahren auch ausgedrückt
werden kann als
```{math}
\langle r_i, d_i \rangle \ = \ \langle r_i, r_i + \beta_i d_{i-1} \rangle \ = \ \langle r_i, r_i \rangle + \beta_i \underbrace{\langle r_i, d_{i-1}\rangle}_{=~0} \ = \ \langle r_i, r_i \rangle.
```
Das Skalarprodukt in obiger Gleichung verschwindet auf Grund von der
Aussage von {prf:ref}`thm:residual_orthogonal` und somit erhalten wir
schlussendlich für den Koeffizienten $\beta_{i+1}$ den simplen Ausdruck:
```{math}
\beta_{i+1} \ = \ \frac{\langle r_{i+1}, r_{i+1}\rangle}{\langle r_i, r_i \rangle}.
```
 ◻

````

Mit der Herleitung von {eq}`eq:conjugated_gradient` können wir nun einen
Algorithmus für das Abstiegsverfahren mit konjugierten Gradienten zum
Lösen eines Gleichungssystems $Ax = b$ angeben.

````{prf:algorithm} Lineares konjugierte Gradientenverfahren
:label: alg:conjugated_gradient

**function** $x^*=$`conjugateGradient`$(A, b, x_0)$  

\# Initialisierung  
$d_0 = r_0 = b - Ax_0$  

\# Führe genau $n$ Schritte durch  
**for** $k = 0,...,n-1$ **do**  
\
    \# Berechne Schrittweite  
    $\alpha_k \ = \ \frac{\langle r_k, d_k\rangle}{\langle d_k, A d_k\rangle}$  
\
    \# Führe Abstiegsschritt durch  
    $x_{k+1} \ = \ x_k + \alpha_k d_k$  
\
    **if** $k < n-1$ **then**  
\
        \# Berechne effizient neues Residuum  
        $r_{k+1} \ = \ r_k - \alpha_k Ad_k$  
\
        \# Berechne Koeffizienten  
        $\beta_{k+1} \ = \ \frac{\langle r_{k+1}, r_{k+1}\rangle}{\langle r_k, r_k \rangle}$  
\
        \# Berechne neue Abstiegsrichung mit Gram-Schmidt  
        $d_{k+1} = r_{k+1} + \beta_{k+1} d_k$  
\
    **end if**  
**end for**  
\
\# Ausgabe des letzten Punktes  
$x^* = x_{k+1}$  

````

(verallgemeinerung-für-nichtlineare-optimierung)=
## Verallgemeinerung für nichtlineare Optimierung

Wir haben festgestellt, dass das konjugierte Gradientenverfahren in
Algorithmus {prf:ref}`alg:conjugated_gradient` als Minimierung der
konvexen quadratischen Funktion
$F(x) \coloneqq \langle x, Ax \rangle - \langle x, b \rangle + c$
konzipiert ist, um das äquivalente lineare Gleichungssystem $Ax = b$ zu
lösen. Es ist natürlich zu fragen, ob wir diesen Ansatz anpassen können,
um allgemeine konvexe Zielfunktionen oder sogar allgemeine nichtlineare
Zielfunktionen $F$ zu minimieren.

Um eine nichtlineare Variante des konjugierte Gradientenverfahrens zu
erhalten, müssen wir einige Anpassungen vornehmen. Fletcher und Reeves
{cite:p}`Fletcher` zeigten, dass man anstelle der explizit berechenbaren
Schrittweite $\alpha_k$ für das quadratische Probleme in
{eq}`eq:optimal_step-size_conjugated` zunächst einmal eine Schrittweite
wählen muss, die ein approximatives Minimum der nichtlinearen
Zielfunktion $F$ entlang der Suchrichtung $d_k \in \R^n$ identifiziert.
Darüber hinaus muss das bisher verwendete Residuum
$r_k = b - Ax_x \in \R^n$ für den Einsatz in der nichtlinearen
Optimierung durch den Gradienten der Zielfunktion $F$ ersetzt werden.

Diese beiden Änderungen führen zu dem folgenden Algorithmus für die
nichtlineare Optimierung.

````{prf:algorithm} Nichtlineares konjugierte Gradientenverfahren
:label: alg:nonlinear_conjugated_gradient

**function** $x^*=$`nonlinearConjugateGradient`$(F, \nabla F, x_0, \alpha_0, \sigma)$  

\# Initialisierung  
$d_0 = - \nabla F(x_0)$  

**while** $||\nabla F(x_k)||^2 > \epsilon$ **do**  
\
    **while** $F(x_k + \alpha_k d_k) > F(x_k)$ **do**  
        \# VerringereSchrittweite um Faktor $\sigma$  
        $\alpha_k = \sigma\alpha_k$  
    **end while**  
\
    \# Führe Abstiegsschritt durch  
    $x_{k+1} \ = \ x_k + \alpha_k d_k$  
\
    \# Berechne Koeffizienten  
    $\beta_{k+1} \ = \ \frac{\langle \nabla F(x_{k+1}), \nabla F(x_{k+1})\rangle}{\langle \nabla F(x_{k}), \nabla F(x_{k}) \rangle}$  
\
    \# Berechne neue Abstiegsrichung mit Gram-Schmidt  
    $d_{k+1} = -\nabla F(x_{k+1}) + \beta_{k+1} d_k$  
\
**end while**

\# Ausgabe des letzten Punktes  
$x^* = x_{k+1}$

````

Man beachte, dass man für die Implementierung von Algorithmus
{prf:ref}`alg:nonlinear_conjugated_gradient` keinerlei Matrix-Vektor
Multiplikationen benötigt und man lediglich Auswertungen der
Zielfunktion $F$ und ihres Gradienten $\nabla F$ braucht. Für den Fall,
dass es sich bei der Zielfunktion $F$ um eine strikt konvexe,
quadratische Funktion handelt und wir in jedem Schritt die optimale
Schrittweite $\alpha_k > 0$ bestimmten können, so entspricht Algorithmus
{prf:ref}`alg:nonlinear_conjugated_gradient` dem linearen konjugierte
Gradientenverfahren in Algorithmus {prf:ref}`alg:conjugated_gradient`.

In der Literatur hat sich unter Anderem eine Modikation des
nichtlinearen konjugierte Gradientenverfahrens von Fletcher und Reeves
etabliert, die sich in numerischen Experimenten häufig als robuster und
effizienter herausgestellt hat.

````{prf:remark} Polak-Ribière Variante
Bei der sogenannten Polak-Ribière Variante des Algorithmus
{prf:ref}`alg:nonlinear_conjugated_gradient` unterscheidet sich
hauptsächlich die Berechnung des Parameters $\beta_{k+1} \in \R$ in
{eq}`eq:conjugated_gradient` für die Anpassung der nächsten
Abstiegsrichtung. Hierbei wird der Faktor $\beta_{k+1}$ nämlich wie
folgt berechnet
```{math}
\beta^{PR}_{k+1} \ \coloneqq \ \frac{\langle \nabla F(x_{k+1}), \nabla F(x_{k+1}) - \nabla F(x_{k}) \rangle}{\langle \nabla F(x_{k}), \nabla F(x_{k}) \rangle}.
```
Man sieht ein, dass im Falle einer strikt konvexen, quadratischen
Zielfunktion $F$ mit optimaler Schrittweitenwahl für die $\alpha_k > 0$
die Orthogonalitätsbedingung für sukzessive Gradienten aus
{prf:ref}`lem:optimale_schrittweite` gilt mit
```{math}
\langle \nabla F(x_{k+1}), \nabla F(x_k) \rangle \ = \ 0.
```
In diesem Fall stimmt der Faktor $\beta^{PR}_{k+1}$ mit dem Faktor
$\beta_{k+1}$ aus Algorithmus
{prf:ref}`alg:nonlinear_conjugated_gradient` überein. In allen anderen
Fällen hingegen unterscheiden sich die Faktoren im Allgemeinen und
führen so zu euben signifikant unterschiedlichen Konvergenzverhalten der
jeweiligen Abstiegsverfahren.

````
