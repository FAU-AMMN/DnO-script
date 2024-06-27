(s:nichtdiffbare_optimierung)=
Nicht-differenzierbare Optimierung
===

Im Folgenden widmen wir uns abschließend der Optimierung von
nicht-differenzierbaren Zielfunktionen, die man häufig in der
Datenanalyse und der mathematischen Bildverarbeitung findet. Hierzu
betrachen wir zunächst das folgende motivierende Beispiel.

````{prf:example} LASSO-Problem
:label: ex:lasso
Wir betrachten das sogenannte LASSO-Problem (engl. für *Least Absolute
Shrinkage and Selection Operator*), das man als Variante eines linearen
Ausgleichsproblem interpretieren kann (siehe {cite:p}`numerik1`) mit
```{math}
F(x) \ = \ \frac{1}2 \Vert b - A x\Vert^2 + \alpha \Vert x \Vert_{\ell^1}.
```
Hierbei ist $A \in \R^{m \times n}, x \in \R^n$ und $b \in \R^m$, wobei
typischerweise $m < n$ gilt.

Mit der Lösung dieses Problems für $\alpha \rightarrow 0$ approximiert
man eine Lösung des linearen Gleichungssystems $Ax=b$ mit minimaler
$\ell^1$-Norm.

Solche Lösungen besitzen in der Regel viele Nulleinträge und sind daher
besonders interessant. Das Problem bei der Minimierung der Zielfunktion
$F$ ist, dass die $\ell^1$-Norm
```{math}
\Vert  x \Vert_{\ell^1} \ = \ \sum_{j=1}^n \vert x_j \vert
```
nicht differenzierbar in allen Punkten $x \in \R^n$ ist, für die
mindestens eine Koordinate Null ist, d.h., für die $x_j = 0$ gilt für
ein $1 \leq j \leq n$.

````

Eine interessante Klasse von Zielfunktion, die wir im Folgenden
betrachten wollen, ist von der Form
```{math}
F(x) = G(x) + H(x),
```
mit einer stetig differenzierbaren Funktion
$G \colon \R^n \rightarrow \R$ und einer konvexen, nicht
notwendigerweise differenzierbaren Funktion
$H \colon \R^n \rightarrow \R$. In diesem Fall können wir keins der in
den vorangegangenen Kapiteln vorgestellten Gradienten-basierten
Optimierungsverfahren verwenden, sondern benötigen einen anderen
methodischen Ansatz.

Zunächst benötigen wir aber noch einige Definition um konvexe Funktionen
analytisch besser beschreiben zu können.

````{prf:definition} Subdifferential und Subgradient
:label: def:subdifferential
Sei $H: \R^n \rightarrow \R$ eine konvexe Funktion. Dann ist das
**Subdifferential** an der Stelle $x \in \R^n$ definiert durch
```{math}
\partial H(x) \ \coloneqq \{ p \in \R^n ~|~ H(x) + \langle p, y-x \rangle \leq H(y) \quad \forall~y \in \R^n \}.
```
Ein Element des Subdifferentials nennen wir **Subgradient**.

````

Man sieht ein, dass ein globales Minimum der konvexen Funktion $H$ durch
die Bedingung $\vec{0}  \in \partial H(x^*)$ charakterisiert wird, denn
mit {prf:ref}`def:subdifferential` erhalten wir in diesem Fall:
```{math}
:label: eq:optimalitaetsbedingung_h
\vec{0} \in \partial H(x^*) \ \Leftrightarrow \ H(x^*) + \underbrace{\langle \vec{0} , y-x^* \rangle}_{=\, 0} \leq H(y) \quad \forall y \in \R^n.
```

Folgendes Beispiel illustriert das Subdifferential einer nicht
differenzierbaren, eindimensionalen Funktion.

````{prf:example} Subdifferential der Betragsfunktion
:label: ex:subdifferential_betrag
Wir betrachten die eindimensionale Betragsfunktion
$H(x) \coloneqq \vert x \vert$. In diesem Fall ist das Subdifferential
von $H$ für alle Punkte $x \neq 0$ gegeben als
```{math}
\partial H(x) \ = \ \{\operatorname{sgn}(x)\}.
```
Im Punkt $x=0$, in dem die Betragsfunktion bekanntlich nicht
differenzierbar ist, gilt für das Subdifferential hingegen
```{math}
\partial H(0) \ = \ [-1,1].
```

````

Basierend auf dieser Beobachtung im Eindimensionalen können wir auch das
Subdifferential der $\ell^1$-Norm im $\R^n$ im folgenden Beispiel
ausrechnen.

````{prf:example} Subdifferential der \\ell\^1-Norm
Wir betrachten die $\ell^1$-Norm, die durch die konvexe Funktion
$H(x) = \Vert x \Vert_{\ell^1}$ gegeben ist. Dann können wir das
Subdifferential der $\ell^1$-Norm für alle Punkte $x \in \R^n$ angeben
als
```{math}
\partial H(x) \ = \ \{ p \in [-1,1]^n ~|~ p_i = \operatorname{sgn}(x_i) \text{ für } x_i \neq 0\}.
```

````

Basierend auf der Optimalitätsbedingung für die konvexe Funktion $H$ in
{eq}`eq:optimalitaetsbedingung_h` können wir eine entsprechende
notwendige Bedingung für ein lokales Minimum für die Summe einer
differenzierbaren Funktion $G$ und einer konvexen Funktion $H$
herleiten, wie das folgende Lemma besagt.

````{prf:lemma} Optimalitätsbedingung
:label: lem:subdifferential_notwendige_bedingung
Sei $F \colon \R^n \rightarrow \R$ eine Zielfunktion, die sich schreiben
lässt als $F=G+H$, wobei $G$ stetig differenzierbar und $H$ konvex ist.
Sei außerdem $x^* \in \R^n$ ein lokaler Minimierer von $F$.

Dann gilt $\vec{0} \in \nabla G(x^*) + \partial H(x^*)$, d.h., es
existiert ein Subgradient $p \in \partial H(x^*)$ mit
$\nabla G(x^*) + p = \vec{0}$.

````

````{prf:proof} 
Sei $x^* \in \R^n$ ein lokaler Minimierer von $F$. Dann gilt für jeden
Punkt $x \in \R^n$ in einer lokalen Umgebung des Minimierers $x^*$ mit
einer Taylor-Entwicklung erster Ordnung der Funktion $G$ schon
```{math}
\begin{split}
G(x^*)  + H(x^*) \ &= \ F(x^*) \ \leq \ F(x) \ = \ G(x) + H(x)\\
&= \ G(x^*) + \langle \nabla G(x^*),x-x^* \rangle + r(x)\Vert x-x^*\Vert + H(x)
\end{split}
```
mit einem Fehlerterm $r \colon \R^n \rightarrow \R$ für den gilt
$r(x) \rightarrow 0$ wenn $x \rightarrow x^*$. Wählen wir nun den
speziellen Vektor $p = -\nabla G(x^*)$ und stellen die Gleichung um, so
gilt also
```{math}
\frac{H(x^*) + \langle p, x - x^* \rangle - H(x)}{\Vert x-x^* \Vert} \ = \ r(x).
```
Dann können wir unter Ausnutzung der Konvexität von $H$ folgende
Abschätzung treffen
```{math}
\begin{split}
0 \ = \ \limsup_{x \rightarrow x^*} r(x) \ &\geq \ \limsup_{x \rightarrow x^*} \frac{H(x^*) + \langle p, x - x^* \rangle - H(x)}{\Vert x-x^* \Vert} \\
&\geq \ \lim_{\lambda \rightarrow 0} \frac{H(x^*) + \langle p, \lambda(x - x^*) \rangle - H((1-\lambda)x^* + \lambda x)}{\lambda \Vert x-x^* \Vert}\\
&\geq \ \lim_{\lambda \rightarrow 0} \frac{H(x^*) + \lambda \langle p, x - x^* \rangle - (1-\lambda)H(x^*) - \lambda H(x)}{\lambda \Vert x-x^* \Vert}\\
&= \ \lim_{\lambda \rightarrow 0} \frac{\langle p, x - x^* \rangle + H(x^*) - H(x)}{\Vert x-x^* \Vert} \ = \ \frac{\langle p, x - x^* \rangle + H(x^*) - H(x)}{\Vert x-x^* \Vert}.
\end{split}
```
Durch Umstellen erhalten wir schließlich
```{math}
H(x^*) + \langle p, x - x^* \rangle \ \leq \ H(x).
```
Dies bedeutet aber schon nach {prf:ref}`def:subdifferential`, dass
$p \in \partial H(x^*)$ gilt. ◻

````

Basierend auf {prf:ref}`lem:subdifferential_notwendige_bedingung`
könnten wir also von der notwendigen Optimalitätsbedingung
$\nabla G(x^*) + p = 0$ ausgehen und damit ein Analogon zum
Gradientenabstiegsverfahren aus {ref}`ss:gradient_descent` aufstellen
mit
```{math}
x_{k+1} = x_k - \alpha_k( \nabla G(x_k) + p_k), \qquad p_k \in \partial H(x_k).
```
Allerdings ist nicht klar welchen Subgradienten
$p_k \in \partial H(x_k)$ wir in jeder Iteration auswählen sollen (falls
mehr als einer existiert) oder ob überhaupt ein Subgradient existiert.
Deshalb werden wir im Folgenden eine Variante betrachten bei der
automatisch Iterationsschritte $x_{k+1} \in \R^n$ mit nichtleerem
Subdifferential und ebenso ein $p_{k+1} \in \partial H(x_{k+1})$
ausgewählt werden.


(ss:proximal_splitting)=
# Proximales Splitting

Die Idee des sogenannten proximalen Splitting, auch *Forward-Backward
Splitting* genannt, ist es den differenzierbaren Teil genauso wie beim
Gradientenabstiegsverfahren auszuwerten, d.h., vorwärts ausgehend vom
Punkt $x_k \in \R^n$, während der Subgradient hingegen bezüglich der
nächsten Iterierten ausgewertet wird, d.h., rückwärts ausgehend vom
Punkt $x_{k+1} \in \R^n$.

Somit lässt sich die Iterationsvorschrift für das proximale Splitting
schreiben als
```{math}
:label: eq:iteration_proximal_splitting
x_{k+1} \ = \ x_k - \alpha_k( \nabla G(x_k) + p_{k+1}), \qquad p_{k+1} \in \partial H(x_{k+1}).
```
Diese Gleichung können wir nun selbstverständlich nicht mehr explizit
auswerten, da der Subgradient $p_{k+1} \in \partial H(x_{k+1})$ vom
bisher unbestimmten Iterationsschritt $x_{k+1} \in \R^n$ abhängt.

Dennoch können wir die Iterationsvorschrift
{eq}`eq:iteration_proximal_splitting` weiter umschreiben zu
```{math}
\frac{1}{\alpha_k} (x_{k+1} - x_k) + \nabla G(x_k) + p_{k+1} \ = \ 0.
```
Die zentrale Idee ist es nun eine Funktion zu finden, deren Ableitung
den Ausdruck auf linken Seite der obigen Gleichung liefert. Man sieht
ein, dass die obige Gleichung die hinreichende Optimalitätsbedingung für
die Minimierung der folgenden strikt konvexen Funktion darstellt
```{math}
F_k(x) \ = \ \frac{1}{2 \alpha_k} \Vert x - x_k + \alpha_k \nabla G(x_k) \Vert^2 + H(x).
```
Hierbei stellt die Iterierte $x_{k+1} \in \R^n$ also einen stationären
Punkt der Zielfunktion $F_k(x)$ dar, während $p_{k+1} \in \R^n$ aus dem
Subdifferential $\partial H(x)$ stammt.

Wir können das Iterationsverfahren in
{eq}`eq:iteration_proximal_splitting` also durchführen, wenn wir strikt
konvexe Funktionen der Form
```{math}
\Phi(x) \ \coloneqq \ \frac{1}{2\alpha} \Vert x - y \Vert^2 + H(x)
```
effizient minimieren können.

Bei der Optimierung der Funktion $\Phi$ versucht man anschaulich einen
Kompromiss einzugehen zwischen dem Ziel die Funktion $H$ zu minimieren
und gleichzeitig nahe dem Punkt $y \in \R^n$ bezüglich der Euklidischen
Norm zu sein. Der Parameter $\alpha > 0$ steuert hierbei die Gewichtung
zwischen diesen beiden Kriterien.

Da es sich bei $\Phi$ um eine strikt konvexe Funktion handelt existiert
ein eindeutiges Minimum. Dieses lässt sich durch den sogenannten
Proximaloperator beschreiben, den wir im Folgenden definieren wollen.

````{prf:definition} Proximaloperator
Sei $H \colon \R^n \rightarrow \R$ eine konvexe, unterhalbstetige
Funktion und sei $\alpha > 0$ ein positiver Parameter. Dann definieren
wir den **Proximaloperator** bezüglich der Funktion $H$ im Punkt
$y\in \R^n$ als
```{math}
\operatorname{prox}_{\alpha H}(y) \ \coloneqq \ \underset{x \in \R^n}{\operatorname{arg\,min}} \left \lbrace \frac{1}{2\alpha} ||x - y||^2 + H(x) \right\rbrace.
```

````

Wir wollen zunächst das Konzept des Proximaloperators im folgenden
Beispiel für die Betragsfunktion bzw. die $\ell^1$-Norm genauer
verstehen.

````{prf:example} Shrinkage-Operator
:label: ex:shrinkage-operator
Sei $H(x) = \vert x \vert$ in diesem Beispiel zunächst die
Betragsfunktion. Dann ist der Proximaloperator
$z =\operatorname{prox}_{\alpha H}(y)$ der eindeutige Minimierer der
strikt konvexen Funktion
```{math}
\Phi(x) \ = \ \frac{1}{2\alpha} (x - y)^2 + \vert x \vert
```
mit der notwendigen Optimalitätsbedingung 1. Ordnung
```{math}
:label: eq:shrinkage_bedingung
\frac{1}\alpha (y-z) \in \partial \vert z\vert.
```
Aus {prf:ref}`ex:subdifferential_betrag` wissen wir bereits, dass
$\partial |z| = \lbrace{ \operatorname{sgn}(z)}\rbrace = \lbrace -1, 1 \rbrace$
für $z \neq 0$ gilt und $\partial |0| = [-1, 1]$ in der Null.

Betrachten wir zunächst den Fall $z > 0$. Dann ist klar, dass
$\partial |z| = 1$ gilt und somit wird aus der notwendigen
Optimalitätsbedingung {eq}`eq:shrinkage_bedingung`
```{math}
\frac{1}\alpha (y-z) \: = \: 1 \quad \Leftrightarrow \quad  z \: = \: y - \alpha.
```
Da wir $z > 0$ angenommen haben, ist dies nur möglich für
$y > \alpha > 0$.

Analog gilt für den Fall $z < 0$, dass $\partial |z| = -1$ ist und somit
erhalten wir aus der notwendigen Optimalitätsbedingung
{eq}`eq:shrinkage_bedingung`
```{math}
\frac{1}\alpha (y-z) \: = \: -1 \quad \Leftrightarrow \quad  z \: = \: \alpha + y.
```
Dies ist aber für $z < 0$ nur möglich, wenn $y < -\alpha < 0$ gilt.

Für den letzten Fall mit $z=0$ lauten die notwendigen
Optimalitätsbedingung {eq}`eq:shrinkage_bedingung`
```{math}
\frac{1}\alpha (y-z) \: = \: \frac{y}{\alpha} \: \in \: [-1, 1]  = \partial |0|.
```
Dies ist äquivalent zu der Bedingung $-1 \leq \frac{y}{\alpha} \leq 1$,
die nur dann erfüllt ist wenn $-\alpha \leq y \leq \alpha$ gilt.

Somit lässt sich der Proximaloperator
$z =\operatorname{prox}_{\alpha H}(y)$ für die Betragsfunktion
$H(x) \coloneqq \vert x \vert$, der auch **Shrinkage-Operator** genannt
wird, wie folgt angeben:
```{math}
\operatorname{prox}_{\alpha \vert \cdot \vert}(y) \ = \ \text{shrink}_\alpha(y) \ \coloneqq \ \left\{ \begin{array}{ll} y - \alpha \quad & \text{falls } y > \alpha \\  0 \quad & \text{falls } y \in [-\alpha, \alpha]  \\ y + \alpha \quad & \text{falls } y < - \alpha. \end{array}\right.
```
Die kontrahierende Wirkung des Shrinkage-Operators wird in
{numref}`fig:shrinkage` illustriert.

Analog können wir auch den Proximaloperator für die $\ell^1$-Norm
$H(x) = \Vert x \Vert_{\ell^1}$ im $\R^n$ angeben als
```{math}
\text{prox}_{\alpha H}(y) = ( \operatorname{shrink}_\alpha(y_i))_{i=1,\ldots,n}.
```

````

Mit Hilfe des Proximaloperators können wir das proximale Splitting in
{eq}`eq:iteration_proximal_splitting` schreiben als
```{math}
x_{k+1} \ = \ \operatorname{prox}_{\alpha_k H}(x_k - \alpha_k  \nabla G(x_k)).
```

Wie wir in {prf:ref}`ex:shrinkage-operator` gesehen haben ist der
Shrinkage-Operator in gewisser Weise kontraktiv. Diese Beobachtung gilt
im Allgemeinen für Proximaloperatoren wie folgendes Lemma feststellt.

````{prf:lemma} Kontraktivität des Proximaloperators
Sei $H$ eine konvexe, unterhalbstetige Funktion. Dann ist der
Proximaloperator $\operatorname{prox}_H$ Lipschitz stetig mit
Lipschitz-Konstante $1$.

````

````{prf:proof} 
Wir betrachten zunächst die beiden Punkte, die für $i = 1,2$ gegeben
sind durch $x_i= \operatorname{prox}_H(y_i)$. Dann gilt wegen der
Optimalitätsbedingung des Proximaloperators $x_i + p_i =  y_i,$ für
einen Subgradienten $p_i \in \partial H(x_i)$. Subtrahieren wir die
beiden Identitäten für $i=1,2$, so erhalten wir
```{math}
x_1 - x_2 + p_1 - p_2 \ = \ y_1 - y_2.
```
Multiplizieren wir diese Gleichung von links mit dem Zeilenvektor
$(x_1-x_2)^T \in \R^n$ so gilt mit Hilfe der *Cauchy-Schwarz
Ungleichung*:
```{math}
:label: eq:lipschitz_cs-ungleichung
\Vert x_1 - x_2 \Vert^2 + \langle x_1-x_2, p_1 - p_2 \rangle = \langle x_1-x_2, y_1 -y_2 \rangle \ \leq \ \Vert y_1-y_2 \Vert \cdot \Vert x_1-x_2 \Vert.
```

Wegen {prf:ref}`def:subdifferential` des Subgradienten können wir
festhalten, dass gilt
```{math}
\begin{aligned}
\langle p_1, x_1- x_2 \rangle \ &\geq \ H(x_1) - H(x_2) \\
\langle p_2, x_2- x_1 \rangle \ &\geq \ H(x_2) - H(x_1).
\end{aligned}
```
Addieren wir diese beiden Ungleichungen, so erhalten wir insgesamt
$\langle p_1 - p_2, x_1-x_2 \rangle \geq 0$. Damit können wir die linke
Seite der Ungleichung {eq}`eq:lipschitz_cs-ungleichung` weiter
abschätzen und erhalten durch Teilen mit dem Wert $||x_1 - x_2||$ auf
beiden Seiten schließlich
```{math}
\Vert x_1 - x_2 \Vert \ \leq \ \Vert  y_1 - y_2 \Vert ,
```
was genau die gewünschte Lipschitz-Stetigkeit des Proximaloperators mit
Lipschitz-Konstante $L=1$ bedeutet. ◻

````

Analog zum Beweis von {prf:ref}`thm:konvergenz_abstieg` können wir auch
vorgehen um die Konvergenz des proximalen Splitting Verfahrens im
folgenden Theorem zu zeigen.

````{prf:theorem} Konvergenz des proximalen Splitting-Verfahrens
:label: thm:konvergenz_splitting
Sei $F \coloneqq \R^n \rightarrow \R$ eine nach unten beschränkte
Zielfunktion, das heißt, dass die Niveaumenge
$K \coloneqq \lbrace x \in \R^n : F(x) \leq F(x_0) \rbrace$ beschränkt
ist. Außerdem lasse sich $F$ als Summe zweier Funktionen schreiben mit
$F(x) = G(x) + H(x)$ für eine zweimal stetig differenzierbare Funktion
$G \colon \R^n \rightarrow \R$ und eine konvexe, unterhalbstetige
Funktion $H \colon \R^n \rightarrow \R$.

Dann konvergiert das proximale Splitting-Verfahren der Form
```{math}
x_{k+1} \ = \ \operatorname{prox}_{\alpha_k H}(x_k - \alpha_k \nabla G(x_k)).
```

````

````{prf:proof} 
Wir nutzen zunächst für einen Subgradienten
$p_{k+1} \in \partial H(x_{k+1})$ die explizite Iterationsvorschrift
{eq}`eq:iteration_proximal_splitting` für das proximale Splitting
```{math}
\frac{1}{\alpha_k} (x_{k+1} - x_k) + \nabla G(x_k) + p_{k+1} \ = \ 0.
```

Multiplizieren wir diese Gleichung von links mit dem Zeilenvektor
$(x_{k+1} - x_k)^T \in \R^n$, dann folgt
```{math}
:label: eq:splitting_skalarprodukt
\frac{1}{\alpha_k} \Vert x_{k+1} - x_k \Vert^2 + \langle x_{k+1} - x_k, \nabla G (x_k) \rangle + \langle x_{k+1} - x_k, p_{k+1} \rangle \ = \ 0.
```
Da die Funktion $G$ nach Voraussetzung zweimal stetig differenzierbar
ist folgt mit einer Taylorapproximation erster Ordnung
```{math}
\langle x_{k+1} - x_k, \nabla G (x_k) \rangle \ = \  G(x_{k+1}) - G(x_k) -r_k.
```
Hierbei lässt sich der Fehlerterm $r_k$ abschätzen durch
```{math}
r_k \ \leq \ \frac{C_k}2 \Vert x_{k+1} - x_k \Vert^2,
```
wobei $C_k$ eine obere Schranke für die Norm der Hesse-Matrix von $G$ in
einer Umgebung von $x_k$ mit Radius $\Vert x_{k+1} - x_k \Vert$ ist.

Aus {prf:ref}`def:subdifferential` folgt für den Subgradienten
$p_{k+1} \in \partial H(x_{k+1})$ außerdem
```{math}
\langle p_{k+1}, x_{k+1} - x_k \rangle \ \geq \ H(x_{k+1}) - H(x_k).
```

Nutzen wir diese Abschätzungen nun für die Identität
{eq}`eq:splitting_skalarprodukt`, so erhalten wir insgesamt
```{math}
\begin{split}
0 \ &= \ \frac{1}{\alpha_k} \Vert x_{k+1} - x_k \Vert^2 + \langle x_{k+1} - x_k, \nabla G (x_k) \rangle + \langle x_{k+1} - x_k, p_{k+1} \rangle \\
&\geq \ \frac{1}{\alpha_k} \Vert x_{k+1} - x_k \Vert^2 + G(x_{k+1}) - G(x_k) - \frac{C_k}2 \Vert x_{k+1} - x_k \Vert^2 + H(x_{k+1}) - H(x_k) \\
&= \ \left(\frac{1}{\alpha_k}-\frac{C_k}{2}\right) \Vert x_{k+1} - x_k \Vert^2 + F(x_{k+1}) - F(x_k)
\end{split}
```
Theoretisch können wir die Schrittweiten $\alpha_k > 0$ so klein wählen,
dass $\frac{1}{\alpha_k}-\frac{C_k}{2} > \epsilon$ für ein beliebiges
$\epsilon > 0$ gilt. Durch Umstellen sehen wir also, dass gilt
```{math}
F(x_k) \ \geq \ \left(\frac{1}{\alpha_k}-\frac{C_k}{2}\right) \Vert x_{k+1} - x_k \Vert^2 + F(x_{k+1})  \ > \ \epsilon \cdot \Vert x_{k+1} - x_k \Vert^2 + F(x_{k+1}).
```
Damit ist offensichtlich, dass $F(x_{k+1}) < F(x_k)$ für alle $k \in \N$
gilt. Damit haben wir gezeigt, dass das proximale Splitting-Verfahren
die Zielfunktion $F$ in jedem Schritt verkleinert.

Analog zum Beweis von {prf:ref}`thm:konvergenz_abstieg` sehen wir ein,
dass für die Folge der Iterationsschritte $(x_k)_{k \in \N} \subset K$
gilt und nach dem *Satz von Bolzano-Weierstrass* eine konvergente
Teilfolge besitzt, deren Grenzwert in der kompakten Menge $K$ liegen
muss.

Es gilt ebenfalls
```{math}
\sum_{k=0}^\infty  \Vert x_{k+1} - x_k \Vert^2 \: < \: \infty.
```
Da die Iterierten auf einer beschränkten Menge bleiben und $G$ zweimal
stetig differenzierbar ist, ist die Norm der Hesse-Matrix von $G$
ebenfalls uniform für alle $k \in \N$ beschränkt für eine Konstante
$C > 0$ mit $0 < C_k < C$ und somit kann die Folge der Schrittweiten
$(\alpha_k)_{k\in\N}$ ebenfalls uniform nach unten beschränkt werden.

Aus der Konvergenz
```{math}
0 \ = \ \lim_{k\rightarrow \infty} \frac{1}{\alpha_k}( x_{k+1} - x_k) \ = \ \lim_{k\rightarrow \infty} \nabla G(x_k) + p_{k+1}, \qquad p_{k+1} \in \partial H(x_{k+1})
```
können wir schließlich folgern, dass jeder Häufungspunkt die
Optimalitätsbedingung aus
{prf:ref}`lem:subdifferential_notwendige_bedingung` erfüllt. ◻

````

(primal-duale-verfahren)=
## Primal-Duale Verfahren

In vielen Fällen lässt sich der Proximaloperator einer beliebigen
konvexen, unterhalbstetigen Funktion $H$ analytisch nicht gut berechnen
und muss gegebenenfalls numerisch approximiert werden. In solchen Fällen
ist das proximale Splitting-Verfahren aus {ref}`ss:proximal_splitting`
häufig weniger effizient anzuwenden.

In vielen interessanten Fällen hat die Funktion $H$ die Form
$H(x) = J(Bx)$, mit einer Matrix $B \in \R^{n\times n}$, die nicht
diagonal ist, und einer konvexen, unterhalbstetigen Funktion $J$, deren
Proximaloperator man analytisch gut berechnen kann. Ein typisches
Beispiel für solch ein Problem ist eine Variante des Lasso-Problems aus
{prf:ref}`ex:lasso`, in dem man folgende Optimierung durchführen will
```{math}
\min_{x \in \R^n} \left\lbrace F(x) \: = \: \frac{1}2 \Vert A x - b\Vert^2 + \alpha \Vert B x \Vert_{\ell^1} \right\rbrace.
```

Obwohl für den Spezialfall der Identität $B = I$ bereits in
{prf:ref}`ex:shrinkage-operator` den Proximaloperator analytisch angeben
konnten, lässt sich dieser für $H(x) = \Vert B x \Vert_{\ell^1}$ im Fall
beliebiger Matrizen $B$ im Allgemeinen nicht explizit angeben.

Um dieses Problem zu umgehen, ist eine Idee, eine Nebenbedingung mit
einer zusätzlichen Variable einzuführen, die wir im Folgenden mit
$y \in \R^n$ bezeichnen. Setzen wir nämlich $y=Bx$, so können wir eine
modifizierte Zielfunktion
```{math}
\tilde F(x,y) \ = \ G(x) + J(y)
```
unter der Nebenbedingung $y=Bx$ minimieren. Um eine Nebenbedingung
dieser Form bei der Minimierung von $\tilde F$ einfach zu
berücksichtigen, können wir die Idee der **Lagrange Multiplikatoren**
verwenden (siehe {cite:p}`tenbrinck2021`). Hierzu definieren wir
zunächst das folgende Lagrange-Funktional
```{math}
L(x,y,z) \ \coloneqq \ \tilde F(x,y) + \langle z, Bx-y \rangle,
```
wobei der Vektor $z \in \R^n$ die Lagrange Multiplikatoren für die
Nebenbedingung $Bx = y$ darstellt.

Man sieht nun ein, dass gilt
```{math}
\inf_{\substack{x,y \in \R^n\\Bx=y}} \tilde F(x,y) \ = \ \inf_{x,y \in \R^n} \sup_{z \in \R^n} L(x,y,z).
```
Dies liegt daran, dass das Supremum über $L(x,y,z)$ unendlich wird,
sobald $Bx \neq y$ ist, das heißt, dass das Infimum von $\tilde F$ in
solchen Fällen sicher nicht angenähert wird. Daher reicht es das
Lagrange-Funktional $L$ auf der Menge
$\lbrace (x,y) \in \R^n \times \R^n : Bx=y \rbrace$ zu betrachten. Auf
dieser Menge wird der zweite Term von $L$ jedoch Null und somit stimmt
das Lagrange-Funktional schon mit der modifizierten Zielfunktion
überein, d.h. es gilt $L(x,y,z) = \tilde F(x,y)$.

Wir sehen also, dass wir im Fall der Optimierung mit Nebenbedingungen
eigentlich ein **Sattelpunktproblem** der folgenden Form lösen wollen
```{math}
\inf_{x,y \in \R^n} \sup_{z \in \R^n} (G(x) + J(y) + \langle z, Bx - y\rangle ).
```

Nehmen wir zunächst an, dass wir das Infimum und das Supremum
vertauschen können (was möglich ist, falls ein Sattelpunkt existiert)
und nutzen die Relationen $\inf_{x,y} - E(x,y) = - \sup_{x,y} E(x,y)$
und $\sup_{z} -F(z) = - \inf_{z} F(z)$ dann lösen wir
```{math}
\begin{aligned}
&\quad \ \sup_{z \in \R^n} \inf_{x,y \in \R^n} (G(x) + J(y) + \langle z, Bx - y\rangle) \\
&= \ \sup_{z \in \R^n} \inf_{x,y \in \R^n} -(\langle z,y \rangle - J(y) + \langle -B^Tz, x \rangle - G(x) )\\
&= \ - \inf_{z \in \R^n}  \sup_{x,y \in \R^n} \ (\langle z,y \rangle - J(y) + \langle -B^Tz, x \rangle - G(x) ) \\
&= \ - \inf_{z \in \R^n} ( \sup_{y \in \R^n} (\langle z,y \rangle - J(y)) + \sup_{x \in \R^n} (\langle -B^Tz, x \rangle - G(x)) ) .
\end{aligned}
```
Bevor wir eine sogenannte primal-duale Formulierungen des ursprünglichen
Problems herleiten, werden wir zunächst die Konvex-Konjugierte
definieren.

````{prf:definition} Legendre-Fenchel-Transformation
:label: def:legendre-transformation
Sei $J \colon \R^n \rightarrow \R$ eine strikt konvexe Funktion. Dann
ist die **Legendre-Fenchel-Transformation** $J^*$ von $J$ (auch
Konvex-Konjugierte genannt) definiert als
```{math}
J^*(z) \ \coloneqq \ \sup_{y\in\R^n} (\langle z,y \rangle - J(y)).
```

````

Mit {prf:ref}`def:legendre-transformation` der Konvex-Konjugierten lässt
sich das ursprünglichen Problem äquivalent umformulieren zu
```{math}
:label: eq:sattelpunkt_aequivalent
\inf_{z\in\R^n} \sup_{x\in\R^n} ( J^*(z) + \langle -B^Tz, x \rangle - G(x) ),
```
wobei wir bemerken, dass wir in dieser Formulierung die explizite
Abhängigkeit von der Hilfsvariable $y = Bx$ vermeiden konnten.

Die neue Formulierung des Sattelproblems in
{eq}`eq:sattelpunkt_aequivalent` dient als Grundlage vieler
primal-dualer Verfahren in der Numerik. Beispielsweise können wir zur
Approximation eines Sattelpunkts abwechselnd einen Abstiegsschritt
mittels Proximaloperator bezüglich der Variable $z\in\R^n$ und einen
Gradientenaufstiegsschritt bezüglich der Variable $x\in\R^n$
durchführen. Dies führt zu einer Variante des sogenannten
**Uzawa-Verfahrens** mit Schrittweiten $\tau_k, \sigma_k \in \R^+$:
```{math}
:label: eq:uzawa
\begin{split}
z_{k+1} \ &= \ \operatorname{prox}_{\tau_k J^*}(z_k + \tau_kB x_k),\\
x_{k+1} \ &= \ x_k - \sigma_k(\nabla G(x_k) + B^T z_k).
\end{split}
```
In diesem Fall müssen wir die Matrix $B$ und $B^T$ nur einmal mit einem
Vektor multiplizieren in jedem Schritt, sowie den Proximaloperator der
Konvex-Konjugierten $J^*$ ausrechnen. Letzteres ist genau dann einfach,
wenn wir den Proximaloperator von $J$ effizient ausrechnen können, wie
folgende Bemerkung festhält.

````{prf:remark} Moreau-Zerlegung für Proximaloperatoren
Sei $J \colon \R^n \rightarrow \R$ eine konvexe, unterhalbstetige
Funktion und $J^* \colon \R^n \rightarrow \R$ die zugehörige
Konvex-Konjugierte von $J$. Sei außerdem $\tau > 0$ ein positiver
Parameter.

Dann kann man die folgende Identität, genannt Moreau-Zerlegung, für die
Proximaloperatoren von $J$ und $J^*$ zeigen:
```{math}
x \ = \ \operatorname{prox}_{\tau J^*}(x) + \tau \operatorname{prox}_{J/\tau}(x/\tau).
```

````

Somit können wir den Proximaloperator $\operatorname{prox}_{\tau_k J^*}$
der Konvex-Konjugierten in {eq}`eq:uzawa` berechnen als:
```{math}
\operatorname{prox}_{\tau_k J^*}(z_k  + \tau_kBx_k) \ = \ z_k + \tau_kBx_k - \tau_k \operatorname{prox}_{\frac{1}{\tau_k}J}\left(\frac{z_k  + \tau_kBx_k}{\tau_k} \right)
```
