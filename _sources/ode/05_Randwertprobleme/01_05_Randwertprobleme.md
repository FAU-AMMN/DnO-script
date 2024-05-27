(s:randwertproblem_existenz_eindeutigkeit)=
# Existenz und Eindeutigkeit von Lösungen

Mit Hilfe der Divergenzform {eq}`eq:divergenzform` können wir zunächst
für den Spezialfall $c(x)\equiv 0$ die sogenannte *Greensche Funktion*
zur analytischen Lösung der Differentialgleichung konstruieren. Es gilt
nämlich in diesem Fall mit einer Integrationskonstante $c_1 \in \R$
```{math}
a(x) \cdot u'(x) \ = \ c_1 - \int_0^x f(z) \, \mathrm{d}z.
```
Da $a(x) = e^{P(x)} > 0$ für alle $x \in \R$ ist, gilt damit auch für
eine weitere Integrationskonstante $c_2 \in \R$
```{math}
u(x) \ = \ c_2 + c_1 \cdot \int_0^x \frac{1}{a(y)} \, \mathrm{d}y - \int_0^x \frac{1}{a(y)}   \int_0^y f(z) \, \mathrm{d}z \, \mathrm{d}y.
```
Definieren wir nun zwei Stammfunktionen
```{math}
A(x) \ \coloneqq \ \int_0^x \frac{1}{a(y)} \, \mathrm{d}y, \qquad F(y) \ \coloneqq \ \int_0^y f(z) \, \mathrm{d}z
```
so erhalten wir durch partielle Integration
```{math}
:label: eq:randwertproblem_loesung_c=0
\begin{split}
u(x) \ &= \ c_2 + c_1 \cdot A(x) - \int_0^x \frac{1}{a(y)} \int_0^y f(z) \, \mathrm{d}z \, \mathrm{d}y \\
&= \  c_2 + c_1 \cdot A(x) - \int_0^x  \frac{1}{a(y)} \cdot F(y) \, \mathrm{d}y \\
&= \ c_2 + c_1 \cdot A(x) - \bigl[A(y) \cdot F(y)\bigr]^x_0 + \int^x_0 A(y) \cdot f(y) \, \mathrm{d}y\\
&= \ c_2 + c_1 \cdot A(x) - \int^x_0 A(x) \cdot f(y) \, \mathrm{d}y + \int^x_0 A(y) \cdot f(y)  \, \mathrm{d}y\\
&= \ c_2 + c_1 \cdot A(x) - \int_0^x (A(x) - A(y)) \cdot f(y)\, \mathrm{d}y.
\end{split}
```
Die unbekannten Integrationskonstanten $c_1, c_2 \in \R$ können wir aus
den Randbedingungen bestimmen, die uns ein Gleichungssystem liefern,
welches es zu lösen gilt. Hierzu betrachten wir im Folgenden die reinen
Dirichlet- oder Neumann-Randbedingungen und vernachlässigen die
gemischten Randbedingungen.

(dirichlet-randbedingungen)=
## Dirichlet-Randbedingungen

Im Fall von Dirichlet-Randbedingungen haben wir $\alpha_i=0$ und
$\beta_i=1$ für $i=0,1$ in {eq}`eq:randwert_rb1` und
{eq}`eq:randwert_rb2` und somit haben wir die Randwerte
```{math}
u(0) \ = \ g_0, \qquad u(1) \: = \: g_1.
```
Dann ist das Gleichungssystem zur Bestimmung der unbestimmten
Integrationskonstanten $c_1, c_2 \in \R$ durch Einsetzen in
{eq}`eq:randwertproblem_loesung_c=0` gegeben durch
```{math}
u(0) \: = \: c_2 \: = \: g_0, \qquad u(1) \ = \, \underbrace{c_2}_{= \: g_0} + \: c_1 \cdot A(1) - \int_0^1 (A(1) - A(y)) \cdot f(y) \, \mathrm{d}y \ = \ g_1.
```
Lösen wir die zweite Randbedingung auf zur Konstanten $c_1 \in \R$ so
erhalten wir:
```{math}
c_1 \ = \ \frac{1}{A(1)} \cdot \left( g_1 - g_0 + \int_0^1 (A(1) - A(y)) \cdot f(y) \, \mathrm{d}y \right).
```
Hieraus können wir für die Lösung der Differentialgleichung in
{eq}`eq:randwertproblem_loesung_c=0` folgern:
```{math}
:label: eq:randwertproblem_loesung_integralform
\begin{split}
u(x) \ &= \ g_0 + \frac{A(x)}{A(1)} \cdot \left( g_1 - g_0 + \int_0^1 (A(1) - A(y)) \cdot f(y) \, \mathrm{d}y\right) - \int_0^x (A(x) - A(y)) \cdot f(y) \, \mathrm{d}y\\
&= \ \left( 1- \frac{A(x)}{A(1)} \right) \cdot g_0 +  \frac{A(x)}{A(1)} \cdot g_1 \\
& \hspace{1.13cm} + \frac{A(x)}{A(1)}  \cdot \int_0^1 (A(1) - A(y)) \cdot f(y) \, \mathrm{d}y - \int_0^x (A(x) - A(y)) \cdot f(y) \, \mathrm{d}y.
\end{split}
```

Wir wollen nun die sogenannte **Greensche Funktion**
$G \colon [0,1]^2 \rightarrow \R$ einführen mit
```{math}
G(x,y) \ \coloneqq \ \left\{ 
\begin{array}{ll}
A(x) \cdot ( 1- \frac{A(y)}{A(1)}) \quad & \text{falls} \ 0 \leq x \leq y \leq 1, \\
A(y) \cdot ( 1-  \frac{A(x)}{A(1)}) \quad & \text{falls} \ 0 \leq y < x \leq 1.
\end{array}\right.
```
Für jede Funktion $a(x) > 0$ gilt für die Stammfunktion
$A(x) > A(y) \geq 0$ für $x > y$ und somit ist die Greensche Funktion
nichtnegativ auf dem Intervall $[0,1]$.

Nun können wir die beiden Integrale in
{eq}`eq:randwertproblem_loesung_integralform` geschickt vereinheitlichen
durch
```{math}
\begin{split}
&\hphantom{=} \ \frac{A(x)}{A(1)}  \cdot \int_0^1 (A(1) - A(y)) \cdot f(y) \, \mathrm{d}y - \int_0^x (A(x) - A(y)) \cdot f(y) \, \mathrm{d}y\\
&= \ \int_0^1 A(x) \cdot \left( 1- \frac{A(y)}{A(1)} \right) \cdot f(y) \, \mathrm{d}y\\
&\hspace{1cm} - \int_0^x \left(A(x) - \frac{A(x)A(y)}{A(1)} - A(y) + \frac{A(y)A(x)}{A(1)} \right) \cdot f(y) \, \mathrm{d}y\\
&= \ \int_0^1 A(x) \cdot \left( 1- \frac{A(y)}{A(1)} \right) \cdot f(y) \, \mathrm{d}y - \int_0^x A(x) \cdot \left(1 - \frac{A(y)}{A(1)} \right) \cdot f(y) \, \mathrm{d}y \\
&\hspace{1cm} + \int_0^x A(y) \cdot \left( 1 - \frac{A(x)}{A(1)} \right) \cdot f(y) \, \mathrm{d}y\\
&= \ \int_0^x A(y) \cdot (1 - \frac{A(x)}{A(1)}) \cdot f(y) \, \mathrm{d}y + \int_x^1 A(x) \cdot (1-\frac{A(y)}{A(1)}) \cdot f(y) \, \mathrm{d}y\\
&= \ \int_0^x G(x,y) \cdot f(y) \, \mathrm{d}y + \int_x^1 G(x,y) \cdot f(y) \, \mathrm{d}y  \ = \ \int_0^1 G(x,y) \cdot f(y) \, \mathrm{d}y.
\end{split}
```

Wir sehen also, dass wir mit Hilfe der Greenschen Funktion $G$ die
analytische Lösung $u$ in {eq}`eq:randwertproblem_loesung_integralform`
des Randwertproblems mit Dirichlet-Bedingungen und $c(x) \equiv 0$
kompakt schreiben können als:
```{math}
:label: eq:randwertproblem_dirichlet_loesung
u(x) \ = \ \left( 1- \frac{A(x)}{A(1)} \right) \cdot g_0 +  \frac{A(x)}{A(1)} \cdot g_1  + \int_0^1 G(x,y) \cdot f(y) \, \mathrm{d}y.
```
Aus dieser Darstellung der analytischen Lösung können wir sofort
folgendes Resultat ableiten.

````{prf:theorem} Existenz und Eindeutigkeitssatz für $c \equiv 0$
Sei $c(x) \equiv 0$ und $a \in C^1([0,1])$ mit $a(x) \geq a_0 > 0$ für
alle $x \in [0,1]$. Sei außerdem $f \in C([0,1])$.

Dann existiert eine eindeutige Lösung $u \in C^2([0,1])$ der
Differentialgleichung {eq}`eq:divergenzform` unter den
Dirichlet-Randwertbedingungen {eq}`eq:randwert_rb1` und
{eq}`eq:randwert_rb2` mit $\alpha_i=0$ und $\beta_i = 1$ für $i=0,1$.
Ist $f(x) \geq 0$ für alle $x \in [0,1]$, so gilt darüber hinaus
folgendes Maximumsprinzip:
```{math}
u(x) \ \geq \ \min\{g_0,g_1\} \qquad \forall x \in [0,1].
```

````

(allgemeiner-fall)=
### Allgemeiner Fall

Wir sehen, dass der *lineare Operator*
```{math}
\begin{split}
K \colon C([0,1]) \ &\rightarrow \ C([0,1]),\\
f  \ &\mapsto \ \int_0^1 G(\cdot,y) \cdot f(y) \, \mathrm{d}y
\end{split}
```
offensichtlich auf dem Raum der stetigen Funktionen $C([0,1])$
wohldefiniert und beschränkt (äquivalent zu stetig) ist. Mit Hilfe
dieser Einsicht können wir auch den allgemeinen Fall
$c(x) \not \equiv 0$ behandeln, da in diesem Fall die Lösung $u$ ein
Fixpunkt der folgenden Gleichung ist:
```{math}
:label: eq:randwertproblem_loesung_fixpunktgleichung
u(x) \ =\  g_0 \cdot w(x) + g_1 \cdot (1-w(x)) + K (f-c\cdot u)(x),
```
wobei wir die Notation $w(x) \coloneqq  1- \frac{A(x)}{A(1)}$ verwendet
haben.

Um die Charakterisierung der Lösung $u(x)$ in Abhängigkeit der Funktion
$c(x)$ besser zu verstehen, betrachten wir erst den Spezialfall einer
**konstanten Funktion** $c(x) \equiv c \in \R$. Da der Operator $K$
linear ist müssen wir die folgende Gleichung lösen:
```{math}
(I+c \cdot K)(u)(x) \ = \ g_0 \cdot w(x) + g_1 \cdot (1-w(x)) + K(f)(x).
```
Bei einer genaueren Analyse lässt sich nun feststellen, dass für die
spezielle Wahl der Konstanten $c = -k^2 \pi^2, k \in \N$ eine
nichttriviale Lösung $u(x) = \sin(k \pi \cdot x)$ des homogenen Systems,
d.h. für die Gleichung $(I - (k \pi)^2 \cdot K)(u)(x) = 0$ existiert. In
diesem Fall ist der Operator $I+c\cdot K$ nicht invertierbar, da dessen
Kern nicht nur die Nullfunktion enthält. Also können wir in diesen
Fällen die inhomogene Gleichung nicht für beliebige rechten Seiten
lösen.

Mit Hilfe der *Spektraltheorie kompakter Operatoren* lässt sich jedoch
zeigen, dass es nur abzählbar viele Werte von $c$ gibt, für die
$I+c \cdot K$ nicht invertierbar ist. Mit dem *Satz von Arzela-Ascoli*
kann man zunächst zeigen, dass $K: C([0,1]) \rightarrow C([0,1])$
kompakt ist, d.h. für jede beschränkte Folge
$(u_n)_{n\in\N} \subset C([0,1])$ hat $K (u_n)$ eine konvergente
Teilfolge. Die Spektraltheorie kompakter Operatoren garantiert nun, dass
es nur eine abzählbare Menge von Eigenwerten $(\lambda_k)_{k\in\N}$
geben kann, so dass der Operator $\lambda_k \cdot I - K$ nicht
invertierbar ist. Es lässt sich zeigen, dass in der Menge der Eigenwerte
Null der einzige Häufungspunkt ist. Setzen wir in diesem Zusammenhang
die Konstante $c \coloneqq - \frac{1}{\lambda_k}$ für $k \in \N$, so
sehen wir, dass es auch nur eine abzählbare Menge von Konstanten $c$
geben kann, für die der Operator $I+c \cdot K$ nicht invertierbar ist.
Diese Menge besitzt entsprechend die (uneigentlichen) Häufungspunkte
$\pm \infty$.

Wenn wir nun eine **allgemeine Funktion** $c(x)$ betrachten, dann können
wir immer noch zeigen, dass der Operator $u \mapsto u + K (c\cdot u)$
kompakt ist. Daraus können wir wieder folgern, dass der Operator genau
dann invertierbar ist, wenn sein Nullraum trivial ist. Bevor wir uns
hierfür eine hinreichende Bedingung erschließen, wollen wir zunächst im
folgenden Lemma eine nützliche Eigenschaft des Integraloperators $K$
nachweisen.

````{prf:lemma} Eigenschaften des Integraloperators
:label: lem:integraloperator_eigenschaften
Für alle stetigen Funktionen $f \in C([0,1])$ gilt die folgende
Abschätzung
```{math}
\int_0^1 f(x) \cdot K(f)(x) \, \mathrm{d}x \ = \ \int_0^1 \int_0^1 f(x) \cdot G(x,y) \cdot f(y) \, \mathrm{d}x \, \mathrm{d}y \ \geq \ 0.
```
Das Integral wird genau dann Null, wenn $f(x) \equiv 0$ gilt.

````

````{prf:proof} 
Um die Eigenschaften des Integraloperators zu zeigen, wählen wir ein
möglichst einfaches Randwertproblem der Form
```{math}
-(a(x) \cdot u'(x))' \ = \ f(x),
```
mit den Dirichlet-Randbedingungen $u(0) = u(1) = 0$. Wir können also die
analytische Lösung $u$ des Randwertproblems entsprechend
{eq}`eq:randwertproblem_loesung_fixpunktgleichung` für $g_0 = g_1 = 0$
und $c \equiv 0$ für $x \in [0,1]$ schreiben als
```{math}
u(x) \ = \ K(f)(x).
```

Somit können wir mit Hilfe von partieller Integration folgern, dass
gilt:
```{math}
\begin{split}
\int_0^1 f(x) \cdot K(f)(x) \, \mathrm{d}x \ &= \ - \int_0^1 (a(x)\cdot u'(x))' \cdot u(x) \, \mathrm{d}x \\
&= \ - \underbrace{\bigl[ a(x) \cdot u'(x) \cdot u(x) \bigr]^1_0}_{= \: 0} \: + \ \int_0^1 a(x) \cdot u'(x) \cdot u'(x) \, \mathrm{d}x\\
&= \ \int_0^1 a(x) \cdot (u'(x))^2 \, \mathrm{d}x \ \geq \ 0.
\end{split}
```
Die Nichtnegativität des letzten Integrals folgt aus dem quadratischen
Term $(u'(x))^2$ und der Eigenschaft, dass $a(x) = e^{P(x)} > 0$ für
alle $x \in [0,1]$ gilt. Damit folgt sofort, dass das Integral genau
dann Null wird, wenn $u'(x) \equiv 0$ gilt. Aus
$-(a(x) \cdot u'(x))' = f(x)$ folgt dann aber auch schon, dass
$f(x) \equiv 0$ gelten muss. ◻

````

Aus dem obigen Spezialfall mit konstanter Funktion $c(x) \equiv c$ sehen
wir, dass ein nichttrivialer Nullraum bei negativem $c$ auftreten kann.
Diese Beobachtung lässt sich nun mit Hilfe der Eigenschaft des
Integraloperators aus
[\[lem:integraloperator_eigenschaften\]](#lem:integraloperator_eigenschaften){reference-type="ref+label"
reference="lem:integraloperator_eigenschaften"} auf den allgemeinen Fall
von $c(x) \not \equiv 0$ übertragen, so dass wir im folgenden ein
allgemeines Resultat zur Existenz und Eindeutigkeit von Lösungen des
Randwertproblems erhalten.

````{prf:theorem} Existenz und Eindeutigkeitssatz für $c \geq 0$
:label: thm:randwertproblem_existenz-eindeutigkeit
Sei $a \in C^1([0,1])$ eine stetig differenzierbare Funktion mit
$a(x) \geq a_0 > 0$ für alle $x \in [0,1]$. Seien außerdem
$c, f  \in C([0,1])$ stetige Funktionen und es gelte $c(x) \geq 0$ für
alle $x\in [0,1]$.

Dann existiert eine eindeutige Lösung $u \in C^2([0,1])$ des
Randwertproblems {eq}`eq:divergenzform` mit Dirichlet-Randbedingungen
{eq}`eq:randwert_rb1` und {eq}`eq:randwert_rb1` mit $\alpha_i = 0$ und
$\beta_i = 1$ für $i=0,1$.

````

````{prf:proof} 
Basierend auf den oben genannten Argumenten aus der Funktionalanalysis
genügt es zu zeigen, dass der lineare Operator
$u\mapsto u+ K(c \cdot u)$ invertierbar ist, was äquivalent dazu ist,
dass er injektiv ist und somit sein Kern nur die Nullfunktion
$u(x) \equiv 0$ enthält.

Sei also $u(x) + K(c \cdot u)(x) = 0$, dann gilt ebenfalls nach
Multiplikation beider Seiten mit $c(x) \cdot u(x)$:
```{math}
c(x) \cdot u^2(x) \ = \ -c(x) \cdot u(x) \cdot K(c \cdot u)(x).
```
Integrieren wir beide Seiten der Gleichung und wenden nun
[\[lem:integraloperator_eigenschaften\]](#lem:integraloperator_eigenschaften){reference-type="ref+label"
reference="lem:integraloperator_eigenschaften"} für die speziell
gewählte stetige Funktion $f(x) \coloneqq c(x) \cdot u(x) \in C([0,1])$
an, so können wir abschätzen:
```{math}
\int_0^1 c(x) \cdot u^2(x) \, \mathrm{d}x \ = \ - \int_0^1 c(x) \cdot u(x) \cdot K(c\cdot u)(x) \, \mathrm{d}x \ \leq \ 0.
```
Da auf der linken Seite ein nichtnegativer Integrand steht muss hier
schon die Gleichheit mit Null vorliegen. Ebenfalls mit
[\[lem:integraloperator_eigenschaften\]](#lem:integraloperator_eigenschaften){reference-type="ref+label"
reference="lem:integraloperator_eigenschaften"} wissen wir, dass das
Integral genau dann Null wird, wenn $f(x) = c(x)\cdot u(x) \equiv 0$
ist. Wenn jedoch $c(x) \cdot u(x) \equiv 0$, so folgt aus der Linearität
von $K$, dass $K(c \cdot u)(x) \equiv 0$ gelten muss. Da wir
$u(x) + K(c \cdot u)(x) = 0$ angenommen haben, muss damit auch
$u(x) \equiv 0$ gelten. Also ist der Nullraum des Operators trivial,
woraus die Invertierbarkeit und damit auch die eindeutige Lösbarkeit des
Randwertproblems folgt. ◻

````

Interessanterweise können wir auch im allgemeinen Fall ein ähnliches
**Maximumsprinzip** zeigen, auch wenn wir nun keine explizite
Darstellung der Lösung als Integral mehr vorliegen haben.

````{prf:corollary} Maximumsprinzip
:label: cor:maximumsprinzip
Sei $a \in C^1([0,1])$ eine stetig differenzierbare Funktion mit
$a(x) \geq a_0 > 0$ für alle $x \in [0,1]$. Seien außerdem
$c, f  \in C([0,1])$ nichtnegative, stetige Funktionen, d.h. es gelte
$c(x) \geq 0$ und $f(x) \geq 0$ für alle $x\in [0,1]$.

Dann gilt für die eindeutige Lösung $u \in C^2([0,1])$ des
Randwertproblems {eq}`eq:divergenzform` mit beliebigen
Dirichlet-Randbedingungen $u(0) = g_0 \in \R$ und $u(1) = g_1 \in \R$,
die folgende Abschätzung:
```{math}
u(x) \ \geq \ \min \lbrace g_0, g_1 \rbrace \qquad \forall x \in [0,1].
```

````

````{prf:proof} 
Die Existenz und Eindeutigkeit einer Lösung $u \in C^2([0,1])$ des
Randwertproblems ist durch
[\[thm:randwertproblem_existenz-eindeutigkeit\]](#thm:randwertproblem_existenz-eindeutigkeit){reference-type="ref+label"
reference="thm:randwertproblem_existenz-eindeutigkeit"} gegeben. Für das
zu zeigende Maximumsprinzip führen wir einen einfachen
Widerspruchsbeweis.

Seien im Folgenden $c, f \in C([0,1])$ nichtnegative Funktionen und wir
nehmen an, dass die analytische Lösung $u(x)$ nicht am Rand des
Intervalls ihr Minimum annimmt, sondern in einem Punkt
$\overline{x} \in (0,1)$. Dann gilt die notwendige
Optimalitätsbedingungen $u'(\overline{x}) = 0$ in $\overline{x}$. Setzen
wir diese in die Differentialgleichung ein, so folgt für alle
$x \in [0,1]$
```{math}
0 \ \leq \ f(\overline{x}) \ = \ -(a(\overline{x}) \cdot \underbrace{u'(\overline{x})}_{=\: 0})' + c(\overline{x}) \cdot u(\overline{x}) \ = \ c(\overline{x}) \cdot u(\overline{x}) \ \leq \ c(\overline{x}) \cdot u(x).
```
Wegen der Annahme, dass $c$ und $f$ nichtnegative Funktionen sind, folgt
damit, dass $u(x)$ ebenfalls nichtnegativ für alle $x \in [0,1]$ ist. Da
wir jedoch beliebige Randwerte $u(0) = g_0$ und $u(1) = g_1$ für das
Randwertproblem erlauben, können wir ebenfalls Randwerte betrachten, so
dass $\min \lbrace g_0, g_1 \rbrace < 0$ gilt. Damit erzeugt die
gefolgerte Nichtnegativität der Lösung $u$ einen Widerspruch zu den
Randwerten und somit nimmt die analytische Lösung $u \in C^2([0,1])$ ihr
Minimum am Rand an und somit haben wir das Maximumsprinzip gezeigt. ◻

````

Wegen der Linearität des Problems kann man aus dem Maximumsprinzip in
[\[cor:maximumsprinzip\]](#cor:maximumsprinzip){reference-type="ref+label"
reference="cor:maximumsprinzip"} auch Stabilitätsaussagen herleiten.
Nehmen wir an $\tilde u \in C^2([0,1])$ sei eine bekannte Lösung für das
Problem
```{math}
-(a(x) \cdot \tilde u'(x))' + c(x) \cdot \tilde u(x) \ = \ \tilde f(x),
```
mit Dirichlet Randwerten $\tilde u(0) = \tilde{g}_0$ und
$\tilde u(1) = \tilde{g}_1$. Betrachten wir nun eine weitere rechte
Seite $f \in C([0,1])$ der Differentialgleichung mit
$\tilde f(x) \geq f(x)$ für alle $x\in [0,1]$, so folgt
```{math}
%TODO Auch hier stand eine 0 im Minimum!
\tilde u(x) - u(x) \ \geq \ \min\{\tilde{g}_0 - g_0, \tilde{g}_1 - g_1 \}.
```
und daraus eine Abschätzung für das Maximum von $u$. Ist $c > 0$, so
können wir etwa sehr einfache konstante Lösungen $\tilde u = \gamma$
konstruieren, indem wir $\tilde f = \gamma c$ setzen. Ist $\gamma$
hinreichend groß, dann ist $\tilde f > f$ und $\tilde{g}_i > g_i$ für
$i=0,1$.

Damit erhalten wir $u(x) \leq \gamma$ für alle $x$. Umgekehrt können wir
die Rollen von $\tilde u$ und $u$ auch vertauschen und
$\tilde u = -\gamma$ wählen, damit erhalten wir eine Abschätzung für
$\vert u (x) \vert.$

(neumann-randbedingungen)=
## Neumann-Randbedingungen

Während allgemeine Randbedingungen sehr ähnlich behandelt werden wie im
oben diskutierten Fall von Dirichlet-Randwertbedingungen, gibt es im
Fall reiner Neumann-Randbedingungen einen Aspekt, den wir besonders
beachten müssen.

Betrachten wir also im Folgenden die Differentialgleichung
{eq}`eq:divergenzform` mit Randwerten {eq}`eq:randwert_rb1` und
{eq}`eq:randwert_rb2` für $\alpha_i =1$ und $\beta_i =0$. Zur
Vereinfachung nehmen wir wieder an, dass $c(x) \equiv 0$ gilt. Dann
können wir die Differentialgleichung zunächst von links und rechts
aufintegrieren zu
```{math}
\begin{split}
&- \int_0^x (a(y) \cdot u'(y))' \, \mathrm{d}y \ = \ -a(x) \cdot u'(x) + a(0) \cdot \underbrace{u'(0)}_{=\: g_0} \ = \ \int_0^x f(y) \, \mathrm{d}y\\
\Rightarrow \quad &-a(x) \cdot u'(x) \ = \ - a(0) \cdot g_0 + \int_0^x f(y) \, \mathrm{d}y.
\end{split}
```
Integrieren wir hingegen mit dem rechten Rand auf so erhalten wir:
```{math}
\begin{split}
&- \int_x^1 (a(y) \cdot u'(y))' \, \mathrm{d}y \ = \ - a(1) \cdot \underbrace{u'(1)}_{=\: g_1} + \, a(x) \cdot u'(x)\ = \ \int_x^1 f(y) \, \mathrm{d}y\\
\Rightarrow \quad &-a(x) \cdot u'(x) \: = \: - a(1) \cdot g_1 - \int_x^1 f(y) \, \mathrm{d}y \: = \: - a(1) \cdot g_1 - \int_0^1 f(y) \, \mathrm{d}y + \int_0^x f(y) \, \mathrm{d}y.
\end{split}
```
Ein direkter Vergleich der beiden Gleichungen zeigt uns, dass die
folgende Bedingung erfüllt sein muss
```{math}
:label: eq:randwertproblem_neumann_konsistenzbedingung
a(0) \cdot g_0 \ = \  a(1) \cdot g_1 + \int_0^1 f(y) \, \mathrm{d}y
```
Unter dieser Bedingung ist die Gleichung prinzipiell lösbar, jedoch ist
eine der Konstanten $a(0), a(1) \in \R$ dann weiterhin unbestimmt. Dies
wird üblicherweise durch eine zusätzliche Normalisierungsbedingung der
folgenden Art erreicht
```{math}
\int_0^1 u(x) \, \mathrm{d}x \ = \ 1.
```
Der Grund für diese zusätzliche Bedingung ist im Gegensatz zu den oben
diskutierten Dirichlet-Randwertbedingungen ist, dass der Nullraum des
Randwertproblems {eq}`eq:divergenzform` mit $c(x) \equiv 0$ unter
Neumann-Randwertbedingungen nichttrival ist, denn jede konstante
Funktion löst bereits das homogene Problem
```{math}
-(a(x) \cdot u'(x))' \ = \ 0
```
mit Randwerten $u'(0) = u'(1) = 0$.

Interessanterweise ist die Theorie für den Fall nichtnegativer
Funktionen $c(x) \geq 0$ hier unterschiedlich. Sobald $c(x) \geq 0$ mit
$c(x) \not \equiv 0$ für alle $x \in [0,1]$ gilt, hat das
Randwertproblem mit Neumann-Randwertbedingungen nur noch einen trivialen
Nullraum. Dies lässt sich durch Multiplikation der homogenen
Differentialgleichung mit $u(x)$ und anschließender Integration sehen.
Es folgt dann nämlich mit partieller Integration und den
Randwertbedingungen $u'(0) = u'(1) = 0$, dass gilt
```{math}
\begin{split}
0 \ &= \ \int_0^1 (-(a(x) \cdot u'(x))' + c(x) \cdot u(x)) \cdot u(x) \, \mathrm{d}x \\
&= \ \underbrace{-\bigl[ a(x) \cdot u'(x) \cdot u(x) \bigr]^1_0}_{= \: 0} + \int_0^1 a(x) \cdot u'(x) \cdot u'(x) \, \mathrm{d}x + \int_0^1 c(x) \cdot u^2(x) \, \mathrm{d}x\\
&= \ \int_0^1 a(x) \cdot (u'(x))^2 + c(x) \cdot u^2(x) \, \mathrm{d}x.
\end{split}
```
Da $a(x) = e^{P(x)} > 0$ und $c(x) \geq 0$ angenommen wurde, folgt
sofort, dass $u'(x) \equiv 0$ für $x\in [0,1]$ gelten muss. Also muss
$u$ eine konstante Funktion sein. Da $c(x) \not \equiv 0$ angenommen
wurde kann
```{math}
\int_0^1 c(x) \cdot u^2(x) \, \mathrm{d}x \ = \ 0
```
für eine konstante Funktion $u$ nur im Fall $u(x) \equiv 0$ gelten.

Der Unterschied dieser Fälle und der potentiell nichttriviale Nullraum
des Randwertproblems mit Neumann-Randwertbedingungen ist natürlich auch
bei der numerischen Lösung zu beachten, welche im nächsten Abschnitt
diskutiert wird. Hier wird sich dieser Nullraum in der
Nichtinvertierbarkeit einer zugehörigen Matrix niederschlagen.

