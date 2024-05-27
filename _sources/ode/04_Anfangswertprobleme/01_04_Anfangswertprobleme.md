(theorie-für-anfangswertprobleme-gewöhnlicher-differentialgleichungen)=
# Theorie für Anfangswertprobleme gewöhnlicher Differentialgleichungen

Im Folgenden beschäftigen wir uns zunächst mit Anfangswertproblemen für
gewöhnliche Differentialgleichungen, im Englischen *Ordinary
differential equations (ODE)* genannt.

````{prf:definition} Anfangswertproblem
Ein mathematischen Problem, bei dem wir eine Lösung
$u \colon \R^+ \rightarrow \R^n$ eines gewöhnlichen
Differentialgleichungssystems $1$. Ordnung mit gegebener
Startwertbedingung der folgenden allgemeinen Form suchen
```{math}
:label: eq:awp
\begin{split}
\frac{du}{dt} \ = \ u'(t) \ &= \ F(t,u(t)), \\ 
u(0) \: &= \: u_0 \in \R^n,
\end{split}
```
wobei $F: \mathbb{R}^+ \times \mathbb{R}^n \rightarrow \mathbb{R}^n$
eine stetige Funktion ist, nennen wir ein **Anfangswertproblem**.

````

Anfangswertprobleme der Form {eq}`eq:awp` treten häufig in der
Modellierung von Zuständen auf, die sich über die Zeit ändern. Nimmt man
an, dass der zu modellierende Zustand zum Zeitpunkt $t > 0$ durch die
Funktion $u \colon \R^+ \rightarrow \R^n$ beschrieben wird, möchte man
verstehen wie sich dieser in hinreichend kleiner Zeit durch die Wirkung
der Funktion $F$ ändern wird. Das heißt wir schreiben
```{math}
u(t+\tau) \ \approx \ u(t) + \tau \cdot F(t,u(t)).
```
Im Grenzwert für immer kleinere Zeitschrittweiten $\tau \rightarrow 0$
erhalten wir dann die Differentialgleichung $u'(t) = F(t, u(t))$.

Solche Differentialgleichungen lassen sich beispielsweise in der Physik
herleiten, wie folgenden Beispiel zeigt.

````{prf:example} Newtonsche Bewegungsgleichung
:label: ex:newton_bewegungsgleichung
Wir betrachten in diesem Beispiel die Bewegung eines Teilchens mit Hilfe
des zweiten Newtonschen Gesetzes. Beschreibt $x(t) \in \R^3$ den Ort des
Teilchens zur Zeit $t > 0$ im dreidimensionalen Raum, $v(t) \in \R^3$
seine Geschwindigkeit, $a(t)$ seine Beschleunigung und $F$ ein Kraftfeld
in dem sich das Teilchen bewegt, so sind folgende Zusammenhänge aus dem
physikalischen Mechanik bekannt:

-   Geschwindigkeit ist Änderung des Orts pro Zeit, d.h.,
    $v(t) = \dot{x}(t) = \frac{dx}{dt}(t)$

-   Beschleunigung ist Änderung der Geschwindigkeit pro Zeit, d.h. es
    gilt $a(t) = \dot{v}(t) = \frac{dv}{dt}(t)$

-   Kraft ist Masse mal Beschleunigung, d.h.,
    $m\cdot a(t) = F(x(t),v(t),t)$

Diese Gesetze können wir auch als Differentialgleichungen für den Ort
$x$ und die Geschwindigkeit $v$ des Teilchens (oder dessen Impuls
$p(t)=m\cdot v(t)$) interpretieren, wenn wir die Beschleunigung $a$
eliminieren. Es gilt dann
```{math}
\frac{d}{dt} (x(t),m\cdot v(t)) = (v(t), F(x(t),v(t),t)).
```
Kennen wir den Anfangsort $x(0) \in \R^3$ des Teilchens und dessen
Anfangsgeschwindigkeit $v(0) \in \R^3$, so haben wir ein
Anfangswertproblem für ein System von Differentialgleichungen mit
insgesamt sechs Gleichungen für sechs Unbekannte.

In großen Systemen mit $n \in \N$ Teilchen hat man dann Gleichungen für
jede Position $x_1,\ldots,x_n \in \R^3$ und die respektiven
Geschwindigkeiten $v_1,\ldots,v_n \in \R^3$ und sucht Lösungen für das
folgende System von Differentialgleichungen:
```{math}
\frac{d}{dt} (x_i(t),v_i(t)) = (v_i(t), F(x_1(t),v_1(t), \ldots , x_n(t),v_n(t),t)), \quad i=1,\ldots,n.
```

Zur Beschreibung realer Vorgänge mit vielen Teilchen (Moleküle, Zellen,
Tiere in Herden, Fußgänger, Autos \...) erhält man also schnell beliebig
komplexe Systeme von Differentialgleichungen.

````

(ss:awp_existenz_eindeutigkeit)=
## Existenz und Eindeutigkeit von Lösungen

Um die in dieser Vorlesung diskutierten Verfahren zur numerischen Lösung
von gewöhnlichen Differentialgleichungen besser zu verstehen müssen wir
im Folgenden zunächst einige theoretische Grundlagen von gewöhnlichen
Differentialgleichungen zusammenfassen. Wir beginnen mit einer
allgemeinen Theorie der Existenz und Eindeutigkeit von Lösungen solcher
Differentialgleichungen. Die wesentliche Grundlage hierfür wird eine
Umformulierung der gewöhnlichen Differentialgleichung in eine
Fixpunktform sein, so dass wir dann einfach einen passenden Fixpunktsatz
anwenden können.

Nehmen wir zunächst an, dass $u \in C^1([0,T])$ eine Lösung des
Anfangswertproblems {eq}`eq:awp` sei. Dann gilt nach dem Hauptsatz der
Differential- und Integralrechnung auch folgende Gleichung
```{math}
u(t) \ = \ u_0 + \int_0^t F(s,u(s))\,\mathrm{d}s , \qquad 0 \leq t \leq T.
```
Diese Darstellung können wir als Fixpunktgleichung $u= {\cal F}(u)$ in
einem Banachraum interpretieren. Um Existenz und Eindeutigkeit eines
Fixpunktes $u$ und somit einer Lösung des Anfangswertproblems zu zeigen,
lassen sich nun zwei unterschiedliche Arten von Fixpunktsätzen anwenden.

Die erste Art (*Satz von Schauder* oder andere Varianten) basiert auf
Kompaktheit, d.h. wenn man den Operator $\mathcal{F}$ auf Funktionen $u$
anwendet, liegt das Bild in einer topologischen Menge mit einem Fixpunkt
des Operators. In unserem Fall erhält man die Kompaktheit aus dem *Satz
von Arzela--Ascoli*, denn man kann zeigen, dass für eine beschränkte
Folge von Funktionen $(u_n)_{n\in\N}$ die Folge
$({\cal F}(u_n))_{n\in\N}$ immer eine konvergente Teilfolge besitzt.
Dies liefert den sogenannten **Satz von Peano**, der die Existenz einer
Lösung für eine stetige Funktion $F$ garantiert. Da man hier jedoch
recht abstrakt und mittels Teilfolgen argumentiert, lässt sich mit
dieser Methode leider nicht die Eindeutigkeit eines Fixpunkts
nachweisen.

Die zweite Art um Existenz und Eindeutigkeit eines Fixpunkts zu zeigen
basiert eigentlich immer auf dem *Banachschen Fixpunktsatz* (siehe
{cite:p}`numerik1`). Diese Herangehensweise wollen wir im Folgenden
etwas näher diskutieren. Dazu beachten wir, dass falls die Funktion $F$
im zweiten Argument Lipschitz-stetig ist (mit Lipschitz-Konstante
$L > 0$), folgendes gilt
```{math}
\Vert {\cal F}(u_1) - {\cal F}(u_2) \Vert_\infty \ &= \ \max_{0 \leq t \leq T} \left\vert \int_0^t F(s,u_1(s)) - F(s,u_2(s)) \, \mathrm{d}s \right\vert \\
&\leq \  \max_{0 \leq t \leq T} \int_0^t \vert  F(s,u_1(s)) - F(s,u_2(s))\vert \, \mathrm{d}s\\
&\leq  \ L \cdot \max_{0 \leq t \leq T} \int_0^t \vert u_1(s) - u_2(s) \vert  \, \mathrm{d}s\\
&\leq  \ L T \cdot \max_{0 \leq s \leq T} \vert u_1(s) - u_2(s) \vert \ = \ L T \cdot \Vert u_1 - u_2 \Vert_\infty. 

```
Wir erkennen an dieser Abschätzung, dass die Abbildung
${\cal F}:C^1([0,T]) \rightarrow C^1([0,T])$ kontraktiv ist, wenn
$T > 0$ klein genug gewählt wird, so dass $L T < 1$ gilt. Mit Hilfe der
Kontraktivität von $\mathcal{F}$ liefert der Banach'sche Fixpunktsatz
die Existenz und Eindeutigkeit eines Fixpunkts $u$ (und damit einer
Lösung der Differentialgleichung) in $C^1([0,T])$ für $T$ hinreichend
klein. Dies ist die erste Version des **Satzes von Picard-Lindelöf,**
die uns die Existenz und Eindeutigkeit für kleine Zeiten $T > 0$
liefert.

In unserem Fall können wir jedoch ein noch besseres Resultat erreichen,
in dem wir einfach die Norm in unserer Abschätzung passend wählen. Die
Idee dazu liefert zunächst eine einfache Differentialgleichung im
folgenden Beispiel.

````{prf:example} Einfache Differentialgleichung mit konstantem Faktor
:label: ex:dgl_einfach
Wir beschäftigen uns in diesem Beispiel mit einer einfachen gewöhnlichen
Differentialgleichung mit einer Konstanten $L > 0$ der Form
```{math}
u'(t) \ = \ L \cdot u(t).
```
Nehmen wir $u \colon \R_+ \rightarrow \R$ als positive Lösung an, dann
folgt mit der Kettenregel der Differentialgleichung
```{math}
\frac{d}{dt} (\log u(t)) \ = \ \frac{u'(t)}{u(t)} \ = \ L.
```
Dies können wir mit Hilfe des Hauptsatzes der Differential- und
Integralrechnung schreiben als
```{math}
\log u(t) - \log u(0) \ = \ \log\left( \frac{u(t)}{u_0} \right) \ = \ L \cdot t.
```
Wenden wir auf beide Seiten die Exponentialfunktion an, erhalten wir die
folgende analytische Lösung der Differentialgleichung:
```{math}
u(t) \ = \ u_0 \cdot e^{L t}.
```

````

Wählen wir nun den Vorfaktor $L$ in {prf:ref}`ex:dgl_einfach` als die
Lipschitzkonstante im zweiten Argument der Funktion $F$, dann erkennen
wir, dass wir ein exponentielles Wachstum mit $e^{L\cdot t}$ erwarten
müssen.

Diese Beobachtung ist auch im Allgemeinen gültig, wie das folgende Lemma
zeigt.

````{prf:lemma} Lemma von Gronwall
:label: lem:gronwall
Sei $v(t)$ eine nichtnegative, stetige Funktion, die die folgende
Ungleichung
```{math}
:label: eq:gronwall_ungleichung_annahme
v(t) \ \leq \ a + \int_0^t b \cdot v(s)\,\mathrm{d}s, \quad \forall \ 0 \leq t \leq T,
```
mit $a > 0$ und $b \in \mathbb{R}$ erfüllt.

Dann lässt sich die Funktion $v$ nach oben abschätzen durch
```{math}
v(t) \ \leq \ a \cdot e^{b t}, \quad \forall \ 0 \leq t \leq T.
```

````

````{prf:proof} 
Wir definieren zunächst die Hilfsfunktion
$w(t) \coloneqq e^{-bt} \cdot v(t) - a$. Damit können wir $v$ auch
schreiben als $v(t) = e^{bt} \cdot(w(t) + a)$. Da die Ungleichung
{eq}`eq:gronwall_ungleichung_annahme` nach Voraussetzung gilt können wir
folgende Abschätzung machen:
```{math}
\begin{split}
w(t) \ &= \ e^{-bt} \cdot v(t) - a \\
&\leq \ e^{-bt} \cdot a + e^{-bt} \cdot \int_0^t b \cdot v(s)\,\mathrm{d}s - a \\
&= \ a \cdot (e^{-bt}-1)+ \int_0^t e^{-bt} \cdot b \cdot v(s)\,\mathrm{d}s\\
&= \ a \cdot (e^{-bt}-1)+ \int_0^t e^{b(s-t)} \cdot b \cdot (w(s)+a)\,\mathrm{d}s \\
&= \ a \cdot (e^{-bt}-1)+ ab\cdot e^{-bt} \cdot \int_0^t e^{bs}\,\mathrm{d}s + b\cdot \int_0^t e^{b(s-t)} \cdot w(s) \, \mathrm{d}s\\
&\leq \ a \cdot (e^{-bt}-1)+ ab \cdot e^{-bt} \cdot \left[ \frac{1}{b} e^{bs} \right]^t_0 +  b \cdot \int_0^t \underbrace{e^{b(s-s)}}_{= \: 1} \cdot \: w(s) \, \mathrm{d}s\\
&= \ \underbrace{a \cdot (e^{-bt}-1)+ a \cdot (1 - e^{-bt})}_{= \: 0} + \: b \int_0^t w(s) \, \mathrm{d}s\\
&= \ \int_0^t b \cdot w(s)\,\mathrm{d}s.
\end{split}
```
Für die letzte Ungleichung haben wir verwendet, dass
$e^{-bs} \geq e^{-bt}$ gilt für alle $0 \leq s \leq t$. Aus dieser
Ungleichung können wir für $t=0$ folgern, dass $w(0) \leq 0$ und somit
$v(0) \leq a \cdot e^0$ gilt. Somit gilt die Abschätzung des Lemmas
bereits in diesem Fall.

Da $v$ als stetige Funktion vorausgesetzt wurde ist $w$ ebenfalls stetig
und somit existiert ein maximaler Zeitpunkt $T \geq T_0 > 0$ für den
$w(t) \leq 0$ für alle $t \leq T_0$ gilt. Ist $T_0=T$, so haben wir die
Abschätzung des Lemmas bereits gezeigt, da aus $w(t) \leq 0$ für alle
$0 \leq t \leq T$ durch Umstellen schon folgt
$v(t) \leq a \cdot e^{bt}$.

Nehmen wir nun an, dass $0 < T_0 < T$ gilt, so gibt es ein hinreichend
kleines Zeitintervall $(T_0,T_0+\delta)$, in dem $w$ positiv ist. Wegen
der Positivität von $w$ folgt dann für ein beliebiges
$t \in (T_0, T_0 + \delta)$ die Ungleichung
```{math}
w(t) \ \leq \ \int_{T_0}^t b \cdot w(s) \,\mathrm{d}s \ \leq \ \delta \cdot \max_{T_0 \leq s \leq T_0 + \delta} b\cdot w(s)
```
und somit auch
```{math}
\max_{T_0 \leq t \leq T_0 + \delta} w(t) \ \leq \ \delta \cdot \max_{T_0 \leq s \leq T_0 + \delta} b \cdot w(s) \ = \ \delta b \cdot \max_{T_0 \leq s \leq T_0 + \delta} w(s).
```
Für $\delta b < 1$ ist dies aber ein Widerspruch zur Positivität von
$w$. Also muss $w(t) \leq 0$ und damit $u(t) \leq a e^{bt}$ für alle
$0 \leq t \leq T$ gelten. ◻

````

Abstrakt gesehen liefert das obige Lemma von Gronwall eine
Stabilitätsaussage für die Differentialgleichung, denn es besagt, dass
die Norm der Lösung in endlicher Zeit nicht beliebig wachsen kann.
Betrachten wir nämlich eine Differentialgleichung mit einer Funktion
$F$, die Lipschitz-stetig im zweiten Argument ist, so folgt mit der Wahl
der Funktion $v(t) = \vert u(t) - u_0 \vert$ nämlich
```{math}
\begin{split}
v(t) \ = \ \vert u(t) -u_0 \vert \ &= \ \left\vert \int_0^t F(s, u(s) \,\mathrm{d}s\right\vert \ \leq \ \int_0^t \vert F(s, u(s)\vert \,\mathrm{d}s\\ 
&= \int_0^t \vert F(s, u(s) - F(s, u_0) + F(s, u_0) \vert \,\mathrm{d}s \\
&\leq \ \int_0^t \vert F(s,u(s)) - F(s,u_0) \vert + \vert F(s,u_0) \vert \,\mathrm{d}s \\
&\leq \ \int_0^t \vert F(s,u(s)) - F(s,u_0) \vert \,\mathrm{d}s + \underbrace{T \cdot \max_{0 \leq t \leq T} \vert F(t,u_0)\vert}_{=: \ a} \\
&\leq \ L \cdot \int_0^t \vert u(s) - u_0\vert \,\mathrm{d}s + a \ = \ b \cdot \int_0^t v(s) \,\mathrm{d}s + a.
\end{split}
```
Wir sehen also, dass die Voraussetzungen für {prf:ref}`lem:gronwall`
erfüllt sind für die Funktion $v(t)$ mit $b \coloneqq L$ und
$a \coloneqq T \cdot \max\limits_{0 \leq t \leq T} \vert F(t,u_0)\vert$.
Somit gilt also nach dem Lemma von Gronwall die Abschätzung
```{math}
\vert u(t) -u_0 \vert \ \leq \ a \cdot e^{bt} \ = \ T \cdot \max_{0 \leq t \leq T} \vert F(t,u_0)\vert \cdot e^{Lt}.
```

Verwenden wir die umgekehrte Dreiecksungleichung
$|u(t)| - |u_0| \leq |u(t) -u_0|$ sehen wir, dass $u(t)$ beschränkt ist
und höchstens wie $e^{Lt}$ wächst:
```{math}
|u(t)| \ \leq \ |u(t) -u_0| + |u_0| \ \leq \ T \cdot \max_{0 \leq t \leq T} \vert F(t,u_0)\vert \cdot e^{Lt} + |u_0|.
```
Dies führt zu der Idee, die folgende gewichtete Norm zu betrachten:
```{math}
\Vert u \Vert_{\infty,L} \ \coloneqq \ \max_{0 \leq t \leq T} e^{-Lt} \cdot \vert u(t) \vert.

```
Da $e^{-Lt}$ nach oben durch Eins und nach unten durch $e^{-LT} \geq 0$
beschränkt ist, ist dies eine äquivalente Norm im Raum der stetig
differenzierbaren Funktionen $C^1([0,T])$.

Wir wiederholen also unsere Abschätzung des Fixpunktoperators in dieser
konstruierten Norm und erhalten somit
```{math}
\Vert {\cal F}(u_1) - {\cal F}(u_2) \Vert_{L,\infty} \ &= \ \max_{0 \leq t \leq T}  e^{-Lt} \cdot \left\vert \int_0^t F(s,u_1(s)) - F(s,u_2(s)) \,\mathrm{d}s \right\vert \\
&\leq \ \max_{0 \leq t \leq T} \int_0^t e^{L(s-t)} \cdot e^{-Ls} \cdot \vert  F(s,u_1(s)) - F(s,u_2(s))\vert\,\mathrm{d}s\\
&\leq  \ \int_0^T e^{-L\tau}\,\mathrm{d}\tau \cdot \max_{0 \leq s \leq T}  e^{-Ls} \cdot L \cdot \vert u_1(s) - u_2(s) \vert \\
&= \ (1-e^{-LT}) \cdot \Vert u_1 - u_2 \Vert_{L,\infty}. 

```
Der Operator ist nun also kontraktiv bezüglich der gewählten Norm für
beliebiges $T > 0$, da stets $1- e^{-LT} < 1$ gilt.

Wir haben mit Hilfe des Banachschen Fixpunktsatzes also die folgende
Variante des Satzes von Picard-Lindelöf bewiesen.

````{prf:theorem} Satz von Picard-Lindelöf
:label: thm:picard-lindelöf
Sei $F$ eine stetige Funktion, die Lipschitz-stetig bezüglich der
zweiten Variable ist.

Dann besitzt das Anfangswertproblem {eq}`eq:awp` eine eindeutige Lösung
in $C^1([0,T])$.

````

Wir werden sehen, dass wir auch bei numerischen Verfahren zur Lösung von
Anfangswertproblemen ähnliche Aussagen und insbesondere eine Variante
des Lemma von Gronwall in {prf:ref}`lem:gronwall` benötigen werden, um
die Stabilität dieser Verfahren garantieren zu können.

(analytische-lösungsverfahren)=
## Analytische Lösungsverfahren

Wir betrachten im Folgenden einige spezielle Fälle von gewöhnlichen
Differentialgleichungen in denen wir analytisch eine geschlossene Form
der Lösung berechnen können.

Andererseits gibt es viele gewöhnliche Differentialgleichungen für die
sich analytisch keine Lösung angeben lässt, wie zum Beispiel die
folgende simple, nichtlineare Gleichung:
```{math}
u'(t) \ = \ t^2 + u^2(t).
```
Ein weiteres bekanntes Phändomen, dass durch ein
Differentialgleichungssystem beschrieben ist, jedoch nicht analytisch
lösbar ist, ist das das *$N$-Körper-Problem*
{cite:p}`N_koerper_problem`, welches auch schon in
{prf:ref}`ex:newton_bewegungsgleichung` durch die geeignete Wahl des
Kraftfeldes beschrieben wird. Aus diesem Grund wollen wir in dieser
Vorlesung numerische Lösungsmethoden für gewöhnliche
Differentialgleichungen diskutieren und untersuchen.

Hierbei werden wir nur kurz die verschiedenen Herangehensweisen zur
Herleitung einer Lösung skizzieren. Für eine formale Beschreibung dieser
Ansätze verweisen wir auf beispielsweise auf {cite:p}`tenbrinck2021`.

(lineare-differentialgleichungen)=
### Lineare Differentialgleichungen

Eine Klasse von gewöhnlichen Differentialgleichungen, die in folgender
Definition erklärt sind, stellen sich als gut verständlich und
analytisch lösbar heraus.

````{prf:definition} Lineare Differentialgleichung
:label: def:lineare_dgl
Wir nennen ein gewöhnliches Differentialgleichungssystem $n$-ter Ordnung
**linear**, wenn es sich in folgender Form schreiben lässt:
```{math}
\sum_{i=0}^n a_i(t)u^{(i)}(t) \ = \ a_0(t) \cdot u(t) + a_1(t) \cdot u'(t) + \ldots + a_n(t) \cdot u^{(n)}(t) \ = \ b(t).
```
Hierbei sind die
$a_i \colon \R_+ \rightarrow \R^{n\times n}, i=0, \ldots,n$ und
$b \colon \R_+ \rightarrow \R^n$ stetige Funktionen, die nicht von $u$
abhängen, und $u^{(i)}$ bezeichnet die $i$-te Ableitung der unbekannten
Funktion $u$.

````

Ein wichtiger Spezialfall von {prf:ref}`def:lineare_dgl`, der uns auch
kanonische Beispiele für numerische Verfahren liefert, ist ein lineares
Differentialgleichungssystem $1$. Ordnung mit konstanten Koeffizienten
der Form
```{math}
:label: eq:lin_dgl
u'(t) \ = \ A \cdot u(t) + b(t),
```
mit einer gegebenen Matrix $A \in \R^{n \times n}$.

(matrixexponential)=
### Matrixexponential

Wir betrachten zunächst den homogenen Fall für $b\equiv\vec{0}$ des
linearen Differentialgleichungssystem in {eq}`eq:lin_dgl`. Ist die
Matrix $A$ diagonalisierbar durch $A = B^{-1}\cdot D\cdot B$ mit einer
Diagonalmatrix $D \in \R^{n \times n}$, so können wir analog eine
unbekannte Hilfsvariable $v(t) = B \cdot u(t)$ betrachten. Damit können
wir nun das lineare Differentialgleichungssystem {eq}`eq:lin_dgl`
schreiben als
```{math}
\begin{split}
u'(t) \ &= \ A \cdot u(t) \ = \ B^{-1}\cdot D\cdot B \cdot u(t) \ = \ B^{-1}\cdot D\cdot v(t)\\
\Rightarrow \ v'(t) \ &= \ (B \cdot u(t))' \ = \ B \cdot u'(t) \ = \ D \cdot v(t).
\end{split}
```
Somit können wir also äquivalent das folgende simple
Differentialgleichungssystem lösen:
```{math}
v'(t) \ = \ D \cdot v(t).
```
Da $D$ eine Diagonalmatrix ist, sind die Gleichungen entkoppelt und das
bedeutet, dass wir für jede Zeile eine gewöhnliche Differentialgleichung
mit konstantem Vorfaktor der Form $v_i'(t) = D_{ii} v_i(t)$ erhalten.
Für diese können wir analog zu {prf:ref}`ex:dgl_einfach` eine explizite
Lösung angeben als
```{math}
v_i(t) \ = \ v_i(0) \cdot e^{D_{ii}\cdot t}.
```
Die unbekannte Lösung $u$ der urspünglichen Differentialgleichung
{eq}`eq:lin_dgl` erhalten wir anschließend durch
```{math}
u(t) \ = \ B^{-1} \cdot v(t).
```

Zur Vereinfachung lässt sich dieser Lösungsansatz mit Hilfe des
sogenannten Matrixexponentials schreiben.

````{prf:definition} Matrixexponential
:label: def:matrixexponential
Sei $A \in \R^{n \times n}$ eine quadratische Matrix. Das
**Matrixexponential** $e^A \in \R^{n \times n}$ ist dann definiert durch
die folgende Potenzreihe
```{math}
e^A \ \coloneqq \ \sum_{k=0}^\infty \frac{A^k}{k!} \ = \  I_n + A + \frac{A^2}{2} +  \ldots
```

````

Für den Fall einer Diagonalmatrix $D \in \R^{n\times n}$ und einem
Skalar $t \in R$ ist das Matrixexponential $e^{Dt}$ in
{prf:ref}`def:matrixexponential` eine Diagonalmatrix mit Einträgen
$(e^{Dt})_{i,i} = e^{D_{i,i}\cdot t}$ für $i=1,\ldots,n$. Mit dieser
Erkenntnis lässt sich die Lösung der Differentialgleichung
{eq}`eq:lin_dgl` kompakt angeben als:
```{math}
:label: eq:lin_dgl_homogene_loesung
u(t) \ = \ B^{-1} \cdot v(t) \ = \ B^{-1} \cdot e^{Dt} \cdot v(0) \ = \  B^{-1} \cdot e^{Dt} \cdot B \cdot u(0) \ =: \ e^{At} \cdot u(0).
```

Wir bekommen durch Diagonalisieren der Matrix $A$ also eine simple Form
des Matrixexponentials, die uns direkt eine Lösung des
Anfangswertproblems liefert. Im Fall einer nicht diagonalisierbaren
Matrix ist ein ähnliches Vorgehen über die Jordan'sche Normalform
möglich, was an dieser Stelle aber über den Rahmen dieser Vorlesung
hinaus gehen würde.

(variation-der-konstanten)=
### Variation der Konstanten

Im inhomogenen Fall $b \neq \vec{0}$ eines linearen
Differentialgleichungssystems erster Ordnung in {eq}`eq:lin_dgl` lässt
sich die sogenannte Variation der Konstanten benutzen um eine Lösung zu
berechnen. Die Idee hierbei ist es einen Produktansatz der unbekannten
Lösung $u$ zu betrachten mit
```{math}
u(t) \ = \ e^{At} \cdot w(t).
```
Anstatt der Konstanten $u(0) = u_0 \in \R^n$ in der homogenen Lösung in
{eq}`eq:lin_dgl_homogene_loesung` betrachten wir also jetzt eine
zeitabhängige Funktion $w(t)$.

Mit Hilfe der Produktregel für Differentiation und der Eigenschaften des
Matrixexponentials lässt sich somit die Ableitung der Funktion $u$
schreiben als
```{math}
u'(t) \ = \ \frac{d}{dt}(e^{At} \cdot w(t)) \ = \ A \cdot \underbrace{e^{At} \cdot w(t)}_{=\, u(t)} + \, e^{At} \cdot w'(t) \ = \ A \cdot u(t) + e^{At}\cdot w'(t),
```
Vergleichen wir nun diesen Ausdruck der Ableitung der unbekannten
Funktion $u$ mit rechten Seite der Differentialgleichung
{eq}`eq:lin_dgl`, so sehen wir ein, dass $b(t) = e^{At} \cdot w'(t)$
gelten muss. Durch Umstellen sehen wir ein, dass für die Ableitung der
unbekannten Hilfsfunktion $w'(t) = e^{-At} \cdot b(t)$ gilt. Integrieren
wir diesen Ausdruck, so lässt sich die unbekannte Lösung $u$ durch die
folgende Identität berechnen:
```{math}
u(t) \ = \ e^{At} \cdot w(t) \ = \ e^{At} \cdot \left( c + \int_0^t w'(s) \,\mathrm{d}s \right) \ = \ e^{At} \cdot \left( c + \int_0^t e^{-As} \cdot b(s) \,\mathrm{d}s\right).
```
Hierbei lässt sich das Matrixexponential
$e^{\pm At} \in \R^{n \times n}$ wie im homogenen Fall beschrieben
bestimmen und die Integrationskonstante $c \in \R$ erhält man durch
Anwendung der Anfangswertbedingung $u(0) = u_0 \in \R^n$.

(trennung-der-variablen)=
### Trennung der Variablen

Ein weiterer Spezialfall für die analytische Bestimmung von Lösungen
gewöhnlicher Differentialgleichungen in Dimension $n=1$ sind sogenannte
separable Gleichungen der folgenden Form:
```{math}
u'(t) \ = \ g(u(t)) \cdot h(t),
```
wobei $g,h \colon \R_+ \rightarrow \R$ zwei stetige Funktionen sind.
Nehmen wir an, dass $g(u(t))$ nicht Null wird, so können wir durch
diesen Term teilen und erhalten
```{math}
\frac{u'(t)}{g(u(t))} \ = \ h(t).
```
Ist $G$ eine Stammfunktion von $\frac{1}g$ und $H$ eine Stammfunktion
von $h$, so schreiben wir diese Gleichung mit Hilfe der Kettenregel der
Differentiation als
```{math}
\frac{d}{dt} G(u(t)) \ = \ G'(u(t)) \cdot u'(t) \ = \ H'(t).
```
Diese Gleichung können wir nun aufintegrieren zu $G(u(t)) = H(t) + c$.
Die Integrationskonstante $c \in \R$ können wir anschließend aus dem
Anfangswert mit $c = G(u_0) - H(0)$ berechnen.

Setzen wir weiter voraus, dass die Stammfunktion $G$ von $g^{-1}$
invertierbar ist, lässt sich die unbekannte Lösung $u$ nun wie folgt
berechnen:
```{math}
u(t) \ = \ G^{-1} (H(t) - H(0) + G(u_0)).
```

(gradientenfluss)=
### Gradientenfluss

Zum Abschluss betrachten wir noch eine Klasse von gewöhnlichen
Differentialgleichungen, die wir aus dem Gradientenabstiegsverfahren der
numerischen Optimierung in {ref}`ss:gradient_descent` erhalten. Hierbei
interpretieren wir die Iterationsschritte $x_k \in \R^n$ als Wert einer
Funktion $u \colon \R_+ \rightarrow \R^n$ zur Zeit $t > 0$ und nehmen
an, dass die Zielfunktion $F$ stetig differenzierbar ist. Damit können
wir das Gradientenabstiegsverfahren schreiben als
```{math}
u(t+\alpha_k ) \ = \ u(t) - \alpha_k \nabla F(u(t)), \qquad \alpha_k > 0.
```
Lassen wir nun die Schrittweiten $\alpha_k > 0$ gegen Null gehen, so
erhalten wir im Grenzwert eine als Gradientenfluss bekannte
Differentialgleichung der Form
```{math}
u'(t) \ = \ - \nabla F(u(t)).
```
Hierbei bleibt die Abstiegseigenschaft erhalten, denn es gilt
```{math}
\frac{d}{dt} F(u(t)) \ = \ \langle \nabla F(u(t)), u'(t) \rangle \ = \ - \Vert \nabla F(u(t)) \Vert^2 \ = \ - \Vert u'(t)\Vert^2 \ \leq \ 0.
```

Ist die Zielfunktion $F$ zusätzlich konvex, d.h. wir wissen
```{math}
\langle v - u, \nabla F(u) \rangle \ \leq \ F(v) - F(u), \qquad \forall \, v,u \in \R^n,
```
dann gilt für einen Minimierer $u^*(t) \equiv u^* \in \R^n$ von $F$ (der
gleichzeitig eine stationäre Lösung des Gradientenflusses darstellt)
sogar
```{math}
\begin{split}
\frac{d}{dt} \left( \frac{1}{2} \Vert u(t) - u^*(t) \Vert^2 \right) \ &= \ \langle u(t) - u^*(t), u'(t) - {u^*}'(t) \rangle\\
&= \ - \langle u(t) - u^*(t), \nabla F(u(t)) - \underbrace{\nabla F(u^*(t))}_{= \, \vec{0}} \rangle \\
&= \ \langle u^*(t) - u(t), \nabla F(u(t)) \rangle \\
&\leq \ F(u^*(t)) - F(u(t)) \ \leq \ 0.
\end{split}
```
Hieran erkennen wir, dass der Abstand der Lösung $u$ zum Minimierer
$u^*$ monoton fällt in der Zeit und wir die Norm der unbekannten Lösung
sogar gleichmäßig durch die Anfangswertbedingung $u(0) = u_0 \in \R^n$
beschränken können.

