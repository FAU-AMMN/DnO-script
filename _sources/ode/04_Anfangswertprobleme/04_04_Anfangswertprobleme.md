(weiterführende-themen)=
# Weiterführende Themen

Im Folgenden diskutieren wir noch einige weiterführende Themen und
Anwendungen zur Numerik von Einschrittverfahren. Dabei beginnen wir
zunächst mit einfachen partiellen Differentialgleichungen und gehen dann
zur Verbindung zwischen Optimierung und Differentialgleichungen über.

(lineare-transportgleichung)=
## Lineare Transportgleichung

Wir betrachten im Folgenden eine *lineare Transportgleichung*, die unter
Anderem verwendet wird um die Ausbreitung eines Stoffes oder eines
Zustands zu modellieren. Wir betrachten dieses Modell der Einfachheit
halber nur in einer räumlichen Dimension mit $n=1$ auf dem
Diskretisierungsgitter $\Omega_h = h \cdot \Z$ und einer örtlichen
Schrittweite $h > 0$. Die Zustandsvariable $u_k(t) \in \R$ beschreibt
einen Zustand im Punkt $k\cdot h \in \Omega_h$ zum Zeitpunkt
$t \in [0, T]$. Wir wählen nun ebenfalls eine äquidistante
Diskretisierung des Zeitintervalls $[0, T]$ mit fester Zeitschrittweite
$\tau > 0$.

Nehmen wir nun an, dass der Zustand mit einer konstanten, positiven
Geschwindigkeit $v \coloneqq \frac{h}{\tau} > 0$, d.h., genau die Breite
einer örtliche Gitterzelle $h$ pro Zeitschritt $\tau$ transportiert
wird, so gilt offensichtlich
```{math}
u_k(t+ \tau) \ = \ u_{k-1}(t).
```
Dies können wir auch durch Erweiterung schreiben als
```{math}
:label: eq:transport_erweiterung
u_k(t+ \tau) \ = \ \underbrace{u_k(t) - u_k(t)}_{= \: 0} + \underbrace{\frac{\tau}{h} v}_{= \: 1 } \cdot \: u_{k-1}(t) \ = \ u_k(t) - \tau \frac{v}h \cdot (u_k(t) - u_{k-1}(t)).
```
Diese Darstellung können wir als *Vorwärts-Euler Verfahren* der
folgenden gewöhnlichen Differentialgleichung mit einer konstanten,
positiven Geschwindigkeit $v \in \R^+$ interpretieren:
```{math}
:label: eq:transport_dgl
u_k'(t) \ = \ - \frac{v}h \cdot (u_k(t) - u_{k-1}(t)).
```
Die rechte Seite hat eine Lipschitz-Konstante der Ordnung
$\mathcal{O}(\frac{1}h)$, welche beliebig groß werden kann für eine
immer feinere örtliche Auflösung, d.h. für $h \rightarrow 0$.

Nach {eq}`eq:transport_erweiterung` gilt
```{math}
u_k(t+\tau) \ = \ (1 - \tau \frac{v}h) \cdot u_k(t) + \tau \frac{v}h \cdot u_{k-1}(t),
```
und wir sehen, dass wir die Stabilität des Einschrittverfahrens zeigen
können solange $\tau \cdot \frac{v}h \leq 1$ gilt. Dann gilt nämlich per
Dreiecksungleichung
```{math}
\vert u_k(t+\tau) \vert \ \leq \ (1 - \tau \frac{v}h) \cdot \vert u_k(t) \vert+ \tau \frac{v}h \cdot \vert u_{k-1}(t) \vert 
\ \leq \ \max\{ \vert u_k(t) \vert,\vert u_{k-1}(t) \vert  \}.
```
Daraus folgt sofort eine Abschätzung für alle örtlichen Gitterpunkte mit
```{math}
\Vert u (t+ \tau) \Vert_\infty \ \leq \ \Vert u (t ) \Vert_\infty.
```

Alternativ können wir ebenfalls das *Rückwärts-Euler Verfahren* für die
gewöhnliche Differentialgleichung {eq}`eq:transport_dgl` betrachten
```{math}
u_k(t+ \tau) \ = \ u_k(t) - \tau \frac{v}h \cdot (u_k(t+\tau) - u_{k-1}(t+\tau)).
```
Wir wollen nun ebenfalls die Stabilität dieses Einschrittverfahrens
näher untersuchen. Der Einfachheit halber betrachten wir nur Lösungen
mit $u_k(0) = 0$ für $k \leq 0$. Es ist klar, dass diese Bedingung dann
auch für alle $t > 0$ gilt. Dies ermöglicht es uns ein gestaffeltes
Gleichungssystem für das implizite Einschrittverfahren herzuleiten und
zu lösen. Es gilt im ersten Ortspunkt zu beliebiger Zeit $t \in [0, T]$
```{math}
\begin{split}
&u_1(t+\tau) \ = \ u_1(t) - \frac{\tau v}{h} \cdot u_1(t+\tau) - \frac{\tau v}{h} \cdot \underbrace{u_0(t+\tau)}_{= \:0}\\
\Leftrightarrow \quad &\left( 1 + \frac{\tau v}{h} \right) \cdot u_1(t+\tau) \ = \ u_1(t)\\
\Leftrightarrow \quad &u_1(t+\tau) \ = \ \frac{h}{h+v\tau} \cdot u_1(t).
\end{split}
```
Für beliebige Punkte $k \cdot h \in \Omega_h$ mit $k > 1$ gilt dann
```{math}
\begin{split}
&u_k(t+\tau) \ = \ \ u_k(t) - \frac{\tau v}{h} \cdot u_k(t+\tau) + \frac{\tau v}{h} \cdot u_{k-1}(t+\tau)\\
\Leftrightarrow \quad &\left( 1 + \frac{\tau v}{h} \right) \cdot u_k(t+\tau) \ = \ u_k(t) + \frac{\tau v}{h} \cdot u_{k-1}(t + \tau)\\
\Leftrightarrow \quad &u_k(t+\tau) \ = \ \frac{h}{h+v\tau} u_k(t) +  \frac{v\tau}{h+v\tau} u_{k-1}(t+\tau).
\end{split}
```
Wir sehen also direkt, dass wir mittels Dreiecksgleichung folgende
Abschätzung treffen können
```{math}
\vert u_k(t+\tau) \vert \ \leq \ \max\{ \vert u_k(t) \vert,  \vert u_{k-1}(t+\tau) \vert\}.
```
Somit folgt induktiv
```{math}
\vert u_k(t+\tau) \vert \ \leq \ \max_{i \leq k} \vert u_i(t) \vert.
```
Die impliziert wieder insbesondere die Stabilitätsabschätzung für alle
örtlichen Gitterpunkte mit
```{math}
\Vert u (t+ \tau) \Vert_\infty \ \leq \ \Vert u (t ) \Vert_\infty.
```
Man beachte, dass in diesem Fall das Einschrittverfahren stabil ist ohne
jegliche Beschränkung an die Zeitschrittweite $\tau > 0$.

Für $h \rightarrow 0$ sehen wir, dass
```{math}
v\cdot \frac{u(k\cdot h,t) - u((k-1)\cdot h,t)}{h} \ \longrightarrow \ v\cdot \partial_x u(x, t)
```
gilt und wir eigentlich eine *partielle Differentialgleichung*
approximiert haben, nämlich die sogenannte **lineare
Transportgleichung** der Form
```{math}
\partial_t u (x,t ) \ = \ - v \cdot \partial_x u(x,t),
```
mit konstanter, positiver Geschwindigkeit $v \in \R^+$.

Im obigen Fall haben wir also die partielle Ableitung in die
Ortskoordinate $x$ durch einen Rückwärtsdifferenzenquotienten
approximiert. Analog könnten wir auch ein Verfahren mit
Vorwärtsdifferenzenquotienten bezüglich der Ortskoordinate $x$
aufschreiben, was zu folgender gewöhnlichen Differentialgleichungen
führt:
```{math}
u_j'(t) \ = \ \frac{v}h \cdot (u_j(t) - u_{j+1}(t)).
```
Unabhängig von der Zeitdiskretisierung ist hier allerdings das
Differentialgleichungssystem schon instabil, wie man formal zeigen kann.
Der Grund hierfür liegt anschaulich betrachtet in der ursprünglichen
Motivation der Transportgleichung. Mit positiver Geschwindigkeit
$v \in \R^+$ beschreibt die Transportgleichung die Ausbreitung eines
Zustands in Richtung steigender Werte der Ortskoordinate $x$. Hierzu
benötigt man ausschließlich Informationen von der linken Seite eines
örtlichen Punkts der Diskretisierung $\Omega_h$. Dies wird durch den
Rückwärtsdifferenzenquotienten korrekt abgebildet, während der
Vorwärtsdifferenzenquotient genau die entgegengesetzte Richtung
verwendet.

(diffusionsgleichung)=
## Diffusionsgleichung

Wir betrachten im Folgenden einen einfachen Diffusionsmodell im
eindimensionalen Raum für $n=1$. Solche Modelle beschreiben
beispielsweise in der Physik die Verteilung eines Stoffes in einem
Medium oder die Ausbreitung von Wärme in einem Material. Wir modellieren
ein Teilchen, dass einen Sprungprozess auf dem Gitter
$\Omega_h \coloneqq h \cdot \Z \cap [0,1]$ mit Schrittweite
$h \coloneqq \frac{1}{N} > 0, N \in \N$ und periodischen Randbedingungen
durchführt, wobei es zu jedem Zeitpunkt $t \in [0,T]$ mit gleicher
Wahrscheinlichkeit $\alpha \in (0, \frac{1}{2})$ nach links oder rechts
springen kann. In diesem Zusammenhang nennt man $\alpha$ einen
Diffusionskoeffizienten für den Sprungprozess. Die Wahrscheinlichkeit
für einen Sprung in einem kleinen Zeitintervall
$[t,t+\tau]\subset [0,T]$ mit Zeitschrittweite $0 < \tau < 1$ ist
dementsprechend $2 \alpha \tau \in (0, 1)$. Dann gilt für die
Wahrscheinlichkeit $p_k(t) \in [0,1]$, dass das Teilchen zur Zeit $t$ im
Gitterpunkt $k\cdot h \in \Omega_h$ ist:
```{math}
:label: eq:diffusion_wahrscheinlichkeit
p_k(t+\tau) \ = \  \alpha \tau \cdot p_{k-1}(t) +  \alpha \tau \cdot p_{k+1}(t) +  (1- 2 \alpha \tau) \cdot p_k(t), \qquad k=0,\ldots,N,
```
wobei wegen der periodischen Randbedingungen $p_{-1}(t) = p_N(t)$ und
$p_{N+1}(t) = p_0(t)$ gilt.

Für immer kleiner werdende Zeitschrittweiten $\tau \rightarrow 0$
erhalten wir entsprechend im Grenzwert das Differentialgleichungssystem
```{math}
:label: eq:diffusion_dgl
p_k'(t) \ = \ \alpha \cdot (p_{k-1}(t) + p_{k+1}(t) - 2 p_k(t)), \qquad k=0,\ldots,N.
```
In diesem Zusammenhang sehen wir ein, dass
{eq}`eq:diffusion_wahrscheinlichkeit` als Vorwärts-Euler Verfahren zur
numerischen Approximation der gewöhnlichen Differentialgleichung
{eq}`eq:diffusion_dgl` interpretiert werden kann. Wir erkennen sofort,
dass dieses Einschrittverfahren stabil ist für $2 \alpha \tau < 1$. Dies
ist potentiell eine starke Einschränkung an die Zeitschrittweite $\tau$
in Fällen in denen der Diffusionskoeffizient $\alpha$ groß ist.

Insbesondere können wir mit der Hilfsvariable
$D \coloneqq \alpha \cdot h^2$ das Einschrittverfahren wieder als
Ortsdiskretisierung einer partiellen Differentialgleichung mittels eins
Differenzenquotient zweiter Ordnung der folgenden Form interpretieren
```{math}
\partial_t p(x,t) \ = \ D \cdot \partial_{xx} p(x,t).
```
Hierbei approximiert man die zweite Ableitung nach der Ortskoordinate
durch
```{math}
\partial_{xx} p(k\cdot h,t) \ \approx \ \frac{p_{k-1}(t) - 2p_k(t) + p_{k+1}(t)}{h^2}.
```
Somit erhalten wir Stabilität des Einschrittverfahrens für
$\tau \sim h^2$, was für sehr kleines $h > 0$ (also für eine sehr feine
Ortsauflösung) einen hohen numerischen Aufwand fordert.

Verwenden wir hingegen das *Rückwarts-Euler Verfahren* zur
Zeitdiskretisierung von {eq}`eq:diffusion_dgl`, so erhalten wir für
$k=0,\ldots,N$
```{math}
:label: eq:diffusion_implizit
p_k(t+\tau ) \ = \ p_k(t) + \alpha \tau \cdot p_{k-1}(t+\tau) + \alpha \tau \cdot p_{k+1}(t+\tau)  -2\alpha \tau \cdot p_k(t + \tau)  .
```
Es stellt sich heraus, dass dieses implizite Einschrittverfahren
wiederum stabil ist ohne Schranke an die Zeitschrittweite $\tau > 0$.
Dies sehen wir folgendermaßen: Sei im Folgenden
$P(t) \coloneqq (p_0(t),\ldots,p_N(t))^T \in \R^{N+1}$, dann können wir
das Einschrittverfahren {eq}`eq:diffusion_implizit` kompakt schreiben
als
```{math}
:label: eq:diffusion_implizit_kompakt
(I + \alpha \tau \cdot B) \cdot P(t+ \tau) \ = \ P(t),
```
mit
```{math}
B = \left( \begin{array}{cccccc}
2 & -1 & 0 & \ldots & 0 & -1\\
-1 &2 & -1 & \ldots & 0 & 0 \\
0 & -1 & 2 & \ldots & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\   
0 & 0 & 0 & \ldots &  2 & -1 \\ 
-1 & 0 & 0 & \ldots &  -1 & 2  
\end{array} \right).
```
Wie man leicht zeigt, gilt für jeden Vektor $x \in \R^{N+1}$
```{math}
\begin{split}
x^T B x \ &= \ x^T \cdot 
\begin{pmatrix} 
2x_0 - x_1 -x_N\\
-x_0 + 2x_1 -x_2\\
\vdots\\
-x_{k-1} + 2x_k - x_{k+1}\\
\vdots\\
-x_0 -x_{N-1} + 2x_N
\end{pmatrix}\\
&= \ x_0 \cdot (2x_0 - x_1 -x_N) \\
& \hspace{0.19cm} + x_1 \cdot (-x_0 + 2x_1 -x_2)\\
& \hspace{0.19cm} + \ldots \\
& \hspace{0.19cm} + x_N \cdot (-x_0 -x_{N-1} + 2 x_N)\\
&= \ \sum_{i=0}^N (x_{i+1} - x_i)^2 \ \geq \ 0.
\end{split}
```
Damit haben wir gezeigt, dass die Matrix $B\in \R^{(N+1) \times (N+1)}$
positiv semidefinit ist. Dies hätten wir auch mit Hilfe der
Gerschgorin-Kreise (vgl. {cite:p}`numerik1`) sehen können, da wir an der
Hauptdiagonalen von $B$ und der Summe der Nebendiagonalen ablesen
können, dass alle Eigenwerte der Matrix $B$ in der Menge
$B_2(2) \coloneqq \lbrace x \in \R : |x - 2| \leq 2 \rbrace$ liegen
müssen. Da $B$ also positiv semidefinit ist, ist
$(I+\tau B) \in \R^{(N+1) \times (N+1)}$ für jedes $\tau > 0$ positiv
definit und invertierbar. Insbesondere erhalten wir folgende
Stabilitätsabschätzung indem wir {eq}`eq:diffusion_implizit_kompakt` von
links mit dem Vektor $P(t + \tau)^T$ multiplizieren:
```{math}
\begin{split}
\Vert P(t+\tau) \Vert_2^2 +  \alpha \tau \cdot P(t+\tau)^T B P(t+\tau) \ &= \ P(t+\tau)^T P(t)\\
&\leq \ \frac{1}2  \Vert P(t+\tau) \Vert_2^2 + \frac{1}2  \Vert P(t) \Vert_2^2.
\end{split}
```
Die Ungleichung folgt aus der Einsicht, dass gilt
```{math}
0 \ \leq \ \frac{1}{2} (P(t+\tau) - P(t))^2 \ = \ \frac{1}{2}||P(t+\tau)||^2 - P(t+\tau)^T P(t) + \frac{1}{2}||P(t)||^2.
```
Wegen der positiven Semidefinitheit von $B$ folgt dann
```{math}
\begin{split}
\Vert P(t+\tau) \Vert_2^2 \ &\leq \ \frac{1}2  \Vert P(t+\tau) \Vert_2^2 + \frac{1}2  \Vert P(t) \Vert_2^2\\
\Leftrightarrow \quad \frac{1}2  \Vert P(t+\tau) \Vert_2^2 \ &\leq \ \frac{1}2  \Vert P(t) \Vert_2^2.
\end{split}
```
Somit können wir nun induktiv folgende Ungleichungskette folgern:
```{math}
\Vert P(t+\tau) \Vert_2 \ \leq \ \Vert P(t ) \Vert_2 \ \leq \ \ldots \ \leq \ \Vert P(0) \Vert_2.
```
Tatsächlich gilt in diesem Fall sogar folgende Abschätzung
```{math}
\min_{k=0,\ldots,N} p_k(0) \ \leq \ \min_{k=0,\ldots,N} p_k(t) \ \leq \ \max_{k=0,\ldots,N} p_k(t) \ \leq \ \max_{k=0,\ldots,N} p_k(0).
```
Abschätzungen dieser Form nennt man im Kontext von partiellen
Differentialgleichungen **Minimums-** und **Maximumsprinzip** und wir
werden ähnliche Argumente im nächsten Kapitel zu Randwertaufgaben näher
diskutieren.

(adjungierte-methode)=
### Adjungierte Methode

Nehmen wir nun an, dass wir uns einen Startpunkt
$x_{k_0} = k_0 \cdot h \in \Omega_h$ des modellierten Teilchens
vorgeben, d.h., die Wahrscheinlichkeit für ein Auftreten des Teilchens
ist zu Anfang $p_{k_0}(0) = 1$ und $p_i(t) = 0$ für alle $i \neq k$.
Dann lässt sich mit Hilfe einer numerischen Lösung der
Differentialgleichung {eq}`eq:diffusion_dgl` im Intervall $t \in (0,T)$
effizient berechnen was die Auftrittswahrscheinlichkeit $p_k(T)$ des
Teilches in allen Punkten $k=0,\ldots,N$ für den Zeitpunkt $T$ ist.

Schwieriger ist jedoch die Beantwortung der umgekehrten Frage. Man
könnte sich die Frage stellen, ob man für einen spezifischen Punkt
$x_{k_0} = k_0 \cdot h \in \Omega_h$ die Wahrscheinlichkeit berechnen
kann, dass das Teilchen zum Zeitpunkt $T$ dort auftritt, d.h., wir
interessieren uns für $p_{k_0}(T)$ für alle möglichen Startpunkte
$k = 0,\ldots,N$ des Teilchens. Hierzu müssten wir eigentlich das
Differentialgleichungssystem mit allen $N+1$ Startpunkten des Teilchens
lösen und die Lösung dann in dem spezifischen Punkt $p_{k_0}(T)$
auswerten.

Eine alternative Berechnungsmöglichkeit für die letzte Fragestellung ist
die sogenannte **adjungierte Methode**, welche aus der Theorie der
optimalen Steuerung stammt. Dazu berechnen wir die Lösung des
adjungierten Problems, welches gegeben ist durch
```{math}
q_k'(t) \ = \ \alpha \cdot (2 q_k(t) - q_{k+1}(t) - q_{k-1}(t)), \qquad k = 0,\ldots,N,
```
mit vorgegebenem Endwert $q_{k_0}(T) = 1$ und $q_i(T) = 0$ für
$i\neq k_0$.

Die Lösung dieses Problems ist genauso zu berechnen wie für das
ursprüngliche System {eq}`eq:diffusion_dgl`. Dies sehen wir mit der
einfachen Variablentransformation $s \coloneqq T-t$ ein, denn dann haben
wir ein Anfangswertproblem mit den gleichen Vorzeichen wie das System
für die Funktionen $p_k(t)$ zu lösen. Hat man die eine numerische Lösung
berechnet, so gilt
```{math}
:label: eq:adjungierte_identität
\begin{split}
p_{k_0}(T) \ &= \ \sum_{i=0}^N p_i(T) \cdot q_i(T) \\
&= \ \sum_{i=0}^N p_i(0) \cdot q_i(0) + \sum_{i=0}^N \int_0^T (p_i(t) \cdot q_i(t))'\,\mathrm{d}t \\
&= \ \sum_{i=0}^N p_i(0) \cdot q_i(0) + \int_0^T \underbrace{\sum_{i=0}^N (p_i'(t) \cdot q_i(t) + p_i(t) \cdot q_i'(t))}_{=\:0}\mathrm{d}t \\
&= \ \sum_{i=0}^N p_i(0) \cdot q_i(0).
\end{split}
```
Bei den obigen Umformungen haben wir ausgenutzt, dass die auftretenden
Terme sich in der Summe wegheben, d.h. es gilt:
```{math}
\begin{split}
&\sum_{i=0}^N (p_i'(t) \cdot q_i(t) + p_i(t) \cdot q_i'(t)) \\
= \ &\sum_{i=0}^N \alpha \cdot (p_{i-1}(t) + p_{i+1}(t) - 2p_i(t)) \cdot q_i(t) + \alpha \cdot p_i(t) \cdot (2q_i(t) - q_{i-1}(t) - q_{i+1}(t)) \\
= \ &0.
\end{split}
```
Nun sind wir in der Lage die obige Frage effizient zu beantworten, denn
wenn wir nun eine Lösung mit Anfangswert $p_\ell(0)=1$ und $p_i(0) = 0$
für $i \neq \ell$ einsetzen, so liefert uns
{eq}`eq:adjungierte_identität` die Identität $p_{k_0}(T) = q_\ell(0)$.
Im Gegensatz zur direkten Berechnung müssen wir nun wieder nur ein
Differentialgleichungssystem für die unbekannten Lösungen
$q_i(t), i=0,\ldots, N$ lösen.

Das zu Grunde liegende allgemeine Prinzip der adjungierten Methode ist
das Folgende. Wir nehmen an, wir haben eine gewöhnliche lineare
Differentialgleichung der Form
```{math}
u'(t) \ = \ A(t) \cdot u(t)
```
gegeben und wir interessieren uns nicht direkt für Werte der Lösung zum
Endzeitpunkt $u(T) \in \R^n$, sondern nur für eine lineare Funktion
$L^T \cdot u(T) \in \R$. Dann lösen wir stattdessen das adjungierte
Problem
```{math}
:label: eq:adjungiertes_problem
v'(t) \ = \ - A(t) \cdot v(t)
```
mit dem Endwert $v(T) = L \in \R^n$. Denn dann gilt
```{math}
L^T \cdot u(T) \ &= \ v(T)^T \cdot u(T) \ = \ v(0)^T \cdot u(0) + \int_0^T (v(t)^T \cdot u(t))'\, \mathrm{d}t \\
&= \ v(0)^T \cdot u(0) + \int_0^T \underbrace{(v'(t)^T \cdot u(t) + v(t)^T \cdot u'(t))}_{=\:0} \, \mathrm{d}t\\
&= \  v(0)^T \cdot u(0).

```
Hierbei nutzt man aus, dass gilt
```{math}
v'(t)^T \cdot u(t) + v(t)^T \cdot u'(t) = - (A(t)^T v(t))^T \cdot u(t) + v(t)^T \cdot A(t) u(t) \ = \ 0.
```
Haben wir die adjungierte Gleichung für $v$ berechnet, können wir die
uns interessierende Größe $L^T \cdot u(T)$ für jeden Anfangswert sofort
durch ein Skalarprodukt $v(0)^T \cdot u(0)$ berechnen. Den Wert
$v(0) \in \R^n$ erhalten wir in dem wir eine Variablentransformation
$s \coloneqq T- t$ in {eq}`eq:adjungiertes_problem` durchführen und das
entstehende Differentialgleichungssystem numerisch lösen.

Bei einer numerischen Lösung müssen wir natürlich die gewöhnliche
Differentialgleichung mit einer geeigneten Methode (Ein- oder
Mehrschrittverfahren) diskretisieren. Es empfiehlt sich in diesem
Kontext die adjungierte Gleichung mit einem passenden Verfahren zu lösen
um die Eigenschaft der Adjungierten auch im Diskreten zu erhalten. Haben
wir für die numerische Approximation $u$ z.B. ein Vorwärts-Euler
Verfahren der Form
```{math}
u(t_{k+1}) \ = \ u(t_k) + \tau \cdot A u(t_k)
```
verwendet, so liefert das Vorwärts-Euler Verfahren in umgekehrter Zeit
```{math}
v(t_k) \ = \ v(t_{k+1}) + \tau A^T v(t_k)
```
genau die richtige Diskretisierung. Es gilt dann nämlich
```{math}
    u(t_N) \cdot v(t_N) \ &= \ u(0) \cdot v(0) + \sum_{k=0}^{N-1} (u(t_{k+1}) \cdot v(t_{k+1}) - u(t_{k }) \cdot v(t_{k })) \\
     &= \ u(0) \cdot v(0) + \sum_{k=0}^{N-1} ((u(t_{k+1})  - u(t_{k }) \cdot v(t_{k+1}) + (v(t_{k+1}) - v(t_{k })) \cdot u(t_k) ) \\
     &= \ u(0) \cdot v(0) + \sum_{k=0}^{N-1} (A u(t_{k })\cdot v(t_{k+1}) - (A^T v(t_{k })) \cdot u(t_k) ) \ = \ u(0)\cdot v(0) .

```
Man sieht ein, dass bei anderen Verfahren für die adjungierte Gleichung,
z.B. einem impliziten Euler-Verfahren, eine solche Identität nicht
gegeben ist. Die Diskretisierung und Adjungierung kommutieren in diesem
Fall nicht. Beim Vorwärts-Euler Verfahren hingegen erhalten wir genau
die Adjungierte der Diskretisierung der Differentialgleichung.

(zusammenhang-zwischen-optimierung-und-differentialgleichungen)=
## Zusammenhang zwischen Optimierung und Differentialgleichungen

Ein häufiges Problem in praktischen Anwendungen ist die Bestimmung von
Parametern in gewöhnlichen Differentialgleichungen. Wir betrachten also
ein Anfangswertproblem
```{math}
u'(t) \ = \ F(t,u(t),w), \qquad u(0) \, = \, u_0(w),
```
bei dem die rechte Seite der gewöhnlichen Differentialgleichung und der
Anfangswert von einem Parametervektor $w \in \R^m$ abhängen. Wir nehmen
im Folgenden an, dass $F$ Lipschitz-stetig ist. Dann existiert für
gegebene Parameter $w \in \R^m$ eine eindeutige Lösung
$u_w \in C^1([0,T])$. Um die Parameter zu einer bestimmten Lösung zu
bestimmen, misst man die Werte einer Funktion $G(u_w) \in R^k$, wobei
typischerweise $k >m$ gilt. Alternativ versucht man die Parameter
$w \in \R^m$ derart zu optimieren, damit ein gewünschter Zustand
$G(u_w) \in \R^k$ erreicht wird. Häufig sind dies Werte der Lösung zu
verschiedenen Zeitpunkten, also beispielweise
$G(u) \coloneqq (u(s_1), \ldots, u(s_k))$.

Für gegebenen Daten $g \in \R^k$ lässt sich dann ein
Optimierungsproblem, etwa das **Kleinste-Quadrate-Problem**
```{math}
\min_{w \in \R^m} \left\lbrace E(w) \, \coloneqq \, \frac{1}2 \Vert H(w) - g \Vert^2  \, = \, \frac{1}2 \Vert G(u_w) - g \Vert^2 \right\rbrace
```
lösen um die Parameter zu bestimmen. Die Frage, die wir uns nun stellen
ist, wie wir in diesem Fall die Lösung des Optimierungsproblems durch
eines der Optimierungsverfahren dieser Vorlesung, wie z.B. das
Gradientenabstiegsverfahrens, bestimmen können. Es ist klar, dass
hierbei die effiziente Berechnung von Gradienten $\nabla_w$ essentiell
ist.

````{prf:example} Parameterbestimmung für lineare DGL
:label: ex:parameterbestimmung_linDGL
Als einfaches Beispiel betrachten wir einen Fall, bei der wir zwar
wissen, dass die rechte Seite zu einer linearen gewöhnlichen
Differentialgleichung erster Ordnung gehört, aber nicht den
Koeffizienten des linearen Terms und auch den Anfangswert nicht kennen.
Somit können wir zwei unbekannte Parameter $w = (w_1, w_2) \in \R^2$
beschreiben durch
```{math}
u_0(w) \: = \: w_1, \qquad F(t,u,w) \: = \: w_2 \cdot u(t).
```
Für Anfangswertprobleme dieser Art können wir eine explizite Lösung
$u_w \in C^1([0,T])$ angeben mit
```{math}
u_w(t) \ = \ w_1 \cdot e^{w_2t}.
```
Geben wir uns nun Werte der Lösung an verschiedenen Zeitpunkten
$G(u) = (u(s_1), \ldots, u(s_k))$ vor, so erhalten wir dann
```{math}
\begin{split}
\partial_{w_1} G(u_w) \ &= \ (e^{w_2 s_1},\ldots,e^{w_2 s_k}), \\
\partial_{w_2} G(u_w) \ &= \ (w_1 s_1 \cdot e^{w_2 s_1},\ldots,w_1 s_k \cdot e^{w_2 s_K}).
\end{split}
```

````

In {prf:ref}`ex:parameterbestimmung_linDGL` konnten wir die Ableitung
der Funktion $G$ nach den unbekannten Parametern $w$ angeben, da wir
eine explizite Lösung der gewöhnlichen Differentialgleichung kannten.
Wie gehen wir aber vor, wenn wir die Differentialgleichung nicht
explizit lösen können? Dazu betrachten wir zunächst die partielle
Ableitung von $u_w$ nach $w_i$, welche wir definieren als
```{math}
u^i_w \ \coloneqq \ \lim_{\delta \rightarrow 0} \frac{u_{w+\delta e_i}(t)- u_w(t)}{\delta},
```
wobei $e_i$ der $i$-te Einheitsvektor ist. Diese Funktion können wir
zwar nicht berechnen, aber wir können ein Anfangswertproblem herleiten,
das von ihr gelöst wird. Unter der Annahme, dass die Abbildung
$w\mapsto u_0(w)$ differenzierbar ist, sehen wir dass gilt
```{math}
u_w^i(0) \ = \ \partial_{w_i} u_0(w).
```
Falls außerdem die rechte Seite $F$ der Differentialgleichung
differenzierbar bezüglich $u$ und $w$ ist, so folgt mit der Kettenregel
```{math}
(u_w^i)'(t) \ = \ \partial_u F(t,u_w(t),w) \cdot u_w^i(t) + \partial_w F(t,u_w(t),w).
```
Wir beachten, dass wir diese lineare Differentialgleichungen für jedes
$u_w^i$ sind, da wir $u_w$ ja schon vorher durch Lösen der
ursprünglichen Anfangswertproblems berechnen können. Ist $G$ ebenfalls
differenzierbar, dann folgt
```{math}
\partial_{w_i} H(w) \ = \ G'(u_w) \cdot u_w^i,
```
daraus bekommen wir also die Jacobi Matrix von $H$ bzw. dann den
Gradienten von $f$ per Kettenregel.

Wir rechnen diesen Zusammenhang nun für das Problem in
{prf:ref}`ex:parameterbestimmung_linDGL` nach. Hier gilt
```{math}
u_w^1(0) \: = \: 1, \qquad u_w^2(0) \: = \: 0,
```
und
```{math}
(u_w^1)'(t) \ = \ w_2 u_w^1(t), \qquad (u_w^2)'(t) \ = \ w_2 u_w^2(t) + u_w,
```
und $\partial_{w_i} G(u) = (u_w^i(s_1),\ldots,u_w^i(s_K))$. Diese können
wir explizit lösen und erhalten $u_w^1(t) = e^{w_2 t}$ und
$u_w^2(t) = w_1 t e^{w_2 t}$. Natürlich stimmen die Ableitungen dann
wieder mit der direkten Differentiation der expliziten Lösung $u_w$ wie
in {prf:ref}`ex:parameterbestimmung_linDGL` berechnet überein.

Ist die Anzahl $M$ der Parameter groß, so ist die Berechnung der
Ableitungen in dieser Form sehr aufwändig, da wir $M$ lineare
Differentialgleichungen lösen müssen. Dies kann aber vermieden werden,
wenn wir uns daran erinnern, dass wir eigentlich
```{math}
\partial_{w_i} f(w) \ = \ (G(u_w) - g) \cdot G'(u_w) u_w^i
```
berechnen wollen, also eine lineare Funktion von $u_w^i$. Es ist
dementsprechend naheliegend wieder eine adjungierte Methode zu
verwenden. Wir betrachten dies wieder näher für
$G(u) = (u_w(s_1),\ldots,u_w(s_K)).$ Wir definieren $v$ als die Lösung
von
```{math}
v'(t) \ = \ - \partial_u F(t,u_w(t),w) \cdot v(t), \qquad t \in (0,T) \setminus \{s_1,\ldots,s_K\},
```
mit $v(T) =0$. An den Messstellen setzen wir
```{math}
v(s_k) \ = \ u_w(s_k) - g_k +  \lim_{t \downarrow s_k} v(t) .
```
Dann gilt mit $s_0=0$, $s_{K+1}=T$,
```{math}
\partial_{w_i} f(w) \ &= \ \sum_k (u_w(s_k) - g_k) \cdot u_w^i(s_k) \\
&= \ \sum_k (\lim_{t \uparrow s_k} v(t) - \lim_{t \downarrow s_k} v(t)) \cdot u_w^i(s_k) \\
&= \ v(0) \cdot u_w^i(0) + \sum_{k=0}^K \int_{s_k}^{s_{k+1}} (v(t) \cdot u_w^i(t))'\,\mathrm{d}t  \\
&= \ v(0) \cdot u_w^i(0) + \int_0^T (v(t) \cdot u_w^i(t))'\,\mathrm{d}t  \\
&= \ v(0) \cdot \partial_{w_i} u_0(w) + \int_0^T v(t) \cdot \partial_{w_i} F(t,u_w(t),w)\,\mathrm{d}t.

```
Damit genügt zur Berechnung des Gradienten die Lösung einer adjungierten
Differentialgleichung, sowie von $M$ Skalarprodukten mit Anfangswerten
und $M$ Integralen mit der Lösung $v$.

(deep-learning)=
## Deep Learning

In modernen Anwendungen des Maschinellen Lernens kommen ähnliche
Techniken wie bei Differentialgleichungen und deren Optimierung zum
Einsatz. Die Idee dabei ist einen parametrisierten Zusammenhang zwischen
Eingangsdaten $x \in \R^n$ und Ausgangsdaten $y \in \R^m$ zu
konstruieren. Dieser Zusammenhang wird beim sogenannten *Deep Learning*
durch ein (künstliches) neuronales Netz mit vielen Schichten (im
Englischen: *Layer*) modelliert, das mathematisch als
Hintereinanderausführung von affinen Abbildungen (Austausch von Impulsen
zwischen Neuronen) und punktweisen Nichtlinearitäten (Aktivierung eines
Neurons durch die eingegangenen Impulse) modelliert.

Ein **neuronales Netzwerk** $f_\Theta \colon \R^n \rightarrow \R^m$ mit
$N \in \N^+$ Schichten $f^k_{\Theta_k}, k=1,\ldots,N$ modelliert die
Relation $x \mapsto y$ dann durch
```{math}
f_\Theta \ \coloneqq \ f_{\Theta_N}^N \circ \ldots \circ f_{\Theta_1}^1,
```
wobei die berechneten Werte der $k$-ten Schicht des neuronalen Netzes
gegeben sind durch
```{math}
u_{k+1} \ = \ f_{\Theta_k}^k(u_k) \ \coloneqq \ \Psi(W_k u_k + b_k) , \qquad k=0,\ldots,L-1,
```
mit frei wählbaren Parametern
$\Theta_k \coloneqq \lbrace W_k \in \R^{n_k \times n_k}, b_k \in \R^{n_k}\rbrace$.
Die **Aktivierungsfunktion** eines Neurons bezeichnen wir mit
$\Psi: \R \rightarrow \R$. Typische Beispiele sind die
*Sigmoid-Funktion*
```{math}
\Psi(x) \: \coloneqq \: \frac{1}{1+e^{-x}}
```
oder die *Rectified Linear Unit (ReLU)*
```{math}
\Psi(x) \: \coloneqq \: \max\{x,0\} .
```
Dazu verwenden wir für Vektoren die Notation
$\Psi(x) = (\Psi(x_i))_{i=1,\ldots,n}$ als punktweise Auswertung. Für
die erste Schicht des neuronalen Netzes setzen wir $u_0=x$ auf die
Eingangsdaten und in der letzten Schicht haben wir für die Antwort des
neuronalen Netzes lediglich eine affine Transformation der Form
```{math}
y \: = \: C \cdot u_N + d,
```
mit $C \in \R^{m \times n_N}$, $d\in \R^m$.

Wir sehen also eine gewisse Analogie zu expliziten Euler-Verfahren für
gewöhnliche Differentialgleichungen. Dies wird noch deutlicher bei
sogenannten residualen Netzwerken von der Form
```{math}
u_{k+1} \ = \ u_k + \tau \cdot \Psi(W_k \cdot u_k + b_k) , \qquad k=0,\ldots,L-1,
```
die wir direkt als Diskretisierung von
```{math}
u'(t) \ = \ \Psi(W(t) \cdot u(t) + b(t))
```
interpretieren können.

Das Training einens neuronalen Netzwerks ist nun die optimale Bestimmung
der freien Parameter $\Theta = ((W_k, b_k)_{k=1,\ldots,N}, C, d)$ aus
einer großen Menge an Trainingsdaten $(x_i,y_i)_{i=1,\ldots,M}$. Dazu
wird ein Minimierungsproblem der Form
```{math}
\min_{\Theta} \left\lbrace E(\Theta) \ \coloneqq \ \frac{1}M \sum_{i=1}^M \mathcal{L}(f_\Theta(x_i),y_i) \right\rbrace,
```
wobei $\mathcal{L}$ eine Metrik ist, die den Abstand zwischen der
Antwort des neuronalen Netzes und den Trainingsdaten misst (im
Englischen: *Loss function*), z.B. einfach
```{math}
\mathcal{L}(f_\Theta(x_i),y_i) \ \coloneqq \ \frac{1}2 \Vert f_\Theta(x_i) - y_i \Vert^2.
```
Dies ist für große $M / N, m / n$ ein riesiges Optimierungsproblem,
dessen approximative Lösung lange Zeit ein großes Hindernis bei der
Umsetzung solcher Lernansätze war. Der heute gängige Ansatz ist die
Berechnung mit einem stochastischen Gradientenverfahren, d.h. durch eine
Iteration
```{math}
w_{j+1} \ = \ w_j - \alpha_j \nabla_w \mathcal{L}(C \cdot u_N(x_{\pi(j)};(W_k,b_k))+d,y_{\pi(j)}),
```
wobei $\pi(j) \in \{1,\ldots,M\}$ zufällig gewählt wird (meist
gleichverteilt). Durch die Auswahl eines einzelnen Datenpaars in jeder
Iteration erspart man sich die $M$-fache Berechnung von
$u_N(x_{i };(A_k,b_k))$, hier muss in jedem Schritt die
Vorwärts-Schleife nur für einen Anfangswert $x_i$ berechnet werden.
Analog zur adjungierten Methode kann man den Gradienten durch Lösung von
```{math}
v_{k-1} \ = \ -(W_k \cdot \Psi'(W_k u_k + b_k))^T \cdot v_k
```
mit
$v_N = \nabla_u \mathcal{L}(C \cdot u_N(x_{i(j)};(W_k,b_k))+d,y_{i(j)})$
berechnen. Dies ist in diesem Zusammenhang als **Backpropagation**
bekannt.
