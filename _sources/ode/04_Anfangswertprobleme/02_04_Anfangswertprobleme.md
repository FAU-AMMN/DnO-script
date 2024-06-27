(s:einschrittverfahren)=
# Einschrittverfahren für Anfangswertprobleme

In diesem Abschnitt wollen wir zunächst eine einfache Möglichkeit zur
numerischen Berechnung von Lösungen gewöhnlicher Differentialgleichungen
der folgenden Form eines Anfangswertproblems behandeln:
```{math}
:label: eq:awp_einschritt
\begin{split}
u'(t) \ &= \ F(t,u(t)), \qquad \forall \ 0 \leq t \leq T,\\ 
u(0) \: &= \: u_0 \in \R^n,
\end{split}
```
Die grundlegende Idee ist es sukzessive (d.h., aufeinander aufbauend)
die Lösung zu verschiedenen Zeitschritten zu approximieren. Hierzu
führen wir eine sogenannte **Zeitdiskretisierung der
Differentialgleichung** durch indem wir das kontinuierliche Intervall
$[0,T] \subset \R^+$ in $N\in\N$ Teilintervalle mit $N+1$ Zeitpunkten
aufteilen. Der Einfachheit halber wählen wir im Folgenden uniforme
Zeitschritte mit einer festgewählten Schrittweite, d.h., wir
diskretisieren das Intervall $\Omega \coloneqq [0,T]$ durch
```{math}
\Omega_\tau \ \coloneqq \ \lbrace t_k \in \Omega \, | \, t_k \coloneqq k \cdot \tau, \ k=0,\ldots,N \rbrace \quad \text{mit} \quad \tau \, \coloneqq \, \frac{T}{N}.
```
Es sei bemerkt, dass ein analoges Vorgehen auch für nicht äquidistante
Schrittweiten möglich ist.

Basierend auf der Diskretisierung der Differentialgleichung können wir
dann sukzessive numerische Approximationen $u_\tau(t_k) \in \R^n$ für
die unbekannte Lösung $u$ des Anfangswertproblems zu diesen diskreten
Zeitschritten berechnen. Hierbei haben wir die Möglichkeit zeitliche
Ableitung $\frac{d}{dt}$ durch geeignete Differenzenquotienten zu den
diskreten Zeitpunkten zu approximieren. Alternativ lässt sich auch eine
numerische Quadraturformel (siehe {cite:p}`numerik1`) für die zugehörige
Integralgleichung in diesen Zeitschritten anwenden. Hierbei macht man
sich die folgende auf dem Hauptsatz der Differential- und
Integralrechnung basierende Beobachtung zu Nutze:
```{math}
:label: eq:anfangswertaufgabe_integralform
u(t_{k+1}) \ = \ u(t_k) + \int_{t_k}^{t_{k+1}} u'(t) \,\mathrm{d}t \ = \ u(t_k) + \int_{t_k}^{t_{k+1}} F(t,u(t))\,\mathrm{d}t , \qquad t_k \in \Omega_\tau.
```

Im Falle von Einschrittverfahren, die wir in diesem Abschnitt zunächst
behandeln wollen, verwenden wir dabei eine Differenzenformel, die nur
einen einzigen Zeitschritt berücksichtigt, d.h. zur Berechnung von
$u_\tau(t_{k+1})$ wird lediglich der Wert $u_\tau(t_k)$ des vorherigen
Zeitschritts verwendet.

Im Folgenden wollen wir zunächst eine allgemeine Form für
Einschrittverfahren zur numerischen Lösung von Anfangswertproblemen
angeben.

````{prf:definition} Einschrittverfahren und Verfahrensfunktion
:label: def:einschrittverfahren
Sei eine Zeitdiskretisierung $\Omega_\tau$ des Anfangswertproblems
{eq}`eq:awp_einschritt` für das Zeitinterval $[0,T] \subset \R^+$ mit
$n+1$ äquidistanten Zeitschritten $t_k = k \cdot \tau$ und
$\tau = T / n$ gegeben.

Dann lassen sich **Einschrittverfahren** zur numerischen Lösung des
Anfangswertproblems in folgender allgemeiner Form angeben:
```{math}
u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \tau \cdot f_\tau(t_{k}, u_\tau(t_{k})),
```
wobei $f_\tau \colon \Omega_\tau \times \R^n \rightarrow \R^n$ eine
numerische Approximation der Funktion $F$ auf der rechten Seite der
gewöhnlichen Differentialgleichung in {eq}`eq:awp_einschritt` darstellt
und allgemein **Verfahrensfunktion** genannt wird.

Hängt die Verfahrensfunktion $f_\tau$ nur vom bekannten Wert
$u_\tau(t_k)$ ab, so nennen wir das Einschrittverfahren **explizit**, da
es eine explizite Vorschrift zur Berechnung des nächsten Zeitschritts
basierend auf den Werten des aktuellen Zeitschritts liefert. Hängt
andererseits die Verfahrensfunktion $f_\tau$ auch vom unbekannten Wert
$u_\tau(t_{k+1})$ ab, so heißt das Verfahren **implizit**, da es nur
eine implizite Bedingung (im Allgemeinen eine nichtlineare Gleichung)
für den neuen Zeitschritt liefert.

````

Ob ein Einschrittverfahren explizit oder implizit ist, hängt häufig von
der Wahl der Approximation der Zeitableitung mittels finiter Differenzen
ab, die in folgender Definition eingeführt werden.

````{prf:definition} Finite Differenzen
:label: def:finite_differenzen
Die Grundidee der sogenannten **finiten Differenzen** ist die
Approximation der Ableitung $u'(t)$ durch Differenzenquotienten mit
Werten auf einem Diskretisierungsgitter $\Omega_\tau$ für eine feste
Schrittweite $\tau \coloneqq \frac{T}{N}, N \in \mathbb{N}$.
Differenzenbildung auf einem Gitter. Hierbei haben wir verschiedene
Möglichkeiten die Zeitableitung zu approximieren und wir unterscheiden
folgende Fälle:
```{math}
\begin{aligned}
u'(t_k) \ &\approx \ \frac{u(t_{k+1}) - u(t_k)}{\tau} \ = \ \frac{u(t_k + \tau) - u(t_k)}{\tau} \ &=: \ D^+ u(t_k),\\
u'(t_k) \ &\approx \ \frac{u(t_k) - u(t_{k-1})}{\tau} \ = \ \frac{u(t_k) - u(t_k - \tau)}{\tau} \ &=: \ D^- u(t_k),\\
u'(t_k) \ &\approx \ \frac{u(t_{k+1}) - u(t_{k-1})}{2\tau} \ = \ \frac{u(t_k + \tau) - u(t_k - \tau)}{2\tau} \hspace{-1cm} &=: \ D^c u(t_k).
\end{aligned}
```
Hierbei werden die obigen Differenzenquotienten $D^+, D^-, D^c$
**Vorwärts-, Rückwärts-** und **zentrale Differenz** genannt.

````

Mittels Taylorformel sieht man ein, dass alle drei finiten Differenzen
in {prf:ref}`def:finite_differenzen` eine gute numerische Approximation
für die erste Ableitung der Funktion $u$ liefern, wenn die Schrittweite
$\tau > 0$ der Diskretisierung hinreichend klein gewählt wird. Da der
Fehler dieser Approximation für $\tau \rightarrow 0$ bei allen drei
Varianten gegen Null konvergiert spricht man auch von *konsistenten
Verfahren* (mehr dazu später). Bei einer quantitativen Analyse des
Fehlerglieds in der Taylorentwicklung lässt sich zeigen, dass die
Vorwärts- und Rückwärtsdifferenz $D^+$ und $D^-$ jeweils die
Konsistenzordnung $1$ besitzen, d.h., der Fehler verschwindet linear
bzgl. der Schrittweite $\tau$, während die zentrale Differenz $D^c$
sogar Konsistenzordnung $2$ besitzt, d.h., der Fehler verschwindet
quadratisch bzgl. der Schrittweite $\tau$.

Wir betrachten nun einige bekannte Einschrittverfahren, die häufig zur
numerischen Lösung von Anfangswertproblemen genutzt werden.

````{prf:example} Einschrittverfahren
Im Folgenden wollen wir drei unterschiedliche Einschrittverfahren
diskutieren, die häufig zur numerischen Lösung von Anfangswertproblemen
genutzt werden.

1.  Das einfachste Beispiel eines Einschrittverfahrens ist das
    **Vorwärts-Euler-Verfahren**
    ```{math}
    u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \tau \cdot F(t_k,u_\tau(t_k)),
    ```
    welches auf einer Approximation der Ableitung $u'(t_k)$ mittels
    Vorwärtsdifferenz aus {prf:ref}`def:finite_differenzen` basiert. Das
    Vorwärts-Euler Verfahren hat die Verfahrensfunktion
    ```{math}
    f_\tau(t_k, u_\tau(t_k)) \ \coloneqq \ F(t_k,u_\tau(t_k)),
    ```
    und ist entsprechend nach {prf:ref}`def:einschrittverfahren` ein
    explizites Verfahren. Daher wird es häufig auch explizites
    Euler-Verfahren genannt.

    In der Integralform der Differentialgleichung bedeutet das, dass wir
    die folgende Approximation benutzen:
    ```{math}
    \begin{aligned}
    u(t_{k+1}) - u(t_k) \ = \ \int_{t_k}^{t_{k+1}} F(t,u(t))\, \mathrm{d}t \ \approx \ \tau \cdot F(t_k, u_\tau(t_k)),
    \end{aligned}
    ```
    d.h., dass wir das Integral durch die Intervallänge mal dem Wert am
    linken Intervallrand annähern. Dies entspricht einer numerischen
    Quadraturformel mit nur einer Stützstelle (siehe
    {cite:p}`numerik1`).

2.  Ein wenig komplizierter ist das **Rückwärts-Euler-Verfahren**
    ```{math}
    u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \tau \cdot F(t_{k+1},u_\tau(t_{k+1})),
    ```
    welches auf einer Approximation der Ableitung $u'(t_k)$ mittels
    Rückwärtsdifferenz aus {prf:ref}`def:finite_differenzen` basiert. In
    diesem Fall müssen wir zunächst eine (möglicherweise nichtlineare)
    Gleichung für den neuen Zeitschritt lösen.

    Das Rückwärts-Euler Verfahren hat die Verfahrensfunktion
    ```{math}
    f_\tau(t_k, u_\tau(t_k)) \ \coloneqq \ F(t_{k+1},u_\tau(t_{k+1})),
    ```
    und ist entsprechend nach {prf:ref}`def:einschrittverfahren` ein
    implizites Verfahren. Daher wird es häufig auch implizites
    Euler-Verfahren genannt.

    In der Integralform der Differentialgleichung bedeutet das, dass wir
    die folgende Approximation benutzen:
    ```{math}
    \begin{aligned}
    u(t_{k+1}) - u(t_k) \ = \ \int_{t_k}^{t_{k+1}} F(t,u(t)) \,\mathrm{d}t \ \approx \ \tau \cdot F(t_{k+1},u_\tau(t_{k+1})),
    \end{aligned}
    ```
    d.h., dass wir das Integral durch die Intervallänge mal dem Wert am
    rechten Intervallrand annähern. Auch dies entspricht einer
    numerischen Quadraturformel mit nur einer Stützstelle wie beim
    expliziten Euler-Verfahren.

3.  Ein Kompromiss aus den beiden obigen Verfahren ist das sogenannte
    **Crank-Nicolson Verfahren** der Form
    ```{math}
    u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \frac{\tau}{2} \cdot \left(F(t_{k },u_\tau(t_{k })) + F(t_{k+1},u_\tau(t_{k+1}))\right).
    ```
    Hierbei ist die Verfahrensfunktion der *Mittelwert* der
    Verfahrensfunktionen des Vorwärts- und Rückwärts-Euler-Verfahrens
    mit
    ```{math}
    \begin{aligned}
    f_\tau(t_k, u_\tau(t_k)) \ \coloneqq \ \frac{1}{2} \left( F(t_{k },u_\tau(t_{k })) + F(t_{k+1},u_\tau(t_{k+1}))\right)
    \end{aligned}
    ```
    und ist entsprechend nach {prf:ref}`def:einschrittverfahren`
    ebenfalls ein implizites Verfahren.

    In der Integralform der Differentialgleichung bedeutet das, dass wir
    die folgende Approximation benutzen:
    ```{math}
    \begin{aligned}
    u(t_{k+1}) - u(t_k) \ &= \ \int_{t_k}^{t_{k+1}} F(t,u(t)) \,\mathrm{d}t \\
    &\approx \ \frac{\tau}{2} \left( F(t_{k },u_\tau(t_{k })) + F(t_{k+1},u_\tau(t_{k+1}))\right),
    \end{aligned}
    ```
    Dies entspricht einer Approximation des Integrals mit Hilfe der
    Trapezregel aus {cite:p}`numerik1`.

````

Wir sehen, dass ein explizites Verfahren sofort wohldefiniert ist, falls
die Verfahrensfunktion $f_\tau$ eine stetige Funktion auf
$\mathbb{R}_+ \times \mathbb{R}^n$ ist, während bei impliziten Verfahren
noch eine Fixpunktgleichung gelöst werden muss. Die Existenz und
Eindeutigkeit dieser Gleichung können wir mit dem **Banachschen
Fixpunktsatz** garantieren, wenn wiederum $\tau$ klein genug ist.

````{prf:theorem} Existenz und Eindeutigkeit von Lösungen
Sei $f_\tau:\R\times\R^n\times\R^n\to\R^n$ stetig und Lipschitz-stetig
bezüglich dem letzten Argument mit Lipschitzkonstante $L > 0$. Dann
existiert für $\tau< \frac{1}{L}$ genau eine Lösung $u_\tau(t_{k+1})$
der Fixpunktgleichung
```{math}
u = u_\tau(t_k) + \tau\ f_\tau(t_k,u_\tau(t_k),u).
```

````

Numerisch müssen wir zur Durchführung eines impliziten Verfahrens immer
noch ein System in $\mathbb{R}^n$ lösen. Ist dieses linear, so können
wir die üblichen Verfahren für lineare Gleichungssysteme anwenden.
Andernfalls bietet sich die Verwendung eines iterativen Verfahrens wie
einer Fixpunktiteration oder des Newton-Verfahrens an (beachte, dass
unter der obigen Bedingung $\mathds{1}-\tau f_\tau'$ invertierbar ist
für $f_\tau \in C^1$). Mit dem Wert $u_\tau(t_k)$ oder einer einfachen
Vorhersage in der Zeit (etwa mit dem expliziten Euler-Verfahren) haben
wir dafür auch einen sehr guten Startwert.

Nachdem wir die Wohldefiniertheit und numerische Umsetzung von
Einschrittverfahren geklärt haben, widmen wir uns nun der Analyse der
Verfahren. Wir wollen dabei den Fehler

```{math}
:label: eq:einschrittfehler 
E_\tau  = \max_{k\in\N} \Vert u_\tau(t_k) - u(t_k) \Vert
```
abschätzen, wobei $u$ die exakte Lösung des Anfangswertproblems ist.
Wollen wir den Fehler an anderen Stellen $t$ abschätzen, so können wir
ein Interpolationsverfahren und die entsprechenden Abschätzungen
anwenden.

Unsere Strategie dabei ist die Folgende: Zunächst schreiben wir eine
Gleichung für den Fehler $e_\tau = u_\tau - u$. Es gilt
```{math}
\begin{aligned}
e_\tau(t_{k+1}) =& e_\tau(t_k) +\\ 
&\tau (f_\tau(t_k,u_\tau(t_k),u_\tau(t_{k+1})) - f_\tau(t_k,u(t_k),u(t_{k+1}))) 
+ \\ & \tau 
\left[ f_\tau(t_k,u(t_k),u(t_{k+1})) - \frac{1}\tau \int_{t_k}^{t_{k+1}} F(t,u(t))~dt
\right]. 
\end{aligned}
```
Nun benötigen wir zwei zentrale Eigenschaften von
Diskretisierungsmethoden:

-   **Konsistenz:** Der Fehler
    ```{math}
    \begin{aligned}
    f_\tau(t_k,u(t_k),u(t_{k+1})) - \frac{1}{\tau}\int_{t_k}^{t_{k+1}} F(t,u(t))~dt,
    \end{aligned}
    ```
    d.h. das Residuum der Lösung des Anfangswertproblems eingesetzt in
    das numerische Verfahren konvergiert gegen Null für
    $\tau \rightarrow 0$.

-   **Stabilität:** Bei der Umsetzung des numerischen Verfahrens mit
    gegebener rechter Seite wird diese nicht beliebig verstärkt,
    insbesondere existiert eine Abschätzung unabhängig von $\tau.$

Zusammen ergeben Konsistenz und Stabilität Konvergenz des Verfahrens,
d.h. $E_\tau \rightarrow 0$. Dies halten wir in einer Definition fest:

````{prf:definition} Konvergenz von Einschrittverfahren
Sei $E_\tau$ definiert durch {eq}`eq:einschrittfehler`, dann heißt das
Verfahren

1.  konvergent, wenn $E_\tau \rightarrow 0$ für $\tau \rightarrow 0$,

2.  konvergent von der Ordnung $p$, wenn $E_\tau = {\cal O}(\tau^p)$ für
    $\tau \rightarrow 0$, d.h. es gibt eine Konstante $C_p$, sodass
    $E_\tau \leq C_p \tau^p$ für $\tau$ hinreichend klein.

````

(konsistenz-von-einschrittverfahren)=
## Konsistenz von Einschrittverfahren

Gemäß der obigen Motivation definieren wir den Konsistenzfehler als
```{math}
:label: eq:konsistenzfehler 
K_\tau = \max_{k\in\N} \Vert g_\tau(t_k) \Vert
```
mit
```{math}
g_\tau(t_k) = f_\tau(t_k,u(t_k),u(t_{k+1})) - \frac{1}\tau \int_{t_k}^{t_{k+1}} F(t,u(t))~dt,
```
wobei $u$ eine Lösung des Anfangswertproblems {eq}`eq:awp` ist.

````{prf:definition} Konsistenzfehler
:label: def:konsistenzfehler
Sei $K_\tau$ definiert durch {eq}`eq:konsistenzfehler`, dann heißt das
Verfahren

1.  konsistent, wenn $K_\tau \rightarrow 0$ für $\tau \rightarrow 0$,

2.  konsistent von der Ordnung $p$, wenn $K_\tau = {\cal O}(\tau^p)$ für
    $\tau \rightarrow 0$ .

````

Die Abschätzung des Konsistenzfehlers erfolgt meist durch
Taylorentwicklung, wir führen dies an zwei Beispielen durch:

````{prf:example} 
Wir betrachten das Vorwärts-Euler Verfahren unter der Annahme, dass $F$
bezüglich beider Variablen Lipschitz-stetig ist. Definieren wir
$\varphi(t) = F(t,u(t))$, dann ist $\varphi$ wegen $u \in C^1$ eine
Lipschitz-stetige Funktion und es gilt
```{math}
\begin{aligned}
\norm{g_\tau(t_k)} &= \norm{\varphi(t_k) - \frac{1}\tau \int_{t_k}^{t_{k+1}} \varphi(t)~dt}\\
&= \norm{\frac{1}\tau \int_{t_k}^{t_{k+1}} (\varphi(t_k)-\varphi(t))~dt}\\
&\leq\frac{1}\tau \int_{t_k}^{t_{k+1}} \norm{\varphi(t_k)-\varphi(t)}~dt   \\
&\leq \frac{1}\tau \int_{t_k}^{t_{k+1}} L_\varphi ( t-t_k)~dt = \frac{L_\varphi}2 \tau.
\end{aligned}
```
Damit das Verfahren die Konsistenzordnung $p=1$, wir sehen im Beispiel
$F(t,u) = t$ auch sofort, dass man im allgemeinen nicht Ordnung zwei
erreichen kann.

````

````{prf:example} 
Wir betrachten das Crank--Nicholson Verfahren unter der Annahme, dass
$F$ bezüglich beider Variablen zweimal stetig differenzierbar ist.
Definieren wir $\varphi(t) = F(t,u(t))$, dann ist $\varphi$ ebenfalls
zweimal stetig differenzierbar, da
```{math}
u''(t) = (F(t,u(t)))' = \partial_t F(t,u(t)) + \partial_u F(t,u(t)) u'(t)
```
stetig ist. Damit gilt
```{math}
\begin{aligned}
\norm{g_\tau(t_k)} &= \norm{\frac{1}2(\varphi(t_k)+\varphi(t_{k+1})) - \frac{1}\tau \int_{t_k}^{t_{k+1}} \varphi(t)~dt}\\
&=    \frac{1}{2\tau}  \Vert \int_{t_k}^{t_{k+1}} (\varphi(t_k)+\varphi(t_{k+1})-2\varphi(t))~dt  \Vert \\
&\leq    \frac{1}{2\tau}   \Vert\int_{t_k}^{t_{k+1}} \varphi'(t_k)(t_k+t_{k+1}-2t) + r_k ~dt \Vert   ,
\end{aligned}
```
mit dem Restglied $r_k={\cal O}(\tau^2)$. Da
$\int_{t_k}^{t_{k+1}} \varphi'(t_k)(t_k+t_{k+1}-2t) ~dt = 0$, folgt
$\Vert g_\tau(t_k) \Vert = {\cal O}(\tau^p).$ Damit das Verfahren die
Konsistenzordnung $p=2$, wir sehen im Beispiel $F(t,u) = t^2$ auch
sofort, dass man im allgemeinen nicht Ordnung zwei erreichen kann.

````

Um eine Konvergenzordnung $p$ zu erhalten, benötigen wir, dass $F$, aber
auch $u$ $p$-mal stetig differenzierbar ist. Aus den Eigenschaften von
$F$ folgt Letzteres aber sofort: wir haben gesehen, dass für $F$ stetig
auch $u$ stetig differenzierbar folgt. Mit dem Argument aus dem letzten
Beispiel sehen wir, dass für $F$ stetig differenzierbar auch $u$ zweimal
stetig differenzierbar ist. Induktiv können wir durch weiteres
differenzieren zeigen, dass $u$ $p$-mal stetig differenzierbar ist, wenn
$F$ $p-1$-mal stetig differenzierbar ist.

(stabilität-und-konvergenz)=
## Stabilität und Konvergenz

Wir widmen uns nun der Frage der Stabilität von Einschrittverfahren.
Hierbei verwenden wir eine diskrete Version des Lemmas von Gronwall.

````{prf:lemma} Diskretes Lemma von Gronwall
:label: lem:gronwall_diskret
Es sei $\beta_j \geq 0, j\in\N_0$ eine Folge nicht-negativer Zahlen und
für die Folge $u_j\in\R, j\in\N_0$ gelte
```{math}
\begin{aligned}
u_0&\leq \alpha\in\R^+_0\\
u_k&\leq \alpha + \sum_{j=0}^{k-1} \beta_j u_j
\end{aligned}
```
für $k\in\N$, dann gilt die Abschätzung
```{math}
\begin{aligned}
u_k \leq \alpha \exp\left(\sum_{j=0}^{k-1} \beta_j\right).
\end{aligned}
```

````

````{prf:proof} 
Übung. ◻

````

Wir zeigen zunächst uniforme Schranken an $u_\tau$.

````{prf:lemma} 
Sei $F_\tau$ stetig und Lipschitz-stetig bezüglich des dritten Arguments
(d.h. $u_\tau(t_{k+1}$) mit Modul unabhängig von $\tau$. Dann existiert
eine Konstante $M(T)$ unabhängig von $\tau$, sodass
```{math}
\max_{t_k \leq T} \Vert u_\tau(t_k) \Vert \leq M
```
gilt für alle $\tau$ hinreichend klein.

````

````{prf:proof} 
Aus der Definition des Verfahrens folgt
```{math}
u_\tau(t_{k+1})-u_0 = u_\tau(t_k) - u_0 + \tau (f_\tau(t_k,u_\tau(t_k),u_\tau(t_{k+1}))-f_\tau(t_k,u_0,u_0)) + \tau f_\tau(t_k,u_0,u_0)
```
und mit der Dreiecksungleichung folgt für
$v_k = \Vert u_\tau(t_k) - u_0 \Vert$

```{math}
\begin{aligned}
 v_{k+1} &\leq v_k + \tau \Vert f_\tau(t_k,u_\tau(t_k),u_\tau(t_{k+1}))-f_\tau(t_k,u_0,u_0) \Vert + \tau \Vert f_\tau(t_k,u_0,u_0) \Vert \\
&\leq  v_k + \tau L (v_k + v_{k+1}) + \tau C. 
\end{aligned}
```
Hier haben wir benutzt, dass $f_\tau$ stetig ist, damit folgt
$f_\tau(t,u_0,u_0)$ ist auf dem kompakten Intervall $[0,T]$ durch eine
Konstante $C$ beschränkt. Dazu bezeichnet $L$ den Lipschitz-Modul von
$f_\tau$ bezüglich zweitem und drittem Argument. Sei nun
$\tau \leq \frac{1}{2L}$, d.h. $1- \tau L \geq \frac{1}2$, dann folgt
```{math}
\begin{aligned}
v_{k+1} \leq 2(1+\tau L) v_k + 2 \tau C.
\end{aligned}
```
Das diskrete Lemma von Gronwall impliziert dann die Beschränktheit von
$v_k$. ◻

````

Mit einem ähnlichen Beweis können wir auch die Stabilität zeigen.

````{prf:theorem} Stabilität von Einschrittverfahren
:label: thm:einschrittverfahren_fehlerabschaetzung
Sei $u_\tau \colon \Omega_\tau \rightarrow \R^n$ die Lösung eines
Einschrittverfahrens mit Lipschitz-stetiger Verfahrensfunktion $f_\tau$
für eine Zeitdiskretisierung
$\Omega_\tau \coloneqq \lbrace t_k \in [0, T] \colon t_k \coloneqq \frac{k\cdot T}{N}, k=0,\ldots,N\rbrace$
mit Schrittweite $\tau \coloneqq \frac{T}{N} > 0$. Sei außerdem
$u \colon [0, T] \rightarrow \R$ die Lösung des Anfangswertproblems
{eq}`eq:awp` mit gleichem Anfangswert $u_0 \in \R^n$. Der lokale
Konsistenzfehler $g_\tau(t_k)$ und der globale Konsistenzfehler $K_\tau$
seien gegeben wie in {prf:ref}`def:konsistenzfehler`.

Dann existiert eine Konstante $C > 0$, so dass für $\tau > 0$
hinreichend klein die folgende Abschätzung für den Fehler des
Einschrittverfahrens gilt:
```{math}
E_\tau \ = \ \max_{t_k} \Vert u_\tau(t_k) - u(t_k) \Vert \ \leq \ C \cdot  \max_{t_k} \Vert g_\tau(t_k) \Vert \ = \ C \cdot K_\tau.
```

````

````{prf:proof} 
TODO
<!---
Wir definieren uns die Hilfsfunktion
$v_k \coloneqq \Vert u_\tau(t_{k}) -u(t_k)) \Vert$. Durch Anwendung der
Dreiecksungleichung und dem Ausnutzen der Lipschitz-Stetigkeit von
$f_\tau$ erhalten wir dann
```{math}
\begin{aligned}
v_{k+1} \ &= \ \Vert u_\tau(t_{k+1}) - u(t_{k+1})\Vert\\
&= \ \left\Vert u_\tau(t_k) + \tau \cdot f_\tau(t_k, t_{k+1}, u_\tau(t_k), u_\tau(t_{k+1})) - u(t_k) -\int_{t_k}^{t_{k+1}} F(t,u(t)) \, \mathrm{d}t \right\Vert\\
&\leq \ v_k + \tau \cdot \norm{f_\tau(t_k, t_{k+1}, u_\tau(t_k), u_\tau(t_{k+1})) - f_\tau(t_k, t_{k+1}, u(t_k), u(t_{k+1}))}\\
& \hspace{1.01cm} + \tau  \cdot \left\Vert f_\tau(t_k, t_{k+1}, u(t_k), u(t_{k+1})) - \frac{1}{\tau} \int_{t_k}^{t_{k+1}} F(t,u(t)) \, \mathrm{d}t \right\Vert\\
&\leq \ v_k + \tau \cdot L \cdot (v_k + v_{k+1}) + \norm{g_\tau(t_k)}.
\end{aligned}
```
Somit erhalten wir durch Umstellen also insgesamt
```{math}
\begin{aligned}
v_{k+1}  \ &\leq \ v_k + L \cdot \tau \cdot (v_k + v_{k+1}) + \Vert g_\tau(t_k) \Vert\\
\Rightarrow \quad v_{k+1} \ &\leq \ \frac{1+\tau L}{1-\tau L} \cdot v_k + \max_{t_k}\norm{g_\tau(t_k)}.
\end{aligned}
```
Für hinreichend kleine Schrittweiten $\tau < \frac{1}{2L}$ erhalten wir
die gewünschte Schranke wieder direkt aus dem diskreten
{prf:ref}`lem:gronwall_diskret` von Gronwall. ◻
-->
````

Eine direkte Folgerung von
{prf:ref}`thm:einschrittverfahren_fehlerabschaetzung` ist die Äquivalenz
von Konsistenz und Konvergenz für Einschrittverfahren, wie wir im
folgenden Korollar feststellen.

````{prf:corollary} Äquivalenz von Konsistenz und Konvergenz
Für ein Einschrittverfahren mit Lipschitz-stetiger Verfahrensfunktion
gilt: ist das Verfahren konsistent von der Ordnung $p$, so ist es auch
konvergent von der Ordnung $p$.

````

Mit Hilfe dieses Korollars und der Abschätzung des Konsistenzfehlers
stellen wir fest, dass das Vorwärts- und Rückwärts-Euler Verfahren
konvergent von der Ordnung eins sind. Das Crank--Nicholson Verfahren
hingegen ist sogar konvergent von der Ordnung zwei.

(rungekutta-verfahren)=
## Runge--Kutta Verfahren

Bisher haben wir Verfahren kennengelernt, die auf einzelne
Funktionsauswertungen von $F$ an den Zeitschritten $t_k$ und $t_{k+1}$
zurückgreifen. Damit haben wir bereits die Konsistenzordnung Eins
erreicht für die beiden Euler-Verfahren und als maximale
Konsistenzordnung Zwei im Falle des Crank--Nicholson Verfahrens. Eine
höhere Konsistenzordnung ist mit so einem Ansatz nicht möglich. Wir
widmen uns im Folgenden also der Frage wie wir Einschrittverfahren
höherer Ordnung konstruieren können. Aus der Analyse der
Konsistenzordnung lässt sich ableiten, dass wir Verfahrensfunktionen der
Art konstruieren müssen, so dass die Taylorentwicklung ein Restglied
höherer Ordnung liefert.

Eine erste Idee zur Steigerung der Konsistenzordnung ist es Ableitungen
von $F$ in den Zeitschritten $t_k$ und $t_{k+1}$ zu berücksichtigen,
womit man offensichtlich die Taylor-Entwicklung besser approximieren und
somit eine höhere Ordnung erreichen kann. Die Berechnung von Ableitungen
von $F$ ist jedoch potentiell numerisch aufwändig und instabil, deswegen
geht alternativ einen anderen Weg und verwendet geschachtelte
Funktionsauswertungen.

Die grundlegende Idee ist es das Integral der Anfangswertaufgabe in
{eq}`eq:anfangswertaufgabe_integralform` durch eine geeignete
Quadraturformel (siehe {cite:p}`numerik1`) der folgenden Form zu
approximieren:
```{math}
\frac{1}{\tau} \int_{t_k}^{t_{k+1}} F(t, u(t)) \,\mathrm{d}t \ \approx \ f_\tau(t_k, u_\tau(t_k)) \ \coloneqq \ \sum_{i=1}^s b_i f_i^k.
```
Hierbei stellen die $b_i \in \R$ Gewichte der Quadraturformel dar und
die Zwischenwerte $f_i^k \in \R^n$ sollen eine Approximation der
Funktionswerte von $F$ liefern für $i=1,\ldots,s$ mit
```{math}
f_i^k \ \approx \ F(t_k+c_i\tau,u(t_k + c_i\tau)).
```
Hierbei sind die $0 \leq c_i \leq 1, i=1,\ldots,s$ aufsteigende
Parameter, die Stützstellen für die jeweiligen Funktionsauswertung
definieren. Da wir die Werte der unbekannten Funktion
$u_\tau(t_k+c_i \tau)$ an den Stützstellen nicht kennen, benötigen wir
ebenfalls Approximationen mit Hilfe von Quadraturformeln der folgenden
Form:
```{math}
u(t_k + c_i\tau) \ = \ u(t_k) + \int_{t_k}^{t_k + c_i\tau} F(t, u(t)) \, \mathrm{d}t \ \approx \  u_\tau(t_k) + \tau \sum_{j=1}^s a_{ij} f_j^k.
```

Diese Idee führt zur Definition von sogenannten Runge-Kutta Verfahren.

````{prf:definition} Runge--Kutta Verfahren der Stufe s
:label: def:runge-kutta
Bei einem **Runge--Kutta Verfahren der Stufe $\mathbf{s}$**,
$s \in \N^+$ berechnet man
```{math}
u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \tau f_\tau(t_k, u_\tau(t_k)) \ = \ u_\tau(t_k) + \tau \sum_{i=1}^s b_i f_i^k.
```

Die Zwischenwerte $f_i^k$ lassen sich als Funktionsauswertungen von $F$
wie folgt berechnen
```{math}
:label: eq:runge_kutta_schachtelung
f_i^k \ \coloneqq \  F\Bigl(t_k + c_i \tau, u_\tau(t_k) + \tau \sum_{j=1}^s a_{ij} f_j^k\Bigr), \qquad i=1,\ldots,s.
```
Um in der Zeit vorwärts zu gehen wählt man die Koeffizienten
$0 \leq c_i \leq 1, i=1,\ldots,s$ als aufsteigende Folge und die Matrix
$(a_{ij})_{i,j=1,\ldots,s}$ als linke, untere Dreiecksmatrix.

Man nennt ein Runge--Kutta Verfahren **explizit** falls für alle
Einträge der Matrix $A\in \R^{s\times s}$ gilt $a_{ij}=0, i=1,\ldots,s$
wenn $j \geq i$.

````

Die zu bestimmenden Koeffizienten $b_i,c_i \in \R$ und $a_{i,j} \in \R$
für $i,j=1,\ldots,s$ in {prf:ref}`def:runge-kutta` erhält man durch
einen Vergleich mit der Taylorentwicklung des zu approximierenden
Integrals $\frac{1}\tau \int_{t_k}^{t_{k+1}} F(t,u(t)) \, \mathrm{d}t$.
Das folgende Beispiel soll die Idee zur Herleitung eines Runge--Kutta
Verfahren für den einfachen Fall der Stufe $s=1$ verdeutlichen.

````{prf:example} Runge--Kutta Verfahren der Stufe 1
:label: ex:runge_kutta_s=1
Für ein Runge--Kutta Verfahren der Stufe $s=1$ ist die
Verfahrensfunktion gegeben durch
```{math}
f_\tau(t_k, u_\tau(t_k)) \ \coloneqq \ b_1 f_1^k,
```
und für den Zwischenwert $f_1^k \in \R^n$ gilt entsprechend:
```{math}
f_1^k \ = \ F(t_k + c_1 \tau, u_\tau(t_k) +  \tau a_{11} f_1^k).
```

Wollen wir nun ein **explizites Verfahren** herleiten, so muss nach
{prf:ref}`def:runge-kutta` schon $a_{11}=0$ gelten. Die
Verfahrensfunktion ist also von der Form
```{math}
f_\tau(t_k, u_\tau(t_k)) \ = \ b_1 f_1^k \ = \ b_1 F(t_k+c_1 \tau,  u_\tau(t_k)).
```
Die einzige sinnvolle Wahl für den Koeffizienten $c_1$ ist Null, da wir
sonst die Funktion $F$ zu einer anderen Zeit auswerten als die
unbekannte Funktion $u$. Um ein konsistentes Verfahren zu erhalten sieht
man außerdem ein, dass $b_1=1$ gelten muss. Also erhalten wir
schlussendlich das Vorwärts-Euler Verfahren.

 \
Im **impliziten Fall** können wir interessanterweise eine höhere Ordnung
erreichen. Wir berechnen hierzu die Taylor-Entwicklung des Zwischenwerts
$f_1^k$ wie folgt
```{math}
:label: eq:runge_kutta_taylor_f1k
\begin{split}
f_1^k \ &= \ F(t_k + c_1 \tau, u_\tau(t_k) +  \tau a_{11} f_1^k) \\
&= \ F(t_k,u(t_k)) + c_1 \tau \partial_t F(t_k,u(t_k)) + \tau a_{11} D_u F(t_k,u(t_k)) f_1^k + {\cal O}(\tau^2).
\end{split}
```
Hierbei bezeichnet $D_u$ die Jacobimatrix bezüglich des zweiten
Arguments $u$ von $F$. Setzen wir auf der rechten Seite nochmal den
ersten Approximationsterm für $f_1^k$ ein, so folgt
```{math}
f_1^k \ = \ F(t_k,u(t_k)) + c_1 \tau \partial_t F(t_k,u(t_k)) + \tau a_{11} D_u F(t_k,u(t_k)) F(t_k,u(t_k))+ {\cal O}(\tau^2).
```

Andererseits gilt mit einer Taylorentwicklung des unbekannten
Funktionswerts $u(t_{k+1})$ im Punkt $t_k$:
```{math}
:label: eq:runge_kutta_taylor_integral
\begin{split}
\frac{1}\tau \int_{t_k}^{t_{k+1}} F(t,u(t))\,\mathrm{d}t \ &= \ \frac{1}{\tau} \left( u(t_{k+1}) - u(t_k) \right) \\
&= \ \frac{1}{\tau} \left( u(t_k) + \tau u'(t_k) + \frac{\tau^2}{2}u''(t_k) + \mathcal{O}(\tau^3) - u(t_k) \right)\\
&= \ u'(t_k) + \frac{\tau}{2} u''(t_k) + \mathcal{O}(\tau^2)\\
&= \ F(t_k,u(t_k)) + \frac{\tau}2 \partial_t F(t_k,u(t_k))\\
& \hspace{2.58cm} +\frac{\tau}2 D_u F(t_k,u(t_k)) F(t_k,u(t_k)) + {\cal O}(\tau^2).
\end{split}
```
Damit gilt
```{math}
f_\tau(t_k, u_\tau(t_k)) \ = \ b_1 f_1^k \ \approx \ \frac{1}{\tau}\int_{t_k}^{t_{k+1}} F(t,u(t))\,\mathrm{d}t,
```
vergleichen wir die beiden Formeln {eq}`eq:runge_kutta_taylor_f1k` und
{eq}`eq:runge_kutta_taylor_integral` und sehen ein, dass wir eine
Konsistenzordnung von Zwei erreichen können, wenn die Koeffizienten als
$b_1 =1$, $c_1=\frac{1}2$ und $a_{11}=\frac{1}2$ gewählt werden.

Die Verfahrensfunktion ist also gegeben durch die Lösung von
```{math}
f_\tau(t_k, u_\tau(t_k)) \ = \ f_1^k \ = \ F(t_k +\frac{\tau}2,u_\tau(t_k) + \frac{\tau}2 f_1^k).
```
Wir können dieses Runge-Kutta Verfahren als *Mittelpunktsregel* im
Intervall $[t_k,t_{k+1}]$ interpretieren, wobei der unbekannte Wert von
$u_\tau$ am Mittelpunkt $t_k +\frac{\tau}2$ des Intervalls durch das
*Rückwärts-Euler-Verfahren* bestimmt wird.

````

Im Folgenden betrachten wir noch $2$-stufige Runge--Kutta Verfahren.

````{prf:example} Runge--Kutta Verfahren der Stufe 2
:label: ex:runge_kutta_s=2
Für ein Runge--Kutta Verfahren der Stufe $s=2$ ist die
Verfahrensfunktion gegeben durch
```{math}
f_\tau(t_k, u_\tau(t_k)) \ \coloneqq \ b_1 f_1^k +b_2 f_2^k,
```
und für die Zwischenwerte $f_1^k, f_2^k \in \R^n$ gilt im **expliziten
Fall** mit $a_{11} = a_{22} = 0$:
```{math}
f_1^k \ = \ F(t_k + c_1 \tau, u_\tau(t_k)), \quad f_2^k \ = \ F(t_k + c_2 \tau, u_\tau(t_k) +  \tau a_{21} f_1^k).
```
Analog zur Argumentation im expliziten Fall von
{prf:ref}`ex:runge_kutta_s=1` sehen wir wieder, dass nur $c_1=0$ eine
sinnvolle Wahl ist. Eine Taylor-Entwicklung des Zwischenwertes $f_2^k$
liefert dann
```{math}
\begin{aligned}
b_1 f_1^k +b_2 f_2^k \ &= \ (b_1 + b_2) F(t_k , u(t_k)) + b_2 c_2 \tau \partial_t F(t_k , u (t_k))\\ 
& \hspace{1cm} + b_2 a_{21} \tau D_u F(t_k , u (t_k)) F(t_k , u (t_k)) + {\cal O}(\tau^2).
\end{aligned}
```
Vergleichen wir dies wieder mit der Taylor-Entwicklung des Integrals in
{eq}`eq:runge_kutta_taylor_integral`, so erhalten wir folgendes
Gleichungssystem
```{math}
b_1+b_2 \, = \, 1, \qquad b_2 c_2 \, = \, \frac{1}2, \qquad b_2 a_{21} \, = \, \frac{1}2.
```
Eine einfache Lösung ist beispielsweise $b_1=0$, $b_2=1$,
$c_2=a_{21} = \frac{1}2$. Diese Wahl der Koeffizienten liefert uns also
ein Verfahren der Konsistenzordnung zwei mit der Verfahrensfunktion
```{math}
f_\tau(t_k, u_\tau(t_k)) \ = \ F(t_k +\frac{\tau}2,u_\tau(t_k) + \frac{\tau}2 F(t_k,u(t_k))).
```
Wir können dieses zweistufige Runge-Kutta Verfahren als eine
*Mittelpunktsregel* im Intervall $[t_k,t_{k+1}]$ interpretieren, wobei
der unbekannte Wert von $u_\tau$ am Mittelpunkt $t_k +\frac{\tau}2$
durch das *Vorwärts-Euler-Verfahren* bestimmt wird.

````

Bei der Konstruktion eines Runge-Kutta Verfahrens in
{prf:ref}`ex:runge_kutta_s=1` und {prf:ref}`ex:runge_kutta_s=2` haben
wir die Koeffizienten $b_i, c_i \in \R$ und $a_{i,j} \in \R$ für
$i,j=1,\ldots,s$ basierend auf vernünftigen Überlegungen zu wählen. Wie
wir speziell im Fall des zweistufigen Runge-Kutta Verfahrens gesehen
haben führt dies im Allgemeinen jedoch nicht zu eindeutigen Lösungen.
Daher ist es sinnvoll Kriterien an die implizierten Runge-Kutta
Verfahren zu stellen, so dass wir eindeutige Koeffizienten ableiten
können.

Am einfachsten ist ein Kriterium für die Konsistenz des Runge-Kutta
Verfahrens für die unbekannten Koeffizienten abzuleiten, wie das
folgende Lemma festhält.

````{prf:lemma} Konsistenzbedingung von Runge-Kutta-Verfahren
:label: lem:runge_kutta_konsistenzbedingung
Ein Runge-Kutta Verfahren der Stufe $s \in \N^+$ ist genau dann
konsistent, wenn die folgende Bedingung erfüllt ist:
```{math}
\sum_{i=1}^s b_i \ = \ 1.
```

````

````{prf:proof} 
Um zu zeigen, dass ein Runge-Kutta Verfahren konsistent ist müssen wir
zeigen, dass der Konsistenzfehler
```{math}
K_\tau \ = \ \max_{k=0,\ldots,N} \left\Vert f_\tau(t_k, u(t_k)) - \frac{1}{\tau} \int_{t_k}^{t_{k+1}} F(t, u(t)) \, \mathrm{d}t \right\Vert
```
mindestens von der Fehlerordnung $\mathcal{O}(\tau)$ ist und für
$\tau \rightarrow 0$ gegen Null geht.

Die Verfahrensfunktion des Runge-Kutta Verfahrens der Stufe $s$ gegeben
ist als
```{math}
f_\tau(t_k, u_\tau(t_k)) \ = \ \sum_{i=1}^s b_i f_i^k
```
und daher können wir die Taylorapproximation erster Ordnung jedes
Zwischenwerts $f_i^k \in \R^n$ für $i=1,\ldots,s$ analog zu
{eq}`eq:runge_kutta_taylor_f1k` betrachten und erhalten somit
```{math}
\sum_{i=1}^s b_i f_i^k \ = \  \sum_{i=1}^s (b_i F(t_k,u(t_k)) + {\cal O}(\tau)) \ = \ \sum_{i=1}^s b_i F(t_k,u(t_k)) + {\cal O}(\tau).
```
Die Taylorapproximation erster Ordnung des Integrals liefert nach
{eq}`eq:runge_kutta_taylor_integral` außerdem
```{math}
\frac{1}\tau \int_{t_k}^{t_{k+1}} F(t,u(t)) \, \mathrm{d}t = F(t_k,u(t_k)) + {\cal O}(\tau).
```

Nun gilt also für den Konsistenzfehler
```{math}
\begin{split}
K_\tau \ &= \ \max_{k=0,\ldots,N} \left\Vert \, f_\tau(t_k, u(t_k))  - \frac{1}{\tau} \int_{t_k}^{t_{k+1}} F(t, u(t)) \, \mathrm{d}t \, \right\Vert \\
&= \ \max_{k=0,\ldots,N} \left\Vert \, \sum_{i=1}^s b_i F(t_k,u(t_k)) - F(t_k,u(t_k)) + {\cal O}(\tau) \, \right\Vert.
\end{split}
```
Wir sehen, dass für den Konsistenzfehler $K_\tau \in \mathcal{O}(\tau)$
genau dann gilt, wenn $\sum_{i=1}^s b_i = 1$ ist. ◻

````

Umgekehrt ist die Konsistenzordnung eines Runge-Kutta Verfahrens nach
oben durch dessen Stufe beschränkt, wie folgendes Lemma zeigt.

````{prf:lemma} Maximale Konsistenzordnung von Runge-Kutta-Verfahren
Für die Konsistenzordnung $p \in \N^+$ eines expliziten Runge-Kutta
Verfahren der Stufe $s \in \N^+$ für ein Anfangswertproblem mit einer
rechten Seite der Differentialgleichung
$F \in C^\infty([0,T] \times \R^n; \R^n)$ gilt, dass die
Konsistenzordnung durch die Stufe des Verfahrens nach oben beschränkt
ist, d.h., es gilt $0 < p \leq s$.

````

````{prf:proof} 
Wir zeigen die Behauptung für den reellen Fall $n=1$ und wählen die
besonders einfache rechte Seite $F(t,u) = u$ mit $u_0=1 \in \R$. Dann
ist die Lösung des Anfangswertproblems
```{math}
u'(t) \ = \ F(t, u(t)) \ = \ u(t), \qquad u(0) \ = \ 1,
```
gegeben durch $u(t)=u_0\cdot e^t = e^t$.

Mit der Taylorentwicklung der Exponentialfunktion
$u(t) = e^t = \sum_{j=0}^\infty \frac{t^j}{j!}$ können wir das Integral
als Polynom in $\tau$ schreiben mit:
```{math}
\begin{aligned}
\frac{1}\tau \int_{t_k}^{t_{k+1}} F(t,u(t))\,\mathrm{d}t \ &= \ \frac{1}\tau \int_{t_k}^{t_{k+1}} u(t)\,\mathrm{d}t \\
&= \frac{1}\tau \int_{t_k}^{t_{k+1}} \underbrace{u(t_k - t_k)}_{=\, 1} \cdot \: u(t)\,\mathrm{d}t\\
&= \frac{1}\tau \int_{t_k}^{t_{k+1}} u(t_k) \cdot u(t  - t_k)\,\mathrm{d}t\\
&= \ \frac{u(t_k)}\tau \int_{t_k}^{t_{k+1}} \sum_{j=0}^p \frac{(t-t_k)^j}{j!} + {\cal O}(\tau^{p+1}) \,\mathrm{d}t \\
&= \ \frac{u(t_k)}\tau \sum_{j=0}^p \frac{(t_{k+1}-t_k)^{j+1}}{(j+1)!} + {\cal O}(\tau^{p}) \\
&= \ \sum_{j=0}^p \tau^{j} \frac{u(t_k)}{(j+1)!} + {\cal O}(\tau^{p}).
\end{aligned}
```
Andererseits sehen wir für ein beliebiges explizites Runge-Kutta
Verfahren der Stufe $s \in \N^+$ mit Konsistenzordnung $p \leq s$, dass
die Zwischenwerte gegeben sind durch:
```{math}
\begin{split}
f_1^k \ &= \ u_\tau(t_k),\\
f_2^k \ &= \ u_\tau(t_k) + \tau a_{21} u_\tau(t_k),\\
f_3^k  \ &= \ u_\tau(t_k) + \tau a_{31} u_\tau(t_k) + \tau a_{32} u_\tau(t_k) + \tau^2 a_{32} a_{21} u_\tau(t_k),\\
&\vdots
\end{split}
```
Wir sehen also, dass jeder Zwischenwert $f_i^k \in \R^n$ für
$i=1,\ldots, s$ ein Polynom in $\tau$ vom Grad kleiner gleich $i-1$ ist
und somit ist die Verfahrensfunktion
$f_\tau(t_k, u_\tau(t_k))  \ = \ \sum_{i=1}^s b_i f_i^k$ ein Polynom in
$\tau$ vom Grad kleiner gleich $s-1$.

Damit können wir keine höhere Ordnung erreichen, da der Fehler ab der
Ordnung $\mathcal{O}(\tau^p)$ im Allgemeinen nicht verschwindet. ◻

````

Man sieht ein, dass sich mittels der Forderung von Konsistenz in
{prf:ref}`lem:runge_kutta_konsistenzbedingung` und dem Bestreben nach
maximaler Konsistenzordnung $p = s$ eines $s$-stufigen Runge-Kutta
Verfahrens, das enstehende Gleichungssystem lösen lässt. Allerdings
bleibt hierbei immer noch eine Uneindeutigkeit der Koeffizienten
$c_i \in \R, i=1,\ldots,s$ übrig. Um hier eine systematische Wahl zu
treffen, fordert man im Allgemeinheit die Invarianz des Verfahrens
gegenüber sogenannter *Autonomisierung*. Hierzu führen wir zunächst den
Begriff einer autonomen Differentialgleichung ein.

````{prf:definition} Autonome gew. Differentialgleichung
Sei $u'(t) = F(t, u(t)), \ u(0) = u_0 \in \R^n$ ein Anfangswertproblem
einer gewöhnlichen Differentialgleichung erster Ordnung.

Wir nennen die Differentialgleichung **autonom**, wenn die Funktion $F$
auf der rechten Seite nicht explizit von $t$ abhängt, d.h., es gilt
$F(t,u(t)) = F(u(t))$.

````

Das allgemeine Anfangswertproblem {eq}`eq:awp` lässt sich
autonomisieren, d.h. in ein äquivalentes System autonomer
Differentialgleichungen für $\tilde u = v,u)$ mit einer Hilfsfunktion
$v$ umschreiben. Hierzu betrachten wir das Differentialgleichungssystem
```{math}
\begin{split}
u'(t) \ &= \ F(v(t),u(t)) \ = \ F(\tilde{u}(t)), \qquad u(0) \, = \, u_0 \in \R^n,\\
v'(t) \ &= \ 1, \hspace{4.7cm} v(0) \, = \, 0.
\end{split}
```
Man erkennt sofort, dass $v(t) =t$ gelten muss, was die Äquivalenz zu
{eq}`eq:awp` impliziert.

Wir nennen die obige Transformation des Anfangswertproblems
Autonomisierung und fordern, dass das Runge-Kutta Verfahren unter der
Autonomisierung invariant sein soll. Dies bedeutet, dass wir fordern,
dass die Anwendung des numerischen Verfahrens auf das neue System die
selbe Lösung $u_\tau$ liefern soll, wie die Anwendung auf das
ursprüngliche Anfangswertproblem {eq}`eq:awp`. Darüber hinaus fordern
wir, dass die Invarianz unter Autonomisierung für jede geeignete
Funktion $F$ auf der rechten Seite gelten soll.

Schreiben wir die Zwischenwerte für beide Formulierungen, so gilt für
das *ursprüngliche Anfangswertproblem*
```{math}
F(t_k + c_i \tau, u(t_k + c_i \tau)) \ \approx \ F(t_k + c_i \tau,u(t_k) + \tau \sum_{j=1}^s a_{ij} f_j^k) \ =: \ f_i^k,
```
und für die *autonome Variante* andererseits
```{math}
F(v(t_k + c_i \tau), u(t_k + c_i \tau)) \ \approx \ F(v(t_k) + \tau \sum_{j=1}^s a_{ij},u(t_k) +  \tau \sum_{j=1}^s a_{ij} f_j^k) \ =: \ f_i^k.
```
Da wegen der exakten numerischen Integration der linearen Funktion
$t \mapsto t$ immer $v(t_k) = t_k$ gilt, sehen wir durch
Koeffizientenvergleich, dass die Invarianz gegenüber Autonomisierung
äquivalent ist zu der Bedingung
```{math}
:label: eq:runge_kutta_ci
c_i \ = \ \sum_{j=1}^s a_{ij}, \qquad i=1,\ldots,s.
```
Aus diesem Grund bestimmt man die Koeffizienten $c_i\in \R$ immer aus
Gleichung {eq}`eq:runge_kutta_ci` und lediglich die Einträge $a_{ij}$
der Matrix $A$ werden aus der Konsistenzbedingung mittels
Taylorentwicklung bestimmt. Damit ist es möglich ein explizites
$s$-stufiges Runge-Kutta Verfahren der Konsistenzordnung $s \in N^+$ zu
konstruieren

(butcher-schema)=
### Butcher-Schema

Allgemein lässt sich ein $s$-stufiges Runge-Kutta Verfahren durch die
Matrix $A$ und die Vektoren $b$, $c$ eindeutig repräsentieren. Diese
lassen sich kompakt in Form eines Schemas angeben.

````{prf:definition} Butcher-Schema
Sei $s \in \N^+$ die Stufe eines Runge-Kutta Verfahrens mit
Koeffizienten, die gegeben sind durch eine Matrix
$A \in \R^{s \times s}$ und den beiden Vektoren
$\vec{c}, \vec{b} \in \R^s$. Dann lässt sich das Runge-Kutta Verfahren
kompakt durch das folgende **Butcher-Schema** repräsentieren:
```{math}
\begin{array}{c|c}
\vec{c} & A \\ \hline & \vec{b}^T
\end{array}
\ = \ \begin{array}{c|cccc}
c_1 & a_{11} & a_{12} & \cdots & a_{1s}\\
c_2 & a_{21} & a_{22} & \cdots & a_{2s}\\
\vdots & \vdots & \vdots & \ddots & \vdots \\
c_s & a_{s1} & a_{s2} & \cdots & a_{ss}\\ \hline
    & b_1 & b_2 & \cdots & b_s 
\end{array}.
```

````

Im folgenden Beispiel geben wir die Butcher-Schemata für vier explizite
Runge-Kutta Verfahren im Vergleich an.

````{prf:example} Explizite Runge--Kutta Verfahren
TODO
<!---

  --------- ---------------------------------------------------------- ------------------------------------------------------------------------------------------------------------------------------------
   $s = 1$                     $\begin{array}{c|c}                     Vorwärts-Euler Verfahren mit: $u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \tau F(t_k, u_\tau(t_k))$. Konsistenzordnung $p=1$.
                                  0 & 0\\ \hline                       
                                       & 1                             
                                   \end{array}$                        
   $s = 2$                     $\begin{array}{c|cc}                    Verbessertes Euler Verfahren mit: $u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \frac{\tau}{2} F(t_k, u_\tau(t_k)) \newline
                                   0 & 0 & 0\\                         \hphantom{aaaaaaaaaas} + \frac{\tau}{2} F(t_k + \tau, u_\tau(t_k) + \tau F(t_k, u_\tau(t_k)))$. Konsistenzordnung $p=2$.
                               1  & 1 & 0\\ \hline                     
                           & \frac{1}{2} & \frac{1}{2}                 
                                   \end{array}$                        
   $s = 3$                    $\begin{array}{c|ccc}                    $u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \frac{\tau}{6} (f_1^k +4f_2^k + f_3^k), \newline
                                 0 & 0 & 0 & 0\\                       \hphantom{aaaaa} f_1^k \ = \ F(t_k, u_\tau(t_k)),\newline
                       \frac{1}{2}  & \frac{1}{2} & 0 & 0\\            \hphantom{aaaaa} f_2^k \ = \ F(t_k + \frac{\tau}{2}, u_\tau(t_k)) + \frac{\tau}{2}f_1^k),\newline
                             1  & -1 & 2 & 0\\ \hline                  \hphantom{aaaaa} f_3^k \ = \ F(t_k + \tau, u_\tau(t_k)) - \tau f_1^k + 2\tau f_2^k)$. Konsistenzordnung $p=3$.
                    & \frac{1}{6} & \frac{4}{6} & \frac{1}{6}          
                                   \end{array}$                        
   $s = 4$                    $\begin{array}{c|cccc}                   Standard Runge-Kutta Verfahren mit: $u_\tau(t_{k+1}) \ = \ u_\tau(t_k) + \frac{\tau}{6} (f_1^k +2f_2^k + 2f_3^k + f_4^k), \newline
                               0 & 0 & 0 & 0 & 0\\                     \hphantom{aaaaa} f_1^k \ = \ F(t_k, u_\tau(t_k)),\newline
                     \frac{1}{2}  & \frac{1}{2} & 0 & 0 & 0\\          \hphantom{aaaaa} f_2^k \ = \ F(t_k + \frac{\tau}{2}, u_\tau(t_k)) + \frac{\tau}{2}f_1^k),\newline
                    \frac{1}{2}  & 0 & \frac{1}{2} & 0 & 0\\           \hphantom{aaaaa} f_3^k \ = \ F(t_k + \frac{\tau}{2}, u_\tau(t_k)) + \frac{\tau}{2}f_2^k),\newline
                           1  & 0 & 0 & 1 & 0\\ \hline                 \hphantom{aaaaa} f_4^k \ = \ F(t_k + \tau, u_\tau(t_k)) - \tau f_3^k)$. Konsistenzordnung $p=4$.
             & \frac{1}{6} & \frac{2}{6} & \frac{2}{6} & \frac{1}{6}   
                                   \end{array}$                        
  --------- ---------------------------------------------------------- ------------------------------------------------------------------------------------------------------------------------------------
-->
````

