(differenzenverfahren-für-randwertprobleme)=
# Differenzenverfahren für Randwertprobleme

Zur Diskretisierung von Randwertproblemen können wir zunächst genauso
vorgehen wie bei der numerischen Lösung von Anfangswertproblemen. Wir
konstruieren also zunächst ein Gitter $0 = x_0 < x_1 < \ldots < x_N = 1$
als Ortsdiskretisierung des Intervalls $[0,1]$. Im einfachsten Fall
wählen wir die Gitterpunkte wieder äquidistant mit Schrittweite
$h\coloneqq \frac{1}{N}$ als
```{math}
\Omega_h \ \coloneqq \ \lbrace x_k \in [0,1] : x_k \coloneqq k\cdot h, k=0,\ldots,N \rbrace.
```
Auf diesem Gitter approximieren wir nun Ableitungen durch finite
Differenzen. Für die zweite Ableitung in $x_k$ verwenden wir in der
Regel ein **Differenzenverfahren zweiter Ordnung** mit
```{math}
u''(x_k) \ \approx \ \frac{u(x_k+h) - 2 u(x_k) + u(x_{k}-h)}{h^2} \ = \ \frac{u(x_{k+1}) - 2 u(x_k) + u(x_{k-1})}{h^2}.
```

Im *nichtäquidistanten Fall* sieht eine entsprechende Approximation
hingegen wie folgt aus:
```{math}
:label: eq:diskretisierung_2.ordnung
u''(x_k) \ \approx \ \frac{1}{h_k} \cdot \left( \frac{u(x_{k+1}) -  u(x_k)}{x_{k+1}-x_k} -\frac{u(x_k)- u(x_{k-1})}{x_k - x_{k-1}} \right).
```
Hier ist zunächst die Frage wie wir den Parameter $h_k > 0$ (abhängig
von den Gitterpunkten $x_{k-1}, x_k, x_{k+1}$) wählen sollen. Dies
können wir als Bedingung an die maximale Konsistenzordnung eines
numerischen Verfahrens formulieren.

Für hinreichend glatte Funktionen $u$ und
```{math}
h \ \coloneqq \!\! \max_{k=0,\ldots,N-1} \vert x_{k+1} -x_k\vert
```
gilt per Taylor-Entwicklung der Funktionswerte $u(x_{k+1})$ und
$u(x_{k-1})$ jeweils im Punkt $x_k \in [0,1]$ dann
```{math}
:label: eq:differenzenverfahren_konsistenz
\begin{split}
    &\hphantom{=} \ \frac{u(x_{k+1}) -  u(x_k)}{x_{k+1}-x_k} -\frac{u(x_k)- u(x_{k-1})}{x_k - x_{k-1}}  \\
    &= \ u'(x_k) + \frac{1}2 u''(x_k) (x_{k+1} -x_k)  +  \frac{1}6 u'''(x_k) (x_{k+1} -x_k)^2  - \\
    & \hspace{0.65cm} u'(x_k) + \frac{1}2 u''(x_k) (x_k-x_{k-1})  -  \frac{1}6 u'''(x_k) (x_k-x_{k-1})^2 + {\cal O}(h^3)\\
    &= \ u''(x_k) \frac{ x_{k+1}-x_{k-1}}2 +   \frac{1}6 u'''(x_k) (x_{k+1} -2 x_k +x_{k-1})( x_{k+1}-x_{k-1} ) + {\cal O}(h^3).
\end{split}
```
Wir sehen, dass wir eine konsistente Approximation der zweiten Ableitung
mit der Wahl von $h_k=\frac{x_{k+1}-x_{k-1}}2$ im Vorfaktor erreichen.

Im äquidistanten Fall erkennen wir an obiger Rechnung wieder, dass die
Konsistenzordnung $\mathcal{O}(h^2)$ erreicht wird, da dann gilt
```{math}
x_{k+1} -2 x_k +x_{k-1} \ = \ \underbrace{(x_{k+1} - x_k)}_{= \:h} - \underbrace{(x_k - x_{k-1})}_{= \:h} \ = \ 0.
```

Die konsistente Wahl des Vorfaktors $h_k=\frac{x_{k+1}-x_{k-1}}2$ lässt
sich ebenfalls anders motivieren: eigentlich berechnen wir numerische
Approximationen erster Ableitungen mit Hilfe des zentralen
Differenzenquotienten als
```{math}
\begin{aligned}
    \frac{u(x_{k+1}) -  u(x_k)}{x_{k+1}-x_k} \ &= \ u'(x_{k+1/2}) + \mathcal{O}(h^2),\\
    \frac{u(x_k)- u(x_{k-1})}{x_k - x_{k-1}} \ &= \ u'(x_{k-1/2}) + \mathcal{O}(h^2)
\end{aligned}
```
mit den beiden (impliziten) Mittelpunkten
```{math}
x_{k-1/2} \ \coloneqq \ \frac{1}2 (x_{k-1} + x_k), \qquad x_{k+1/2} \ \coloneqq \ \frac{1}2 (x_{k+1} + x_k).
```
Diese beiden Mittelpunkte liegen also auf einem versetzten Gitter.

Nun können wir die zweite Ableitung wiederum als numerische
Approximation der Ableitung der Diskretisierungen der ersten Ableitungen
betrachten. Eine weitere Anwendung des zentralen Differenzenquotienten
auf die Werte der versetzten Gitterpunkte $x_{k\pm1/2}$ liefert nämlich
genau eine Approximation der zweiten Ableitung im Punkt
$x_k \in \Omega_h$. Die Schrittweite für diesen zweiten
Diskretisierungsschritt ist genau der von uns berechnete Vorfaktor
```{math}
h_k \ = \ \frac{x_{k+1} - x_{k-1}}2 \ = \ \frac{x_{k+1}+x_k}2 - \frac{x_k+x_{k-1}}2 \ = \ x_{k+1/2} - x_{k-1/2}.
```

(dirichlet-randbedingungen-1)=
### Dirichlet Randbedingungen

Mit den oben eingeführten Differenzenverfahren können wir nun bereits
Gleichungen der Form
```{math}
:label: eq:randwertproblem_numerik_a=0
- u''(x) + c(x) \cdot u(x) \ = \ f(x)
```
für $x \in [0,1]$ direkt diskretisieren. In diesem Kontext betrachten
wir im Folgenden das Randwertproblem mit Dirichlet-Randwertbedingungen,
welches wir als lineares Gleichungssystem für den unbekannten
Lösungsvektor
```{math}
U_h \ = \ (u_1,\ldots,u_{N-1})^T \in \R^{N-1}
```
schreiben können. Die beiden Werte $u_0 = u(0) = g_0$ und
$u_N = u(1) = g_1$ ergeben sich ganz natürlich aus den
Dirichlet-Randbedingungen und können direkt in einen Lösungsvektor für
das gesamte numerische Gitter $\Omega_h$ eingesetzt werden.

Verwenden wir nun das Differenzenschema zweiter Ordnung in
{eq}`eq:diskretisierung_2.ordnung` zur numerischen Diskretisierung der
Differentialgleichung, so können wir die Koeffizienten des
Differenzenschemas für jeden Punkt in eine Zeile der Systemmatrix
$A \in \R^{N-1 \times N-1}$ des Problems schreiben. Hierzu können wir
genauer für die Matrix-Einträge von $A$ definieren als:
```{math}
\begin{aligned}
A_{k,k} \ &\coloneqq \ \frac{1}{h_k} \left( \frac{1}{x_{k+1}-x_k} +  \frac{1}{x_k+x_{k-1}} \right) + c(x_k), \qquad k=1,\ldots,N-1\\
A_{k,k-1} \ \coloneqq \ - \, &\frac{1}{h_k(x_k - x_{k-1})}, \qquad  A_{k-1,k} \ \coloneqq \ - \, \frac{1}{h_k(x_{k+1} - x_{k})} \qquad k=2,\ldots,N-1\\
& \hspace{2cm} A_{k,j} \ \coloneqq \ 0 \qquad \text{sonst, }
\end{aligned}
```
Wir erhalten damit also eine Tridiagonalmatrix
$A \in \R^{N-1 \times N-1}$, die sich im äquidistanten Fall vereinfacht
zu
```{math}
:label: eq:randwertproblem_systemmatrix
A \ \coloneqq \ \frac{1}{h^2} \left( \begin{array}{ccccc} 2 + h^2 c(x_1) & - 1 & 0 & \ldots & 0 \\
-1 & 2 + h^2 c(x_2)  & -1 & \ldots & 0 \\ 0 & -1 & 2 + h^2 c(x_3) & \ldots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 &0 & \ldots & 2 + h^2 c(x_{N-1})\end{array} \right).
```

Für die rechte Seite des linearen Gleichungssystems müssen wir zum einen
die rechte Seite $f(x)$ der Differentialgleichung abbilden, als auch die
Randwertbedingungen $u(0) = g_0$ und $u(1) = g_1$ berücksichtigen.
Hierzu definieren wir uns entsprechend den Vektor $F \in \R^{N-1}$ mit:
```{math}
\begin{aligned}
F_k \ &\coloneqq \  f(x_k),  \qquad k \in \{2, \ldots,N-2\} \\
F_1 \ &\coloneqq \ f(x_1) + \frac{g_0}{h_1(x_1-x_0)}, \\
F_{N-1} \ &\coloneqq \ f(x_{N-1}) + \frac{g_1}{h_{N-1}(x_N-x_{N-1})}.
\end{aligned}
```
Insgesamt müssen wir also das lineare $(N-1) \times (N-1)$
Gleichungssystem $A U_h  \ = \ F$ numerisch lösen um eine diskrete
Lösung des Randwertproblems {eq}`eq:randwertproblem_numerik_a=0` mit
Dirichlet-Randbedingungen $u(0) = g_0$ und $u(1) = g_1$ zu erhalten.

Im allgemeinen Fall einer nichtkonstanten Funktion
$a(x) \not \equiv c \in \R$ können wir die selbe Einsicht über
verschobene Gitter wie oben benutzen. Die erste Ableitung
$(a(x_k)\cdot u'(x_k)'$ interpretieren wir diskret als Wert am
Mittelpunkt der Gitterzellen $x_{k\pm1/2}$. Da wir diese mit der
Funktion $a$ multiplizieren, verwenden wir auch den Wert von $a$ an
diesen Gitterpunkten. Dementsprechend approximieren wir
```{math}
:label: eq:diskretisierung_2.ordnung_mit_a
(a u')'(x_k) \ \approx \ \frac{1}{h_k} \left( a(x_{k+1/2}) \frac{u_{k+1} -  u_k}{x_{k+1}-x_k} - a(x_{k-1/2}) \frac{u_k- u_{k-1}}{x_k - x_{k-1}} \right).
```

In der folgenden Bemerkung wollen wir auf Eigenschaften der Systemmatrix
$A$ eingehen, die allgemein bei der Lösung von Randwertproblemen eine
wichtige Rolle spielen.

````{prf:remark} Eigenschaften der Systemmatrix A
:label: rem:systemmatrix_eigenschaften
Wir können folgende Beobachtungen zu den Eigenschaften der Matrix $A$ in
{eq}`eq:randwertproblem_systemmatrix` machen.

-   Die Matrix $A$ ist **dünnbesetzt** (im Englischen: *sparse*), d.h.
    die meisten Einträge von $A$ sind gleich null. Genauer noch
    existieren von den $(N-1)^2$ möglichen Matrixeinträgen von $A$ nur
    $3\cdot (N-1)-2 = 3N - 5$ Einträge, die nicht verschwinden. Dies hat
    einige Konsequenzen für die effiziente Speicherung von $A$, sowie
    bei der Lösung des Systems $AU_h =F$. Wir können die Matrix im
    sogenannten **Sparse-Format** speichern. Hierbei handelt es sich in
    der Regel um eine lineare Liste für alle Nichtnulleinträge der Form
    $(i,j,A_{ij})$. Statt $(n-1)^2$ reeller Zahlen speichern wir hier
    nur $6N-10$ ganze und $3N-5$ reelle Zahlen. Gegenüber potentiell
    $(N-1)^2$ zu speichernden reellen Zahlen ist dies eine enorme
    Einsparung für großes $N \in \N$.

    Bei der numerischen Lösung des linearen Systems ist es ebenfalls
    vorteilhaft diese Struktur zu benutzen. Bei direkten Verfahren wie
    der LR-Zerlegung jedoch ist dies nicht der Fall. Dort kann es zum
    sogenannten **Fill-In** kommen, d.h. $L$ und $R$ sind nicht mehr
    dünnbesetzt. Grob können wir uns den Effekt bei der LR-Zerlegung wie
    bei der Integration von $-u''=f$, $u(0) = g_0$, $u(1) =g_1$
    vorstellen. Hier erhalten wir durch die Integration von $0$ bis $x$
    eine Integraldarstellung für $u'(x)$, die im diskreten einer linken
    unteren Dreiecksmatrix enstpricht. Integrieren wir wiederum von $1$
    bis $x$ um $u(x)$ zu erhalten, entspricht dies einer rechten oberen
    Dreiecksmatrix. Beide Matrizen haben dann keine dünnbesetzte
    Struktur mehr, da das diskrete Integral alle vorhandenen
    Gitterpunkte benutzt.

-   Für nichtnegative Funktionen $c(x) \geq 0$ ist die Matrix $A$
    **schwach diagonaldominant**, d.h. es gilt
    $A_{kk} \geq \sum_{j \neq k} \vert A_{kj} \vert$ (bzw. auch in den
    Spalten $A_{kk} \geq \sum_{j \neq k} \vert A_{jk} \vert$). In der
    ersten Zeile für $k=1$ und in der letzten Zeile für $k=N-1$ der
    Matrix $A$ sind diese Ungleichungen sogar strikt.

    Darüber hinaus hat $A$ nur positive Diagonaleinträge und
    nichtpositive Nebendiagonaleinträge. Wie wir später sehen werden ist
    die Matrix $A$ dann invertierbar und die Inverse enthält nur
    nichtnegative Einträge. Dies ist das diskrete Äquivalent zum
    Maximumsprinzip, denn hat $F$ nur nichtnegative Einträge, so gilt
    auch $U_h = A^{-1} F \geq 0$.

-   Die Matrix $A$ ist **symmetrisch**, d.h. es gilt $A_{jk}  = A_{kj}$
    für alle $k,j$. Dies ist eine Konsequenz aus der Divergenzform
    {eq}`eq:divergenzform`, da der Operator $L: u \mapsto - (au')' + cu$
    formal selbstadjungiert ist. Es gilt mit partieller Integration für
    $u$ und $v$ und angenommenen Dirichlet-Nullrandwerten:
    ```{math}
    \begin{aligned}
    \langle Lu, v \rangle_{L^2} \ &= \ \int_0^1 (Lu)(x)~v(x)\, \mathrm{d}x \ = \ \int_0^1 (-(a(x) u'(x))v(x) + c(x) u(x) v(x)\, \mathrm{d}x \\
    &= \ \int_0^1 a(x) u'(x) v'(x)  + c(x) u(x) v(x)\, \mathrm{d}x \\
    &= \ \int_0^1 (-(a(x) v'(x)) u(x) + c(x)  v(x)u(x))\, \mathrm{d}x \ = \ \langle v, Lu \rangle_{L^2} .
    \end{aligned}
    ```

    Ist $c(x) \geq 0$, dann ist die Matrix $A$ sogar symmetrisch positiv
    definit. Dies ist für $u(x) \not \equiv 0$ eine Konsequenz aus
    ```{math}
    \langle L u, u \rangle \ = \ \int_0^1 a(x) \cdot |u'(x)|^2   + c(x) \cdot u(x)^2 \mathrm{d}x \ > \ 0.
    ```

````

Mit den oben diskutierten Eigenschaften der Systemmatrix $A$ erhalten
wir insbesondere ihre Invertierbarkeit aus der positiven Definitheit.
Also existiert eine eindeutige Lösung $U_h \in \R^{N-1}$ des linearen
Gleichungssystems $AU_h=F$.

(konvergenz-von-differenzenverfahren)=
## Konvergenz von Differenzenverfahren

In diesem Abschnitt wollen wir uns mit der Frage beschäftigen welche
Bedingungen gelten müssen, damit ein Differenzenverfahren zur
numerischen Lösung eines Randwertproblems gegen die echte Lösung der
Differentialgleichung konvergiert. Wie wir bereits gesehen haben führt
das Differenzenverfahren zu einem linearen Gleichungssystem der Form
$A U_h = F$, wobei $U_h \in \R^{N-1}$ die numerische Approximation der
echten Lösung $u \in C^2([0,1])$ in den Gitterpunkten darstellt, d.h.,
es gilt $(U_h) _k = u_k \approx u(x_k)$ für $k=1,\ldots,N-1$.

(m-matrizen)=
### M-Matrizen

Wie wir im Folgenden feststellen werden hängt die Stabilität und
folglich damit die Konvergenz eines Differenzenverfahrens maßgeblich von
den Eigenschaften der Systemmatrix $A \in \R^{N-1 \times N-1}$ ab. Die
in {prf:ref}`rem:systemmatrix_eigenschaften` diskutierten Eigenschaften
der speziellen Systemmatrix $A$ motivieren die folgende Definition.

````{prf:definition} M-Matrix
:label: def:m-matrix
Eine Matrix $A \in \R^{n \times n}$ heißt **M-Matrix**, wenn sie
folgende Eigenschaften erfüllt:

-   $A$ hat nur positive Diagonaleinträge und nichtpositive
    Nebendiagonaleinträge.

-   $A$ ist schwach diagonaldominant, d.h. es gilt
    $A_{kk} \: \geq \: \sum_{j \neq k} \vert A_{jk} \vert$.

-   Für mindestens ein $j \in \{1,\ldots,n\}$ gilt
    $A_{kk} \: > \: \sum_{j \neq k} \vert A_{jk} \vert$.

````

Wie wir leicht nachrechnen können erfüllt unsere Diskretisierung diese
Eigenschaft:

````{prf:corollary} 
Sei $A$ die Matrix aus dem obigen Differenzenverfahren in
{eq}`eq:randwertproblem_systemmatrix` für eine nichtnegative Funktion
$c(x) \geq 0$ für $x \in [0,1]$. Dann ist $A$ eine M-Matrix.

````

Die letzte Eigenschaft einer M-Matrix in {prf:ref}`def:m-matrix` ist
entscheidend für deren Invertierbarkeit, da sonst beispielsweise auch
die folgende Matrix zulässig wäre:
```{math}
A \ \coloneqq \ \left( \begin{array}{cc} 1 & -1 \\ -1 & 1 \end{array} \right)
```
Diese Eigenschaft können wir allgemein für M-Matrizen beweisen und
darüber hinaus etwas über die Einträge der Inversen aussagen.

````{prf:lemma} Inverse einer M-Matrix
:label: lem:m-matrix_inverse
Sei $A \in \R^{n \times n}$ für $n \in \N$ eine M-Matrix. Dann ist $A$
invertierbar und $A^{-1}$ hat nur nicht-negative Einträge.

````

````{prf:proof} 
Wir nehmen zunächst an, dass $A$ irreduzibel ist, d.h. es existiert eine
Permutation $\pi \in \Pi_n$ von $\{1,\ldots,n\}$, sodass
$A_{\pi(i)\pi(i+1)} \neq 0$ gilt. Andernfalls können wir analog wie im
Folgenden für alle irreduziblen Blöcke von $A$ vorgehen.

Zunächst zeigen wir, dass $A$ invertierbar ist, d.h. wir müssen zeigen,
dass der Nullraum trivial ist. Sei $z \in {\cal N}(A)$ ein Vektor des
Nullraums und es sei ein Vektoreintrag $z_m$ für
$m \in \lbrace 1, \ldots, n\rbrace$ so, dass $z_m \leq z_k$ für alle $k$
gilt. Nehmen wir nun an, dass $z_m < 0$ ist, dann folgt
```{math}
A_{mm} z_m \ = \ - \sum_{j \neq m} A_{mj} z_j \ =\ \sum_{j \neq m} \vert A_{mj} \vert  z_j \geq \sum_{j \neq m} \vert A_{mj} \vert  z_m .
```
Ist $z_m <0$, so folgt wegen der Diagonaldominanz
$\sum_{j \neq m} \vert A_{mj} \vert  z_m \geq A_{mm} z_m$. Dies ist aber
nur möglich, wenn $z_j = z_m$ für alle $j$ mit $A_{mj} \neq 0$ gilt.
Damit ist auch $z_j$ minimal und wir können das gleiche Argument auf
$z_j$ anwenden, wegen der Irreduzibilität erreichen wir dann
schrittweise alle Einträge von $z$. Dies bedeutet der konstante Vektor
$z_j = 1$ für alle $j$ ist eine Lösung, was aber der Bedingung
```{math}
A_{kk} \ > \ \sum_{j \neq k} \vert A_{jk} \vert \ = \ -\sum_{j \neq k}  A_{jk},
```
für ein $k$ gilt. Also muss $\min_k z_k\geq 0$ gelten. Analog können wir
aber auch $\max_k z_k \leq 0$ zeigen, also bleibt nur $z \equiv 0$ im
Nullraum.

Um zu zeigen, dass die Inverse von $A$ nur nicht-negative Einträge hat,
können wir zeigen, dass die Lösung $z$ von $Az=y$, mit $y$ gleich einem
Einheitsvektor, nur nicht-negative Einträge hat. Dies ist natürlich
äquivalent dazu, dass die Lösung von $Az = y$ mit nicht-negativem $y$
nur nicht-negative Einträge hat. Dazu benutzen wir das selbe Argument
wie oben: Sei $z_k = \min_j z_j < 0$. Dann ist
```{math}
A_{kk} z_k \ = \ - \sum_{j\neq k} A_{kj} z_j + y_k \geq -  \sum_{j\neq k} A_{kj} z_k \ \geq \ A_{kk} z_k
```
mit Gleichheit wenn $y_k = 0$ und $x_j = x_k$ für alle $j$. Dann folgt
aber auch $y_j =0$ für alle $j$ und deshalb $z = 0$, ein Widerspruch zu
$z_k < 0$. ◻

````

(konvergenz-von-differenzenverfahren-1)=
### Konvergenz von Differenzenverfahren

Um die Konvergenz eines numerischen Differenzenverfahrens für ein
Randwertproblem zu untersuchen, benötigen wir wieder die üblichen
Begriffe der Konsistenz und Stabilität. Diese wollen wir im Folgenden
einführen bevor wir hinreichende Bedingungen angeben. Sei
$u \in C^2([0,1])$ die Lösung des Randwertproblems und dementsprechend
sei $U \coloneqq (u(x_k)) \in \R^{N-1}$ die Auswertung der Lösung auf
dem Diskretisierungsgitter. Sei außerdem $U_h \in \R^{N-1}$ die
numerische Lösung des linearen Gleichungssystems $A U_h = F$ des
Randwertproblems. Dann können wir den Unterschied zwischen der
numerischen Lösung und der echten Lösung der Differentialgleichung
bezüglich der rechten Seite $F$ betrachten mit
```{math}
:label: eq:randwertproblem_konsistenz
A (U_h - U) \ = \ F - A U \ =: \ G \in \R^{N-1}.
```
Die rechte Seite $G$ misst also die Auswirkung des
Diskretisierungsfehlers, der durch das Aufstellen der Systemmatrix $A$
entsteht, bei Anwendung auf die echte Lösung der Differentialgleichung.
Für die in {eq}`eq:diskretisierung_2.ordnung_mit_a` betrachtete
Diskretisierung ergibt sich somit für die rechte Seite $G$ unter
Verwendung von $F_k = f(x_k) = (a\cdot u')'(x_k)$ für $k=1,\ldots,N-1$
```{math}
G_k \ = \ (a u')'(x_k) - \frac{1}{h_k} \left( a(x_{k+1/2}) \frac{u_{k+1} -  u_k}{x_{k+1}-x_k} - a(x_{k-1/2}) \frac{u_k- u_{k-1}}{x_k - x_{k-1}} \right).
```
Wir beachten dabei, dass der Term $c(x)\cdot u(x)$ als
Funktionsauswertung exakt diskretisiert wird und daher nicht weiter zum
Fehler beiträgt. Hierauf aufbauend können wir also den Konsistenzfehler
definieren.

````{prf:definition} Konsistenzfehler und -ordnung
:label: def:randwertproblem_konsistenzfehler
Sei $A \in \R^{N-1 \times N-1}$ eine Systemmatrix, die ein
Differenzenverfahren für ein Randwertproblem in Divergenzform
{eq}`eq:divergenzform` realisiert. Dann definieren wir den **globalen
Konsistenzfehler** des Verfahrens als
```{math}
K_h \ \coloneqq \ \Vert G\Vert_\infty \ = \! \max_{k=1,\ldots,N-1} \vert G_k \vert
```
Ein Differenzenverfahren heißt **konsistent von der Ordnung**
$\mathbf{p \in \N}$, wenn für den Konsistenzfehler $K_h = {\cal O}(h^p)$
gilt.

````

Wie wir bereits in {eq}`eq:differenzenverfahren_konsistenz`
nachgerechnet haben für den Spezialfall einer konstanten Funktion
$a(x) \equiv 1$ für $x \in [0,1]$ ist das diskutierte
Differenzenverfahren konsistent von Ordnung $2$. Dieses Ergebnis lässt
sich leicht mittels Taylorentwicklung für allgemeine Funktionen $a(x)$
erweitern.

Um eine Konvergenz des Differenzenverfahrens zu erhalten benötigen wir
eine Abschätzung an die Inverse von $A$, denn ist $A$ invertierbar so
gilt für den Konvergenzfehler $E_h$ mit der Submultiplikativität der
Matrixnorm:
```{math}
E_h \ := \ \Vert U_h - U\Vert_\infty \ = \ \Vert A^{-1} \cdot \underbrace{A(U_h - U)}_{= \:G} \Vert_\infty \ \leq \ \Vert A^{-1}  \Vert_\infty \cdot \Vert G \Vert_\infty \ = \ \Vert A^{-1}  \Vert_\infty \cdot K_h.
```
Wir sehen aus dieser Abschätzung sofort, dass aus Konsistenzordnung
$p\in \N$ auch schon direkt Konvergenzordnung $p$ folgt, wenn
$\Vert A^{-1}\Vert_\infty$ beschränkt ist für $N\rightarrow \infty$ bzw.
$h \rightarrow 0$ und somit $E_h \in \mathcal{O}(h^p)$ gilt. Dies
bezeichnen wir folglich als Stabilität und halten dies entsprechend in
einer Definition fest.

````{prf:definition} Stabilität von Differenzenverfahren
Wir nennen ein Differenzenverfahren eines Randwertproblems **stabil**,
falls für den Konvergenzfehler folgende Abschätzung gilt:
```{math}
E_h \ = \ \Vert U_h - U\Vert_\infty \ \leq \ C \cdot K_h,
```
wobei $K_h$ der globale Konsistenzfehler aus
{prf:ref}`def:randwertproblem_konsistenzfehler` ist und $C > 0$ eine von
der Diskretisierungsschrittweite
$h \coloneqq \max_{k=1,\ldots,N-1} |x_{k-1} - x_k|$ unabhängige
Konstante.

````

Wie wir sehen werden ist die M-Matrix Eigenschaft eine hinreichende
Bedingung für die Stabilität eines Differenzenverfahrens. Wie beim
Maximumsprinzip für die Differentialgleichung in
{ref}`s:randwertproblem_existenz_eindeutigkeit` werden wir
Vergleichslösungen suchen, um Fehlerschranken zu erhalten. Dazu ist es
wichtig, dass die Vergleichslösung unabhängig von der Schrittweite
$h > 0$ (bzw. der Anzahl der Gitterpunkte $N \in \N$) ist.

Im Folgenden sei $v \in C^2([0,1])$ die Lösung des Randwertproblems mit
*homogenen Dirichlet-Randwertbedingungen*
```{math}
\begin{split}
 - (a(x)\cdot v'(x))' +c(x) \cdot v(x) \ &= \ f(x), \qquad x \in (0,1),\\
 v(0) \ = \ v(1) \ &= \ 0.
 \end{split}
```
Außerdem sei $u \in C^2([0,1])$ die Lösung des Randwertproblems mit
*allgemeinen Dirichlet-Randwertbedingungen*
```{math}
\begin{split}
 - (a(x)\cdot u'(x))' +c(x) \cdot u(x) \ &= \ f(x), \qquad x \in (0,1),\\
 u(0) \ = \ g_0, \quad u(1) \ &= \ g_1.
 \end{split}
```
Die numerische Lösung von $A U_h = F$ bezeichnen wir mit
$U_h = (u_k)_{k=1,\ldots,N-1}$. Mit dem üblichen Konsistenzfehler $G$
gilt
```{math}
A (U - U_h) \ = \ G,
```
wobei $U = (u(x_k))_{k=1,\ldots,N-1}$ die Auswertung der Lösung des
allgemeinen Dirichletproblems an den Gitterpunkten ist. Nun betrachten
wir die beiden Vektoren
```{math}
H_\pm \ \coloneqq \ U - U_h \ \pm \ 2 \cdot K_h \cdot V
```
mit dem globalen Konsistenzfehler $K_h = \Vert G \Vert_\infty$ und der
Lösung des homogenen Dirichletproblems ausgewertet an den Gitterpunkten
$V = (v(x_k))_{k=1,\ldots,N-1}$. Das folgende Lemma liefert uns
praktische Abschätzungen.

````{prf:lemma} Abschätzungen durch Vergleichslösungen
:label: lem:abschaetzungen_vergleichsloesungen
Mit den obigen Definitionen von $U,U_h,V \in \R^{N-1}$ und dem globalen
Konsistenzfehler $K_h$ betrachten wir die beiden Vektoren
$H_\pm \coloneqq U - U_h \ \pm \ 2 \cdot K_h \cdot V \in \R^{N-1}$. Sei
außerdem $A \in \R^{N-1 \times N-1}$ eine M-Matrix.

Dann gelten für hinreichend kleine Schrittweiten $h > 0$ die
Ungleichungen
```{math}
:label: eq:randwertproblem_fehler_ungleichungen
\begin{split}
    A \cdot H_+ \ = \ A \cdot (U - U_h+2 \cdot K_h \cdot V) \ &\geq \ 0, \\
    A \cdot H_- \ = \ A \cdot (U - U_h-2 \cdot K_h \cdot V) \ &\leq \ 0.
\end{split}
```
Die Ungleichungen {eq}`eq:randwertproblem_fehler_ungleichungen` sind
hierbei komponentenweise gemeint.

````

````{prf:proof} 
Wegen der Konsistenz des Verfahrens wird für hinreichend kleine
Schrittweiten $h > 0$ die Differenz
```{math}
A \cdot V - (-(av')'(x_k) + c(x_k))_{k=1,\ldots,N-1} \ = \ A \cdot V - (1)_{k=1,\ldots,N-1}
```
beliebig klein, insbesondere gilt dann
```{math}
A \cdot V \ \geq \ (\frac{1}2)_{k=1,\ldots,N-1}.
```
Setzen wir nun ein, so folgt
```{math}
A \cdot  (U - U_h+2 \cdot K_h \cdot V) \ = \ G_h + 2 \cdot K_h \cdot A \cdot V \ \geq \ G + K_h (1)_{k=1,\ldots,N-1} \ \geq \ 0
```
und da $K_h = \Vert G \Vert_\infty$ gilt. Analog folgt die zweite
Ungleichung. ◻

````

Nun können wir eine hinreichende Bedingung für die Stabilität eines
Differenzenverfahrens formulieren.

````{prf:theorem} Stabilität eines Differenzenverfahrens
Sei $A \in \R^{N-1 \times N-1}$ eine Systemmatrix, die ein
Differenzenverfahren für ein Randwertproblem in Divergenzform
{eq}`eq:divergenzform` realisiert. Seien $U,U_h,V \in \R^{N-1}$
definiert wie oben und $K_h$ sei der globale Konsistenzfehler. Außerdem
sei $A$ eine M-Matrix.

Dann ist das Differenzenverfahren stabil.

````

````{prf:proof} 
Um die Stabilität des Differenzenverfahrens zu zeigen betrachten wir
zunächst wieder die Vektoren
```{math}
H_\pm \ \coloneqq \ U - U_h \, \pm \, 2 \cdot K_h \cdot V.
```
Da $A$ eine M-Matrix ist nach Voraussetzung können wir
{prf:ref}`lem:abschaetzungen_vergleichsloesungen` anwenden und es gelten
somit die Abschätzungen {eq}`eq:randwertproblem_fehler_ungleichungen`.
Wegen der Nichtnegativität der Matrixeinträge der Inversen $A^{-1}$
können wir dann komponentenweise folgern
```{math}
\begin{split}
H_+ \ = \ U - U_h + 2 \cdot K_h \cdot V \ \geq \ A^{-1} \cdot 0 \ = \ 0,\\
H_- \ = \ U - U_h - 2 \cdot K_h \cdot V \ \leq \ A^{-1} \cdot 0 \ = \ 0,\\
\end{split}
```

Somit können wir also den Fehler zwischen der Lösung der
Differentialgleichung $U$ und der numerischen Lösung des
Differenzenverfahren $U_h$ nach oben und unten abschätzen durch:
```{math}
- 2 \cdot K_h \cdot V \ \leq \ U - U_h \ \leq \ 2 \cdot K_h \cdot V.
```

Da die Lösung $v \in C^2([0,1])$ insbesondere stetig ist, nimmt sie ihr
betragliches Maximum auf dem kompakten Intervall $[0,1] \subset \R$ an
und wir können somit ebenfalls unabhängig von der Schrittweite $h > 0$
des Gitters abschätzen:
```{math}
\Vert V \Vert_\infty \ \leq \ \max_{x \in [0,1]} \vert v(x) \vert.
```

Mit diesen Abschätzungen erhalten wir schließlich eine gleichmäßige
Schranke für den Konvergenzfehler $E_h$ durch:
```{math}
E_h \ = \ ||U_h - U||_\infty \ \leq \ 2 \cdot ||V||_\infty \cdot K_h \ \leq \ \underbrace{2 \cdot \max_{x \in [0,1]} \vert v(x) \vert}_{=: \: C \: > \: 0} \cdot \: K_h.
```
 ◻

````

Mit der Stabilität eines Differenzenverfahrens erhalten wir als direkte
Folgerung das folgende Konvergenzresultat.

````{prf:corollary} Konvergenz eines Differenzenverfahrens
Seien $U,U_h,V \in \R^{N-1}$ definiert wie oben und $K_h$ sei der
globale Konsistenzfehler. Dann gilt für $h$ hinreichend klein
```{math}
E_h \ = \ \Vert U_h - U \Vert_\infty \ \leq \ 2 \cdot \Vert v \Vert_\infty \cdot K_h.
```
Insbesondere stimmen die Konsistenz- und Konvergenzordnung überein.

````

Wir sehen, dass die obige Theorie auch leicht auf andere
Differenzenverfahren anwendbar ist, solange Konsistenz vorliegt und die
entsprechende M-Matrix Eigenschaft erfüllt ist. So kann man die Aussagen
auch auf partielle Differentialgleichungen erweitern, etwa wenn wir
```{math}
- \nabla \cdot (a(x) \nabla u(x) ) + c(x) \cdot u(x)  \ = \ f(x), \qquad x \in [0,1]^2
```
auf einem rechteckigen Gebiet lösen wollen, mit
$\nabla \cdot (a(x) \nabla u(x) ) = \sum_{i=1}^n \partial_{x_i} (a(x) \partial_{x_i} u(x)).$
Legen wir darüber ein Gitter und diskretisieren die zweiten Ableitungen
in jeder Richtung analog, so erhalten wir wieder ein konsistentes
Verfahren, das durch ein lineares System mit einer M-Matrix beschrieben
wird. Damit können wir eine völlig analoge Theorie durchführen und
Konvergenz beweisen. Der einzige kleine Unterschied ist, dass wir keine
Tridiagonalmatrix erhalten, sondern eine etwas allgemeinere dünnbesetzte
Matrix. Dies ändert aber wenig an der Struktur, nur die numerische
Lösung des Gleichungssystems $AU=F$ wird deutlich aufwändiger, da wir
viel mehr Gitterpunkte benötigen. Wir werden uns deshalb später noch mit
der numerischen Lösung dünnbesetzter linearer Systeme beschäftigen. Der
einzige Nachteil in mehreren Dimensionen ist, dass finite Differenzen
kanonisch für rechteckige Gitter verwendbar sind. Hat man andere Gebiete
auf denen man partielle Differentialgleichungen lösen will, kommen eher
Verfahren wie finite Elemente zum Einsatz, deren Idee wir im Folgenden
noch kurz diskutieren wollen.
