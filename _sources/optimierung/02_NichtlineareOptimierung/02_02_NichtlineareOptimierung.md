(s:abstiegsverfahren)=
# Abstiegsverfahren

Zu Anfang dieser Vorlesung möchten wir eine Klasse von Algorithmen zur
Optimierung von Funktionen besprechen, die einer simplen und
anschaulichen Idee folgen: die sogenannten *Abstiegsverfahren* oder auch
*Liniensuchverfahren* (im Englischen: *line search methods*). Diese
wurden bereits kurz im Zusammenhang mit dem Gauss-Newton Verfahren in
der Vorlesung Einführung in die Numerik in {cite:p}`numerik1` erwähnt,
jedoch nicht ausführlich diskutiert. Dies wollen wir im Folgenden
nachholen.

Sei im Folgenden $\Omega \subset \mathbb{R}^n$ ein offenes,
zusammenhängendes Gebiet und sei
$F \colon \Omega \rightarrow \mathbb{R}$ eine differenzierbare,
reellwertige Zielfunktion. Die allgemeine Idee der Abstiegsverfahren ist
es nun ausgehend von einem aktuellen Punkt $x_k \in \Omega$ einen
Schritt in eine Richtung $p_k \in \mathbb{R}^n$ zu machen, so dass der
Funktionswert von $F$ in dem neuen Punkt $x_{k+1}$ abnimmt. Das bedeutet
wir versuchen ein **allgemeines Abstiegsverfahren** der Form
```{math}
:label: eq:abstiegsverfahren
x_{k+1} \ = \ x_k + \alpha_k p_k, \quad \alpha_k \in \R^+,
```
zu konstruieren, das in jedem Schritt den Funktionswert $F(x_k)$
verringert bis keine Verbesserung mehr möglich ist. Die Schrittweite
$\alpha_k > 0$ in Richtung $p_k \in \mathbb{R}^n$ wird häufig in jedem
Schritt des Abstiegsverfahrens neu gewählt.

(ss:gradient_descent)=
## Gradientenabstiegsverfahren

Zuerst beschäftigen wir uns mit dem wohl bekanntesten Abstiegsverfahren,
dem *Gradientenabstiegsverfahren* (im Englischen: *gradient descent
(GD)*). Dieses wird wegen seiner Einfachheit für viele numerische
Optimierungsprobleme verwendet.

Da $F$ differenzierbar ist können wir dessen Gradienten in jedem Punkt
$x \in \Omega$ betrachten. Wir gehen im Folgenden immer davon aus, dass
wir den Gradienten $\nabla F$ der Zielfunktion $F$ analytisch bestimmen
können und eine Auswertung $\nabla F(x_k)$ für eine Folge von Punkten
$(x_i)_{i\in\N} \subset \Omega$ numerisch durchführbar ist. Wir werden
also nicht versuchen den Gradienten mit Hilfe von *finiten
Differenzenverfahren* numerisch zu approximieren.

Der Gradient $\nabla F(x)$ beschreibt bekanntlich eine Richtung des
stärksten Anstiegs der Funktion $F$ im Punkt $x$ und dementsprechend
zeigt der negative Gradient $-\nabla F(x)$ in die Richtung des stärksten
Abstiegs. Das lässt sich auch formal zeigen indem wir uns die
Taylorapproximation erster Ordnung der Funktion $F$ in allgemeine
Richtung $p$ mit Schrittweite $\alpha > 0$ für den $k$-ten Schritt des
Abstiegsverfahren näher anschauen. Hier gilt nämlich
```{math}
F(x_k + \alpha p) \ = \ F(x_k) + \alpha \langle \nabla F(x_k), p \rangle + \mathcal{O}( \nabla^2 F(x_k)),
```
wobei $\nabla^2 F$ die Hesse-Matrix der Zielfunktion $F$ bezeichnet. Für
kleine Schrittweiten $\alpha$ wird klar, dass die Änderung der
Funktionswerte von $F$ im Wesentlichen von der Größe des Terms
$\langle \nabla F(x_k), p \rangle$ bestimmt wird. Wir können also für
eine maximale Verringerung der Funktion $F$ nach derjenigen Richtung mit
Einheitslänge suchen, die das folgende Minimierungsproblem löst:
```{math}
:label: eq:optimization_direction
\min_{p \in \mathbb{R}^n} \langle \nabla F(x_k), p \rangle \quad \text{ mit } \quad ||p|| = 1.
```
Da außerdem gilt
```{math}
:label: eq:optimization_direction_cosine
 \langle \nabla F(x_k), p \rangle \ = \ ||\nabla F(x_k)|| \cdot ||p|| \cdot \cos(\theta),
```
wobei $\theta$ der Winkel zwischen den Vektoren $\nabla F(x_k)$ und $p$
bildet, können wir das Problem {eq}`eq:optimization_direction`
umschreiben zu
```{math}
\min_{\theta \in [0,2\pi)} ||\nabla F(x_k)|| \cdot \cos(\theta).
```
Der Kosinus nimmt sein Minimum von $\cos(\theta)=-1$ in $\theta = \pi$
an. Eingesetzt in {eq}`eq:optimization_direction_cosine` erhalten wir
damit, dass die optimale Richtung $p \in \mathbb{R}^n$ folgenden
Zusammenhang erfüllen muss:
```{math}
\langle p, \frac{\nabla F(x_k)}{||\nabla F(x_k)||}\rangle \ = \ -1.
```
Da die unbekannte Richtung $p \in \R^n$ aber normiert sein soll, folgt
damit aber auch schon, dass gilt
```{math}
:label: eq:optimal_direction
p \ = \ -\frac{\nabla F(x_k)}{||\nabla F(x_k)||}.
```
Das bedeutet, dass der Funktionswert von $F$ am stärksten in Richtung
des negativen Gradienten abnimmt.

Da wir daran interessiert sind die Zielfunktion $F$ schnellstmöglich zu
minimieren, macht es Sinn in eben dieser Richtung nach einem lokalen
Minimum zu suchen. Aus dieser Idee heraus lässt sich bereits ein sehr
simpler Algorithmus zur Minimierung von $F$ formulieren. Sei
$x_0 \in \Omega$ ein beliebiger Startwert. Dann können wir iterativ eine
Folge von Punkten $x_1,x_2,\ldots$ in $\Omega$ bestimmen, so dass die
entsprechenden Funktionswerte von $F$ monoton fallen sollten:
```{math}
:label: eq:gradient_descent_simple
x_{k+1} \ = \ x_k - \nabla F(x_k).
```
Intuitiv stoppt man das Iterationsschema
{eq}`eq:gradient_descent_simple` sobald die Folge der Funktionswerte
$F(x_k)$ nicht mehr kleiner wird. Man sieht leicht ein, dass das simple
Iterationsschema {eq}`eq:gradient_descent_simple` ein Spezialfall des
allgemeinen Abstiegsverfahrens {eq}`eq:abstiegsverfahren` mit einer
festen Schrittweite $\alpha_k = 1$ und Richtungsvektor
$p_k = -\nabla F(x_k)$ ist. Diese einfache Methode lässt sich durch
folgenden Algorithmus implementieren.

````{prf:algorithm} Simples Gradientenabstiegsverfahren
:label: alg:gradient_descent_simple

**function** $[x^*, F(x^*)]=$`gradientDescentSimple`$(F,\nabla F, x_0)$  

\# Initialisierung  
$x_k = x_0$  
$F(x_{k+1}) = -\infty$  

**while** $F(x_{k+1})-F(x_k) < 0$ **do**  
    \# Update in Richtung des größten Gradientenabstiegs  
    $x_{k+1} = x_k - \nabla F(x_k)$  
**end while**

\# Ausgabe des letzten Punktes  
$x^* = x_k$  
$F(x^*) = F(x_k)$
````

Unglücklicherweise ist Algorithmus
{prf:ref}`alg:gradient_descent_simple` in dieser Form praktisch nicht
anwendbar. Warum dies so ist sieht man leicht an folgendem Beispiel.

````{prf:example} Simpler Gradientenabstieg
:label: ex:gradient_descent_simple
Sei $\Omega = \mathbb{R}$ und sei $F(x) \coloneqq |x|$. Man beachte,
dass $F$ überall differenzierbar ist, außer an der Stelle $x=0$. Das
eindeutig bestimmte, globale Minimum der konvexen, nichtlinearen
Funktion $F$ wird ebenfalls in $x^* = 0$ angenommen.

Definiert man $\nabla F(0) \coloneqq 0$ in der Singularität, so erhält
man das Minimum $x^* = 0$ durch Algorithmus
{prf:ref}`alg:gradient_descent_simple` nur für Startwerte
$x_0 \in \mathbb{Z}$. Wähle zum Beispiel den Startwert $x_0 = 0.5$, so
terminiert das Gradientenabstiegsverfahren bereits nach dem ersten
Schritt ohne eine gute Näherung an $x^* = 0$ zu liefern.

````

Wie das {prf:ref}`ex:gradient_descent_simple` zeigt, besteht bei dem
Gradientenabstiegsverfahren in Algorithmus
{prf:ref}`alg:gradient_descent_simple` die Gefahr einen stationären
Punkt und somit ein potentielles Minimum zu überspringen. Aus diesem
Grund kommt man auf die Idee die Schrittweite $\alpha_k > 0$ in
{eq}`eq:abstiegsverfahren` **genügend klein** zu wählen. Damit diese
feste Schrittweite unabhängig von der Magnitude des Gradienten
$\nabla F$ ist, normiert man in der Regel die Richtung des steilsten
Gradientenabstiegs durch die Norm des Gradienten, d.h., wir erhalten
eine **steuerbare Version des Gradientenabstiegsverfahrens** in
{eq}`eq:gradient_descent_simple` durch:
```{math}
:label: eq:gradient_descent_control
x_{k+1} \ = \ x_k - \tau \frac{\nabla F(x_k)}{||\nabla F(x_k)||}, \quad \tau > 0 .
```
Das Iterationsschema {eq}`eq:gradient_descent_adaptive` ist wiederum ein
Spezialfall des allgemeinen Abstiegsverfahrens in
{eq}`eq:abstiegsverfahren` für eine feste Schritteweite
$\alpha_k = \tau > 0$ und den Richtungsvektor
$p_k = -\nabla F(x_k) / ||\nabla F(x_k)||$.

Es ist klar, dass durch eine kleinere Schrittweite $\tau > 0$ ein
stationärer Punkt $x^* \in \Omega$ der Zielfunktion $F$ immer besser
angenähert werden kann. Leider erhöht sich aber gleichzeit die benötigte
Iterationszahl zur Erreichung einer erwünschten Genauigkeit
$|F(x^*) - F(x_k)| < \epsilon$ je kleiner man die Schrittweite $\tau$
wählt. Man muss also bei der Wahl der Schrittweite einen Kompromiss
zwischen Genauigkeit der numerischen Approximation und der Laufzeit des
Verfahrens eingehen.

Eine weiterführende Idee ist es die Schrittweiten **adaptiv** zu wählen,
dass heißt man passt sie innerhalb des Iterationsschemas an die
Funktionswerte von $F$ geeignet an. Ein Iterationsschema, dass eine
immer kleiner werdende Schrittweite $\alpha_k > 0$ verwendet, lässt sich
für einen fest gewählten Reduktionsfaktor $0 < \sigma < 1$ wie folgt
angeben:
```{math}
:label: eq:gradient_descent_adaptive
\begin{split}
\alpha_{k+1} \ &= \
\begin{cases}
\ \alpha_k, \quad &\text{ falls } F\left(x_k - \alpha_k \frac{\nabla F(x_k)}{||\nabla F(x_k)||}\right) < F(x_k) , \\
\ \sigma \alpha_k, \quad &\text{ sonst}.
\end{cases},\\
x_{k+1} \ &= \ x_k - \alpha_{k+1} \frac{\nabla F(x_k)}{||\nabla F(x_k)||}.
\end{split}
```
Das adaptive Gradientenabstiegsverfahren in
{eq}`eq:gradient_descent_adaptive` lässt sich mit folgendem Algorithmus
umsetzen.

 \

````{prf:algorithm} Adaptives Gradientenabstiegsverfahren
:label: alg:gradient_descent_adaptive

**function** $[x^*, F(x^*)]=$`gradientDescentAdaptive`$(F,\nabla F, x_0, \alpha_0, \sigma, \epsilon$)  

\# Initialisierung  
$\alpha_k = \alpha_0$  
$x_k = x_0$  
$F(x_{k+1}) = +\infty$  

\# Iteriere bis gewünschte Genauigkeit erreicht ist  
**while** $|F(x_{k+1}) - F(x_k)| > \epsilon$ **do**  
    **while** $F(x_k) - \alpha_k \nabla F(x_k)/|| \nabla F(x_k)||) > F(x_k)$ **do**  
        \# Verringere Schrittweite um Faktor $\sigma$  
        $\alpha_k = \sigma\alpha_k$  
    **end while**  
    \# Update in Richtung des größten Gradientenabstiegs  
    $x_{k+1} = x_k - \alpha_k\nabla F(x_k)/||\nabla F(x_k)||$  
**end while**  

\# Ausgabe des letzten Punktes  
$x^* = x_k$  
$F(x^*) = F(x_k)$
````

````{prf:remark} Zurücksetzen der Schrittweite
In manchen Anwendungsfällen macht es Sinn die Schrittweite
$\alpha_{k+1}$ in {eq}`eq:gradient_descent_adaptive` in jedem Schritt
wieder auf den initialen Wert $\alpha_0 > 0$ zurückzusetzen, um eine
verbesserte Konvergenzgeschwindigkeit zu erhalten. Dies macht vor allem
dann SInn, wenn eine Auswertung der Zielfunktion $F$ und ihres
Gradienten $\nabla F$ numerisch günstig zu realisieren ist.

````
```{figure} ../atelier/img/gradient_descent.png
---
name: "fig:gradient_descent_adaptive"
---
Approximation des Minimierers einer Funktion $F$ in zwei Variablen mit
Hilfe des adaptiven Gradientenverfahrens
{eq}`eq:gradient_descent_adaptive`.
```

In {numref}`fig:gradient_descent_adaptive` ist ein typischer Verlauf des
adaptiven Gradientenverfahrens in Algorithmus
{prf:ref}`alg:gradient_descent_adaptive` bei der Minimierung einer
konvexen Zielfunktion $F \colon \mathbb{R}^2 \rightarrow \mathbb{R}$ zu
sehen. Man erkennt, dass die Schrittweiten immer kleiner werden, je
näher man sich dem lokalen Minimum $x_*$ nähert. Außerdem sieht man,
dass die Richtung des steilsten Gradientenabstiegs immer orthogonal zu
den Niveaulinien der zu minimierenden Funktion steht.

````{prf:remark} 
Das in diesem Abschnitt beschriebene Gradientenabstiegsverfahren mit
adaptiver Schrittweite $\alpha_k > 0$ ist ein gängiger Algorithmus zur
Minimierung einer Funktion $F$, wenn deren Ableitung $\nabla F$ bekannt
und numerisch günstig zu berechnen ist. Dennoch gibt es Situationen in
denen es ratsam ist alternative Optimierungsalgorithmen zu verwenden.
Zum Beispiel ist ein häufiges Problem des Gradientenverfahrens die
starke Verlangsamung in der Nähe eines Sattelpunktes, was zu sehr langen
Laufzeiten des Algorithmus führt. Außerdem passiert die Minimierung
einer Funktion $F$ mit Hilfe des Gradientenabstiegsverfahrens in der
Regel entlang eines Zickzack-Pfades (siehe
{numref}`fig:gradient_descent_adaptive`), welcher in den meisten Fällen
offensichtlich suboptimal ist. Aus diesen Gründen wollen wir uns in den
nächsten Abschnitten mit alternativen Minimierungsmethoden beschäftigen.

````

(ss:coordinate_descent)=
## Koordinatenabstiegsverfahren

Eine weitere Variante des in {ref}`ss:gradient_descent` behandelten
Gradientenabstiegsverfahrens ist das *Koordinatenabstiegsverfahren* (im
Englischen: *coordinate descent* (CD)).

Die grundlegende Idee des Koordinatenabstiegsverfahrens ist es in jedem
Schritt des Iterationsschemas eine *Koordinatenrichtung* auszuwählen und
einen Abstieg in diese Richtung durchzuführen. Damit lässt sich ein
möglicherweise kompliziertes multivariates Optimierungsproblem durch
eine Reihe von einfachen univariaten Optimierungsproblemen behandeln.
Die Auswahl der Koordinatenrichtung kann entweder mit Hilfe einer
Auswahlregel, z.B. mit einem *Rundlaufverfahren*, oder aber *zufällig*
geschehen.

Wir wollen im Folgenden den Fall einer **zufälligen Wahl der
Koordinatenrichtung** diskutieren. Für einen zufälligen Index
$j \in \lbrace 1,\ldots,n\rbrace$ und den entsprechenden zufälligen
Einheitsvektor $e_j \in \mathbb{R}^n$ lässt sich das
Koordinatenabstiegsverfahren schreiben als:
```{math}
:label: eq:coordinate_descent_stochastic
x_{k+1} \ = \ x_k - \alpha_k \langle \nabla F(x_k), e_j \rangle e_j \ = \ x_k - \alpha_k \frac{\partial F}{\partial x^i}(x_k) e_j.
```
Das bedeutet, dass man in jeder Iteration nur eine Koordinate des
aktuellen Parametervektors $x_k \in \Omega$ verändern muss. Die
Schrittweite $\alpha_k > 0$ in {eq}`eq:coordinate_descent_stochastic`
kann hierbei ähnlich wie in {ref}`ss:gradient_descent` *fest* oder
*adaptiv* gewählt werden. Das Koordinatenabstiegsverfahren in
{eq}`eq:coordinate_descent_stochastic` mit adaptiver Schrittweite lässt
sich mit folgendem Algorithmus umsetzen.

 \

````{prf:algorithm} Koordinatenabstiegsverfahren
label="alg:coordinate_descent_stochastic"}

**function** $[x^*, F(x^*)]=$`coordinateDescentStochastic`$(F,\nabla F, x_0, \alpha_0, \sigma, \epsilon$)  

\# Initialisierung  
$\alpha_k = \alpha_0$  
$x_k = x_0$  
$F(x_{k+1}) = +\infty$  

\# Iteriere bis gewünschte Genauigkeit erreicht ist  
**while** $|F(x_{k+1} - F(x_k)|>\epsilon$ **do**  
\
    \# Wähle zufällige Koordinatenrichtung  
    $i = \texttt{randomDraw}([1:n])$  
\
    \# Berechne Ableitung in Koordinatenrichtung  
    $p_k = \frac{\partial}{\partial x_k^i}F(x_k) \cdot e_i$  
\
    **while** $F(x_k - \alpha_k p_k) > F(x_k)$ **do**  
        \# Verringere Schrittweite um Faktor $\sigma$  
        $\alpha_k = \sigma\alpha_k$  
    **end while**  
\
    \# Update in Richtung des größten Gradientenabstiegs  
    $x_{k+1} = x_k - \alpha_k p_k$  
**end while**  
\
\# Ausgabe des letzten Punktes  
$x^* = x_k$  
$F(x^*) = F(x_k)$

````

Das Koordinatenabstiegsverfahren benötigt in der Regel deutlich mehr
Iterationen als das normale Gradientenabstiegsverfahren und beschreibt
häufig noch mehr einen Zickzack-Pfad bei der Minimierung. Dennoch bietet
das Verfahren Vorteile gegenüber dem Gradientenabstiegsverfahren gerade
in Optimierungsproblemen mit vielen Variablen (z.B. für das Training
eines künstlichen neuronalen Netzes), da jedes eindimensionale
Optimierungsproblem wesentlich leichter zu lösen ist als die Berechnung
des gesamten Gradienten in jedem Schritt.

````{prf:remark} Blockkoordinatenabstiegsverfahren
Um den Zufallseffekten und den damit verbundenen Zickzack Pfad bei der
Minimierung der Funktion $F$ durch das Koordinatenabstiegsverfahren
entgegen zu wirken, kann man einen Kompromiss zwischen der Verwendung
einer einzelnen Koordinatenrichtung und dem gesamten Gradienten
eingehen. Hierbei spricht man von den sogenannten
*Blockkoordinatenabtiegsverfahren* (im Englischen *block coordinate
descent* (BCD)). Hierbei wählt man zuerst die Größe
$s \in \lbrace 1,\ldots,n\rbrace$ der Koordinatenblöcke, d.h., die Größe
der Teilmenge der verwendeten Richtungsableitungen des Gradienten
$\nabla F$. Anschließend wird in jedem Schritt ein Block an
Koordinatenrichtungen deterministisch oder zufällig ausgewählt und in
dessen Richtung minimiert. Für eine zufällige Wahl der Koordinatenblöcke
ergibt sich somit:
```{math}
:label: eq:block_coordinate_descent
x_{k+1} \ = \ x_k - \alpha_k \sum_{i=1}^s \langle \nabla F(x_k), e_{\sigma(i)} \rangle e_{\sigma(i)}.
```
Hierbei ist
$\sigma \colon \lbrace1,\ldots,n\rbrace \rightarrow \ \lbrace1,\ldots,n\rbrace$
eine zufällige Permutation der Indizes $1,\ldots,n$. Es ist klar, dass
das Koordinatenabstiegsverfahren in
{eq}`eq:coordinate_descent_stochastic` und das normale
Gradientenabstiegsverfahren in {eq}`eq:gradient_descent_adaptive`
Spezialfälle des Blockkoordinatenabstiegsverfahrens in
{eq}`eq:block_coordinate_descent` für Blockgrößen $s = 1$ und $s=n$
sind.

````

(ss:stochastic_gradient_descent)=
## Stochastisches Gradientenabstiegsverfahren

Eine aktuell weit verbreitete Variante des Gradientenabstiegsverfahrens
in {eq}`eq:gradient_descent_adaptive` ist das *stochastische
Gradientenabstiegsverfahren* (im Englischen: *stochastic gradient
descent* (SGD)). Wie der Name schon verrät handelt es sich hierbei nicht
um einen deterministischen Algorithmus. Das bedeutet, dass man bei
mehrmaliger Anwendung des Verfahrens bei gleichbleibenden
Startbedingungen in der Regel unterschiedliche Ergebnisse in
unterschiedlichen Laufzeiten erhält. Was auf den ersten Blick wie ein
Nachteil wirkt, kann in manchen Fällen jedoch praktische Eigenschaften
mit sich bringen. So kann die Zufallsnatur des stochastischen
Gradientenverfahrens dazu führen, dass Sattelpunkte und ungewollte,
lokale Minima der Funktion durch die Folge der Punkte vermieden werden.
Das Verfahren findet aktuell vor allem beim Training von neuronalen
Netzen bei der sogenannten *Backpropagation* in verschiedenen
Variationen Anwendung, da man hierdurch dem bekannten Problem des
*Übertrainierens* des neuronalen Netzes entgegen wirken kann.

Beim stochastischen Gradientenverfahren geht man davon aus, dass sich
die zu minimierende Zielfunktion
$F \colon \Omega \rightarrow \mathbb{R}$ als eine Summe der folgenden
Gestalt schreiben lässt:
```{math}
:label: eq:objective_function_sum
F(x) \ = \ \sum_{i=1}^m F_i(x), \quad \text{ für alle } x \in \Omega.
```
Solche Zielfunktionen treten natürlicherweise in vielen
Problemstellungen auf, zum Beispiel bei Maximimum-Likelihood Ansätzen
oder der Methode der kleinsten Quadrate. Im Bereich des maschinellen
Lernens lässt sich der Trainingsfehler über alle Trainingsdaten in der
Regel als eine solche Summe schreiben. In diesem Fall lässt sich das
normale Gradientenabstiegsverfahren in
{eq}`eq:gradient_descent_adaptive` umschreiben zu:
```{math}
:label: eq:gradient_descent_sum
x_{k+1} \ = \ x_k - \alpha_k \frac{\nabla F(x_k)}{||\nabla F(x_k)||} \ = \ x_k - \alpha_k \frac{\sum_{i=1}^m\nabla F_i(x_k)}{||\sum_{i=1}^m\nabla F_i(x_k)||}, \quad \alpha_k > 0 .
```

Die Idee des stochastischen Gradientenverfahrens ist es nun einen
zufälligen Summanden aus {eq}`eq:objective_function_sum` zu wählen und
nur den Gradienten bezüglich dieses Summanden zu betrachten. Durch diese
starke Vereinfachung von {eq}`eq:gradient_descent_sum` führt man mit
einem zufällig ausgewählten Index $j \in \lbrace 1,\ldots,m\rbrace$ nun
einen Gradientenabstieg der Form
```{math}
:label: eq:stochastic_gradient_descent
x_{k+1} \ = \ x_k - \alpha_k \frac{\nabla F_j(x_k)}{||\nabla F_j(x_k)||}, \quad \alpha_k > 0
```
durch.

````{prf:remark} 
Ähnlich wie im Fall des Koordinatenabstiegsverfahrens in Kapitel
{ref}`ss:coordinate_descent`, gibt es auch beim stochastischen
Gradientenverfahren die Möglichkeit einen Kompromiss zwischen dem
normalen Gradientenabstieg in {eq}`eq:gradient_descent_control` und dem
auf einen Summanden beschränkten Gradientenabstieg
{eq}`eq:stochastic_gradient_descent` einzugehen. Indem man eine
zufällige Untermenge von fester Größe $s \in \lbrace 1,\ldots,m\rbrace$
von Summanden von $F$ auswählt, lässt sich das sogenannte *stochastische
Minibatch-Gradientenabstiegsverfahren* formulieren:
```{math}
x_{k+1} = x_k - \alpha_k \frac{\sum_{i=1}^s\nabla F_{\sigma(i)}(x_k)}{||\sum_{i=1}^s\nabla F_{\sigma(i)}(x_k)||}, \quad \alpha_k > 0 .
```
Hierbei ist
$\sigma \colon \lbrace 1,\ldots,m\rbrace \rightarrow \lbrace 1,\ldots,m\rbrace$
eine zufällige Permutation der Indizes $1,\ldots,m$.

````

(ss:newton)=
## Newton Verfahren

In diesem Abschnitt wollen wir uns das bereits bekannte Newton-Verfahren
in Erinnerung rufen und dieses geeignet zur Optimierung von
nichtlinearen Funktionen verallgemeinern. In {cite:p}`numerik1` haben
wir das *Newton Verfahren* zur Approximation von Nullstellen
nichtlinearer Gleichungssysteme hergeleitet. Wir haben zunächst die
Taylorapproximation einer nichtlinearen Nullstellengleichung
$F(x^*) = 0$ von der folgenden Form betrachtet
```{math}
0 \ = \ F(x^*) \ \approx \ F(x) + F'(x)(x^* - x),
```
wobei $F'$ die als regulär angenommene Jacobi-Matrix der
differenzierbaren Funktion $F \colon \R^n \rightarrow \R^n$ bezeichnet.
Hierauf basierend haben wir die folgende **Fixpunktfunktion** als
Approximation erster Ordnung angegeben:
```{math}
:label: eq:newton_fixpunkt
G(x) \ = \ x - (F'(x))^{-1} F(x), \quad \text{ für } F'(x) \text{ regulär}.
```
Hierbei haben wir die Fixpunktgleichung als erfüllt gesehen, wenn wir
ein $x^* \in \Omega$ gefunden haben, so dass für die Fixpunktgleichung
{eq}`eq:newton_fixpunkt` gilt $x^* = G(x^*)$. Unter dieser Beobachtung
haben wir das *Newton-Verfahren* als iteratives Schema zur Bestimmung
eines solchen Fixpunktes $x^* \in \Omega$ hergeleitet:
```{math}
x_{k+1} \ = \ x_k - (F'(x_k))^{-1} F(x_k),, \quad \text{ für } F'(x_k) \text{ regulär}.
```
Hierfür benötigten wir einen geeigneten Startwert $x_0 \in \Omega$ in
einer lokalen Umgebung $U \subset \Omega$ des Fixpunktes $x^* \in U$.

Der folgende Satz formuliert Bedingungen für die lokale Konvergenz des
Newton-Verfahrens.

````{prf:theorem} Lokale Konvergenz des Newton Verfahrens
Sei $F: \R^n \rightarrow \R^n$ in einer Umgebung von
$\overline{x} \in \R^n$ stetig differenzierbar und $\overline{x}$ sei
eine Nullstelle von $F$ mit $F(\overline{x}) = 0$. Sei außerdem die
Jacobi-Matrix $F'$ lokal Lipschitz-stetig und $F'(\overline{x})$ regulär
in der Nullstelle.

Dann existiert eine lokale Umgebung $B_R(\overline{x})$, so dass das
Newton-Verfahren für jeden Startwert $x_0 \in B_R(\overline{x})$ gegen
die Nullstelle $\overline{x}$ konvergiert, d.h. es gilt
$\lim_{x\rightarrow \infty} x_k =  \overline{x}$.

````

````{prf:proof} 
Siehe {cite:p}`numerik1`. ◻

````

Anstatt nun eine Nullstelle der Funktion $F$ zu suchen, wollen wir das
Newton-Verfahren nutzen, um eine Nullstelle des Gradienten $\nabla F$
(d.h einen stationären Punkt von $F$) zu approximieren und damit die
notwendigen Optimalitätsbedingungen in {prf:ref}`thm:minimum_notwendig`
zu erfüllen. Im Folgenden sei $\Omega \subset \mathbb{R}^n$ ein offenes,
zusammenhängendes Gebiet und $F \colon \Omega \rightarrow \mathbb{R}$
eine differenzierbare, reellwertige Funktion. Wir betrachten wieder die
Taylorapproximation der Funktion $F$ in eine Abstiegsrichtung
$x_k + p \in \Omega$ des allgemeinen Iterationsschemas
{eq}`eq:abstiegsverfahren`, aber berücksichtigen diesmal auch Terme von
zweiter Ordnung:
```{math}
:label: eq:modellfunktion
F(x_k + p) \ \approx \ F(x_k) + \langle p, \nabla F(x_k) \rangle + \frac{1}{2} \langle p, \nabla^2F(x_k)p \rangle \ \eqqcolon \ m_k(p).
```

Unter gewissen Bedingungen an die Hessematrix $\nabla^2 F(x_k)$, lässt
sich ein eindeutiges Minimum der Modellfunktion $m_k(p)$ in
{eq}`eq:modellfunktion` bestimmen, wie folgendes Theorem besagt.

````{prf:theorem} Newton-Abstiegsrichtung
Sei $\Omega \subset \R^n$ ein offenes, zusammenhängendes Gebiet. Sei
außerdem $F \colon \Omega \rightarrow \R$ eine in einer lokalen Umgebung
eines Punktes $x_k \in \Omega$ zweimal stetig differenzierbare
Zielfunktion, deren Hessematrix $\nabla^2 F(x_k)$ im Punkt $x_k$ positiv
definit ist.

Dann ist der **Newton-Abstiegsrichtung** benannte Vektor
$p_k^N \in \R^n$ mit
```{math}
:label: eq:newton_richtung
p_k^N \ = \ -(\nabla^2 F(x_k))^{-1} \nabla F(x_k)
```
das eindeutige Minimum der Modellfunktion $m_k(p)$ in
{eq}`eq:modellfunktion`.

````

````{prf:proof} 
In den Übungsaufgaben zu zeigen. ◻

````

Mit der Newton-Abstiegsrichtung in {eq}`eq:newton_richtung` lässt sich
ein iteratives Abstiegsverfahren für einen initialen Punkt
$x_0 \in \Omega$, welcher geeignet in der Nähe des stationären Punktes
$x^* \in \Omega$ gewählt wird, wie folgt konstruieren:
```{math}
:label: eq:newton_abstieg
x_{k+1} \ = \ x_k + p_k^N \ = \ x_k - (\nabla^2 F(x_k))^{-1} \nabla F(x_k).
```
Damit das Newton-Abstiegsverfahren in {eq}`eq:newton_abstieg` überhaupt
sinnvoll ist, müssen wir fordern, dass die Hessematrix in jedem Punkt
$x_k \in \Omega$ der Iterationsfolge *regulär* und somit invertierbar
ist. Um sicher zu gehen, dass es sich tatsächlich um eine
Abstiegsrichtung handelt müssen wir fordern, dass die Hessematrix
$\nabla^2 F(x_k)$ nicht nur invertierbar für alle $x_k \in \Omega$ der
Iterationsfolge ist, sondern auch *positiv definit* in jedem Punkt $x_k$
ist. Denn dann ergibt eine Taylorapproximation zweiter Ordnung die
folgende Abschätzung:
```{math}
\begin{split}
F(x_{k+1}) \ &= \ F(x_k + p_k^N) \ \approx \ F(x_k) + \langle p_k^N, \nabla F(x_k) \rangle + \frac{1}{2} \langle p_k^N, \nabla^2F(x_k) p_k^N \rangle \\
\ &= \ F(x_k) - \langle p_k^N, \nabla^2 F(x_k) p_k^N \rangle + \frac{1}{2} \langle p_k^N, \nabla^2F(x_k) p_k^N \rangle \\
\ &= \ F(x_k) - \frac{1}{2} \underbrace{\langle p_k^N, \nabla^2 F(x_k) p_k^N}_{>~0} \rangle.
\end{split}
```
Wir sehen also, dass wir einen echten Abstieg der Funktionswerte
erhalten, wenn die Hessematrix $\nabla^2 F(x_k)$ positiv definit ist für
alle $x_k \in \Omega$ der Iterationsfolge. Sollte die Hessematrix nicht
positiv definit in einem Punkt $x_k$ der Iterationsfolge sein, so muss
zumindest eine Abnahme der Funktionswerte vorliegen, d.h., es muss für
die Newton-Abstiegsrichtung gelten:
```{math}
\langle (\nabla^2F(x_k))^{-1} \nabla F(x_k), \nabla F(x_k) \rangle \ > \ 0.
```
Sollte dies nicht der Fall sein, so existieren Methoden um dennoch einen
Abstieg zu erzwingen, siehe zum Beispiel {cite:p}`nocedal_1999`. Auf
diese werden wir jedoch im weiteren Verlauf der Vorlesung nicht näher
eingehen.

````{prf:remark} Schrittweite und Konvergenz
Das Newton-Abstiegsverfahren in {eq}`eq:newton_abstieg` ist ein
Abstiegsverfahren der Art {eq}`eq:abstiegsverfahren` dessen Schrittweite
$\alpha_k > 0$ implizit durch die lokale Krümmung und die Ableitung der
Funktion $F$ bestimmt ist. In diesem Fall können wir $\alpha_k \equiv 1$
für alle $k \in \mathbb{N}$ setzen. Das Newton-Abstiegsverfahren
konvergiert in der Regel *quadratisch* gegen einen stationären Punkt
$x^* \in \Omega$ mit $\nabla F(x^*)=0$, d.h. man erreicht sehr schnell
eine hohe Genauigkeit bei der Approximation von $x^*$.

````

(ss:quasi-newton)=
## Quasi-Newton Verfahren

Im {ref}`ss:newton` haben wir das Newton Verfahren zur iterativen
Approximation eines stationären Punktes $x^* \in \Omega$ einer Funktion
$F$ mit $\nabla F(x^*) = 0$ hergeleitet. Hierbei haben wir im Gegensatz
zu den vorherigen numerischen Verfahren auch Ableitungen höherer Ordnung
hinzugezogen. Dies führt in der Regel zu einem verbesserten
Konvergenzverhalten im Vergleich zu den Verfahren, die nur die lokale
Ableitung $\nabla F$ der Zielfunktion $F$ verwenden.

Dennoch ist das Newton Verfahren aus numerischer Sicht noch nicht ideal,
da es einige Probleme mit sich bringt. Zuerst mussten wir fordern, dass
die Hessematrix $\nabla^2 F(x_k)$ in jedem Punkt des
Iterationsverfahrens positiv definit ist, da ansonsten kein Abstieg der
Funktionswerte garantiert werden kann. Zweitens muss für die Berechnung
der Newton-Richtung in {eq}`eq:newton_richtung` zuerst die Hessematrix
bestimmt und anschließend invertiert werden. Dies ist aus
Effizienzgründen unerwünscht, da die Inversion einer $n \times n$-Matrix
bekanntlich in $\mathcal{O}(n^3)$ Rechenoperationen liegt (siehe
{cite:p}`numerik1`. Da die Bestimmung und die Inversion der Hessematrix
in jedem Iterationsschritt passieren müssen, ist das Newton Verfahren
nur eingeschränkt empfehlenswert für die numerische Optimierung.

Eine naheliegende Idee ist es nun die echte Hessematrix in jedem
Iterationsschritt durch eine geeignete Matrix zu approximieren, so dass
der numerische Aufwand geringer wird, d.h., wir suchen nach einer Matrix
$B_k \in \R^{n \times n}$
```{math}
:label: eq:matrix_bk
B_k \approx \nabla^2 F(x_k).
```
Damit können wir die Modellfunktion $m_k(p)$ in {eq}`eq:modellfunktion`
schreiben als:
```{math}
m_k(p) \ = \ F(x_k) + \langle p, \nabla F(x_k)\rangle + \frac{1}{2} \langle p, B_k p \rangle,
```
das heißt, wir approximieren die Zielfunktion $F$ im $k$-ten
Iterationsschritt entlang der Richtung $p \in \mathbb{R}^n$ lokal durch
eine quadratische Funktion. Für sehr kleine Schrittweiten können wir
davon ausgehen, dass der Fehler dieser Approximation gering ist, da wir
davon ausgehen, dass $F$ stetig differenzierbar in einer lokalen
Umgebung $U \subset \Omega$ des stationären Punktes $x^* \in \Omega$ ist
und für $p = 0$ die Approximation exakt ist, da gilt
```{math}
m_k(0) \ = \ F(x_k).
```
Wenn wir fordern, dass $B_k$ in {eq}`eq:matrix_bk` eine positiv definite
Matrix ist, so lässt sich ein Abstiegsschritt des Iterationsverfahrens
{eq}`eq:abstiegsverfahren` analog zur Herleitung des Newton
Abstiegsverfahrens in {ref}`ss:newton` angeben als:
```{math}
:label: eq:quasi-newton-abstieg
x_{k+1} \ = \ x_k + \alpha_k p_k, \qquad p_k \ = \ -B_k^{-1} \nabla F(x_k).
```
Die sogenannten *Quasi-Newton Verfahren* verfolgen diesen Ansatz.

````{prf:remark} Konvergenzgeschwindigkeit Quasi-Newton Verfahren
Durch die Approximation der echten Hessematrix verlieren Quasi-Newton
Verfahren an Genauigkeit, wodurch ihre Konvergenzgeschwindigkeit
*superlinear* anstatt *quadratisch* ist. Dafür gewinnen sie zusätzliche
Geschwindigkeit durch die Vermeidung der Bestimmung und Inversion von
$\nabla^2 F(x_k)$. Der Vorteil der Quasi-Newton Methoden ist es, dass
man nur den Gradienten $\nabla F$ für einen Schritt des numerischen
Optimierungsverfahrens benötigt und keine expliziten Informationen über
die zweiten Ableitungen. Dadurch werden sie in bestimmten Problemen
sogar effizienter bei der Approximation eines stationären Punktes als
das Newton Abstiegsverfahren in {ref}`ss:newton`.

````

(sekantengleichung-und-krümmungsbedingung)=
### Sekantengleichung und Krümmungsbedingung

Die entscheidende Frage bei der Konstruktion eines Quasi-Newton
Abstiegsverfahrens der Form {eq}`eq:quasi-newton-abstieg` ist es, wie
die positiv definite Matrix $B_k$ in jedem Schritt möglichst effizient
bestimmt werden kann. Anstatt die Näherung $B_k$ der Hessematrix
$\nabla^2 F(x_k)$ in jedem Schritt von Grund auf neu zu berechnen, wäre
es wünschenswert ein initiales $B_0$ zu bestimmen, das in jedem Schritt
des Iterationsverfahrens nur aktualisiert werden muss. Hierbei ist es
möglich die durch den Iterationsschritt erhaltenen Informationen über
den Gradienten $\nabla F$ zu Hilfe zu nehmen.

Wir nehmen an, wir haben bereits einen Abstiegsschritt durchgeführt und
so einen neuen Punkt $x_{k+1} = x_k + \alpha_kp$ erhalten. Unsere
quadratische Approximation in diesem neuen Punkt für eine neue Richtung
$p \in \mathbb{R}^n$ sieht dementsprechend wie folgt aus:
```{math}
:label: eq:quasi-newton_modellfunktion
m_{k+1}(p) \ = \ F(x_{k+1}) + \langle \nabla F(x_{k+1}), p \rangle + \frac{1}{2}\langle p, B_{k+1}p \rangle.
```
Es ist leicht einzusehen, dass die Modellfunktion $m_{k+1}$ im Punkt
$x_{k+1} \in \R^n$ zentriert ist und für $p=0$ mit dem Funktionswert der
Zielfunktion im Punkt $x_{k+1}$ übereinstimmt, d.h., es gilt
$m_{k+1}(0) = F(x_{k+1})$.

Da wir uns bei der Wahl der positiv definiten Matrix $B_{k+1}$ noch
nicht festgelegt haben, wird durch {eq}`eq:quasi-newton_modellfunktion`
eine Funktionenschar beschrieben. Eine Forderung, die man nun die
Modellfunktion $m_{k+1}$ stellen kann, um eine sinnvolle Matrix
$B_{k+1}$ zu bestimmen, ist, dass ihre Ableitung $\nabla m_{k+1}$ mit
der Ableitung der Zielfunktion $F$ in den letzten beiden Punkten $x_k$
und $x_{k+1}$ übereinstimmt. Dies bedeutet, dass man die Matrix
$B_{k+1}$ versucht so zu bestimmen, dass die Modellfunktion $m_{k+1}$
die Krümmung der Zielfunktion $F$ gut approximiert. Da bereits gilt
```{math}
\nabla m_{k+1}(0) \ = \ \nabla F(x_{k+1}),
```
ist eine der beiden Forderungen automatisch erfüllt. Für die zweite
Forderung können wir nutzen, dass $x_k = x_{k+1} - \alpha_k p_k$ gilt
und wir erhalten somit:
```{math}
:label: eq:quasi-newton_forderung
\nabla F(x_k) \ \overset{!}{=} \ \nabla m_{k+1}(-\alpha_kp_k) \ = \ \nabla F(x_{k+1}) - B_{k+1}\alpha_kp_k.
```
Durch Umstellen von {eq}`eq:quasi-newton_forderung` erhalten wir die
Bedingung
```{math}
\nabla F(x_{k+1}) - \nabla F(x_k)  \ \overset{!}{=} \ B_{k+1}\alpha_k p_k \ = \ B_{k+1}(x_{k+1} - x_k).
```
Eine vernünftige Wahl der Matrix $B_{k+1}$ in {eq}`eq:matrix_bk` sollte
diese Eigenschaft, auch bekannt als **Sekantengleichung**, versuchen zu
imitieren. Im eindimensionalen Fall mit
$F \colon \Omega \subset \mathbb{R} \rightarrow \mathbb{R}$ bedeutet die
Sekantengleichung nichts anderes, als dass der Faktor $B_{k+1}$ eine
Approximation des zweiten Ableitung von $F$ im Sinne eines
Differenzenquotienten ist, d.h., im Fall $n=1$ soll gelten:
```{math}
B_{k+1} \ \overset{!}{=} \ \frac{F'(x_{k+1})-F'(x_k)}{x_{k+1} - x_k}.
```

Für unser allgemeines Quasi-Newton Verfahren in
{eq}`eq:quasi-newton-abstieg` suchen wir also einen Weg die bereits
bekannte Approximation der Hessematrix $B_k \approx \nabla^2F(x_k)$ zu
einer Matrix $B_{k+1}$ zu aktualisieren, so dass der folgende
Zusammenhang für den nächsten Punkt $x_{k+1} \in \Omega$ erfüllt wird:
```{math}
:label: eq:sekantengleichung
B_{k+1} s_k \ = \ y_k,
```
wobei
```{math}
s_k \ = \ x_{k+1} - x_k, \qquad y_k = \nabla F(x_{k+1}) - \nabla F(x_k).
```
Es wird klar, dass diese Forderung alleine nicht genügt für die
Konstruktion eines Abstiegsverfahrens, da die Sekantengleichung in
{eq}`eq:sekantengleichung` für $n > 1$ unterbestimmt ist, d.h., dass es
mehr unbekannte Einträge der Matrix $B_{k+1} \in \R^{n \times n}$ gibt
als durch die $n$ Gleichungen festgelegt werden. Daher versuchen wir im
Folgenden weitere Forderungen an die Matrix $B_{k+1}$ zu stellen.

Um die positive Definitheit der Matrix $B_{k+1}$ in Schrittrichtung
$x_{k+1} - x_k = \alpha p_k \in \R^n$ zu gewährleisten müssen wir
fordern, dass die Vektoren $y_k$ und $s_k$ die sogenannte
**Krümmungsbedingung** erfüllen:
```{math}
:label: eq:kruemmungsbedingung
\langle s_k, y_k \rangle \ > \ 0.
```
Dies ist eine hinreichende Bedingung für die positive Definitheit von
$B_{k+1}$ bezüglich der Richtung $\alpha_k p_k$, da wir einfach die
Sekantengleichung {eq}`eq:sekantengleichung` von links mit dem Vektor
$s_k^T$ multiplizieren können und so erhalten wir mit der Forderung
{eq}`eq:kruemmungsbedingung` schon:
```{math}
\langle s_k, B_{k+1}s_k \rangle \ = \ \langle s_k, y_k \rangle \ > \ 0.
```

````{prf:remark} Krümmungsbedingung und Konvexität
Falls die Zielunktion $F$ strikt konvex ist, so ist die
Krümmungsbedingung {eq}`eq:kruemmungsbedingung` für alle Punktepaare
$x_k, x_{k+1} \in \Omega$ erfüllt und die Matrix $B_{k+1}$ wird damit
positiv definit. Für nichtkonvexe Funktionen hingegen muss man die
Krümmungsbedingung explizit forcieren, um ein Abstiegsverfahren zu
erhalten.

````

Falls die Krümmungsbedingung {eq}`eq:kruemmungsbedingung` erfüllt ist,
so existiert mindestens eine Lösung $B_{k+1}$ der Sekantengleichung
{eq}`eq:sekantengleichung`. Man sieht ein, dass es in der Tat sogar
unendlich viele Lösungen $B_{k+1}$ gibt, da eine symmetrische
$n \times n$ Matrix $n(n+1)/2$ Freiheitsgrade besitzt und die
Sekantengleichung {eq}`eq:sekantengleichung` nur $n$ Bedingungen an
$B_{k+1}$ stellt. Zusätzlich erhält man $n$ Bedingungen an $B_{k+1}$
durch die Forderung von positiver Definitheit, da alle $n$ Hauptminoren
von $B_{k+1}$ positiv sein müssen. Dies reicht jedoch nicht für die
eindeutige Bestimmung der Matrix $B_{k+1}$. Hierfür müssen wir
zusätzlich fordern, dass die Matrix $B_{k+1}$ diejenige Matrix unter
allen möglichen Lösungen ist, die der vorherigen Matrix $B_k$ am
nächsten bezüglich eines geeigneten Maßes ist. Das heißt wir suchen eine
Lösung des folgenden Optimierungsproblems:
```{math}
:label: eq:quasi-newton_optimization-problem
\begin{split}
\min_{B \in \R^{n\times n}} || B - B_k ||, \quad \text{ unter den Nebenbedingungen: } \\
B \ = \ B^T, \qquad B s_k \ = \ y_k, \qquad \langle p, Bp \rangle > 0, \forall p \in \mathbb{R}^n / \lbrace 0 \rbrace,
\end{split}
```
wobei $s_k$ und $y_k$ definiert sind wie in der Sekantengleichung
{eq}`eq:sekantengleichung`. Man beachte, dass man eine unterschiedliche
Lösung $B_{k+1}$ des Optimierungsproblems
{eq}`eq:quasi-newton_optimization-problem` in Abhängigkeit der gewählten
Matrixnorm erhält und somit auch ein unterschiedliches Quasi-Newton
Verfahren herleiten kann.

(das-davidon-fletcher-powell-verfahren)=
### Das Davidon-Fletcher-Powell Verfahren

Im ursprünglich im Jahr $1959$ von Davidon vorgeschlagenen Verfahren
{cite:p}`davidon_1959` (das im Übrigen bei der Erstbegutachtung
abgelehnt wurde) wählt man für die Norm im Optimierungsproblem
{eq}`eq:quasi-newton_optimization-problem` eine gewichtete Frobeniusnorm
der Form
```{math}
||A||_W \ \coloneqq \ || W^\frac{1}{2}A W^{-\frac{1}{2}} ||_F.
```
Die Gewichtungsmatrix $W$ dient dazu, dass das implizierte Quasi-Newton
Verfahren zur Approximation eines stationären Punktes $x^* \in \Omega$
skalierungs-invariant wird. Hierzu wählt man eine beliebige Matrix für
die die Relation $W y_k = s_k$ gilt, d.h., eine Matrix $W$, die sich wie
die Inverse der Matrix $B$ in {eq}`eq:quasi-newton_optimization-problem`
verhält. Ein konkretes Beispiel für solch eine Gewichtungsmatrix wäre
$W \coloneqq G_k^{-1}$, wobei $G_k$ die *durchschnittliche Hessematrix*
von $F$ entlang des letzten Abstiegsschritts von
$x_k \rightarrow x_{k+1}$ ist mit
```{math}
G_k \ \coloneqq \ \int_0^1 \nabla^2 F(x_k + t \alpha_k p_k)~\mathrm{d}t.
```
Mit der konkreten Wahl dieser Gewichtungsmatrix $W = G_k^{-1}$ wird die
gewichtete Frobeniusnorm dimensionslos und man erhält eine eindeutige
Lösung des Optimierungsproblems
{eq}`eq:quasi-newton_optimization-problem` wie folgt:
```{math}
:label: eq:DFP
B_{k+1} \ = \ (I-\gamma_k y_ks_k^T) B_k (I - \gamma_k s_k y_k^T) + \gamma_k y_k y_k^T, \quad \text{ mit } \gamma_k \ \coloneqq \ \frac{1}{\langle y_k, s_k \rangle}.
```
Die Gleichung {eq}`eq:DFP` wird auch **DFP-Schritt** genannt, da sie
zuerst von Davidon vorgeschlagen und später von Fletcher und Powell
untersucht und verbreitet wurde.

Obwohl wir die explizite Berechnung der Hessematrix $\nabla^2 F(x_k)$
vermieden haben und die Aktualisierung der Matrix $B_k$ zu $B_{k+1}$
lediglich auf den Gradienteninformationen von $F$ basiert ist der
numerische Aufwand bei direkter Verwendung von $B_{k+1}$ in {eq}`eq:DFP`
noch zu hoch. Das liegt daran, dass wir für einen Schritt des
Quasi-Newton Verfahrens in {eq}`eq:quasi-newton-abstieg` die Inverse der
Matrix $B_k$ benötigen und die Inversion einen numerischen Aufwand von
$\mathcal{O}(n^3)$ besitzt. Glücklicherweise gibt es einen Trick, wie
wir die Inverse von $B_k$ in jedem Schritt des Iterationsverfahrens
numerisch günstig erhalten können. Sei $H_k \coloneqq B_k^{-1}$, dann
können wir die sogenannte Sherman-Morrison-Woodbury Formel (siehe
{cite:p}`SMW_formel`) auf Gleichung {eq}`eq:DFP` anwenden um die neue
Inverse $H_{k+1}$ durch eine Aktualisierung der Matrix $H_k$ zu
berechnen:
```{math}
:label: eq:smw-formel
H_{k+1} \ = \ H_k - \frac{H_k y_k y_k^T H_k}{\langle y_k, H_k y_k\rangle} + \frac{s_k s_k^T}{\langle y_k, s_k \rangle} .
```
Wie man einssieht liegt der numerische Rechenaufwand für das Update von
$H_{k+1}$ in {eq}`eq:smw-formel` in $\mathcal{O}(n^2)$. Es fällt
außerdem auf, dass $H_k$ nur durch die Addition zweier Matrizen mit Rang
$1$ verändert wird, also insgesamt eine Änderung von höchstes Rang $2$
erfährt. Das passt gut zu der Forderung, dass wir erwarten, dass sich
die Approximation der Hessematrix $\nabla^2 F$ in einer lokalen Umgebung
nur wenig ändert.

(das-broydenfletchergoldfarbshanno-verfahren)=
### Das Broyden--Fletcher--Goldfarb--Shanno Verfahren

Das Davidon--Fletcher--Powell Verfahren wurde trotz seiner Effektivität
bald schon durch ein Verfahren abgelöst, das noch besser war und bis
heute zu den effizientesten Quasi-Newton Verfahren gehört: das
Broyden--Fletcher--Goldfard--Shanno (BFGS) Verfahren in
{cite:p}`broyden_1970`. Die Idee des BFGS Verfahrens leitet sich
unmittelbar aus der Idee des DFP Verfahren ab. Anstatt das
Optimierungsproblem {eq}`eq:quasi-newton_optimization-problem` mit
bestimmten Bedingungen an die Approximation $B_{k+1}$ der Hessematrix
$\nabla^2 F(x_k)$ zu stellen, versucht man direkt die Inverse der
Hessematrix $(\nabla^2 F(x_k))^{-1}$ geeignet zu approximieren. Hierfür
nehmen wir an, dass wir eine Matrix $H_{k+1}$ als geringfügige
Aktualisierung einer bereits vorher bestimmten Matrix $H_k$ suchen, die
gleichzeitig symmetrisch und positiv definit ist und zusätzlich die
Sekantenbedingung in umgeschriebener Form erfüllt:
```{math}
H_{k+1}y_k \ = \ s_k.
```
Hierzu formuliert man ein analoges Optimierungsproblem zu
{eq}`eq:quasi-newton_optimization-problem` von der Form:
```{math}
:label: eq:bfgs_optimierungsproblem
\begin{split}
\min_{H} || H - H_k ||, \quad \text{ unter den Nebenbedingungen: } \\
H \ = \ H^T, \qquad H y_k \ = \ s_k, \qquad \langle p, Hp \rangle > 0, \forall p \in \mathbb{R}^n / \lbrace 0 \rbrace.
\end{split}
```
Unter der Verwendung der gewichteten Frobeniusnorm und einer beliebigen
Gewichtsfunktion, die die Sekantengleichung $Ws_k = y_k$ erfüllt, erhält
man wiederum die eindeutige Lösung des Minimierungsproblems
{eq}`eq:bfgs_optimierungsproblem` als:
```{math}
:label: eq:bfgs
H_{k+1} \ = \ (I - \rho_k s_k y_k^T) H_k (I - \rho_k y_k s_k^T) + \rho_k s_k s_k^T, \quad \text{ mit } \rho_k \coloneqq \frac{1}{\langle y_k, s_k \rangle}.
```
Das Update der Matrix $H_k$ in {eq}`eq:bfgs` kann numerisch in
$\mathcal{O}(n^2)$ durchgeführt werden, was man schnell einsieht, wenn
man das Produkt ausschreibt:
```{math}
H_{k+1} \ = \ H_k - H_k\rho_k y_k s_k^T - \rho_k s_k y_k^T H_k + \rho_k s_k y_k^T H_k \rho_k y_k s_k^T + \rho_k s_k s_k^T.
```
In dieser Schreibweise sieht man gut, dass man lediglich Skalarprodukte
in $\mathcal{O}(n)$, Matrix-Vektor Multiplikationen in
$\mathcal{O}(n^2)$ und dyadische Produkte in $\mathcal{O}(n^2)$
berechnen muss. Im Gegensatz hierzu würde eine naive Implementierung des
BFGS-Updates in {eq}`eq:bfgs` zu einem numerischen Aufwand von
$\mathcal{O}(n^3)$ führen.

Abschließend bleibt die Frage was eine gute Initialisierung der Matrix
$H_0$ ist. Idealerweise hat man bereits Informationen über die Inverse
der Hessematrix $(\nabla^2 F(x_0))^{-1}$ im Initialisierungspunkt
$x_0 \in \Omega$, zum Beispiel durch eine numerische Approximation
mittels finiter Differenzen (später in der Vorlesung!). Andererseits
erwarten wir, dass die Aktualisierung von $H_k$ im $k$-ten Schritt des
Iterationsverfahrens {eq}`eq:bfgs` zu $H_{k+1}$ die aktuellen
Informationen über den Verlauf der Gradienten $\nabla F(x_k)$ und
$\nabla F(x_{k+1})$ berücksichtigt. Darum ist eine häufige Wahl von
$H_0$ die Initialisierung als Einheitsmatrix $I_n$ oder ein Vielfaches
der Einheitsmatrix, wobei die Vorfaktoren der Diagonaleinträge
entsprechend der Skalierung der Variablen gewählt werden.

