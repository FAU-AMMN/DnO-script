(wahl-der-schrittweite)=
Wahl der Schrittweite
===

Wir haben in den vorausgegangenen Abschnitten bereits verschiedene
Abstiegsverfahren der Form
```{math}
x_{k+1} \ = \ x_k + \alpha_k p_k, \qquad \alpha_k > 0, \quad p_k \in \R^n,
```
wie in {eq}`eq:abstiegsverfahren` kennen gelernt, die basierend auf
verschiedenen Annahmen unterschiedliche Abstiegsrichtungen
$p_k \in \R^n$ realisiert haben. Bis auf eine heuristische Wahl von
adaptiven Schrittweiten in {eq}`eq:gradient_descent_adaptive` und der
Wahl einer optimalen Schrittweite im Fall von strikt konvexen,
quadratischen Zielfunktionen in {eq}`eq:optimal_step-size_conjugated`
haben wir uns bisher noch nicht weiter mit der Frage einer geeigneten
Wahl der Schrittweiten $\alpha_k > 0$ beschäftigt. Dies wollen wir im
Folgenden nachholen.

Wir werden ab jetzt immer annehmen, dass $p_k \in \R^n$ eine fest
gewählte Abstiegsrichtung ist, d.h. für die Richtungsableitung in
Richtung $p_k$ in jedem Iterationsschritt $k=0,1,\ldots$ gilt
```{math}
\langle \nabla F(x_k), p_k \rangle \, < \, 0.
```
Nehmen wir an, dass die Zielfunktion $F$ stetig differenzierbar ist und
Ist die gewählte $\alpha_k > 0$ sehr klein, so bleiben wir sicher in
einer lokalen Umgebung von $x_k$ in der die folgende
Taylor-Approximation erster Ordnung gilt:
```{math}
F(x_k +\alpha_k p_k) \ \approx \ F(x_k) + \alpha_k \langle \nabla F(x_k), p_k \rangle  \: < \: F(x_k).
```
Allerdings würde das Iterationsverfahren dann in der Regel sehr langsam
konvergieren. Ist andererseits $\alpha_k$ zu groß, dann ist die
Abstiegsbedingung nicht mehr garantiert. Es könnte zum Beispiel
passieren, dass man beim Iterationsschritt
$x_{k+1} = x_k + \alpha_k p_k$ zu weit über ein lokales Minimum springt.
Deshalb benötigen wir intuitiv zwei Bedingungen an die Schrittweite
$\alpha_k > 0$, die zu kleine und zu große Schritte verhindern sollen.

Zunächst wollen wir analog zu {eq}`eq:optimal_step-size_conjugated` eine
theoretische Möglichkeit der optimalen Wahl der Schrittweite
$\alpha_k > 0$ untersuchen, nämlich jene Schrittweite, die zum
größtmöglichen Abstieg führt:
```{math}
\alpha_k  \ \coloneqq \ \text{arg}\min_{\alpha \in \R^+} F(x_k + \alpha p_k).
```
Um eine optimale Schrittweite $\alpha_k$ ausrechnen zu können, müssen
wir beliebige eindimensionale Probleme analytisch lösen können, was
jedoch im Allgemeinen schwierig ist. Deshalb fordern wir nicht die
analytische Optimalität der Schrittweite $\alpha_k >0$, sondern
versuchen lediglich Bedingungen zu finden, die wir numerisch leicht
überprüfen können und für die wir Konvergenz des Abstiegverfahrens
garantieren können.

Die Idee hierbei ist es für eine vorgegebene Schrittweite $\alpha > 0$
eine Linearisierung des Problems zu betrachten und den linearisierten
Abstieg mit dem echten Abstieg zu vergleichen. Hierzu machen wir
zunächst folgende Definitionen.

````{prf:definition} Erwarteter und tatsächlicher Abstieg
:label: def:erwarteter_abstieg
Sei $F \colon \Omega \rightarrow \R$ eine stetig differenzierbare
Funktion für die wir das Abstiegsverfahren $x_{k
+1} = x_k + \alpha p_k$ für $x_k \in \Omega$,
$p_k \in \R^n \setminus \lbrace 0\rbrace$ und $\alpha > 0$ betrachten.
Wir definieren basierend auf der Taylorapproximation erster Ordnung den
**erwarteten Abstieg** im Punkt $x_k$ in Richtung $p_k$ mit Schrittweite
$\alpha$ als
```{math}
:label: eq:abstieg_erwartet
E_k(\alpha) \ \coloneqq \ F(x_k) + \alpha  \langle \nabla F(x_k), p_k \rangle - F(x_k) \ = \ \alpha \langle \nabla F(x_k), p_k \rangle.
```
Darüber hinaus definieren wir für die gleichen Größen den
**tatsächlichen Abstieg** als
```{math}
:label: eq:abstieg_tatsaechlich
D_k(\alpha) \ \coloneqq \ F(x_k  + \alpha p_k) - F(x_k) .
```

````

Unsere beiden Bedingungen an eine geeignete Schrittweite $\alpha > 0$
können wir nun über die Abweichung des erwarteten und tatsächlichen
Abstiegs $D_k(\alpha)$ und $E_k(\alpha)$ in
{prf:ref}`def:erwarteter_abstieg` formulieren.

````{prf:definition} Armijo-Goldstein Bedingungen
:label: def:armijo-goldstein
Sei ein allgemeines Abstiegsverfahren der Form
$x_{k+1} = x_k + \alpha_k p_k$ für eine vorgegebene Abstiegsrichtung
$p_k \in \R^n$ gegeben und seien $c_1, c_2 \in \R^+$ Konstanten mit
$0 < c_1 < c_2 < 1$.

Dann formuliert man die sogenannten **Armijo-Goldstein Bedingungen** für
die Wahl einer geeigneten Schrittweite $\alpha_k > 0$ des
Abstiegsverfahrens als
```{math}
:label: eq:armijo-goldstein
c_1  E_k(\alpha_k) \: > \: D_k(\alpha_k) \: > \: c_2 E_k(\alpha_k),
```
wobei $E_k$ und $D_k$ den erwarteten und tatsächlichen Abstieg aus
{eq}`eq:abstieg_erwartet` und {eq}`eq:abstieg_tatsaechlich` definieren.
Da wir den erwarteten Abstieg für kleine Schrittweiten $\alpha_k > 0$
als negativ annehmen, d.h., es gilt $E_k(\alpha_k) < 0$, ist
{eq}`eq:armijo-goldstein` sinnvoll definiert.

````

Die erste Bedingung auf der linken Seite der Armijo-Goldstein
Bedingungen in {eq}`eq:armijo-goldstein` garantiert, dass zumindest ein
gewisser Teil des Abstiegs erreicht wird. Die zweite Bedingung auf der
rechten Seite verhindert, dass wir uns zu stark dem Fall $\alpha_k =0$
annähern, in dem die rechte Seite zu einer Gleichheit mit Konstante
gleich eins wird. Eine typische Wahl der Parameter in
{prf:ref}`def:armijo-goldstein` ist $c_1 =0.1$ und $c_2 = 0.9$.

In der Praxis lassen sich die Armijo-Goldstein Bedingungen wie folgt
einsetzen. Man beginnt mit einer Schrittweite von
$\alpha_{k} = \alpha_{k-1} > 0$, die im letzten Iterationsschritt $k-1$
zu einem Abstieg geführt hat und testet mit dieser die Armijo-Goldstein
Bedingungen aus {prf:ref}`def:armijo-goldstein`. Ist die erste
Ungleichung nicht erfüllt, d.h., für den tatsächlichen Abstieg gilt
$c_1 E_k(\alpha_k) \leq D_k(\alpha_k)$, so verkleinert man die
Schrittweite (zum Beispiel durch Halbierung). Ist andererseits die
zweite Ungleichung nicht erfüllt, d.h., für den tatsächlichen Abstieg
gilt $D_k(\alpha_k) \leq c_2 E_k(\alpha_k)$, so vergrößert man die
Schrittweite entsprechend. Um nicht in einen periodischen Zyklus zu
geraten, sollte man zur Vergrößerung der Schrittweite einen anderen
Faktor als zur Verkleinerung wählen, etwa $\sigma = 1.5$. Die Wahl der
Schrittweite nach den Armijo-Goldstein Regeln ist also relativ einfach
durchführbar und führt zu einer beweisbaren Konvergenz eines
Abstiegsverfahrens, wie das folgende Theorem zeigt.

````{prf:theorem} Konvergenz von Abstiegsverfahren
:label: thm:konvergenz_abstieg
Sei ein Abstiegsverfahren der Form $x_{k+1} = x_k + \alpha_k p_k$
gegeben, mit einer Menge von Vektoren $p_k \in \R^n$, die für jeden
Punkt $x_k \in \Omega$, der kein stationärer Punkt ist, eine uniforme
Abstiegsrichung liefern, d.h., es existieren fixe Konstanten
$\beta, \gamma > 0$, so dass gilt
```{math}
\langle \nabla F(x_k), p_k \rangle \ < \ - \gamma \cdot \Vert \nabla F(x_k)  \Vert^{\beta+1}.
```
Die Folge der Schrittweiten $(\alpha_k)_{k\in\N}$ erfülle die
Armijo-Goldstein Bedingungen. Außerdem sei
$F: \mathbb{R}^n \rightarrow \mathbb{R}$ eine nach unten beschränkte,
stetig differenzierbare Zielfunktion, für die somit die Niveaumenge
$K \coloneqq \{ x \in \R^n ~|~ F(x) \leq F(x_0) \}$ beschränkt ist.

Ist darüber hinaus die Folge der Abstiegsrichtungen $(p_k)_{k\in\N}$
beschränkt, dann bestitzt die Folge der Iterationsschritte
$(x_k)_{k\in\N}$ eine konvergente Teilfolge und jeder Häufungspunkt der
Folge ist ein stationärer Punkt der Zielfunktion $F$.

````

````{prf:proof} 
Falls für ein $k \in \N$ gilt, dass der Vektor $p_k  = \vec{0}$ ist, so
haben wir bereits einen stationären Punkt erreicht und die Aussage des
Theorems ist trivialerweise erfüllt.

Nehmen wir also im Folgenden an, dass $p_k \neq \vec{0}$ für alle
$k \in \N$ gilt und wir damit einen echten Abstieg vorliegen haben. In
diesem Fall impliziert die erste Armijo-Goldstein Bedingung
$c_1 E_k(\alpha_k) > D_k(\alpha_k)$, dass gilt
```{math}
F(x_{k+1}) - F(x_k) \: < \: c_1 \cdot \alpha_k \langle \nabla F(x_k) ,p_k \rangle \: < \: 0.
```
Induktiv gilt somit ebenfalls
```{math}
:label: eq:monotoner_abstieg
F(x_{k+1}) < F(x_k) < \ldots < F(x_0).
```

Also liegt die gesamte Folge $(x_k)_{k\in\N}$ in der nach Voraussetzung
beschränkten Menge $K \coloneqq \{ x \in \R^n ~|~ F(x) \leq F(x_0) \}$
und hat nach dem *Satz von Bolzano-Weierstrass* somit eine konvergente
Teilfolge $(x_{k_\ell})_{\ell \in \N}$. Darüber hinaus sehen wir ein,
dass die Menge $K$ ebenfalls abgeschlossen und somit nach dem *Satz von
Heine-Borell* sogar kompakt ist. Daher liegt der Grenzwert $x^* \in K$
ebenfalls in der Menge $K$.

Durch Einsetzen in die folgende Teleskopsumme
```{math}
F(x_k) - F(x_0) \ = \ \sum_{j=0}^{k-1} F(x_{j+1}) - F(x_j)
```
erhalten wir die stärkere Bedingung
```{math}
:label: eq:starke_bedingung
F(x_k) + c_1 \sum_{j=0}^{k-1} - \alpha_j \langle \nabla F(x_j), p_j \rangle \ < \ F(x_0).
```
Mit der nach Voraussetzung geltenden uniformen Schranke folgt dann für
fixe Konstanten $\beta, \gamma > 0$ und für alle $k \in \N$ die
Ungleichung
```{math}
|| \nabla F(x_k) ||^{\beta + 1} \: < \: -\frac{1}{\gamma} \langle \nabla F(x_k), p_k \rangle.
```
Da die Niveaumenge $K$ kompakt ist wissen wir, dass ein Punkt
$x^* \in K$ existiert in der die Zielfunktion $F$ ihr Minimum auf $K$
annimmt. Dann gilt wegen {eq}`eq:monotoner_abstieg` offensichtlich für
alle $k \in \N$
```{math}
0 \: < \: F(x_0) - F(x_k) \: < \: F(x_0) - F(x^*) \ =: \ M.
```
Zusammen mit der Ungleichung {eq}`eq:starke_bedingung` können wir somit
folgern, dass gilt
```{math}
\sum_{j=0}^{k-1}  \alpha_j \Vert \nabla F(x_j) \Vert^{\beta+1}  \ \leq \ \frac{1}\gamma \sum_{j=0}^{k-1} - \alpha_j \langle \nabla F(x_j), p_j \rangle \ \leq \ \frac{1}{\gamma c_1} (F(x_0) - F(x_k)) \: < \: \frac{M}{\gamma c_1}.
```
Damit gilt offensichtlich
$\alpha_k ||\nabla F(x_k)||^{\beta + 1} \rightarrow 0$ und somit gilt
ebenfalls $\alpha_k \nabla F(x_k) \rightarrow \vec{0}$ für
$k \rightarrow \infty$.

Nun müssen wir abschließend noch zeigen, dass die Folge der
Schrittweiten $(\alpha_k)_{k \in \N}$ selbst nicht gegen Null
konvergiert, damit $\nabla F(x_k) \rightarrow \vec{0}$ gilt und somit
die Folge der Iterationsschritte $(x_k)_{k\in\N}$ gegen einen
stationären Punkt der Zielfunktion $F$ konvergiert. Nehmen wir also das
Gegenteil für eine Teilfolge $(\alpha_{k_\ell})_{\ell \in \N}$ an, die
gegen Null konvergiert, so gilt auch
$\alpha_{k_\ell} p_{k_\ell} \rightarrow \vec{0}$. Damit gilt aber schon
für beliebiges $\epsilon > 0$ mit $(1-\epsilon) > c_2$, dass für
hinreichend große $k_\ell \in \N$ die folgende Ungleichung erfüllt ist:
```{math}
F(x_{k_\ell} + \alpha_{k_\ell} p_{k_\ell}) - F(x_{k_\ell}) \: \leq \: (1-\epsilon) \cdot \alpha_{k_\ell} \langle \nabla F(x_{k_\ell}), p_{k_\ell} \rangle \: < \: c_2 \cdot \alpha_{k_\ell} \langle \nabla F(x_{k_\ell}), p_{k_\ell}\rangle.
```
Dies ist jedoch nicht möglich, da die Armijo-Goldstein Bedingungen nach
Voraussetzung erfüllt sind. ◻

````

Für das Gradientenabstiegsverfahren bzw. dessen Varianten in
{ref}`ss:gradient_descent` können wir im folgenden Korollar noch mehr
zeigen,, da die Abstiegsvektoren $p_k \in \R^n$ in direkter Verbindung
zum Gradienten der Zielfunktion $-\nabla F(x_k)$ stehen.

````{prf:corollary} Konvergenz des Gradientenabstiegverfahrens
Sei $F: \mathbb{R}^n \rightarrow \mathbb{R}$ eine nach unten
beschränkte, stetig differenzierbare Funktion, so dass die Niveaumenge
$K \coloneqq \{ x \in \R^n ~|~ F(x) \leq F(x_0) \}$ beschränkt ist.
Gegeben sei außerdem eine Wahl an Abstiegsrichtungen der Form
```{math}
p_k \ = \ -A_k \nabla F(x_k),
```
wobei für jedes $k \in \N$ die Matrix $A_k \in \mathbb{R}^{n \times n}$
symmetrisch positiv definit ist. Der kleinste und größte Eigenwert der
Matrizen $A_k$ seien darüber hinaus für jedes $k \in \N$ uniform durch
$\lambda_{max} \geq \lambda_{min} > 0$ nach unten bzw. nach oben
beschränkt.

Dann hat die Folge $(x_k)_{k \in \N}$ eine konvergente Teilfolge und
jeder Häufungspunkt ist ein stationärer Punkt von $F$.

````

````{prf:proof} 
In den Übungsaufgaben zu zeigen. ◻

````
