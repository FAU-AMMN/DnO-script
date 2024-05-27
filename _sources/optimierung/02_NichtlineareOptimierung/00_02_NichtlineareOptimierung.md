(ch:optimierung)=
Numerische Optimierung
===

Optimierung ist ein omnipräsentes Phänomen, dass nicht nur in der
abstrakten Welt der Mathematik existiert. Viel mehr stellt es ein
naturgegebenes Prinzip dar, welches überall um uns herum zur Anwendung
kommt. In der Physik beispielsweise spielt Optimierung eine wesentliche
Rolle bei der Modellierung von Energiezuständen auf unterschiedlichen
Skalen: Moleküle formieren sich in einer Art, die die Gesamtenergie des
Teilchensystems unter Berücksichtigung aller wechselseitigen Kräfte
minimiert. Gleichzeitig strebt das Universum mit all seinen Planeten,
Sternen und Galaxien nach einem Zustand von maximaler Verteilung,
beschrieben durch die thermodynamische Größe der Entropie. Auch hier
folgt die Zunahme der Entropie dem Prinzip der Energieminimierung des
Gesamtsystems. Menschen betreiben seit jeher Optimierung in den
verschiedensten Anwendungen, oft mit unterschiedlichen Motivationen.
Flugzeuge werden von Ingenieuren so entworfen und gebaut, dass sie
möglichst stromlinienförmig aussehen, um damit den Reibungswiderstand in
der Luft zu minimieren und gleichzeitig den nötigen Auftrieb für einen
sicheren Flug zu erzeugen. Fondmanager streben danach Portfolios zu
erstellen, deren Gewinn möglichst maximal ist und dennoch
Spekulationsrisiken vermeiden. Die gesamte Automatisierung der
Industrie, angefangen bei den ersten Manufakturen hin zu modernen
vollautomatischen Roboterfabriken, dient lediglich dem Prinzip der
Gewinnmaximierung durch Minimierung der Produktionskosten.

Es ist also nicht wirklich überraschend, dass sich ein gesamtes
Teilgebiet der Angewandten Mathematik mit der Theorie der Optimierung
befasst und somit dazu beiträgt viele Optimierungsprobleme besser zu
verstehen und zu lösen. Im folgenden Abschnitt wollen wir uns
hauptsächlich mit der unbeschränkten (oder: unrestringierten)
Optimierung beschäftigen und uns nützliche Werkzeuge zum Lösen von
numerischen Problemstellungen herleiten. Nach einer *allgemeinen
mathematischen Einführung* in {ref}`s:opt_grundlagen` beginnen wir in
Kapitel {ref}`s:abstiegsverfahren` mit einer Klasse von Algorithmen, die
einer einfachen Idee folgen: die *numerischen Abstiegsverfahren*. In
{ref}`s:cg_verfahren` behandeln wir insbesondere das Verfahren der
*konjugierten Gradienten*, welches zur iterativen Lösung von großen,
linearen Gleichungssystemen mit besonderen Eigenschaften genutzt werden
kann. Zum Schluss untersuchen wir in {ref}`s:nichtdiffbare_optimierung`
zwei moderne Optimierungsalgorithmen zur Lösung von *konvexen,
nicht-differenzierbaren Problemen*.

