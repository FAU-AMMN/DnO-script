(numerische-lösungsverfahren-für-anfangswertprobleme)=
Numerische Lösungsverfahren für Anfangswertprobleme
===

Differentialgleichungen spielen eine wichtige Rolle in der Modellierung
vieler naturwissenschaftlicher Phänomene. Insbesondere in der Physik und
den Ingenieurswissenschaften ist das Verständnis und das Lösen von
Differentialgleichungen essentiell, da sich viele physikalische Gesetze
und Zusammenhänge durch diese beschreiben lassen. Leider existieren für
viele Differentialgleichungen keine geschlossen darstellbaren
analytischen Lösungen, so dass man Lösungen mit Hilfe numerischer
Verfahren approximieren muss. In diesem und nächstem Kapitel widmen wir
uns daher der numerischen Lösung von verschiedenen Problemen, die
Differentialgleichungen beinhalten. Dies sind insbesondere Anfangswert-
und Randwertprobleme.

Mathematisch gesehen ist eine Differentialgleichung eine Gleichung, in
der eine unbekannte Funktion und ihre Ableitungen auftreten. Hierbei
werden folgende Kriterien zur Unterscheidung von Differentialgleichungen
genutzt.

1.  Die größte auftretende Ableitung der unbekannten Funktion bestimmt
    die **Ordnung der Differentialgleichung**.

2.  Je nachdem ob Ableitungen der unbekannten Funktion bezüglich einer
    einzigen Variablen oder unterschiedlicher Variablen auftreten
    sprechen wir von einer **gewöhnlichen Differentialgleichung** oder
    einer **partiellen Differentialgleichung**.

3.  Werden gleich mehrere solche Funktionen durch mehrere Gleichungen
    beschrieben, so spricht man von einem
    **Differentialgleichungssystem**.

Wir werden in den folgenden Abschnitten zunächst die wichtigsten
theoretischen Erkenntnisse für gewöhnliche Differentialgleichungen
wiederholen und insbesondere Existenz- und Eindeutigkeitsaussagen
diskutieren. Diese werden für die spätere Konstruktion und Analyse von
numerischen Lösungsverfahren hilfreich sein.

