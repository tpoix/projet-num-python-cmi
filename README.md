# Projet numérique Python - CMI M1

### Groupe : B
### Sujet : Spectrométrie de masse
Membres :
- Sarah DELORD
- Axelle MOREAU
- Thomas POIX

## Objectif du projet
L'objectif de ce projet est de pouvoir trouver les compositions possibles d'une molécule à partir de sa masse. Cela permet d'aider à trouver une molécule inconnue dont on vient de mesurer la masse expérimentale.

## Fonctionnement du projet
Le fonctionnement du projet est simple dans la théorie : il faut tester toutes les combinaisons possibles de groupements d'atomes et calculer la masse théorique totale de ces combinaisons, à partir des données fournies dans le fichier Excel.<br>
Si une combinaison a une masse théorique très proche de la masse expérimentale trouvée, alors on considère que la molécule cherchée peut être cette combinaison. On présente donc au final toutes les molécules candidates à cette masse expérimentale.<br>
Etant donné qu'il faut tester un nombre conséquent de combinaisons, l'algorithme prend du temps à s'exécuter dans son ensemble.
