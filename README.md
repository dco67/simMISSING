# simMISSING — L'analyse des données manquantes

**Comprendre MCAR, MAR et MNAR, voir ce que chaque traitement fait aux résultats, et vérifier par simulation lequel résiste.**

simMISSING est un laboratoire pédagogique **R/Shiny**, en français, sur les données manquantes. On efface des valeurs selon un mécanisme choisi, puis on compare les traitements : cas complets, imputation par la moyenne, par régression (déterministe ou stochastique), imputation multiple et mise en commun selon les règles de Rubin. Les vraies valeurs effacées restent visibles, pour juger chaque méthode.

- [Utiliser l'application en ligne](https://shiny.dcstats.org/simMISSING/)
- [Consulter le manuel pédagogique](https://dcstats.org/simMISSING/)

## Objectifs pédagogiques

- distinguer les trois mécanismes de Rubin (1976) : MCAR, MAR et MNAR;
- comprendre ce que les données peuvent diagnostiquer, et ce qu'elles ne peuvent pas trancher (MAR ou MNAR);
- voir comment l'imputation par la moyenne écrase la variabilité et biaise les relations;
- comprendre l'imputation multiple : chaînes, nombre de jeux, mise en commun;
- évaluer la sensibilité des conclusions à une hypothèse MNAR.

## Fonctionnalités

| Volet | Contenu |
| --- | --- |
| Introduction | Concepts, bac à sable (un petit jeu de 20 observations), quiz, bibliographie. |
| Données et visualisation | Génération des données (N, variables, corrélations), structure du manque. |
| Mécanisme et traitement | Traitements comparés, mécanisme du manque, diagnostic, sensibilité MNAR. |
| Évaluer et valider | Étude de simulation, qualité de l'imputation (convergence des chaînes, nombre de jeux), données réelles. |

## Exécuter l'application localement

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "MASS", "mice"))
shiny::runApp("simMISSING")
```

## Licence

GPL-3.0. Voir le fichier `LICENSE`.

## Auteur

Daniel Coulombe, Ph.D. — <https://dcstats.org>
