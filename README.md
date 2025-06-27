**Analyse numérique d’équations aux dérivées partielles par différences finies et implémentation optimisée pour le calcul haute performance**

*par Jean-Baptiste Gaillot*

# Informations et introduction

Le but de ce projet est de concevoir des schémas d’approximation numérique pour résoudre des problèmes d’équations aux dérivées partielles sur un maillage fini, en utilisant la méthode des différences finies, puis de les implémenter en différentes versions afin de comparer leurs performances.

Des remerciements sont adressés à Francesco Bonaldi, enseignant-chercheur, pour la confiance accordée en vue de la réalisation autonome de ce projet.

Pour chaque problème, l'approche sera la suivante :

**Partie mathématiques**
- concevoir un ou plusieurs schémas numériques pour obtenir une solution approchée du problème,  
- s'assurer de l'existence et de l'unicité de la solution approchée,  
- s'assurer des bonnes propriétés du schéma (consistance, convergence et erreur locale),  
- concevoir un ou plusieurs schémas de résolution de l'éventuel système linéaire associé à cette méthode.  
  
**Partie informatique**
- implémenter des fonctions de résolutions du problème et un programme principal,  
- implémenter le calcul d'une solution exacte connue dans le but de calculer l'erreur entre la solution approchée et la solution exacte,  
- implémenter la résolution du problème en différentes version dans le but de comparer les performances (différents schémas en versions naïves,   séquentielles, parallèles et utilisation de bibliothèque).  

Les différents résultats (erreurs et temps d'exécutions) seront présentés sous forme de tableaux et de graphiques. Le langage de programmation utilisé est ```C```.

# Partie mathématique — Conventions

> **Notations**  
> - $N + 1$ est le nombre de nœuds dans une direction spatiale, $N_t + 1$ est le nombre de nœuds en temps,  
> - $h$ est le pas de discrétisation en espace, $h_t$ est le pas de discrétisation en temps,  
> - $x_i := ih, y_j := jh$ et $t_k := kh_t$ pour $i, j \in$ { $0, ..., N$ }, $k \in$ { $0, ..., N_t$ },
> - $u\left(x_i, y_j, t_k\right) :\approx u_{i,j}^k$,  
> - En 1D : $u := \left(u_1 ~ \cdots ~ u_{N-1}\right)^T$ est le vecteur de la solution approchée,  
> - En 2D : $u_j := \left(u_{1,j} ~ \cdots ~ u_{N-1, j} \right)^T$ est le vecteur d'une ligne de la solution approchée,  
> - $E_h$ est l'erreur de troncature en espace, $E_{h_t}$ est l'erreur de troncature en temps,  
> - $\lVert e_{h,h_t} \rVert_{\infty}$ est l'erreur locale en fonction de $h$ et $h_t$.

> **Définitions**  
> - Un schéma numérique est *consistant en espace* lorsque $\lim_{h \to 0} \lVert E_h \rVert = 0$ et est *consistant en temps* lorsque $\lim_{h_t \to 0} \lVert E_{h_t} \rVert = 0$.  
> - Un schéma numérique est *convergent* lorsque $\lim_{h,h_t \to 0} \lVert e_{h,h_t} \rVert_{\infty} = 0$.

# Partie informatique — Organisation du projet

## Structure du projet

- **Équation de Poisson en dimension 1** : `Probleme-1D`
  - Version 0 (`base`) : Base (Gauss)
  - Version 1 (`sequentiel-1`) : Méthode itérative, séquentiel (Jacobi)
  - Version 2 (`parallele-1`) : Méthode itérative, parallèle OpenMP (Jacobi)
  - Version 3 (`parallele-2`) : Méthode itérative, parallèle MPI (Jacobi)
  - Version 4 (`sequentiel-2`) : Méthode directe (Cholesky)

- **Équation de Poisson en dimension 2** : `Probleme-2D`
  - Version 0 (`base`) : Base (Gauss)
  - Version 1 (`sequentiel-1`) : Méthode itérative, séquentiel (Jacobi)
  - Version 2 (`parallele-1`) : Méthode itérative, parallèle OpenMP (Jacobi)
  - Version 3 (`parallele-2`) : Méthode itérative, parallèle MPI bloquant (Jacobi)
  - Version 4 (`parallele-3`) : Méthode itérative, parallèle MPI non bloquant (Jacobi)
  - Version 5 (`sequentiel-2`) : Méthode directe (Cholesky)
  - Version 6 (`sequentiel-3`) : Méthode directe, bibliothèque ```cholmod``` (Cholesky)

- **Équation des ondes en dimension 1** : `Probleme-Ondes`
  - Version $1$ (```sequentiel-1```) : Schéma explicite, séquentiel

- **Équation de la chaleur en dimension 2** : `Probleme-Chaleur`
  - Version 1 (`sequentiel-1`) : Schéma explicite, séquentiel
  - Version 2 (`parallele-1`) : Schéma explicite, parallèle OpenMP
  - Version 3 (`parallele-2`) : Schéma explicite, parallèle MPI
  - Version 4 (`sequentiel-2`) : Schéma implicite, bibliothèque `cholmod`
  - Version 5 (`sequentiel-3`) : Schéma semi-implicite, bibliothèque `cholmod`

L'arborescence complète est visualisable dans le fichier ```Arborescence.pdf```.

## Informations

- Le répertoire ```Fonctions-communes``` contient des fichiers de fonctions qui sont appelées pour chaque problème (affichages, opérations sur des tableaux, sauvegardes de résultats dans un fichier, ...)
- Pour chaque problème, il y a les répertoires ```Sources```, ```Librairies``` (qui contient les déclarations de fonctions), ```Objets```, ```Binaires``` et ```Textes``` (qui contient les résultats). À l'intérieur de chaque répertoire, il y a un sous-répertoire pour chaque version du problème.
- Pour chaque problème, il y a un Makefile incluant les règles nécessaires pour compiler et/ou nettoyer chaque version et un script pour exécuter chaque version pour les paramètres voulus. Pour chaque exécution, des informations de debug sur les paramètres du problème, l'erreur et le temps d'exécution sont affichées dans le terminal.
- Pour chaque version d'un problème, il y a principalement les fichiers suivants :
  - ```main.c``` : contient le programme principal,
  - ```resolution.c``` : contient les fonctions de résolution,
  - ```parallele.c``` : contient les fonctions pour préparer les données MPI (informations sur les nœuds à traiter, topologie cartésienne, échange des halos, ...) (uniquement pour les versions MPI).
- Des variables définies dans ```main.c``` ou ```parallele.c``` (informations sur le maillage, informations MPI, ...) sont souvent globales et externes pour être utilisées par le fichier ```resolution.c```.
- Certaines fonctions qui sont appelées de nombreuses fois sont mises inline.
- Les flags de compilations ```-O3``` et ```-Wall``` sont utilisés pour chaque version.
- Les résultats qui seront présentés ont été exécutés sur une machine ordinaire (8 CPU et 16 Go de mémoire) et non un cluster de calcul.

