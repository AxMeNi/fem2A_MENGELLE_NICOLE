# README.txt
Utilisation de programmes de simulations de processus : fem2A_MENGELLE_NICOLE

### Initialisation
Pour lancer les programmes, il faut : 
1) Tout d'abord clone le repo en local
2) Se placer dans le repo
3) Compiler le code avec make
4) Executer le code avec ./build/fem2a suivi d'une commande qui peut être -t (pour les tests), -s (pour les simulations), -h (pour obtenir de l'aide), -v (pour print les versions verbosées avec beaucoup de détails, mais celle-ci n'a pas eu le temps d'être mise en place)

### Tests
Plusieurs fonctions peuvent être tout d'abord testées. La liste des fonctions testables est accessible dans le main.cpp.
Pour tester une fonction il faut
1) Choisir quelle fonction doit être testée en définissant les booléens de la fonction run_tests dans le main.cpp sur true
2) Exécuter les tests avec le ./build/fem2a -t

/!\Attention, certains tests nécessitent de décommenter des lignes présents dans des fonctions du fichier fem.cpp afin de permettre les affichages. C'est le cas de la fonction local_to_global_matrix et de assemble_elementary matrix.

/!\Attention, pour le test de la fonction apply dirichlet_boundary_condition, il serait bien de définir un chemin de sortie pour le fichier de solution afin d'afficher autre chose que les valeurs x qui sont illisibles. Cela a été fait au cours des simulations, raison pour laquelle ce test n'a pas produit de fichier mesh et n'a pas de solution visualisable.

### Simulations
Afin de réaliser une simulation il faut :
1) Choisir quelle simulation effectuée dans le main.cpp en passant le booléen qui lui est associé dans run_simu sur true.
2) Choisir le maillage sur lequel sera effectué la simulation dans le main.cpp à l'intérieur du if associé au test (les maillages 
3) Exécuter avec ./build/fem2a -m

Le fichier solution.bb et le fichier solution.mesh associé au problème sont stockés dans data/output. Le nom des solutions sont modifiables dans le fichier simu.h

### Affichage simulations
Pour afficher les solutions, éxecuter (via le projet Medit) la commande : path_to_medit data/output/nom_de_la_solution.mesh

