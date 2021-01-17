phyloseq
================

  - [Importer dans phyloseq :](#importer-dans-phyloseq)
  - [Prokaryote Diversity](#prokaryote-diversity)
  - [Bar plot:](#bar-plot)
      - [Bar plot :](#bar-plot-1)
      - [Ordinate:](#ordinate)

## Importer dans phyloseq :

Ici on va appeler le package phyloseq, et ensuite savoir sa version. Ce
package permet d’importer, stocker, analyser et afficher graphiquement
des données de séquençage phylogénétique complexes qui ont déjà été
regroupées en unités taxonomiques opérationnelles (OTU):

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.34.0'

Ici on va Rapporter le numéro de version du Biostrings, ce package
permet une manipulation efficace des chaînes biologiques:

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## [1] '2.58.0'

Ici on va Rapporter le numéro de version du package ggplot2, ce package
permet la création déclarative de graphiques:

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.2'

La fonction theme\_set() permet de remplacer le thème actuel,et choisir
le thème classique dark-on light de ggplot2 qui pourrait mieux
fonctionner pour les présentations graphiques.

``` r
theme_set(theme_bw())
```

``` r
load("02_data-analysis")
```

# Prokaryote Diversity

Le but c’est de construire une matrice à partir des informations
encodées dans le fichier seqtab.mochim. Pour celà je vais créer une
variable sample.out contenant les noms des lignes se trouvant dans le
ficher seqtab.mochib en utilisant la fonction rawnames().

j’ai utilisé ensuite la fonction strsplit() qui va diviser le vecteur de
chaîne d’entrée samples.out en sous-chaînes. et puis sapply () qui va
prendre le résultat de la fonction strsplit() comme entrée et donne une
sortie une matrice. la variable profondeur va contenir toutes les noms
de lignes de fichiers seqtab.mochim. la variable date va contenir toute
les noms de lignes se trouvant dans Depht. après je vais créer une
variable samdf qui va être une tableau contenant les variables Depht,
Month et Zone comme colonnes.

``` r
samples.out <- rownames(seqtab.nochim)
Depht <- sapply(strsplit(samples.out, "D"), `[`, 1)
Month <- substr(Depht,0,71)
Zone <- substr(Depht,0,71)
samdf <- data.frame(Depht=Depht, Month=Month, Zone=Zone)
```

Ici je voulais structurer mais donner dans le tableau, et comme je
savais pas comment faire cela de manière plus élégante avec moins de
ligne de code, j’ai fait ça comme ça, l’essentiel que ça marche :D Ces
lignes de code vont permettre de mettre de remplir les lignes de la
colonne Depht avec les valeurs de profondeur correspondant aux séquences
appropriées.

``` r
samdf$Depht[c(17:22,30,31)] <- ("0m")
samdf$Depht[c(1,2,7,12,23:25,27,43,46,50,53,58)] <- ("1m")
samdf$Depht[c(67)] <- ("5m")
samdf$Depht[c(40,41,44)] <- ("10m")
samdf$Depht[c(3,8,11,26,28,29,32,37,47,51,62,64,68)] <- ("20m")
samdf$Depht[c(54,59)] <- ("25m")
samdf$Depht[33] <- ("120m")
samdf$Depht[52] <- ("200m")
samdf$Depht[c(65,66)] <- ("300m")
samdf$Depht[42] <- ("365m")
samdf$Depht[45] <- ("375m")
samdf$Depht[c(4,9,48,55,56,60,69)] <- ("500m")
samdf$Depht[63] <- ("75m")
samdf$Depht[c(5,6,10,16,35,36,38,39,49,57,61,70,71)] <- ("1000m")
samdf$Depht[c(13:15,34)] <- ("320m")
```

les zones putatives épipélagiques (\<200 m de profondeur) et
mésopélagiques (200-1000 m de profondeur) de l’océan Arctique .

``` r
samdf$Zone[c(1:3,7,8,11,12,17:33,37,40:41,43:44,46:51,53:54,58:59,62:64,67:68)] <- ("Photic Zone")
samdf$Zone[c(4:6,9:10,13:16,34:36,38:39,42,45,48:49,52,55:57,60:61,65,66,69:71)] <- ("Mesopelagic Zone")
```

Ces lignes de code vont permettre de mettre les noms des mois appropriés
à leurs séquences dans les lignes de la colonnes Month.

``` r
samdf$Month[c(1:10)] <- c("january")
samdf$Month[c(11:35)] <- c("mars")
samdf$Month[c(36:49)] <- c("may")
samdf$Month[c(50:61)] <- c("august")
samdf$Month[c(62:71)] <- c("november")
```

Avec les lignes précédentes j’ai pu faire un tableau avec trois
colonnes: Month et Depht, après je voulais mettre les nons de lignes se
trouvant dans la variable samples.out dans mon tableau.

``` r
rownames(samdf) <- samples.out
```

# Bar plot:

On crée d’abbord une variable ps et on assigne à cette variable les
valeurs des résultas de la fonction phyloseq, qui va permettre de créer
des objets de classe phyloseq. phyloseq() est une méthode de
constructive, C’est la principale méthode suggérée pour construire un
objet de niveau expérience (classe phyloseq) à partir de ses données de
composant (classes de données de composant: classe otu\_table,
classe\_données\_échantillon, classe-taxonomyTable, classe-phylo ).
sample\_data c’est la méthode suggérée pour à la fois construire et
accéder à une table de variables de niveau échantillon (samdf), si le
premier argument est un objet de niveau expérience (classe phyloseq),
alors le sample\_data correspondant est renvoyé. sample\_names() va
permettre d’obtenir les noms de l’échantillion (ps) et puis
prune\_samples() va permettre de filtrer les échantillons indésirables
en définissant ceux que vous souhaitez conserver.

``` r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
```

avec la commande suivante, on va faire apparaître les nouveaux taxons
courts dans des tableaux et des graphiques. Et on va pouvoir récupérer
les séquences d’ADN correspondant à chaque ASV (c’est plus pratique
d’utiliser les noms cours pour nos ASV que la séquence d’ADN
compléte), et donc on va stocker les séquences d’ADN dans les ASV
d’abbord pour pouvoir faire ça. Le package Biostrings contient des
classes et des fonctions pour représenter des chaînes biologiques telles
que l’ADN, l’ARN et les acides aminés.

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 11957 taxa and 71 samples ]
    ## sample_data() Sample Data:       [ 71 samples by 3 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 11957 taxa by 6 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 11957 reference sequences ]

## Bar plot :

Alors j’ai essayé d’avoir le même bar plot que celui de l’article, mais
j’y arrivais pas, le problème c’est qu’il faut pas afficher tous les
classes, mais qlq un, en effet la matrice taxa cpntient 71742 elements
ce qui fait que c’est quasi un possible de faire le tri sauf si y a
moyen de faire cela de manière automatisé et c’est sûr qu’on peut faire
ça mais je n’y arrivais pas, j’ai cherché sur internet mais comme je
sais pas ce que je suis sensée trouver la recherche n’a pas abouti à un
résultat.

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:37]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, y="Depht", x="Abundance", fill="Class") + facet_wrap(~Month, scales="free_x")
```

![](03_phyloseq_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Ordinate:

ici on va Transformer les données en proportions appropriées pour les
distances de Bray-Curtis, en utilisant l’indice Unifrac.

``` r
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.pcoA.bray <- ordinate(ps.prop,"unifrac", method="PCoA", distance="bray")
```

``` r
evals <- ord.pcoA.bray$values$Eigenvalues
plot_ordination(ps.prop, ord.pcoA.bray, color="Zone" ,title="pcoA",axes = c(1,2)) + 
  labs(col="Zone") + scale_colour_manual(values = c("red", "blue"))
```

![](03_phyloseq_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

plot\_ordination va permettre de tracer l’ordination PCoA qui est une
approche d’analyse de représenter la (dis) similitude inter-objets dans
un espace euclidien de faible dimension.
