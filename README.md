Is Sharing Caring?
==================

# Investigating the Relationship Between HGT Rates and CRISPR-Cas Presence Using Network Analysis

---
abstract: |
    [hgt]{acronym-label="hgt" acronym-form="singular+short"} is a mechanism
    by which organisms (mainly prokaryotes) can share genetic material
    outside of inheritance. [hgt]{acronym-label="hgt"
    acronym-form="singular+short"} has proven to have significant effects on
    bacterial genome evolution, allowing for increased genetic diversity and
    advanced niche adaptation. [crspc]{acronym-label="crspc"
    acronym-form="singular+short"} is an adaptive immune system in
    prokaryotes that has garnered a lot of research attention recently,
    largely due to it's applications in gene editing. Due to the nature of
    how it works, using guide RNA to cut DNA, [crspc]{acronym-label="crspc"
    acronym-form="singular+short"} systems have been thought to affect rates
    of [hgt]{acronym-label="hgt" acronym-form="singular+short"}. Effort has
    mostly been focused on how [crspc]{acronym-label="crspc"
    acronym-form="singular+short"} systems affect the mechanisms of
    [hgt]{acronym-label="hgt" acronym-form="singular+short"} and thus little
    is known about its effects on [hgt]{acronym-label="hgt"
    acronym-form="singular+short"} rates. This work uses is a
    network-theoretic approach to better characterize the effects of the
    presence of [crspc]{acronym-label="crspc" acronym-form="singular+short"}
    systems on [hgt]{acronym-label="hgt" acronym-form="singular+short"}
    rates within bacterial populations. This approaches makes use of
    phylogenetic methods for estimating [hgt]{acronym-label="hgt"
    acronym-form="singular+short"} rates, improving on estimates using
    methods based on genome composition used previously. Understanding the
    effects of [crspc]{acronym-label="crspc" acronym-form="singular+short"}
    on [hgt]{acronym-label="hgt" acronym-form="singular+short"} may help
    uncover potential targets for curbing the spread antibiotic resistance
    genes.
bibliography: 'writing/refs.bib'
---

Background
==========

What is CRISPR-Cas?
-------------------

[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems are
sets of nucleotide motifs (spacers) interspaced with nucleotide repeats
([crsp]{acronym-label="crsp" acronym-form="singular+short"}s) and
[cas]{acronym-label="cas" acronym-form="singular+short"} proteins (with
sequences usually adjacent to the [crsp]{acronym-label="crsp"
acronym-form="singular+short"} motifs) that have an adaptive immune
function in many bacteria and archaea (Rath et al. 2015). Each
nucleotide motif is the result of a DNA molecule that was previously
taken up by the host, serving as a marker for the
[cas]{acronym-label="cas" acronym-form="singular+short"} proteins to
degrade any DNA matching the motif (Rath et al. 2015). If a bacterium
possessing a [crspc]{acronym-label="crspc"
acronym-form="singular+short"} system is infected with a phage and
survives, a motif representative of that phage can be integrated as a
spacer. If the bacterium is reinfected with the same phage strain the
spacers will guide [cas]{acronym-label="cas"
acronym-form="singular+short"} proteins to the invading phage DNA and
degrade it, hence adaptive immunity. Although
[crspc]{acronym-label="crspc" acronym-form="singular+short"} appears to
have evolved to degrade viral DNA system, the majority of
[crsp]{acronym-label="crsp" acronym-form="singular+short"} spacers have
been found to match bacterial [mge]{acronym-label="mge"
acronym-form="singular+short"}s with no known viral match (Shmakov et
al. 2017).

As of 2017, over $45\%$ of bacterial genomes analyzed ($n=6782$) appear
to contain [crsp]{acronym-label="crsp" acronym-form="singular+short"}
motifs (Grissa, I. and Drevet, C. and Couvin, D. 2017). Moreover,
[crsp]{acronym-label="crsp" acronym-form="singular+short"} motifs show
significant diversity between individual bacteriums as they represent a
chronological history of spacer acquisition (usually via viral infection
or [mge]{acronym-label="mge" acronym-form="singular+short"} "infection")
for that specific bacterium(Rath et al. 2015). There still exist many
bacterial strains, and even entire genera with no *known*
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems,
although they may have simply not been discovered yet (Zhang and Ye
2017; Haft et al. 2005). Furthermore, the diversification of CRISPR-Cas
systems is driven further by [hgt]{acronym-label="hgt"
acronym-form="singular+short"} acting on [crsp]{acronym-label="crsp"
acronym-form="singular+short"} and [cas]{acronym-label="cas"
acronym-form="singular+short"} components independently, adding another
level of complexity to the propagation of CRISPR-Cas systems (Rath et
al. 2015).

Horizontal Gene Transfer
------------------------

### Mechanisms

[hgt]{acronym-label="hgt" acronym-form="singular+short"} can be defined
as the exchange of genetic information across lineages (Zhaxybayeva and
Doolittle 2011), as opposed to vertical gene transfer between parents
and offspring (Ravenhall et al. 2015). It is a source of genetic
variation, allowing organisms to adapt quickly by copying a gene with a
specific function, rather than evolving it themselves (Ravenhall et al.
2015; Marri, Hao, and Golding 2007). There are 3 main mechanisms of
[hgt]{acronym-label="hgt" acronym-form="singular+short"}\
**Transformation** Free floating DNA is taken up by a bacterium and
incorporated into the genome (Zhaxybayeva and Doolittle 2011). DNA is
not always incorporated successfully even if it taken up.\
**Conjugation** The sharing of genetic material through cell-to-cell
bridges, the genes for which are usually carried on a plasmid (Davison
1999).\
**Transduction** Genes or DNA fragments can be transferred through
either lytic or lysogenic bacteriophages (Griffiths et al. 2000).
Bacterial DNA can be accidentally packaged into the lysogenic phage head
during cell lysis and integrate into the next infected host. (Griffiths
et al. 2000). Lysogenic phages can take up bacterial DNA flanking the
viral sequence and bring it with them to the next host(Griffiths et al.
2000).

It should be noted that *successful* [hgt]{acronym-label="hgt"
acronym-form="singular+short"} requires that a gene be maintained,
either by genomic integration or plasmid replication. Frequently,
putatively transferred genes are either lost quickly or diverge quickly
due to minimal selective pressure maintaining them (Hao and Golding
2006).

### Rate Influencing Factors

The rate of [hgt]{acronym-label="hgt" acronym-form="singular+short"} in
bacteria is constantly in flux(Popa and Dagan 2011). The more exogenous
DNA, higher population density or higher phage density means more DNA is
available for transfer (Zhaxybayeva and Doolittle 2011). Just like
mutation rates, [hgt]{acronym-label="hgt" acronym-form="singular+short"}
rates are also thought to evolve in response to environmental factors or
selective pressure (Wielgoss et al. 2013; Mozhayskiy and Tagkopoulos
2012). The metaboliccost or the possibility of receiving toxic or
incompatible genes can dis incentivize a cell to produce the machinery
required for [hgt]{acronym-label="hgt" acronym-form="singular+short"}
(Baltrus 2013). But for bacteria in hospitals, the potential benefit of
receiving antibiotic resistance genes can outweigh any potential danger
or metabolic cost, inducing increased bacterial competence (Dzidic and
Bedeković 2003) It has been suggested that genes acquired via
[hgt]{acronym-label="hgt" acronym-form="singular+short"} are often
quickly lost since they often confer no advantage to a cell's current
selective pressures (Hao and Golding 2006). Ultimately
[hgt]{acronym-label="hgt" acronym-form="singular+short"} rates are
influenced by a variety of factors related balancing the potential
fitness costs and benefits.

Phylogenomic Networks
---------------------

[hgt]{acronym-label="hgt" acronym-form="singular+short"} is an important
factor in understanding evolution in prokaryotes. In graph theory a tree
is defined as a graph where there is only one path between every pair of
nodes. In phylogenetics this implies there is only one path for genetic
material to transfer between organisms, that path being vertical
inheritance. The existence of [hgt]{acronym-label="hgt"
acronym-form="singular+short"} demonstrates that the tree model is
clearly an incomplete representation of genetic relationships between
bacterial [otu]{acronym-label="otu" acronym-form="singular+short"}s.
Genetic material can be transferred outside of reproduction, creating
multiple paths through which a gene can exist in two different
[otu]{acronym-label="otu" acronym-form="singular+short"}s (Zhaxybayeva
and Doolittle 2011). The frequency of [hgt]{acronym-label="hgt"
acronym-form="singular+short"} among prokaryotes has lead many to
re-evaluate the concept of a "prokaryotic tree of life", which ignores
these horizontal interactions (Kunin et al. 2005). This prompted the
idea of a prokaryotic network of life (as opposed to a tree), with edges
indicating both vertical and horizontal transfers of genetic material
(Kunin et al. 2005). Edges can now connect closely or distantly related
[otu]{acronym-label="otu" acronym-form="singular+short"}s if
[hgt]{acronym-label="hgt" acronym-form="singular+short"} has occurred
between them.

Detection
---------

While understanding that [hgt]{acronym-label="hgt"
acronym-form="singular+short"} is important to bacterial evolution and
networks provide a useful theoretic framework to study it, constructing
such networks is not trivial. Previously researches have used methods
based on gene composition (GC content, codon usage, di/tri-nucleotide
frequency) to detect transferred genes. These methods are often
inadequate as there are several reasons why gene composition may differ
in one region vs the rest of the genome that are not related to
transfer. Further they often cannot discriminate
[mge]{acronym-label="mge" acronym-form="singular+short"}s which may be
fasely annotated as transfers. Phylogenetic methods, which rely on
recognizing discordance between gene trees and species trees, are
currently considered the best for inferring [hgt]{acronym-label="hgt"
acronym-form="singular+short"} events. If a gene tree is found to have a
significantly different topology from a species tree, this difference
may be the result of an [hgt]{acronym-label="hgt"
acronym-form="singular+short"} event (Than et al. 2007). Generally
phylogenetic methods are preferred for multiple reasons:

-   Can make use of multiple genomes at once (Ravenhall et al. 2015)

-   Require explicit evolutionary models, which come with their own
    framework for hypothesis testing and model selection (Ravenhall et
    al. 2015).

-   [hgt]{acronym-label="hgt" acronym-form="singular+short"} events
    identified by parametric methods are often found by phylogenetic
    methods as well (Ravenhall et al. 2015).

-   In recent years, the requirements of computing power and multiple
    well sequenced genomes for phylogenetic methods have become easier
    and easier to meet (Ravenhall et al. 2015).

While detecting [hgt]{acronym-label="hgt" acronym-form="singular+short"}
events with high degrees of certainty is still difficult, much progress
has been made in recent years,especially using phylogenetic methods
(Ravenhall et al. 2015). Events that may lead to false diagnosis of
[hgt]{acronym-label="hgt" acronym-form="singular+short"} are: incomplete
lineage sorting, gene duplication followed by loss in one of the
descendant lineages or homologous recombination (Ravenhall et al. 2015;
Than et al. 2007).

How Does CRISPR Affect HGT?
---------------------------

### Interference Mechanisms

[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems
have also been found to interfere with transformation-mediated
[hgt]{acronym-label="hgt" acronym-form="singular+short"}, by degrading
foreign DNA taken up by a cell (Zhang et al. 2013). They have been shown
to interfere with conjugation by targeting genes on conjugative plasmid
(Marraffini and Sontheimer 2008). They have also been shown to interfere
with transduction by creating immunity to phage infection (Marraffini
and Sontheimer 2008). Thus it been hypothesized that lower rates of
[hgt]{acronym-label="hgt" acronym-form="singular+short"} will be
observed in bacterial strains with [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems than without (Marraffini and
Sontheimer 2008).

### Complexities And Costs Of CRISPR-Cas Systems

Since antibiotic resistance genes are often transferred on plasmids
maintaining a [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems can present a large opportunity
cost, especially in environments like hospitals or trees (Dzidic and
Bedeković 2003). [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems incur a metabolic cost, as
[cas]{acronym-label="cas" acronym-form="singular+short"} proteins, guide
RNAs and spacer acquisition proteins must all be expressed
consistently(Rath et al. 2015). Much like the human immune system,
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems can
have off-target effects, sometimes affecting [hgt]{acronym-label="hgt"
acronym-form="singular+short"} (Bondy-Denomy and Davidson 2014). While
resisting lytic phage infection clearly provides some fitness benefit,
[crspc]{acronym-label="crspc" acronym-form="singular+short"} has also
been shown to help resist prophages which can provide super-infection
immunity or reduce competitor populations (Bondy-Denomy and Davidson
2014; Watson, Staals, and Fineran 2018). It has also been shown that
spacers targeting a bacterium's own DNA can be acquired, leading to an
auto-immune response (Stern et al. 2010). As
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems
prevent phage infection, mechanisms, denoted as
anti-[crsp]{acronym-label="crsp" acronym-form="singular+short"}s, have
evolved in certain phages making them immune to
[crsp]{acronym-label="crsp" acronym-form="singular+short"} (Bondy-Denomy
and Davidson 2014). This has a two-fold effect, as it can increase the
susceptibility of the host to infection reducing the fitness benefit of
[crsp]{acronym-label="crsp" acronym-form="singular+short"}, but it can
also increase transduction-mediated [hgt]{acronym-label="hgt"
acronym-form="singular+short"} (Bondy-Denomy and Davidson 2014).

### Balancing the Cost of [crsp]{acronym-label="crsp" acronym-form="singular+short"} and [hgt]{acronym-label="hgt" acronym-form="singular+short"}

Due to the myriad of fitness costs associated with consistently
expressing [crspc]{acronym-label="crspc" acronym-form="singular+short"}
systems, bacteria have appeared to develop strategies to mitigate these
costs. It has been posited that [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems need only be present in some
proportion of a population as they can be horizontally transferred
themselves between members(Godde and Bickerton 2006). This allows
populations to maintain phage immunity while isolating the metabolic
cost to only a few organisms (Bondy-Denomy and Davidson 2014). It should
also be noted that [crspc]{acronym-label="crspc"
acronym-form="singular+short"} genes are not necessarily constitutively
transcribed, potentially allowing bacteria to tune their
[crspc]{acronym-label="crspc" acronym-form="singular+short"} to suit
their selective pressures (Bondy-Denomy and Davidson 2014). The presence
of [crspc]{acronym-label="crspc" acronym-form="singular+short"} systems
have also been shown to actually enhance [hgt]{acronym-label="hgt"
acronym-form="singular+short"} at a population level via transduction by
reducing total phage abundance (Watson, Staals, and Fineran 2018).

A bioinformatic analysis has shown increased levels of gene insertion
and deletion as in Firmicutes with [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems compared to closely related
outgroups without (Zambelis, Dang, and Golding 2015). The effects of
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems on
[hgt]{acronym-label="hgt" acronym-form="singular+short"} rates are
highly complex, owning in no small part to the broad range of these
effects, how [crsp]{acronym-label="crsp" acronym-form="singular+short"}
activity can be modulated and the transfer of
[crsp]{acronym-label="crsp" acronym-form="singular+short"} systems
within a population (Bondy-Denomy and Davidson 2014). Taking a
systematic approach may help elucidate the dynamics between
[crsp]{acronym-label="crsp" acronym-form="singular+short"} system
presence and [hgt]{acronym-label="hgt" acronym-form="singular+short"}
rate.

Discussion
==========

Gene Indel Rates are Different for [crsp]{acronym-label="crsp" acronym-form="singular+short"} and Non-[crsp]{acronym-label="crsp" acronym-form="singular+short"} [otu]{acronym-label="otu" acronym-form="singular+short"}s
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

For most genera, the gene indel rate for non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} genera is larger than for
[crsp]{acronym-label="crsp" acronym-form="singular+short"} genera. This
is in-line with literature surrounding the mechanisms of
[hgt]{acronym-label="hgt" acronym-form="singular+short"} as
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems are
meant to stop the integration of foreign DNA into the bacterial genome.
Despite this, the mean node degree of [crsp]{acronym-label="crsp"
acronym-form="singular+short"} and non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s is relatively similar across genera. One
reason for this discrepancy may be that there are often more
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"}s than
[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"}s, thus more
genes are exchanged between all non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s, but each
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"} transfers genes
at a similar rate to [crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s. One possible explanation for certain
genera having very high gene indel rates for [crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s may be that it is an efficient way to
acquire new spacers. [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems may enhance
[hgt]{acronym-label="hgt" acronym-form="singular+short"} to preemptively
acquire new spacers from the environment in response to environmental
phage density.

Phylogenomic Networks Have Low Assortativity
--------------------------------------------

There seems to be significant [hgt]{acronym-label="hgt"
acronym-form="singular+short"} between [crsp]{acronym-label="crsp"
acronym-form="singular+short"} and non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s, with no clear clustering, assortativity
or modularity among most of the networks examined.
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems do
not appear to have a segregating effect on the network, but do appear to
have a population level effect of decreasing the rate of
[hgt]{acronym-label="hgt" acronym-form="singular+short"}. As suggested
by (Watson, Staals, and Fineran 2018) [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems may have a population level
effect on [hgt]{acronym-label="hgt" acronym-form="singular+short"} rate,
but it appears to be suppressive in this case, as opposed to theirs.

[hgt]{acronym-label="hgt" acronym-form="singular+short"} Dynamics Vary Across Bacterial Genera
----------------------------------------------------------------------------------------------

Despite some trends being observable, the one constant is that there is
significant variability between genera with regards to
[hgt]{acronym-label="hgt" acronym-form="singular+short"}.
[hgt]{acronym-label="hgt" acronym-form="singular+short"} rates can be
similar or significantly different between [crsp]{acronym-label="crsp"
acronym-form="singular+short"} and non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s, either can be larger, both can be
similar and either large or small. While the means are similar, the
shapes of the distributions of network assortativity and modularity are
not homogeneous.

Methods
=======

Summary
-------

The goal of the project is to create a phylogenetic network from a set
of protein fasta files and the corresponding nucleotide sequences. In
this case all full genomes for a given bacterial genus, for analysis of
[hgt]{acronym-label="hgt" acronym-form="singular+short"}. The workflow
is as follows for a single genus:

1.  Download fasta files

2.  Filter mobile genetic elements from genomes

3.  Cluster all genes into families using Diamond (% identity $> 80\%$)

4.  Construct a presence/absence matrix of gene families for organisms

5.  Estimate gene family indel rates separately for the CRISPR and
    non-CRISPR containing genomes using the R package markophylo

6.  Construct a species tree using sections from a whole genome
    alignment of all genomes of the genus

    1.  Create whole genome alignment with parsnp

    2.  Filter all aligned regions that are uninformative (identical
        sequence for every genome)

    3.  Concatenate all aligned regions into 1 alignment as a nexus file

    4.  Build the tree using Mr Bayes (10000 generations, 25% burn in
        rate)

7.  Construct the gene trees ($\leq 1500$)

    1.  Only consider families with a gene belonging in at least 40% of
        the genomes analyzed (ex: a family with 6 genes in 6 of 15
        genomes)

    2.  Align each family using mafft with default settings

    3.  Build a tree for each alignment using Mr Bayes (10000
        generations, 25% burn in)

8.  Create 1000 subsets of 50 gene trees through bootstrap sampling

9.  For each subset, use the program HiDe to infer a phylogenetic
    network from the species tree and the 50 gene trees.

10. Annotate each network with CRISPR data scraped from the CRISPR-one
    database.

11. Using the gene indel rates estimated and the annotated networks
    examine if there are any trends or effects on the dynamics of
    [hgt]{acronym-label="hgt" acronym-form="singular+short"} between
    organisms with and without CRISPR-Cas systems.

Data Collection
---------------

Protein and nucleotide fastas of all CDS sequences and whole genome
fastas for all bacterial strains included in this analysis were
downloaded from NCBI RefSeq. Only bacterial data from genomes annotated
as complete by RefSeq were used. Data was downloaded using custom
scripts available in the github repo
`./scripts/scrapingData/{NCBINameScraper.py,get_gbffs.sh}`.
[crsp]{acronym-label="crsp" acronym-form="singular+short"} annotations
of [cas]{acronym-label="cas" acronym-form="singular+short"} and Cfp
proteins from the [crsp]{acronym-label="crsp"
acronym-form="singular+short"}one tool from Zhang and Ye were used to
determine the presence of [crsp]{acronym-label="crsp"
acronym-form="singular+short"} systems (Zhang and Ye 2017). Data from
their online results were parsed and downloaded using a custom python
script available in the github repo
`./scripts/scrapingData/crispr_one_scraper.py` An
[otu]{acronym-label="otu" acronym-form="singular+short"} was annotated
as having a [crsp]{acronym-label="crsp" acronym-form="singular+short"}
system if it was annotated as having a [crsp]{acronym-label="crsp"
acronym-form="singular+short"} system and at least 1
[crsp]{acronym-label="crsp" acronym-form="singular+short"} array
sequence according to the CRISPROne database (Zhang and Ye 2017).

Gene Presence/Absence Matrix
----------------------------

In order to use the program markophylo to estimate indel rates, a
[pa]{acronym-label="pa" acronym-form="singular+short"} matrix of gene
families and organisms and a species tree are required. First any genes
classified as [mge]{acronym-label="mge" acronym-form="singular+short"}s
(from NCBI annotations) are removed. Next genes are grouped into
families by reciprocal BLAST hits and single link clustering. The
remaining unclassified genes are compared to the NCBI non-redundant
database with BLAST to check if they are genes, and if they are then
they are considered their own family with one member. The
[pa]{acronym-label="pa" acronym-form="singular+short"} matrix is
constructed as follows, for each [otu]{acronym-label="otu"
acronym-form="singular+short"} a binary vector is created, where each
entry represents a gene family and a 1 indicates that that
[otu]{acronym-label="otu" acronym-form="singular+short"} contains at
least 1 gene in that family. This is repeated for all
[otu]{acronym-label="otu" acronym-form="singular+short"}s, creating a
$G \times O$ binary matrix, where $G$ is the total number of gene
families and $O$ is the number [otu]{acronym-label="otu"
acronym-form="singular+short"}s.

There are many ways to construct a species tree, but for this project
the tree will be constructed with 16S rRNA genes, using Bayesian
methods, as implemented in the program MrBayes.

Makophylo Rate Estimations
--------------------------

Given a species tree and a gene family [pa]{acronym-label="pa"
acronym-form="singular+short"} matrix for the [otu]{acronym-label="otu"
acronym-form="singular+short"}s of the species tree the R package
*markophylo* can provide gene insertion and gene deletion rate estimates
(Dang and Golding 2016). The presence or absence of gene families are
considered 2 discrete states, for which a $(2\times 2)$ transition rate
matrix (of a [mc]{acronym-label="mc" acronym-form="singular+short"}
model) can be estimated using maximum likelihood techniques. The values
in this estimated transition matrix are the insertion rate (transition
probability of gene absence $\to$ presence) and deletion rate
(transition probability of gene presence $\to$ absence) (Dang and
Golding 2016).

Network Construction
--------------------

Quartet decomposition is method by which [hgt]{acronym-label="hgt"
acronym-form="singular+short"} events can be identified using a set of
gene trees and a species tree. Given a tree $T$ a quartet is a subtree
contain 4 of the leaf nodes in $T$, meaning that for a tree with $N$
leaf nodes (or [otu]{acronym-label="otu"
acronym-form="singular+short"}s) there are $\binom{N}{4}$ unique
quartets in that tree. A quartet $Q$ is considered consistent with a
tree if $Q = T|Le(Q)$ where $T|Le(Q)$ is the tree obtained by
suppressing all degree-two nodes in $T[X]$ and $T[X]$ is the minimal
subtree of T with all nodes in $X$, which is a leaf set of $T$ (Bansal
et al. 2013). To calculate the weight of an edge for the network, given
a species tree $S$ and a set of gene trees $G$ (Bansal et al. 2013):

1.  Pick a horizontal edge $H = ((u,v),(v,u))$ from $S$

2.  Pick a gene tree $G_i$ in $G$

3.  Decompose $G_i$ into it's set of quartets $\phi_i$

4.  Remove all quartets from $\phi_i$ either consistent with $S$ or
    discordant with $S$ but accounted for previously by a quartet set
    from another tree $G_j \in G$

5.  Set $RS((u,v),\phi_i)$ to be the number of quartets in $\phi_i$ that
    support the existence of edge $(u,v)$

6.  Set $NS((u,v),\phi_i) = \frac{RS((u,v),\phi_i)}{\lambda}$, where
    $\lambda$ is the total number of quartets in $S$ that are consistent
    with the existence of edge $(u,v)$.

7.  The score for the edge $H$ for tree $G_i$ is
    $max\{NS((u,v),\phi_i),NS((v,u),\phi_i)\} \in [0,1]$

8.  The total score for the edge $H$ is the sum of scores for each tree
    $G_i \in G$

This total score calculation is repeated for each horizontal edge $H_i$
in S, resulting in a list of edges, which is a complete description of
the network. This is further explained in the original work, (Bansal et
al. 2013).

Network Statistics
------------------

All networks will be comprised of nodes representing
[otu]{acronym-label="otu" acronym-form="singular+short"}s and weighted
edges represent the estimated amount of [hgt]{acronym-label="hgt"
acronym-form="singular+short"} events between the two incident nodes.
Bootstrap replicates of a network for each genus can be computed by
creating a set of networks from a set of random samples of gene trees
and the same species tree. For each genus 1000 bootstrap replicate
networks were produced, each with 50 gene trees sampled randomly from
all gene trees produced for that genus. The population of gene trees for
a genus is all gene families that contained genes such that at least
$40\%$ of all [otu]{acronym-label="otu" acronym-form="singular+short"}s
in the genus had 1 gene member of the family. If the genus Proteus
contains 20 [otu]{acronym-label="otu" acronym-form="singular+short"}s
(genomes) then all gene trees (up to a max of 3000) of Proteus genes
with at least $0.4*20=8$ different [otu]{acronym-label="otu"
acronym-form="singular+short"}s are sampled from. In addition extra
samples are computed beyond the 1000 replicates if there are gene trees
that were not sampled at least once in any sample i.e sampling is
continued until the union of all samples is equal to the set of all gene
trees produced, althought there were no cases of additional sampling
required. Given a network, with a set of nodes $V = \{V_0$ ...$V_i\}$ of
cardinality $N$ and a set of weighted edges (an unordered 2-tuple and
weight) $T = \{((V_1,V_2),W_{1,2})$ ...$((V_i,V_j),W_{i,j})\}$ with
cardinality $E$ descriptive statistics can be computed as follows
(Newman 2003):

-   **Average Node Degree**: $\frac{1}{|N_u|}\sum_{uv}^{N_u} w_{uv}$
    where $N_u$ is the set of nodes incidenent to $u$

-   **Average Edge Weight**: $\frac{1}{N_c}\sum_i w_i$, The average edge
    weight for all nodes with [crsp]{acronym-label="crsp"
    acronym-form="singular+short"} or without
    [crsp]{acronym-label="crsp" acronym-form="singular+short"}

-   **Node Clustering
    Coefficient**:$\frac{1}{k_u(k_u-1)} \sum_{vw}^{T(u)} (\hat{w}_{uw} \hat{w}_{vw} \hat{w}_{uv})^{\frac{1}{3}}$
    where $T(u)$ is the set of traingles containing $u$ (Onnela et
    al. 2005)

-   **Node Assortativity**: $A = \frac{Tr(M)-||M^2||}{1-||M^2||}$ Where
    $M$ is the mixing matrix of a given attribute and $||M||$ is the sum
    of all elements of $M$. $A \in [-1,1]$.(Newman 2002)

-   **Network Modularity**:
    $Q=\frac{1}{2m}\sum_{uv}^W [W_{uv} - \frac{k_u k_v}{2m}]\delta(u,v)$
    where $m$ is the total weight of alledges, $k_u$ is the degree of
    $u$ and $\delta(u,v)$ is 1 if $u$ and $v$ both have or do not have
    [crsp]{acronym-label="crsp" acronym-form="singular+short"} systems
    and 0 otherwise. $Q \in [-1,1]$ (Newman 2004)

For the statistics estimated separately for [crsp]{acronym-label="crsp"
acronym-form="singular+short"} and Non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} nodes, the average and standard error was
calculated across all replicates. All network statistics were calculated
using the networkx python package, except for modularity which was
calculated using my own implementation as no modularity function was
available. All subsequent statistical analysis was done using base R and
all figures were produced using ggplot2.

Results
=======

**Note**: The term indel refers to gene insertion/deletion events. It is
impossible to tell between two [otu]{acronym-label="otu"
acronym-form="singular+short"}s if a gene was deleted from one or
inserted in the other, leading to an ambiguity. Also a
[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"} is an
[otu]{acronym-label="otu" acronym-form="singular+short"} that I
annotated as having a [crsp]{acronym-label="crsp"
acronym-form="singular+short"} system as described in the methods
section.

In Figure [\[net\]](#net){reference-type="ref" reference="net"} it
appears that several nodes have weak connections with most other nodes
but strong connections with a few nodes. Further both
[crsp]{acronym-label="crsp" acronym-form="singular+short"} and
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"} nodes
both show distributions of strong and weak connections with other
[crsp]{acronym-label="crsp" acronym-form="singular+short"} and
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"} nodes
both. Also the standard error of each edge appears proportional to it's
weight. This is likely due to the sampling, as if more genes were
transferred along an edge, the more likely some of those genes were left
out of any individual bootstrap sample, as the size of each bootstrap
sample was $\frac{50}{376}$ of the total number of gene trees.

What is immediately clear from Figure
[\[degplot\]](#degplot){reference-type="ref" reference="degplot"} is
that the node degrees are much more similar than the mean edge weights
for the [crsp]{acronym-label="crsp" acronym-form="singular+short"} and
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"}s. This is also
supported by preforming a Wilcoxon Rank Sign test comparing the
[crsp]{acronym-label="crsp" acronym-form="singular+short"} and
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"} means.
For the mean node degree the Wilcoxon p value is 0.07412 (n=197 genera).
For the mean edge weight the Wilcoxon p value is 0.001981 (n=197
genera). So there appears to be a significant difference between the
mean edge weight for all edges but not for the degree.

Figure [\[is\]](#is){reference-type="ref" reference="is"} shows the
estimated gene indel rates for the [crsp]{acronym-label="crsp"
acronym-form="singular+short"} and non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} branch partitions for 128 genera. The
Wilcoxon signed rank test gives a p-value of $0.594$, pointing to a lack
of a significant difference in indel rates for
[crsp]{acronym-label="crsp" acronym-form="singular+short"} and
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"}s. This is
despite the fact that the slope of this fitted trend line is only
0.39464, but this may be because the trend line is influenced by
specific outliers.

Figure [\[cfrd\]](#cfrd){reference-type="ref" reference="cfrd"} show
that the fraction of all [otu]{acronym-label="otu"
acronym-form="singular+short"}s in a genus with a
[crsp]{acronym-label="crsp" acronym-form="singular+short"} system
appears to be independent of the gene indel rate. Neither the
[crsp]{acronym-label="crsp" acronym-form="singular+short"} or
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"} show a
relationship with [crsp]{acronym-label="crsp"
acronym-form="singular+short"} fraction, having regression line slope
and p values of 6.099065e-06, 0.9956342 and -0.0003156696, 0.7553688
respectively.

[\[asso\_mod\]]{#asso_mod label="asso_mod"}

Figure [\[asso\_mod\]](#asso_mod){reference-type="ref"
reference="asso_mod"} **A)** shows that the network modularity is
centered near $0$ for most networks, implying that
[crsp]{acronym-label="crsp" acronym-form="singular+short"} and
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"}s do not form
distinct communities. There does appear to be some relationship between
[crsp]{acronym-label="crsp" acronym-form="singular+short"} fraction and
modularity, where the closer the [crsp]{acronym-label="crsp"
acronym-form="singular+short"} fraction is to 0.5 the lower the
modularity is. This mkes intuitive sense, as if more members of the
network share an attribute, the more terms are included in the sum in
the definition of modularity. Figure
[\[asso\_mod\]](#asso_mod){reference-type="ref" reference="asso_mod"}
**B)** shows distribution of network assortativity is centered near
$0.10$ for most networks, implying non-assortative mixing between
[crsp]{acronym-label="crsp" acronym-form="singular+short"} and
non-[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"}s. Assortativity
is usually defined in terms of a node's degree but here node similarity
is defined as whether 2 nodes are both [crsp]{acronym-label="crsp"
acronym-form="singular+short"} or both non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"}. This wold suggest that there is no clear
separation in how [crsp]{acronym-label="crsp"
acronym-form="singular+short"} and non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s edges are organized in a network, with
some specific exceptions.

Conclusion
==========

This work highlights the large degree of variability in
[hgt]{acronym-label="hgt" acronym-form="singular+short"} rate between
bacterial genera. While there seems to be a broad association of
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems
with decreased rates of [hgt]{acronym-label="hgt"
acronym-form="singular+short"} compared to [otu]{acronym-label="otu"
acronym-form="singular+short"}s without such systems, prominent
exceptions exist. Further [crsp]{acronym-label="crsp"
acronym-form="singular+short"} or non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s do not appear to transfer genes
preferentially to either [crsp]{acronym-label="crsp"
acronym-form="singular+short"} or non-[crsp]{acronym-label="crsp"
acronym-form="singular+short"} [otu]{acronym-label="otu"
acronym-form="singular+short"}s, respectively. This lack of clear,
significant effects for individual [otu]{acronym-label="otu"
acronym-form="singular+short"}s may suggest population level effects of
[crspc]{acronym-label="crspc" acronym-form="singular+short"} systems on
[hgt]{acronym-label="hgt" acronym-form="singular+short"}. This is
further supported by the negative relationship between indel rates for
Non-[crsp]{acronym-label="crsp" acronym-form="singular+short"}
[otu]{acronym-label="otu" acronym-form="singular+short"}s and the
fraction of [otu]{acronym-label="otu" acronym-form="singular+short"}s in
a genus with a [crspc]{acronym-label="crspc"
acronym-form="singular+short"} system. One of the simplest such
explanations is that if there are more [crspc]{acronym-label="crspc"
acronym-form="singular+short"} systems in a population, they are
decreasing the overall level of phage and free-floating DNA affecting
all [otu]{acronym-label="otu" acronym-form="singular+short"}s.
Ultimately, this pipeline provides a fairly straightforward way to study
trends in [hgt]{acronym-label="hgt" acronym-form="singular+short"} for a
set of bacterial [otu]{acronym-label="otu"
acronym-form="singular+short"}s. Clearly the dynamics of
[crsp]{acronym-label="crsp" acronym-form="singular+short"} and
[hgt]{acronym-label="hgt" acronym-form="singular+short"} warrent further
investigation and should be studied within individual genera such as
Streptomyces, which have been model systems for studying
[crspc]{acronym-label="crspc" acronym-form="singular+short"} previously.
There are multiple ways to expand this analysis to answer other
questions related to the transfer of genes. Ways to expand the work
shown here include inferring transfer direction, using continuous
estimation of [crspc]{acronym-label="crspc"
acronym-form="singular+short"} activity, incorportaing ecological data
or comparing transfer for functional categories of genes.

::: {#refs .references}
::: {#ref-hgtcost}
Baltrus, David A. 2013. "Exploring the Costs of Horizontal Gene
Transfer." *Trends in Ecology and Evolution* 28 (8): 489--95.
<https://doi.org/https://doi.org/10.1016/j.tree.2013.04.002>.
:::

::: {#ref-hide}
Bansal, Mukul S., Guy Banay, Timothy J. Harlow, J. Peter Gogarten, and
Ron Shamir. 2013. "Systematic Inference of Highways of Horizontal Gene
Transfer in Prokaryotes." *Bioinformatics* 29 (5): 571--79.
<https://doi.org/10.1093/bioinformatics/btt021>.
:::

::: {#ref-acqorres}
Bondy-Denomy, J., and A. R. Davidson. 2014. "To Acquire or Resist:The
Complex Biological Effects of Crispr-Cas Systems." *Trends Microbio.* 22
(4): 218--25. <https://doi.org/10.1016/j.tim.2014.01.007>.
:::

::: {#ref-marko}
Dang, Utkarsh J., and G. Brian Golding. 2016. "Markophylo: Markov Chain
Analysis on Phylogenetic Trees." *Bioinformatics* 32 (1): 130--32.
<https://doi.org/10.1093/bioinformatics/btv541>.
:::

::: {#ref-conjug}
Davison, J. 1999. "Genetic Exchange Between Bacteria in the
Environment." *Plasmid* 42: 73--91.
:::

::: {#ref-hospital}
Dzidic, Senka, and Vladimir Bedeković. 2003. "Horizontal Gene
Transfer-Emerging Multidrug Resistance in Hospital Bacteria." *Acta
Pharmacologica Sinica* 24 (6): 519---526.
<http://www.chinaphar.com/1671-4083/24/519.htm>.
:::

::: {#ref-crisprlgt}
Godde, James S., and Amanda Bickerton. 2006. "The Repetitive Dna
Elements Called Crisprs and Their Associated Genes: Evidence of
Horizontal Transfer Among Prokaryotes." *Journal of Molecular Evolution*
62 (6): 718--29. <https://doi.org/10.1007/s00239-005-0223-z>.
:::

::: {#ref-transd}
Griffiths, A. J. F., S. R. Wessler, R. C. Lewontin, W. M. Gelbart, D. T.
Suzuki, and J. H. Miller. 2000. *An Introduction to Genetic Analysis
$7^{th}$ Edition*. W.H. Freeman.
:::

::: {#ref-crispdb}
Grissa, I. and Drevet, C. and Couvin, D. 2017. "CRISPRdb."
<http://crispr.i2bc.paris-saclay.fr/>.
:::

::: {#ref-casguild}
Haft, D. H., J. Selengut, E. F. Mongodin, and K. E. Nelson. 2005. "A
guild of 45 CRISPR-associated (Cas) protein families and multiple
CRISPR/Cas subtypes exist in prokaryotic genomes." *PLoS Comput. Biol.*
1 (6): e60.
[https://doi.org/ https://doi.org/10.1371/journal.pcbi.0010060](https://doi.org/ https://doi.org/10.1371/journal.pcbi.0010060).
:::

::: {#ref-fastlane}
Hao, W., and G. B. Golding. 2006. "The fate of laterally transferred
genes: life in the fast lane to adaptation or death." *Genome Res.* 16
(5): 636--43.
:::

::: {#ref-netoflife}
Kunin, V., L. Goldovsky, N. Darzentas, and C. A. Ouzounis. 2005. "The
net of life: reconstructing the microbial phylogenetic network." *Genome
Res.* 15 (7): 954--59.
:::

::: {#ref-staphlim}
Marraffini, Luciano A., and Erik J. Sontheimer. 2008. "CRISPR
Interference Limits Horizontal Gene Transfer in Staphylococci by
Targeting Dna." *Science* 322 (5909): 1843--5.
<https://doi.org/10.1126/science.1165771>.
:::

::: {#ref-adaevo}
Marri, P. R., W. Hao, and G. B. Golding. 2007. "The role of laterally
transferred genes in adaptive evolution." *BMC Evol. Biol.* 7 Suppl 1:
S8.
:::

::: {#ref-hgtrate}
Mozhayskiy, Vadim, and Ilias Tagkopoulos. 2012. "Horizontal Gene
Transfer Dynamics and Distribution of Fitness Effects During Microbial
in Silico Evolution." *BMC Bioinformatics* 13 (10): S13.
<https://doi.org/10.1186/1471-2105-13-S10-S13>.
:::

::: {#ref-netstat}
Newman, M. 2003. "The Structure and Function of Complex Networks." *SIAM
Review* 45 (2): 167--256. <https://doi.org/10.1137/S003614450342480>.
:::

::: {#ref-newmanmix}
Newman, M. E. 2002. "Assortative mixing in networks." *Phys. Rev. Lett.*
89 (20): 208701.
:::

::: {#ref-modularity}
---------. 2004. "Analysis of weighted networks." *Phys Rev E Stat
Nonlin Soft Matter Phys* 70 (5 Pt 2): 056131.
:::

::: {#ref-clustering}
Onnela, J. P., J. Saramaki, J. Kertesz, and K. Kaski. 2005. "Intensity
and coherence of motifs in weighted complex networks." *Phys Rev E Stat
Nonlin Soft Matter Phys* 71 (6 Pt 2): 065103.
:::

::: {#ref-trendbs}
Popa, Ovidiu, and Tal Dagan. 2011. "Trends and Barriers to Lateral Gene
Transfer in Prokaryotes." *Current Opinion in Microbiology* 14 (5):
615--23. <https://doi.org/https://doi.org/10.1016/j.mib.2011.07.027>.
:::

::: {#ref-crispgen}
Rath, Devashish, Lina Amlinger, Archana Rath, and Magnus Lundgren. 2015.
"The Crispr-Cas Immune System: Biology, Mechanisms and Applications."
*Biochimie* 117: 119--28.
<https://doi.org/https://doi.org/10.1016/j.biochi.2015.03.025>.
:::

::: {#ref-ihgt}
Ravenhall, Matt, Nives Škunca, Florent Lassalle, and Christophe
Dessimoz. 2015. "Inferring Horizontal Gene Transfer." *PLoS
Computational Biology* 11 (5): 1--16.
<https://doi.org/10.1371/journal.pcbi.1004095>.
:::

::: {#ref-nonvspacer}
Shmakov, Sergey A., Vassilii Sitnik, Kira S. Makarova, Yuri I. Wolf,
Konstantin V. Severinov, and Eugene V. Koonin. 2017. "The Crispr Spacer
Space Is Dominated by Sequences from Species-Specific Mobilomes." Edited
by Michael S. Gilmore, Rotem Sorek, and Rodolphe Barrangou. *mBio* 8
(5). <https://doi.org/10.1128/mBio.01397-17>.
:::

::: {#ref-selfcrisp}
Stern, Adi, Leeat Keren, Omri Wurtzel, Gil Amitai, and Rotem Sorek.
2010. "Self-Targeting by Crispr: Gene Regulation or Autoimmunity?"
*Trends in Genetics* 26 (8): 335--40.
<https://doi.org/https://doi.org/10.1016/j.tig.2010.05.008>.
:::

::: {#ref-hgterr}
Than, C., D. Ruths, H. Innan, and L. Nakhleh. 2007. "Confounding factors
in HGT detection: statistical error, coalescent effects, and multiple
solutions." *J. Comput. Biol.* 14 (4): 517--35.
:::

::: {#ref-transhgt}
Watson, Bridget N. J., Raymond H. J. Staals, and Peter C. Fineran. 2018.
"CRISPR-Cas-Mediated Phage Resistance Enhances Horizontal Gene Transfer
by Transduction." Edited by Joseph Bondy-Denomy and Michael S. Gilmore.
*mBio* 9 (1). <https://doi.org/10.1128/mBio.02406-17>.
:::

::: {#ref-mtrate}
Wielgoss, Sebastien, Jeffrey E. Barrick, Olivier Tenaillon, Michael J.
Wiser, W. James Dittmar, Stephane Cruveiller, Beatrice Chane-Woon-Ming,
Claudine Medigue, Richard E. Lenski, and Dominique Schneider. 2013.
"Mutation Rate Dynamics in a Bacterial Population Reflect Tension
Between Adaptation and Genetic Load." *Proceedings of the National
Academy of Sciences* 110 (1): 222--27.
<https://doi.org/10.1073/pnas.1219574110>.
:::

::: {#ref-athena}
Zambelis, A., U. J. Dang, and G. B. Golding. 2015. "Effects of
Crispr-Cas System Presence on Lateral Gene Transfer Rates in Bacteria."
:::

::: {#ref-ineqcas}
Zhang, Quan, and Yuzhen Ye. 2017. "Not All Predicted Crispr--Cas Systems
Are Equal: Isolated Cas Genes and Classes of Crispr Like Elements." *BMC
Bioinformatics* 18 (1): 92. <https://doi.org/10.1186/s12859-017-1512-4>.
:::

::: {#ref-climtrans}
Zhang, Yan, Nadja Heidrich, Biju Joseph Ampattu, Carl W. Gunderson, H.
Steven Seifert, Christoph Schoen, Jörg Vogel, and Erik J. Sontheimer.
2013. "Processing-Independent Crispr Rnas Limit Natural Transformation
in Neisseria Meningitidis." *Molecular Cell* 50 (4): 488--503.
<https://doi.org/https://doi.org/10.1016/j.molcel.2013.05.001>.
:::

::: {#ref-lgt}
Zhaxybayeva, Olga, and W. Ford Doolittle. 2011. "Lateral Gene Transfer."
*Current Biology* 21 (7): R242--R246.
<https://doi.org/https://doi.org/10.1016/j.cub.2011.01.045>.
:::
:::

