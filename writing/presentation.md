---
nocite: '@*'
theme: Ilmenau
colortheme: dolphin
linkcolor: blue
header-includes:
- \setbeamercolor{block title}{bg=blue!50}
- \setbeamercolor{block body}{bg=blue!10}
- \setbeamertemplate{itemize item}{\scriptsize$\blacktriangleright$}
- \setbeamertemplate{itemize subitem}{\scriptsize$\diamond$}
- \setbeamertemplate{enumerate items}{\insertenumlabel.}
- \setbeamertemplate{section page}{\begin{centering} \usebeamerfont{section title}\insertsection\par\end{centering}}
- \usepackage{tikz-network}
- \usepackage{hyperref}

title: Is Sharing Caring? 
subtitle: How CRISPR-Cas systems affect rates of Horizonal Gene Transfer
author: Siddharth Reed
date: \today
institute: Golding Lab, McMaster University
---

# Background

## What is Horizontal Gene Transfer?

:::: {.columns}
::: {.column width="45%"}
![ @trendslgt](hgt_mechanisims_trendslgt.png){height="70%"}
:::
::: {.column width="60%"}
\vspace{0.25cm}
- Transformation
  - Incorporation of free-floating DNA into the genome
- Conjugation
  - Transfer of DNA through cell-cell connections
- Transduction
  - Transfer of DNA through phage
- Rates of HGT can be impacted by
  - Amount of exogenous DNA/cell density/phage density
  - Selective pressures
  - Metabolic costs
  - Sequence compatibility
:::
::::

## What is A Network?
:::: {.columns}
::: {.column width="45%"}
![ @bondy](graphExample_bondy.png)
:::
::: {.column width="55%"}
- Useful mathematical abstraction of real world system
- Nodes represent objects
- Edges represent relationships
- Nodes and edges can have attributes
- __Nodes are OTUs, edges are inferred HGT rates__
:::
::::

## What are CRISPR-Cas systems?
:::: {.columns}
::: {.column width="45%"}
![ @crispgen](CRISPR-immunity_crispgen.jpg){height=80%}
:::
::: {.column width="65%"}
\vspace{0.25cm}
- Adaptive immune system in bacteria
- Failed “infection” $\to$ spacer acquisition $\to$ targeted degradation for next “infection”
- Protects against foreign DNA absorption/integration 
- Requires CRISPR loci + Cas proteins 
- $45\%$ of bacteria have CRISPR loci $(n=6782)$ [@crispdb]
:::
::::

## Do CRISPR Systems Affect Horizontal Gene Transfer?

\centering
\huge Yes

---

- Gophna et al. 2015 found no relation between the presence of CRISPR systems and HGT over short evolutionary timescales
  - Assume all singletons arose from HGT
  - Used GC% to identify HGT
- Contradicted by a former undergraduate thesis student [@athena]
  - Can see inhibitory effects of CRISPR on HGT over short evolutionary time scales
  - Higher gene indel rates for CRISPR containing OTUs than non-CRISPR containing outgroups

## CRISPR Cost Complexity and Curbing It 

### Fitness Cost Factors
- Metabolic maintenance [@crispgen]
- Off-target effects [@selfcrisp]
- Environmental pressures [@hospital]
- Phage density [@acqorres]
- Anti-CRISPR systems [@acqorres]
- Prophage abundance [@transhgt]

\onslide<6->
### Cost Reduction Strategies
- Selective CRISPR inactivation [@crispgen]
- CRISPRs get transferred $\implies$ population level immunity [@crisprlgt]
- CRISPR can enhance transduction-mediated HGT [@transhgt]

# My Project

## Goals

### Within Network Comparisons
For genera with CRISPR containing OTUs, compare the node statistics of CRISPR containing OTUs to non-CRISPR containing OTUs.

\onslide<2->
### Gene Indel Rates vs. Network Statistics
Compare gene Indel rates to node/network statistics for CRISPR containing and non-CRISPR containing OTUs

## InDel Rate Estimates

:::: {.columns}
::: {.column width="45%"}
![](partition_example.png)
:::
::: {.column width="60%"}

### Markophylo [@marko]
- Input: a presence/absence matrix of gene families + species tree
- Branches are partitioned by the OTU having a CRISPR system
- Internal branches ignored
- Gene birth-death rates are estimated for each branch partition 
:::
::::

## Network Sampling
![](netsample.png)

## Network Statistics

- __Mean Node Degree:__ $\frac{1}{|N_u|}\sum_{uv}^{N_u} w_{uv}$ where $N_u$ is the set of nodes incident to $u$
- __Node Clustering Coefficient:__$\frac{1}{k_u(k_u-1)} \sum_{vw}^{T(u)} (\hat{w}_{uw} \hat{w}_{vw} \hat{w}_{uv})^{\frac{1}{3}}$ where $T(u)$ is the set of triangles containing $u$ [ @clustering]
- __Network Assortativity:__ $A = \frac{Tr(M)-||M^2||}{1-||M^2||}$ Where $M$ is the mixing matrix of a given attribute and $||M||$ is the sum of all elements of $M$. $A \in [-1,1]$. [ @newmanmix]
- __Network Modularity:__ $Q=\frac{1}{2m}\sum_{uv}^W [W_{uv} - \frac{k_u k_v}{2m}]\delta(u,v)$ where $m$ is the total edge weight, $k_u$ is the degree of $u$ and $\delta(u,v)$ is 1 if $u$ and $v$ both have or do not have CRISPR systems and 0 otherwise. $Q \in [-1,1]$ [ @modularity]
<!-- -  __Average Edge Weight:__ $\frac{1}{N_c}\sum_i w_i$, The average edge weight for all nodes with or without CRISPRs -->
<!-- -  __Node Eigenvector Centrality:__ $\frac{N-1}{\sum_v d(u,v)}$ where $d(x,y)$ is the length of the shortest path $v \to u$. [ @egcen] -->

## Workflow (per genus)

\centering
\begin{tikzpicture}
    \SetDistanceScale{2.5}
    \SetVertexStyle[Shape=rectangle,LineOpacity=0.0,FillOpacity=0.8,MinSize=1.75\DefaultUnit]
    \onslide<2->  \Vertex[x=0,y=2,color=green]{ns}
                  \Text[x=0,y=2]{Get Fastas}
    \onslide<3->  \Vertex[x=0,y=1]{n2}
                  \Text[x=0,y=1,width=1.75cm]{Filter MGEs}
                  \only<3->{\Edge[Direct](ns)(n2)}
    \onslide<4->  \Vertex[x=0,y=0]{n3}
                  \Text[x=0,y=0,width=1.75cm]{Cluster Genes}
                  \only<4->{\Edge[Direct](n2)(n3)}
    \onslide<5->  \Vertex[x=1.5,y=1]{n41}
                  \Text[x=1.5,y=1,width=1.75cm]{Species Tree}
                  \only<5->{\Edge[Direct](n3)(n41)}
    \onslide<6->  \Vertex[x=1,y=0]{n40}
                  \Text[x=1,y=0,width=1.75cm]{P/A Matrix}
                  \only<6->{\Edge[Direct](n3)(n40)}
    \onslide<7->  \Vertex[x=2,y=0]{n50}
                  \Text[x=2,y=0,width=1.75cm]{Gene InDel Rates}
                  \only<7->{\Edge[Direct](n40)(n50)}
                  \only<7->{\Edge[Direct](n41)(n50)}
    \onslide<8->  \Vertex[x=1,y=2]{n42}
                  \Text[x=1,y=2,width=1.75cm]{Gene Trees}
                  \only<8->{\Edge[Direct](n3)(n42)}
    \onslide<9->  \Vertex[x=2,y=2]{n52}
                  \Text[x=2,y=2,width=1.75cm]{Bootstrap Network Replicates}
                  \only<9->{\Edge[Direct](n42)(n52)}
                  \only<9->{\Edge[Direct](n41)(n52)}
    \onslide<10-> \Vertex[x=3,y=2]{n62}
                  \Text[x=3,y=2,width=1.75cm]{Annotate with CRISPR}
                  \only<10->{ \Edge[Direct](n52)(n62)}
    \onslide<11-> \Vertex[x=3,y=1]{n72}
                  \Text[x=3,y=1,width=1.75cm]{Compte Network Stats}
                  \only<11->{\Edge[Direct](n62)(n72)}
    \onslide<12-> \Vertex[x=3,y=0,color=red]{ne}
                  \Text[x=3,y=0,width=1.75cm]{Analyze Results}
                  \only<12->{\Edge[Direct](n72)(ne)}
                  \only<12->{\Edge[Direct](n50)(ne)}

\end{tikzpicture}

# Results

## Example “Consensus” Network
![](network.png){height=95%}

## Genus Size Distribution
![](genus_size_dist.png)

## Mean Node Degree
![](c_nc_deg_bar.png)

## Gene Indel Rates
![](c_nc_indel_bar.png)

---

![](c_nc_rate_scatter.png)

## Gene Indel Rate Vs. Fraction of CRISPR OTUs
![](crate_cfrac_scatter.png)

## Gene Indel Rate Vs. Fraction of CRISPR OTUs
![](ncrate_cfrac_scatter.png)

## Mean Node Weighted Clustering Coefficient
![](c_nc_clust_scatter.png)

## Assortativity Distributions
![](asst_violin.png)

## Modularity Distributions
![](mod_violin.png)

## Indel Rate Pair Plot
\centering
![](pairplot.png){height="75%"}

# Moral of the Study

## Findings
- Large variation in HGT rate between genera.
- CRISPR systems $\sim$ associated with lower HGT rates 
  - Prominent exceptions exist
  - High mixing between CRISPR and non-CRISPR OTUs
- Population level effects of CRISPR systems may decrease HGT rates
- Interplay of CRISPR systems and HGT is complex and warrants further study

## Possible Future Directions
- __Intergenic comparisons:__ What if we analyze networks with multiple genera e.g. ones that share a microbiome?
- __Inferring direction:__ Inferring direction of transfer $\to$ more analytic tools available
- __CRISPR Label:__ binary (presence) $\to$ continuous (activity) e.g. array length, transciptomic data, etc.
- __Gene function analysis:__ How are HGT dynamics different for different functional groups?
- __Transfer of CRISPR systems:__ How do CRISPR systems get transfered around?

## Conclusion
\centering

\Huge Is Sharing Caring?

\vspace{0.2in}

\Large
- \onslide<2-> Yes, for researchers
- \onslide<3-> Jury's still out for bacteria

## Thanks

:::: {.columns}
::: {.column width="45%"}
Thank you to

> - Dr. G. Brian Golding
> - Dr. Ben Evans
> - The Golding lab
>   - Caitlin Simopoulos
>   - Daniella Lato
>   - Zachery Dickson
>   - Sam Long
>   - George Long
>   - Lucy Zhang
>   - Brianne Laverty
>   - Nicole Zhang
> - Everyone here for listening
:::
::: {.column width="60%"}
![](mcmaster_logo.png)

### Code Availability

All code written written for this project is available at
[https://github.com/DJSiddharthVader/Undergrad_Thesis](https://github.com/DJSiddharthVader/Undergrad_Thesis)

:::
::::

## Bibliography {.allowframebreaks}
