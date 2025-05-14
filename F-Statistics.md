# Introduction to F3- and F4-Statistics

f-statistics (Patterson et al., 2012; Reich, Thangaraj, Patterson, Price, & Singh, 2009) are a widely used 
toolkit for making inferences about phylogeny and admixture from population genetic data, particularly in humans. 
The statistics measure correlations in allele frequencies among sets of two, three or four populations. 
Observed values reflect degrees of shared ancestry and can serve as a means for testing hypotheses regarding 
population split orders and past gene flow events under historical models. As compared to some other common methods 
in population genetics, f-statistics are quite simple and flexible, but interpreting them is not always 
straightforward.

## F3 Statistics
F3 statistics are a useful analytical tool to understand population relationships. F3 statistics, just as F4 
and F2 statistics measure allele frequency correlations between populations and were introduced by 
Nick Patterson ([Patterson et al. 2012](https://doi.org/10.1534/genetics.112.145037)), but see also 
([Peter 2016](https://doi.org/10.1534/genetics.115.183913)) for another introduction.

F3 statistics are used for two purposes: 
i) as a test whether a target population (C) is admixed between two source populations (A and B), and 
ii) to measure shared drift between two test populations (A and B) from an outgroup (C).

F3 statistics are in both cases defined as the product of allele frequency differences between population C to A 
and B, respectively:

_F3(A,B;C)=⟨(c−a)(c−b)⟩_

Here,⟨.⟩ denotes the average over all genotyped sites, and a, b and c denote the allele frequency for a given 
site in the three populations A, B and C.

### Admixture - F3 Statistics
It can be shown that if F3(A,B;C)  is negative, it provides unambiguous proof that population C is admixed between 
populations A and B, as in the following phylogeny (taken from Figure 1 from (Patterson et al. 2012):

![](aDNA_popgen_analyses/images/f3_phylogeny.png)

Intuitively, an F3 statistics becomes negative if the allele frequency of the target population C is on average intermediate between the allele frequencies of A and B. 
Consider as an extreme example a genomic site where α = 0, b = 1 and c = 0.5. Then we have (c−a)(c−b)=−0.25, which is negative. So if the entire statistics is negative, 
it suggests that in many positions, the allele frequency c is indeed intermediate, suggesting admixture between the two sources.
NOTE: If an F3 statistics is not negative, it does not proof that there is no admixture!
NOTE: Admixture F3 cannot be calculated if the target (C) is represented by a single individual with pseudohaploid data (no variation across individuals within a population can be detected)

One way to understand this is by looking what happens to a list of SNPs and allele frequencies for groups A, B and C:

| SNP  | A | B | C | (c−a)(c−b) |
|------|---|---|---|------------|
| #1   | 1 | 0 | 0.5 | -0.25          |
| #2   | 0.8 | 0 | 0 | 0          |
| #3   | 0 | 0.7 | 0 | 0         |
| #4   | 0.1 | 0.5 | 0.3 | -0.04          |
| #5   | 0 | 0.1 | 0.2 | 0.02         |
| #6   | 1 | 0.2 | 0.9 | -0.07          |
| **F3(A,B;C)** |   |   |   |   | **-0.057** |

Every SNP where C has an allele frequency intermediate between A and B contributes negatively. Here, the average is also negative, providing evidence for admixture. 
For statistical certainty, an error bar for this estimate is needed, which is typically computed via Jackknife.

### Outgroup - F3 Statistics
Outgroup F3 statistics are a special case how to use F3 statistics. The definition is the same as for Admixture F3 statistics, but instead of a target C and two source 
populations A and B, one now gives an outgroup C and two test populations A and B. To get an intuition for this statistics, consider the following tree:

FIGURE

In this scenario, the statistic F3(A, B; C) measures the branch length from C to the common ancestor of A and B, coloured red. So this statistic is simply a measure of how 
closely two population A and B are related with each other, as measured from a distant outgroup. It is thus a similarity measure: The higher the statistic, the more genetically 
similar A and B are to one another.

It is helpful to consider an example using a SNP list, this time assuming that every population is just a single (haploid) individual, so each allele frequency can just be 0 or 1. 
For example:

| SNP  | A | B | C | (c−a)(c−b) |
|------|---|---|---|------------|
| #1   | 1 | 0 | 0 | 0          |
| #2   | 1 | 0 | 0 | 0          |
| #3   | 0 | 0 | 1 | 1         |
| #4   | 1 | 0 | 1 | 0          |
| #5   | 1 | 1 | 0 | 1         |
| #6   | 0 | 0 | 1 | 1          |
| **F3(A,B;C)** |   |   |   |   | **0.5** |

You can see that each position which is similar between A and B, but different to C contributes 1, all other SNPs 0. So it directly measures similarity between A and B on 
alleles that differ from the outgroup C.

### Computing F3-statistics

You will need to install admixtools as a conda environment:
```
conda activate
conda create -n admixtools 
conda install bioconda::admixtools #v7.0.2 or v5.1, preferebly the second one for qpAdm
conda activate admixtools
```
Prepare the poplist following this order:

```
Source1 Source2 Target
```
or
```
PopA PopB Outgroup
```

Run the following script:
```
in1=$1 #geno and snp files prefix
in2=$2 #ind file with the pop label to use
in3=$3 #poplist

#run qp3
qp3Pop -p <(echo "genotypename:        ${in1}.geno
snpname:        ${in1}.snp
indivname:      ${in2}
popfilename: ${in3}") \
> ${in3}.qp3Pop.out

#print results for R pltting
grep 'result:' ${in3}.qp3Pop.out | awk '{print $2, $3, $4, $5, $6, $7, $8, $9}' > ${in3}.qp3Pop.R.out
```

The results (before cleanup) have the following format:

```
result:   Source1  Source2   Target f_3  std.err Z 
```

Plot using RStudio:
```
library("tidyverse")
library(gridExtra)

setwd("PATH/TO/YOUR/DIRECTORY")

#read in data
f3dat = read.table("f3.R.qp3Pop.out",
                   col.names=c("PopA", "PopB", "PopC", "F3", "StdErr", "Z", "SNPs"))

#order by F3 depending on the Pop you are most interested
Ordered <- arrange(f3dat,F3) 

#set pop order as factor so ggplot won't re-order automatically
Ordered$PopB <- factor(Ordered$PopB, levels = Ordered$PopB) 

# ggplot
ggplot(data = Ordered) + 
  geom_vline(xintercept = 0, linetype = 3) + #use if admixture F3
  geom_errorbar(aes(y=PopB, xmin=F3-StdErr, xmax=F3+StdErr, #plot error bars
                    width=0)) + # change whisker length
  geom_point(aes(x=F3, y=PopB, fill = abs(Z)), #plot F3 points
             stroke=1, size = 4, shape=21) + #set point size
  theme_light() + #set theme
  labs(x = "f3",
       y = "test", title = "f3(PopA, test; Outgroup)") + 
  geom_text(aes(x=F3, y=PopB, label=SNPs, hjust=0.1, vjust=-0.95), size = 3) +
  scale_fill_gradient(low = "blue", high = "yellow") 
```

For more information and parameters: https://github.com/DReichLab/AdmixTools

## F4 Statistics
A different way to test for admixture is by “F4 statistics” (or “D statistics” which is very similar), also introduced in (Patterson et al. 2012).
F4 statistics are also defined in terms of correlations of allele frequency differences, similarly to F3 statistics (see above), but involving four 
different populations, not just three. Specifically we define:

_F4(A,B;C,D)=⟨(a−b)(c−d)⟩_

Consider the following tree:

FIGURE

In this tree, without any additional admixture, the allele frequency difference between A and B should be completely independent from the allele 
frequency difference between C and D. In that case, F4(A, B; C, D) should be zero, or at least not statistically different from zero. 
However, if there was gene flow from C or D into A or B, the statistic should be different from zero. Specifically, if the statistic is significantly negative, 
it implies gene flow between either C and B, or D and A. If it is significantly positive, it implies gene flow between A and C, or B and D.

The way this statistic is often used, is to put a divergent outgroup as population A, for which we know for sure that there was no admixture into either C or D. 
With this setup, we can then test for gene flow between B and D (if the statistic is positive), or B and C (if it is negative).

### ABBA-BABA sites
It is helpful to consider an example using a SNP list, this time assuming that every population is just a single (haploid) individual, so each allele frequency can just be 0 or 1. 
For example:

| SNP  | A | B | C | D | (a−b)(c−d) |
|------|---|---|---|---|------------|
| #1   | 1 | 0 | 0 | 0 | 0          |
| #2   | 1 | 0 | 1 | 1 | 0          |
| #3   | 0 | 1 | 1 | 0 | -1         |
| #4   | 0 | 1 | 0 | 1 | 1          |
| #5   | 1 | 0 | 0 | 1 | -1         |
| #6   | 1 | 0 | 0 | 0 | 0          |
| **F4(A,B;C,D)** |   |   |   |   | **-0.0167** |

the only SNPs that contribute positively to this statistics are SNPs where the alleles are distributed as 1010 and 0101, and the only SNPs that contribute negatively are 1001 and 0110. 
In the literature, the two patterns have been dubbed “ABBA” and “BABA”, which is why the statistical test behind this statistic (see below) was sometimes called the ABBA-BABA test 
(see for example ([Martin, Davey, and Jiggins 2015](https://academic.oup.com/mbe/article/32/1/244/2925550)).

In positions that are polymorphic in both (A,B) and (C,D), this statistic asks whether B is genetically more similar to C than it is to D. This is most useful as a test for “treeness”: 
If A, B, C, D are related to each other as indicated in the above tree, then B should be equally closely related to C as to D. But if we actually find evidence that B is closer to C than to D, or vice versa, 
then this means that the tree above cannot be correct, but that there must be a closer connection between B and C or B and D, depending on the sign of the statistic.

So the ABBA- and BABA-categories of SNPs help shape intuition for how this statistic behaves for single haploid genomes. But what about population allele frequencies? 
Looking back at the formula ⟨(a−b)(c−d)⟩ this doesn’t help very much with intuition how this behaves with frequencies. Well, a nice feature of F4-Statistics is that averages factor out. 
This means, that if you have multiple samples in one or multiple slots A, B, C or D, the total F4-statistic of the groups is exactly equal to the average of F4-Statistics of the individuals.
Here is a more mathematical definition.

Let’s say we have 2 individuals in each of A and B, so we may perhaps write A = {A1,A2} and B = {B1,B2}. Then one can show to have

F4(A,B;C,D) =  Average of [F4(A1,B1;C,D), F4(A1,B2;C,D), F4(A2,B1;C,D), F4(A2,B2;C,D)]

so just the average over all individual-based F4-statistics. And this can be shown to be true for arbitrary numbers of samples. 
So in other words: An F4-Statistic always measures the average excess of pairwise BABA SNPs over ABBA SNPs. 

NOTE: F4-statistics have been famously used to show that Neanderthals are more closely related to Non-African populations than to Africans, suggesting gene-flow between Neanderthals and Non-Africans 
(shown in (Green et al. 2010)). You can reproduce this famous result with F4(Chimp.REF,Altai_published.DG;Yoruba,French) and F4(Chimp.REF,Altai_published.DG;Sardinian,French) which shows that the first 
statistic is significantly positive with a Z-score of 7.99, while the second one is insignificantly different from zero (Z=1.01).

### Computing F3-statistics

```
conda activate admixtools

in1=$1
in2=$2
in3=$3

#run
qpDstat -p <(echo "genotypename: ${in1}.geno
snpname:      ${in1}.snp
indivname:    ${in2}
popfilename: ${in3}
f4mode:   YES") \
> ${in3}.qpDstat.f4.out

#clean up output file for R
grep "result" ${in3}.qpDstat.f4.out | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' > ${in3}.qpDstat.f4.R.out
```

The results (before cleanup) have the following format:

```
result:   Pop1 (W)  Pop2 (X) : Pop3 (Y)  Pop4 (Z)  D-stat	Z	BABA	ABBA	#_of_SNPs_that_all_populations_have_data
```
Plot using RStudio:
```
library("tidyverse")
library(gridExtra)

setwd("PATH/TO/YOUR/DIRECTORY")

#read in data
f4dat = read.table("qpDstat.f4.R.out",
                   col.names=c("Pop1", "Pop2", "Pop3", "Pop4", "Dstat", "Z", "BABA", "ABBA", "SNPs"))

#order by Dstat
Ordered <- arrange(f4dat,Dstat) 

#set pop order as factor so ggplot won't re-order automatically
Ordered$Pop3 <- factor(Ordered$Pop3, levels = Ordered$Pop3) 

# ggplot
ggplot(data = Ordered) + 
  geom_vline(xintercept = 0, linetype = 3) +
  geom_errorbar(aes(y=Pop3, xmin = Dstat - (Dstat / Z), xmax = Dstat + (Dstat / Z),
                    width=0)) + # change whisker length
  geom_point(aes(x=Dstat, y=Pop3, fill = abs(Z)), shape=21, #plot Dstat points
             size = 4) + #set point size
  theme_light() + #set theme
  labs(x = "f4",
       y = "test", title = "f4(Pop1, Pop2; Pop3, Pop4)") + 
  geom_text(aes(x=Dstat, y=Pop3, label=SNPs, hjust=0.2, vjust=-1.7), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high="red", midpoint = 3) 
```
	
## REFERENCES
- https://compvar-workshop.readthedocs.io/en/latest/contents/03_f3stats/f3stats.html
