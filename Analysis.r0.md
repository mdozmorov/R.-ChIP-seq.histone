





Analysis of histone marks after filtering out reads overlapping low complexity regions
=======================================================================================

About 70% of such reads were removed.

Differentially expressed peaks, as identified by MACS2
------------------------------------------------------
MACS2 identified differential peaks by comparing KO and WT conditions directly, without considering input signal.




* First, we check the number of differentially bound peaks.


```
## [1] "============================================"
## [1] "Condition"         "KO_vs_WT_H3K27me3"
## [1] "Number of Peaks:" "931"             
## [1] "============================================"
## [1] "Condition"        "KO_vs_WT_H3K4me2"
## [1] "Number of Peaks:" "2"               
## [1] "============================================"
## [1] "Condition"     "KO_vs_WT_H3K9"
## [1] "Number of Peaks:" "599"
```


Conclusion: Main differences are observed in H3K27me3 and H3K9 conditions.

* Second, we check if these peaks. overlap.

<img src="img/venn1.png" title="plot of chunk venn1" alt="plot of chunk venn1" width="700" />


Conclusion: Differentially bound peaks are largely non-overlapping.

* Third, we annotate the peaks by their relation to genes, and save the results in files.


```
## Adding symbol ... done
## prepare output ... done
```

```
## Adding symbol ... done
## prepare output ... done
```

```
## Adding symbol ... done
## prepare output ... done
```


* Fourth, we check where the peaks occur most frequently, in relation to genes.

<img src="img/pie1.png" title="plot of chunk pie1" alt="plot of chunk pie1" width="700" />


<img src="img/pie3.png" title="plot of chunk pie3" alt="plot of chunk pie3" width="700" />


Conclusion: Majority of the peaks are located upstream of genes, some are insight, some are downstream

* Fifth, we look at the enriched gene ontologies and KEGG pathways that may be affected by genes where differential peaks are located

Test files contain more detailed information about the enrichments, including which genes belong to which enriched category.

No multiple testing is applied, to be more inclusive. Caveat emptor.

KO_vs_WT_H3K27me3.GO
---------------------
<img src="img/enrichedGO1.png" title="plot of chunk enrichedGO1" alt="plot of chunk enrichedGO1" width="700" />


KO_vs_WT_H3K27me3.KEGG
----------------------
<img src="img/enrichKEGG1.png" title="plot of chunk enrichKEGG1" alt="plot of chunk enrichKEGG1" width="700" />


KO_vs_WT_H3K9.GO
-----------------
<img src="img/enrichedGO3.png" title="plot of chunk enrichedGO3" alt="plot of chunk enrichedGO3" width="700" />


KO_vs_WT_H3K9.KEGG
--------------------

```
## Error: No enriched pathway can be found.
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'as.matrix': Error in eval(expr, envir, enclos) : 
##   object 'KO_vs_WT_H3K9.KEGG' not found
## Calls: cbind -> standardGeneric -> eval -> eval -> eval
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unique': Error: object 'KO_vs_WT_H3K9.KEGG' not found
```








