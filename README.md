This repository contains a collection of function that makes it easier to use mixed models in DMRcate.

There are two functions in this repository: `cpg.annotate_mixed_model` and `cpg.annotate_dupp_cor`.

# `cpg.annotate_mixed_model`

This function extends `cpg.annotate` and allows for mixed model. Most of the parameters have identical effect as those of `cpg.annotate`.
the new parameter `formula` must be a mixed model formula (e.g. `~x+(1|id)`).

`design` must be a design matrix however without the random effect part. 

For example: `designMatrix <- model.matrix(~x, data=dataset)`

`meta_data` must be a data frame having all the variables that are included in the model (e.g. `x` and `id`). 

So far, the `differential` and `ANOVA`  analysis type has been implemented. 

# `cpg.annotate_dupp_cor`

This is based on `duplicateCorrelation` in `limma`. 
The only new argument is `block` which is the variable showing the hierarchical structure of the data. For example `id` in the example above.
Similar to above, the `design` must be a design matrix however without the block part. 


# Example
```R
formula_fix<-~x+(1|id)
designMatrix <- model.matrix(~x, data=sample.sheet) 
contMatrix <- makeContrasts(AvsB=A-B, levels=designMatrix)
	
myAnnotation <- annotate_mix_model(object = GRset,formula = ~x+(1|id),
meta_data =dataset , datatype = "array", analysis.type="differential", design=designMatrix,
arraytype="450K", contrasts = TRUE, cont.matrix = contMatrix,
coef="AvsB", fdr = 0.05)
                               
myAnnotation_2 <- cpg.annotate_dupp_cor(object = GRset, datatype = "array", 
analysis.type="differential", design=designMatrix,
arraytype="450K", contrasts = TRUE, cont.matrix = contMatrix,
coef="AvsB", fdr = 0.05,block=dataset[,"id"])                            

```
# Notes

`cpg.annotate_mixed_model` is slow when running on a large data set. It can take days to finished. If mixed models are not strickly needed, `cpg.annotate_dupp_cor`.
