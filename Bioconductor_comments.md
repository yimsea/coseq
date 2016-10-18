### Grammar

1. Avoid using 1:number for creating a sequence, use the more robust method of seq_len(number)

2. Don’t overwrite / assign object names that are already function (i.e., I <- which(label == k)

3. Use descriptive argument name (rather than x use coseq)

4. Use S4Methods rather than S3

5. Use is(object, “class”) to test for class.

6. Some functions look like they’re meant to dispatch using S4 methods

7. Nested lists may necessitate an S4 structure to avoid using chain $ calls
 
8. Avoid using periods or under_scores in function names

9. if(is.null(arg.user$Kmin.init)) arg.user$Kmin.init <- “small-em” -- Seems like you should have a named list where you can input values 
based on is.null results.

10. No need to use ; semi-colons

11. Minor: use which.max, paste0 exists, if (length(conds)) vs if (length(conds) > 0)

12. Avoid the use of @ in extractor functions -- create accessor functions/methods

13. Don’t need logic like if (transpose == TRUE) or if (is.numeric(K) == TRUE)

14. Why so many options for plotting? If you do need all of those consider the data structure and how it is supposed to facilitate your plots.

15. Don’t use functions as variables (i.e., “c” in {for c in sequence})

### Vignette

1. Could you provide a compilable vignette for the package? It would be better for users to be able to follow and run the functions that 
you provide in the package interactively.

### Structure

1. It would be better to use available S4 structures for your class (i.e., using DataFrame instead of data.frame)

2. You can combine related functions into one file.
    I can't imagine that there are no Bioconductor classes/structures available for your use. We would prefer that you re-use any existing functionality.
