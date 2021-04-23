# syndRomics 0.0.6

**Bug fix:** 

* variable order for plotting CI error vars did not match variable reordering based on the var_order argument in the barmap communality plots.

* add new catch error in permutation functions for repeated run failure. Now if there are 3 consecutive failures at running the pca would stop.

* there was a bug preventing some call of princals to finish as it will get stuck in an loop during permutation. It has been fixed.

# syndRomics 0.0.5

- Bug fix of resampling functions for objects of class "princals". Data value was not pass properly to Gifi::princals()

# syndRomics 0.0.4

- Catching errors for passed number of dimensions bigger than present in dataset. Now the funciton for component similarity only returns the minimal number of comparisons possible, with a warning if this number is lower than the ndim specified by the user

- add examples

# syndRomics 0.0.3

- Bug fix of resampling functions for objects of class "princals"
