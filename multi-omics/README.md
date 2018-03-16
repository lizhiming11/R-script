The project is a multi-omics analysis of the relevant scripts.

    The script adonis_each.r and adonis_each.r are calculated to affect size.
    
    The operation is as follows:
        1: run adonis_each.r  -> Calculate the impact of each variable on the matrix.
            The following parameters:
                -i : input affect file
                -d : input affected file
                -m : Select the test distance matrix for the calculation method.
                -p : permutations number
                -o : output file
        2: run effect_size.r ->  For the effects previously obtained, get rid of redundant,to calculate the total impact size.
            The following parmeters:
                -i : input affect file
                -d : input affected file
                -m : Select the test distance matrix for the calculation method.
                -p : permutations number
                -o : output file
                -a : input each variable adonis result(Results of script 1)
                -p_v: adonis screening conditions(Results of script 1),default = 0.05
                -c_v:  get rid of redundant screening conditions(pearson correlation coefficient),default = 0.5



