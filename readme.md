
### these scripts is used to extract the prognostic features from the features and samples matrix based on the survival infomation.

###  R packages

library(ggpubr)

library(reshape2)

library("survival")

library("survminer")

library(dplyr)

library(purrr)

library(optparse)

library(data.table)


### the main steps of the pipeline:

#### train the models repeat hundrands times
Usage: .\RadioML_main.R [options]
Options:
        -f FILE, --file=FILE
                Choose the features data file

        -p PREFIX, --prefix=PREFIX
                make the prefix of the result, default is prefix

        -t PARTIAL, --partial=PARTIAL
                the ratio of split the samples into train and test (integer: 1-9, default 5)

        -n NUMPERM, --numperm=NUMPERM
                the repeat times  (integer, default 100)
##### examples
Rscript .\RadioML_main.R -f demo_feature_data.txt -p prefix -t 5 -n 1000

#### get the high frequency featurs and related coefficient
Usage: .\RadioML_features.R [options]
Options:
        -p PREFIX, --prefix=PREFIX
                make the prefix of the result, default is prefix

        -t PARTIAL, --partial=PARTIAL
                the ratio of split the samples into train and test (integer: 1-9, default 5)

        -n NUMPERM, --numperm=NUMPERM
                the repeat times  (integer, default 100, it should be same with RadioML_main.R option)

        -R RDATA, --RData=RDATA
                the previous main outout RData
##### examples
Rscript .\RadioML_features.R -p TIANTAN -t 5 -n 1000 -R TIANTAN_1000_5_training.RData


#### predict the Risk score of samples
Usage: .\Risk_score.R [options]
Options:
        -f FILE, --file=FILE
                Choose the feature and samples dataframe

        -p PREFIX, --prefix=PREFIX
                make the prefix of the result, default is prefix

        -c FEATURES, --features=FEATURES
                Choose the feature and coe data

        -m METHODS, --methods=METHODS
                marker the previous methods
##### examples
Rscript .\Risk_score.R -f demo_feature_data.txt -p prefix -c perm5.830.feature_coe.txt -m perm5.830
