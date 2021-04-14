
### This model is to identify prognostic radiomics features using clinical information and radiomics features extracted from T2-weighted MRI imagings.

###  The required R packages

library(ggpubr)

library(reshape2)

library("survival")

library("survminer")

library(dplyr)

library(purrr)

library(optparse)

library(data.table)


### The main steps of the pipeline:

#### Train RadioML
Usage: .\RadioML_main.R [options]
Options:
        -f FILE, --file=FILE
                Choose the features data file

        -p PREFIX, --prefix=PREFIX
                Make the prefix of the result, default is prefix

        -t PARTIAL, --partial=PARTIAL
                The ratio of split the samples into train and test (integer: 1-9, default 5)

        -n NUMPERM, --numperm=NUMPERM
                The number of random permutations  (integer, default 100)
##### Examples
Rscript .\RadioML_main.R -f demo_feature_data.txt -p prefix -t 5 -n 1000

#### Extract robust radiomics features and the corresponding coefficients
Usage: .\RadioML_features.R [options]
Options:
        -p PREFIX, --prefix=PREFIX
                make the prefix of the result, default is prefix

        -t PARTIAL, --partial=PARTIAL
                The ratio of split the samples into train and test (integer: 1-9, default 5)

        -n NUMPERM, --numperm=NUMPERM
                The number of random permutations  (integer, default 100, it should be the same with RadioML_main.R option -n)

        -R RDATA, --RData=RDATA
                the previous main outout RData
##### Examples
Rscript .\RadioML_features.R -p TIANTAN -t 5 -n 1000 -R TIANTAN_1000_5_training.RData


#### Calculate the Risk score of the samples
Usage: .\Risk_score.R [options]
Options:
        -f FILE, --file=FILE
                Choose the feature and samples dataframe

        -p PREFIX, --prefix=PREFIX
                Make the prefix of the result, default is prefix

        -c FEATURES, --features=FEATURES
                Choose the feature and coe data

        -m METHODS, --methods=METHODS
                Marke the previous methods
##### Examples
Rscript .\Risk_score.R -f demo_feature_data.txt -p prefix -c perm5.830.feature_coe.txt -m perm5.830
