.. |eufund| image:: _static/eu_co-funded.png
  :width: 200
  :alt: Alternative text

.. |chfund| image:: _static/ch-logo-200x50.png
  :width: 210
  :alt: Alternative text

.. |ukrifund| image:: _static/ukri-logo-200x59.png
  :width: 150
  :alt: Alternative text

.. |hmsc1| image:: _static/hmsc1.png
  :width: 650
  :alt: Alternative text

.. |output_icon| image:: _static/output_icon.png
  :width: 50
  :alt: Alternative text

.. |logo_BGE_small| image:: _static/logo_BGE_alpha.png
  :width: 120
  :alt: Alternative text
  :target: https://biodiversitygenomics.eu/

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


Hierarchical Modelling of Species Communities
*********************************************

Hierarchical Modelling of Species Communities (HMSC) is a framework for Joint Species Distribution Modelling; 
a model-based approach for analyzing community 
ecological data (`Ovaskainen et al.2017 <https://doi.org/10.1111/ele.12757>`_).


The **obligatory data** for HMSC-analyses includes a matrix of **species occurrences or abundances** and a 
matrix of **environmental covariates** (sampling units as **rows**). Optionally additional input data include species 
**traits** and **phylogeny**, and information about the spatiotemporal context of the 
sampling design. HMSC yields inference both at species and community levels. 

|hmsc1|


Starting point 
~~~~~~~~~~~~~~

The input for HMSC analysis consists of species/OTU community matrix (Y matrix) accompanied with 
the matrix of environmental covariates (X matrix), and optionally species 
traits (T) matrix and phylogeny (C) matrix.

.. admonition:: Community and Environment (metadata)

  *example of input Y (community) matrix*:

  +-------------+--------+--------+--------+--------+-----+
  |             | OTU_01 | OTU_02 | OTU_03 | OTU_04 | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample1** | 579    | 405    | 0      | 0      | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample2** | 0      | 345    | 0      | 62     | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample3** | 0      | 449    | 231    | 345    | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample4** | 0      | 430    | 69     | 0      | ... |
  +-------------+--------+--------+--------+--------+-----+

  *example of input X (environmental metadata) matrix*:

  +-------------+---------+-----------------+----------+-----------+-----+
  |             | site    | collection_date | latitude | longitude | ... |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample1** | site_01 | 2022_07_15      | 0.1      | 0.2       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample2** | site_01 | 2022_07_21      | 0.1      | 0.2       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample3** | site_02 | 2022_07_28      | 1.0      | 1.0       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample4** | site_02 | 2022_08_05      | 1.0      | 1.0       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+


.. admonition:: Traits

  Including the traits matrix helps to understand how species-specific traits influence the community assebly processes. 
  Traits matrix may include data for example about body size, shape, feeding type, etc. 

  +-------------+-----------+----------+---------------+
  |             | body_size | shape    | trophic_guild |
  +-------------+-----------+----------+---------------+
  | **OTU_01**  | 2         | square   | 1             |
  +-------------+-----------+----------+---------------+
  | **OTU_02**  | 2         | round    | 1             |
  +-------------+-----------+----------+---------------+
  | **OTU_03**  | 0.1       | round    | 2             |
  +-------------+-----------+----------+---------------+
  | **OTU_04**  | 0.2       | variable | 3             |
  +-------------+-----------+----------+---------------+


.. admonition:: Phylogeny

  Including the phylogeny matrix helps to understand if the species responses to the environmental covariates are phylogenetically structured, i.e., do similar species  respond similarly.

  The phylogeny matrix may be presented as **Newick tree** file. 
  Alternatively, data on **taxonomic ranks** may used as a proxy of phylogenetic relatedness. 
  In the latter case, the UNCLASSIFIED taxonomic ranks *(that are common in the metabarcoding data)* 
  should be informative   in a sense that 
  **not all e.g. family level unclassified OTUs would be considered as closely related** because 
  they have the 'same label'.   For example, the distance between unclassified OTUs could 
  be calculated and then OTUs falling within the user defined distance threshold could be 
  classified as various levels of *pseudotaxa*. 

  *Example of the taxonomy table:*

  +------------+-----+------------+------------------+-----------------+---------------+
  |            | ... | Class      | Order            | Family          | Genus         |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_01** | ... | Collembola | Entomobryomorpha | Entomobryidae   | Entomobrya    |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_02** | ... | Collembola | Entomobryomorpha | **Family_0032** | **Genus_001** |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_03** | ... | Collembola | Entomobryomorpha | Isotomidae      | Parisotoma    |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_04** | ... | Collembola | Entomobryomorpha | **Family_0032** | **Genus_022** |
  +------------+-----+------------+------------------+-----------------+---------------+

  Here, for example the sequences of **OTU_02 and OTU_04 differ 9.8%**; so we consider those are originating from different genera but likely from the same family. 

___________________________________________________

Install HMSC
~~~~~~~~~~~~

Install hmsc R package (if already not installed).

.. code-block:: R
   :caption: install hmsc R package  
   :linenos:

   #!/usr/bin/Rscript

   # install 'devtools'; if not yet installed
   install.packages("devtools") 
   library(devtools)
   
   # install hmsc package
   install_github("hmsc-r/HMSC")

   # load hmsc
   library(Hmsc)

   # check the version
   packageVersion("Hmsc")

___________________________________________________

.. _select_data:

Select data
~~~~~~~~~~~

In this 'select data' section, we are assuming that 
our input Y and X matrixes are tab delimited text files where 
samples are in rows; and taxonomy table format follows the above example.

Here, we are also deciding if we want to proceed with the **presence-absence or abundance** (read count)
community matrix.

.. note:: 

  Before using to the full dataset, try fitting the model with a **small subset**
  for faster model testing and validation.


.. code-block:: R
   :caption: select data and prep. data
   :linenos:

   #!/usr/bin/Rscript

   # input matrices file names; according to the above HMSC figure.
      # here, expecting all files to be tab-delimited.
   Y_file = "OTU_table.txt" # Community; samples are rows
   X_file = "env_meta.txt"  # Environment; samples are rows
   C_file = "taxonomy.txt"  # Phylogeny (optional) 
                            # (herein a 'pseudo-phylogeny' based on the 
                            # assigned taxonomic ranks; species are rows)
   T_file = ""              # Traits (optional)
   sp_prevalence = 100      # set a species prevalence threshold; 
                            # meaning that perform HMSC for species 
                            # in a Y-matrix, that occur >= 100 samples in this case. 
    #-------------------------------------------------------------------------------------------#

    # load community matrix Y
    Y = read.table(Y_file, sep = "\t", 
                      check.names = F, 
                      header = T, 
                      row.names = 1)
    # load metadata matrix X
    X = read.table(X_file, sep = "\t", 
                        check.names = F, 
                        header = T, 
                        row.names = 1)
    # load taxonomy matrix C; 
      # so we can use it as a proxy for phylogeny 
    taxonomy = read.table(C_file, sep = "\t", 
                          check.names = F, 
                          header = T, 
                          row.names = 1)

    # herein, converting Y matrix to presence-absence
    Y = 1*(Y>0)

    # select OTUs/species that are present at least in $sp_prevalence (specified above) samples
    prevalence = colSums(Y != 0)
    select.sp = prevalence >= sp_prevalence
    Y = Y[, select.sp]
    taxonomy = taxonomy[select.sp, ]

    ### creating a phylogenetic tree from a set of nested taxonomic ranks in the taxonomy table
     # ranks as.factors; assuming that ranks start from the 2nd column in the taxonomy table 
        # and we have 7 ranks
    for(i in 2:8) taxonomy[,i] = as.factor(taxonomy[,i])
     # convert tax to phylo tree
    phy.tree = as.phylo(~Phylum/Class/Order/Family/Genus/Species, 
                    data = taxonomy, collapse = FALSE)
     
    # this "pseudo-tree" does not have any branch lengths, which are needed for the model;
    # assign arbitrary branch lengths
      if (is.null(phy.tree$edge.length)) {
        phy.tree$edge.length = rep(1, nrow(phy.tree$edge))
      }
     # rename tree tip lables according to the labels in the Y matrix.
     phy.tree$tip.label = colnames(Y)

     # check the tree 
     plot(phy.tree, cex=0.6)

___________________________________________________

.. _define_model: 

Define model
~~~~~~~~~~~~

According to our dataset, we are defining the model for HMSC.
That is, we specify the structure, including the response variable (community data), 
covariates (environmental predictors), random effects, and phylogenetic relationships. 

.. code-block:: R
   :caption: define model
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)

    # defining our study design; structure of the data
    studyDesign = data.frame(
                    sample = as.factor(rownames(X)),    # rownames(X) = sample names
                    site = as.factor(X$Site)
                    )

    ### incorporating random effects into the  HMSC model. 
     # (to capture the influence of unmeasured factors that vary across different 
     # levels of the data; e.g., among sites, samples). 
    # sampling units
    rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
    # sampling sites
    rL.site = HmscRandomLevel(units = levels(studyDesign$site)) 

    # convert sample collection dates into Julian days relative to a specific start date 
    da = as.Date(meta$CollectionDate)
    jday = 1 + julian(da) - julian(as.Date("2024-01-01"))

    # not covered here, but DOWNLOAD RELEVANT COVARIATES (E.G., CLIMATE, WEATHER, LANDCOVER) 
    # FROM DATABASES BASED ON COORDINATES AND SAMPLING TIMES

    # create a data frame of covariates (predictor variables) that will be included in the model
    XData = data.frame(seqdepth = log(meta$seq_depth),  # number of sequences per sample
                        jday)                           # Julian days

    ### specify the formula for the fixed effects
    # using 3.141593 instead of pi to prevent issues when 'pi' is considered as a variable; 
    XFormula = ~cos(2*3.141593*jday/365) +  # model seasonal effects; annual cycles
                sin(2*3.141593*jday/365) +  # model seasonal effects; annual cycles
                cos(4*3.141593*jday/365) +  # model seasonal effects; semiannual cycles
                sin(4*3.141593*jday/365) +  # model seasonal effects; semiannual cycles
                seqdepth                    # number of sequences per sample

    ### define a model 
    m = Hmsc(Y = Y,             # response matrix
          distr = "probit",      # distribution model for the response variable ('probit' for PA)
          XData = XData,          # predictor variables 
          XFormula = XFormula,     # fixed effects in the model
          phyloTree = phy.tree,     # phylogenetic tree object
          studyDesign = studyDesign, # study design object
          ranLevels = list(sample = rL.sample, site = rL.site)) # random level objects
    
    # organize, name, and save your HMSC models (to easily manage multiple models if needed)
    models = list(m)
    names(models) = c("model_1")
    save(models, file = paste0("models/unfitted_models.RData"))

    # check models
    models

___________________________________________________

.. _fit_model:

Fit model
~~~~~~~~~

The following HMSC pipeline is a modified version of the pipeline presented at the `ISEC 2024 Hmsc workshop <https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc>`_.

.. admonition:: input:
  
  Unfitted_models file (**unfitted_models.RData**)

.. admonition:: output:
  
  | Fitted models, with fitting done for multiple RUNs.
  | First short MCMC chains are for **quick, preliminary checks** to ensure the model is running correctly. Then the fitting is done with increasingly longer MCMC chains (until MCMC convergence or computational limit is reached).
  | **Output files**:
  | *models**/models_thin_1_samples_5_chains_4.Rdata* (RUN 0),
  | *models/models_thin_1_samples_250_chains_4.Rdata* (RUN 1),
  | *models/models_thin_10_samples_250_chains_4.Rdata* (RUN 2), ...


.. code-block:: R
   :caption: fit models
   :linenos:

   #!/usr/bin/Rscript

   # set working directory 
   localDir = "."
   # specify 'models' dir (input/output dir) 
   modelDir = file.path(localDir, "models")

   # load input (unfitted_models)
   load(file=file.path(modelDir,"unfitted_models.RData"))

   # number of parallel chains (CPUs) for MCMC sampling
   nParallel = NULL   # if NULL, then nParallel = nChains

   # number of samples and thinning intervals for each iteration
   samples_list = c(5, 250, 250, 250, 250, 250)   # number of iterations of the MCMC
   thin_list = c(1, 1, 10, 100, 1000, 10000)      # thinning interval;  keep every k-th sample 
   # number of MCMC chains
   nChains = 4
   #-------------------------------------------------------------------------------------------#

   library(Hmsc)
   
   # iterate over sample and thin lists to fit models with different configurations
   nm = length(models)  # number of models (for loop)
   if(is.null(nParallel)) nParallel = nChains
   Lst = 1
   while(Lst <= length(samples_list)){
     thin = thin_list[Lst]
     samples = samples_list[Lst]
     print(paste0("thin = ",as.character(thin),"; samples = ", as.character(samples)))
     filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                               "_samples_", as.character(samples),
                                               "_chains_", as.character(nChains),
                                               ".Rdata", sep = ""))
     if(file.exists(filename)){
       print("model had been fitted already")
     } else {
       print(date())
       for (mi in 1:nm) {
         print(paste0("model = ",names(models)[mi]))
         m = models[[mi]]
         m = sampleMcmc(m, samples = samples, thin=thin,
                       adaptNf=rep(ceiling(0.4*samples*thin), m$nr), 
                       transient = ceiling(0.5*samples*thin),
                       nChains = nChains,
                       nParallel = nParallel) 
         models[[mi]] = m
       }
       save(models, file=filename)
     }
     Lst = Lst + 1
   }

____________________________________________________

.. _evaluate_convergence:

Evaluate convergence
~~~~~~~~~~~~~~~~~~~~

Evaluating convergence ensures the accuracy and reliability of the model's results by 
confirming that the MCMC algorithm has adequately explored the parameter space 
and that the chains have stabilized. 

.. admonition:: input:
  
  Fitted models directory ``models``. 
  **Check the input path**, and typically other parts of the scripts do not need modifications.

.. admonition:: output:
  
  MCMC convergence statistics for selected model parameters, 
  illustrated (for all runs) in the file "results/MCMC_convergence.pdf", 
  and the text file "results/MCMC_convergence.txt".


.. code-block:: R
   :caption: evaluate convergence
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)
   library(colorspace)
   library(vioplot)

   # make the script reproducible
   set.seed(1)

   # set working directory 
   localDir = "."
   # specify 'models' dir (input dir) 
   modelDir = file.path(localDir, "models")
   # specify 'results' dir (output dir) 
   resultDir = file.path(localDir, "results")
   if (!dir.exists(resultDir)) dir.create(resultDir)

   # number of samples and thinning intervals for each iteration (AS ABOVE)
   samples_list = c(5, 250, 250, 250, 250, 250)   # number of iterations of the MCMC
   thin_list = c(1, 1, 10, 100, 1000, 10000)      # thinning interval;  keep every k-th sample 
   # number of MCMC chains
   nChains = 4
   #-------------------------------------------------------------------------------------------#

   # setting commonly adjusted parameters
   showBeta = TRUE  # default = TRUE, converg. shown for beta-parameters
   showGamma = TRUE # default = FALSE, converg. not shown for gamma-parameters
   showOmega = TRUE # default = FALSE, converg. not shown for Omega-parameters
   maxOmega = 100   # default = 50, converg. of Omega shown for 50 randomly selected species pairs
   showRho = TRUE   # default = FALSE, converg. not shown for rho-parameters
   showAlpha = TRUE # default = FALSE, converg. not shown for alpha-parameters

   # set up a text file to store MCMC convergence statistics
   text.file = file.path(resultDir,"/MCMC_convergence.txt")
   cat("MCMC Convergennce statistics\n\n", file = text.file, sep="")

   ma.beta = NULL
   na.beta = NULL
   ma.gamma = NULL
   na.gamma = NULL
   ma.omega= NULL
   na.omega = NULL
   ma.alpha = NULL
   na.alpha = NULL  
   ma.rho = NULL
   na.rho = NULL
   Lst = 1
   while(Lst <= nst){
     thin = thin_list[Lst]
     samples = samples_list[Lst]
      
     filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
     if(file.exists(filename)){
       load(filename)
       cat(c("\n",filename,"\n\n"),file=text.file,sep="",append=TRUE)
       nm = length(models)
       for(j in 1:nm){
         mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
         nr = models[[j]]$nr
         cat(c("\n",names(models)[j],"\n\n"),file=text.file,sep="",append=TRUE)
         if(showBeta){
           psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
           tmp = summary(psrf)
           cat("\nbeta\n\n",file=text.file,sep="",append=TRUE)
           cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
           if(is.null(ma.beta)){
             ma.beta = psrf[,1]
             na.beta = paste0(as.character(thin),",",as.character(samples))
           } else {
             ma.beta = cbind(ma.beta,psrf[,1])
             if(j==1){
               na.beta = c(na.beta,paste0(as.character(thin),",",as.character(samples)))
             } else {
               na.beta = c(na.beta,"")
             }
           }
         }
         if(showGamma){
           psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
           tmp = summary(psrf)
           cat("\ngamma\n\n",file=text.file,sep="",append=TRUE)
           cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
           if(is.null(ma.gamma)){
             ma.gamma = psrf[,1]
             na.gamma = paste0(as.character(thin),",",as.character(samples))
           } else {
             ma.gamma = cbind(ma.gamma,psrf[,1])
             if(j==1){
               na.gamma = c(na.gamma,paste0(as.character(thin),",",as.character(samples)))
             } else {
               na.gamma = c(na.gamma,"")
             }
           }
         }
         if(showRho & !is.null(mpost$Rho)){
           psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
           cat("\nrho\n\n",file=text.file,sep="",append=TRUE)
           cat(psrf[1],file=text.file,sep="\n",append=TRUE)
         }
         if(showOmega & nr>0){
           cat("\nomega\n\n",file=text.file,sep="",append=TRUE)
           for(k in 1:nr){
             cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
             tmp = mpost$Omega[[k]]
             z = dim(tmp[[1]])[2]
             if(z > maxOmega){
               sel = sample(1:z, size = maxOmega)
               for(i in 1:length(tmp)){
                 tmp[[i]] = tmp[[i]][,sel]
               }
             }
             psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
             tmp = summary(psrf)
             cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
             if(is.null(ma.omega)){
               ma.omega = psrf[,1]
               na.omega = paste0(as.character(thin),",",as.character(samples))
             } else {
               ma.omega = cbind(ma.omega,psrf[,1])
               if(j==1){
                 na.omega = c(na.omega,paste0(as.character(thin),",",as.character(samples)))
               } else {
                 na.omega = c(na.omega,"")
               }
             }
           }
         }
         if(showAlpha & nr>0){
           for(k in 1:nr){
             if(models[[j]]$ranLevels[[k]]$sDim>0){
               cat("\nalpha\n\n",file=text.file,sep="\n",append=TRUE)
               cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
               psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
               cat(psrf[,1],file=text.file,sep="\n",append=TRUE)            
             }
           }
         }
       }
     }
     Lst = Lst + 1
   }

   pdf(file= file.path(resultDir,"/MCMC_convergence.pdf"))
   if(showBeta){
     par(mfrow=c(2,1))
     vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
     legend("topright",legend = names(models), fill=rainbow_hcl(nm))
     vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
   }
   if(showGamma){
     par(mfrow=c(2,1))
     vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
     legend("topright",legend = names(models), fill=rainbow_hcl(nm))
     vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
   }
   if(showOmega & !is.null(ma.omega)){
     par(mfrow=c(2,1))
     vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
     legend("topright",legend = names(models), fill=rainbow_hcl(nm))
     vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
   }
   dev.off()


.. admonition:: check the diagnostics

  We need to confirm that the MCMC chain has converged. If the chains have not converged, the 
  samples privided by the MCMC chain can yield biased parameter estimates. 

___________________________________________________

.. _compute_model_fit:

Compute model fit
~~~~~~~~~~~~~~~~~

Assessing the model's performance and validating its accuracy. ... 

.. admonition:: input:
  
  Fitted models directory ``models``. 
  **Check the input path**, and typically other parts of the scripts do not need modifications.

.. admonition:: output:
  
  | Model fits computed by the cross-validation, with fitting (which is part of cross-validation) done for multiple RUNs.
  | First short MCMC chains (to provide some results fast), and then with increasingly long MCMC chains (up to the longest run performed in :ref:`Fit models <fit_model>`). 
  | **Output files**:
  | *models/MF_thin_1_samples_5_chains_4.Rdata* (RUN 0),
  | *models/MF_thin_1_samples_250_chains_4.Rdata* (RUN 1),
  | *models/MF_thin_10_samples_250_chains_4.Rdata* (RUN 2), 
  | *models/MF_thin_100_samples_250_chains_4.Rdata* (RUN 3), and so on.


.. code-block:: R
   :caption: evaluate convergence
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)

   # make the script reproducible
   set.seed(1)

   # set working directory 
   localDir = "."
   # specify 'models' dir (input/output dir) 
   modelDir = file.path(localDir, "models")

   # number of samples and thinning intervals for each iteration (AS ABOVE)
   samples_list = c(5, 250, 250, 250, 250, 250)   # number of iterations of the MCMC
   thin_list = c(1, 1, 10, 100, 1000, 10000)      # thinning interval;  keep every k-th sample 
   # number of MCMC chains
   nChains = 4

   # number of parallel chains (CPUs) for MCMC sampling
   nParallel = NULL   # if NULL, then nParallel = nChains
   # number of cross-validations
   nfolds = 2         # default: (2) two-fold cross-validation
   #-------------------------------------------------------------------------------------------#

   if(is.null(nParallel)) nParallel = nChains
   Lst = 1
   while(Lst <= length(samples_list)){
     thin = thin_list[Lst]
     samples = samples_list[Lst]
     filename.in = file.path(modelDir,paste("models_thin_", as.character(thin),
                                           "_samples_", as.character(samples),
                                           "_chains_",as.character(nChains),
                                           ".Rdata",sep = ""))
     filename.out = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                             "_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             "_nfolds_", as.character(nfolds),
                                             ".Rdata",sep = ""))
     if(file.exists(filename.out)){
       print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
       print("model fit had been computed already")
     } else {
       if(file.exists(filename.in)){
         print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
         print(date())
         load(file = filename.in) #models
         nm = length(models)
        
         MF = list()
         MFCV = list()
         WAIC = list()
        
         for(mi in 1:nm){
           print(paste0("model = ",names(models)[mi]))
           m = models[[mi]]
           preds = computePredictedValues(m)
           MF[[mi]] = evaluateModelFit(hM=m, predY=preds)
           partition = createPartition(m, nfolds = nfolds) #USE column = ...
           preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
           MFCV[[mi]] = evaluateModelFit(hM=m, predY=preds)
           WAIC[[mi]] = computeWAIC(m)
         }
         names(MF)=names(models)
         names(MFCV)=names(models)
         names(WAIC)=names(models)
        
         save(MF,MFCV,WAIC,file = filename.out)
       }
     }
     Lst = Lst + 1
   }

___________________________________________________

.. _show_model_fit: 

Show model fit
~~~~~~~~~~~~~~

Evaluate and visualize the model fit for different thinning and sampling configurations.

.. admonition:: input:
  
  Fitted models directory ``models``. 
  **Check the input path**, and typically other parts of the scripts do not need modifications.

.. admonition:: output:
  
  Model fits illustrated (for highest RUN of :ref:`compute model fit <compute_model_fit>`) in the file *results/model_fit.pdf*.


.. code-block:: R
   :caption: show model fit
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)

   # make the script reproducible
   set.seed(1)

   # set working directory 
   localDir = "."
   # specify 'models' dir (input dir) 
   modelDir = file.path(localDir, "models")
   # specify 'results' dir (output dir) 
   resultDir = file.path(localDir, "results")
   if (!dir.exists(resultDir)) dir.create(resultDir)

   # number of samples and thinning intervals for each iteration (AS ABOVE)
   samples_list = c(5, 250, 250, 250, 250, 250)   # number of iterations of the MCMC
   thin_list = c(1, 1, 10, 100, 1000, 10000)      # thinning interval;  keep every k-th sample 
   # number of MCMC chains
   nChains = 4

   # number of cross-validations
   nfolds = 2         # default: (2) two-fold cross-validation
   #-------------------------------------------------------------------------------------------#

   for (Lst in nst:1) {
     thin = thin_list[Lst]
     samples = samples_list[Lst]
     
     filename = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         "_nfolds_", as.character(nfolds),
                                         ".Rdata",sep = ""))
     if(file.exists(filename)){break}
   }
   if(file.exists(filename)){
     load(filename)
     
     nm = length(MF)
     modelnames = names(MF)
     pdf(file = file.path(resultDir,paste0("/model_fit_nfolds_",nfolds,".pdf")))
     for(j in 1:nm){
       cMF = MF[[j]]
       cMFCV = MFCV[[j]]
       if(!is.null(cMF$TjurR2)){
         plot(cMF$TjurR2,cMFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
             xlab = "explanatory power",
             ylab = "predictive power",
             main=paste0(modelnames[j],", thin = ",
                         as.character(thin),
                         ", samples = ",as.character(samples),
                         ": Tjur R2.\n",
                         "mean(MF) = ",as.character(mean(cMF$TjurR2,na.rm=TRUE)),
                         ", mean(MFCV) = ",as.character(mean(cMFCV$TjurR2,na.rm=TRUE))))
         abline(0,1)
         abline(v=0)
         abline(h=0)
       }
       if(!is.null(cMF$R2)){
         plot(cMF$R2,cMFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
             xlab = "explanatory power",
             ylab = "predictive power",
             main=paste0(modelnames[[j]],", thin = ",as.character(thin),
                         ", samples = ",as.character(samples),
                         ": R2. \n",
                         "mean(MF) = ",as.character(mean(cMF$R2,na.rm=TRUE)),
                         ", mean(MFCV) = ",as.character(mean(cMFCV$R2,na.rm=TRUE))))
         abline(0,1)
         abline(v=0)
         abline(h=0)
       }
       if(!is.null(cMF$AUC)){
         plot(cMF$AUC,cMFCV$AUC,xlim=c(0,1),ylim=c(0,1),
             xlab = "explanatory power",
             ylab = "predictive power",
             main=paste0(modelnames[[j]],", thin = ",as.character(thin),
                         ", samples = ",as.character(samples),
                         ": AUC. \n",
                         "mean(MF) = ",as.character(mean(cMF$AUC,na.rm=TRUE)),
                         ", mean(MFCV) = ",as.character(mean(cMFCV$AUC,na.rm=TRUE))))
         abline(0,1)
         abline(v=0.5)
         abline(h=0.5)
       }
       if(FALSE && !is.null(cMF$O.TjurR2)){
         plot(cMF$O.TjurR2,cMFCV$O.TjurR2,xlim=c(-1,1),ylim=c(-1,1),
             xlab = "explanatory power",
             ylab = "predictive power",
             main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": O.Tjur R2"))
         abline(0,1)
         abline(v=0)
         abline(h=0)
       }
       if(FALSE && !is.null(cMF$O.AUC)){
         plot(cMF$O.AUC,cMFCV$O.AUC,xlim=c(0,1),ylim=c(0,1),
             xlab = "explanatory power",
             ylab = "predictive power",
             main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": O.AUC"))
         abline(0,1)
         abline(v=0.5)
         abline(h=0.5)
       }      
       if(!is.null(cMF$SR2)){
         plot(cMF$SR2,cMFCV$SR2,xlim=c(-1,1),ylim=c(-1,1),
             xlab = "explanatory power",
             ylab = "predictive power",
             main=paste0(modelnames[[j]],", thin = ",as.character(thin),
                         ", samples = ",as.character(samples),
                         ": SR2. \n",
                         "mean(MF) = ",as.character(mean(cMF$SR2,na.rm=TRUE)),
                         ", mean(MFCV) = ",as.character(mean(cMFCV$SR2,na.rm=TRUE))))
         abline(0,1)
         abline(v=0)
         abline(h=0)
       }    
       if(FALSE && !is.null(cMF$C.SR2)){
         plot(cMF$C.SR2,cMFCV$C.SR2,xlim=c(-1,1),ylim=c(-1,1),
             xlab = "explanatory power",
             ylab = "predictive power",
             main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": C.SR2"))
         abline(0,1)
         abline(v=0)
         abline(h=0)
       }  
     }
     dev.off()
   }

___________________________________________________

.. _show_parameter_estimates:

Show parameter estimates
~~~~~~~~~~~~~~~~~~~~~~~~

Assessing the model's performance and validating its accuracy. ... 

.. admonition:: input:
  
  Fitted models directory ``models``. 
  **Check the input path**, and typically other parts of the scripts do not need modifications.

.. admonition:: output:
  
  | Parameter estimates illustrated (for highest RUN of :ref:`Fit models <fit_model>`) in the file
  | *results/parameter_estimates.pdf*, the text file *results/parameter_estimates.txt*, 
  | as well as given numerically in multiple csv files (one per parameter type) named 
  | *results/parameter_estimates_[parameter_name].csv*.

.. code-block:: R
   :caption: show parameter estimates
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)
   library(colorspace)
   library(corrplot)
   library(writexl)

   # make the script reproducible
   set.seed(1)

   # set working directory 
   localDir = "."
   # specify 'models' dir (input dir) 
   modelDir = file.path(localDir, "models")
   # specify 'results' dir (output dir) 
   resultDir = file.path(localDir, "results")
   if (!dir.exists(resultDir)) dir.create(resultDir)

   # number of samples and thinning intervals for each iteration (AS ABOVE)
   samples_list = c(5, 250, 250, 250, 250, 250)   # number of iterations of the MCMC
   thin_list = c(1, 1, 10, 100, 1000, 10000)      # thinning interval;  keep every k-th sample 
   # number of MCMC chains
   nChains = 4

   # thresholds for the posterior support levels
   support.level.beta = 0.95       # default: 0.95
   support.level.gamma = 0.95      # default: 0.95
   support.level.omega = 0.9       # default: 0.9


   var.part.order.explained = NULL # default: in variance partitioning of explained variance, 
                                   #     species are shown in the order they are in the model.
     # -> use var.part.order.explained to select which order species are shown in the 
     #                                                       raw variance partitioning.
     # var.part.order.raw should be a list of length the number of models. 
     # for each element provide either 0 (use default);
     # or a vector of species indices;
     # or "decreasing" if you wish to order according to explanatory power
     # var.part.order.explained = list()
     # var.part.order.explained[[1]] = 0
     # var.part.order.explained[[2]] = c(2,1)

   var.part.order.raw = NULL       # default: in variance partitioning of raw variance, species 
                                   #              are shown in the order they are in the model.
     # -> use var.part.order.raw to select which order species are shown in the 
     #                                           explained variance partitioning.
     # same options apply as for var.part.order.explained
     # var.part.order.raw = list()
     # var.part.order.raw[[1]] = "decreasing"
     # var.part.order.raw[[2]] = c(1,2)

   
   show.sp.names.beta = NULL    # default: species names shown in beta plot if there are at 
                                #                          most 30 species and no phylogeny
     # -> use show.sp.names.beta to choose to show / not show species names in betaPlot
     # if given, show.sp.names.beta should be a vector with length equalling number of models
     # show.sp.names.beta = c(TRUE,FALSE) 

   plotTree = NULL              # default: tree is plotted in Beta plot if the model includes it.
     # -> use plotTree to choose to plot / not plot the tree in betaPlot.
     # if given, plotTree should be a vector with length equalling number of models
     # plotTree = c(FALSE,FALSE)

   omega.order = NULL           # default: species shown in the order they are in the model.
     # -> use omega.order to select which order species are shown in omega plots
     # omega.order should be a list of length the number of models. 
     # for each element provide either 0 (use default);
     # or a vector of species indices;
     # or "AOE" if you wish to use the angular order of the eigenvectors.
     # omega.order = list()
     # omega.order[[1]] = "AOE"
     # omega.order[[2]] = c(2,1)
   
   show.sp.names.omega = NULL   # default: species names shown in beta plot if there are at 
                                #                                           most 30 species.
     # show.sp.names.omega = c(TRUE,FALSE)
   #-------------------------------------------------------------------------------------------#


   text.file = file.path(resultDir,"/parameter_estimates.txt")
   cat(c("This file contains additional information regarding parameter estimates.","\n","\n", 
                                                                      sep=""), file=text.file)

   for (Lst in nst:1) {
     thin = thin_list[Lst]
     samples = samples_list[Lst]
     filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
     if(file.exists(filename)){break}
   }
   if(file.exists(filename)){
     load(filename)
     cat(c("\n",filename,"\n","\n"),file=text.file,sep="",append=TRUE)
     nm = length(models)
     if(is.null(var.part.order.explained)){
       var.part.order.explained = list()
       for(j in 1:nm) var.part.order.explained[[j]] = 0
     }
     if(is.null(var.part.order.raw)){
       var.part.order.raw = list()
       for(j in 1:nm) var.part.order.raw[[j]] = 0
     }
     if(is.null(omega.order)){
       omega.order = list()
       for(j in 1:nm) omega.order[[j]] = 0
     }
     
     modelnames = names(models)
     
     pdf(file= file.path(resultDir,"parameter_estimates.pdf"))
     for(j in 1:nm){
       cat(c("\n",names(models)[j],"\n","\n"),file=text.file,sep="",append=TRUE)
       m = models[[j]]
       if(m$XFormula=="~."){
         covariates = colnames(m$X)[-1]
       } else {
         covariates = attr(terms(m$XFormula),"term.labels")
       }
       if(m$nr+length(covariates)>1 & m$ns>1){
         preds = computePredictedValues(m)
         VP = computeVariancePartitioning(m)
         vals = VP$vals
         mycols = rainbow(nrow(VP$vals))
         MF = evaluateModelFit(hM=m, predY=preds)
         R2 = NULL
         if(!is.null(MF$TjurR2)){
           TjurR2 = MF$TjurR2
           vals = rbind(vals,TjurR2)
           R2=TjurR2
         }
         if(!is.null(MF$R2)){
           R2=MF$R2
           vals = rbind(vals,R2)
         }
         if(!is.null(MF$SR2)){
           R2=MF$SR2
           vals = rbind(vals,R2)
         }
         filename = file.path(resultDir, paste("parameter_estimates_VP_",modelnames[j],".csv"))
         write.csv(vals,file=filename)
         if(!is.null(VP$R2T$Beta)){
           filename = file.path(resultDir,paste("parameter_estimates_VP_R2T_Beta",modelnames[j],
                                                                                       ".csv"))
           write.csv(VP$R2T$Beta,file=filename)
         }
         if(!is.null(VP$R2T$Y)){
           filename = file.path(resultDir, paste("parameter_estimates_VP_R2T_Y",modelnames[j],
                                                                                       ".csv"))
           write.csv(VP$R2T$Y,file=filename)
         }
         if(all(var.part.order.explained[[j]]==0)){
           c.var.part.order.explained = 1:m$ns
         } else {
           if(all(var.part.order.explained[[j]]=="decreasing")){
             c.var.part.order.explained = order(R2, decreasing = TRUE)
           }
           else {
             c.var.part.order.explained  = var.part.order.explained[[j]]
           }
         }
         VP$vals = VP$vals[,c.var.part.order.explained]
         cat(c("\n","var.part.order.explained","\n","\n"),file=text.file,sep="",append=TRUE)
         cat(m$spNames[c.var.part.order.explained],file=text.file,sep="\n",append=TRUE)
         plotVariancePartitioning(hM=m, VP=VP, main = paste0("Proportion of explained variance, ",
                      modelnames[j]), cex.main=0.8, cols = mycols, 
                                      args.leg=list(bg="white",cex=0.7))
         if(all(var.part.order.raw[[j]]==0)){
           c.var.part.order.raw = 1:m$ns
         } else {
           if(all(var.part.order.raw[[j]]=="decreasing")){
             c.var.part.order.raw = order(R2, decreasing = TRUE)
           }
           else {
             c.var.part.order.raw  = var.part.order.raw[[j]]
           }
         }
         if(!is.null(R2)){
           VPr = VP
           for(k in 1:m$ns){
             VPr$vals[,k] = R2[k]*VPr$vals[,k]
           }
           VPr$vals = VPr$vals[,c.var.part.order.raw]
           cat(c("\n","var.part.order.raw","\n","\n"),file=text.file,sep="",append=TRUE)
           cat(m$spNames[c.var.part.order.raw],file=text.file,sep="\n",append=TRUE)
           plotVariancePartitioning(hM=m, VP=VPr,main=paste0("Proportion of raw variance, ",
                modelnames[j]),cex.main=0.8, cols = mycols, args.leg=list(bg="white",cex=0.7),
                                                                                  ylim=c(0,1))
         }
       }
     }
     for(j in 1:nm){
       m = models[[j]]
       if(m$nc>1){
         postBeta = getPostEstimate(m, parName="Beta")
         filename = file.path(resultDir, 
                            paste("parameter_estimates_Beta_",modelnames[j],".xlsx"))
         me = as.data.frame(t(postBeta$mean))
         me = cbind(m$spNames,me)
         colnames(me) = c("Species",m$covNames)
         po = as.data.frame(t(postBeta$support))
         po = cbind(m$spNames,po)
         colnames(po) = c("Species",m$covNames)
         ne = as.data.frame(t(postBeta$supportNeg))
         ne = cbind(m$spNames,ne)
         colnames(ne) = c("Species",m$covNames)
         vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
         writexl::write_xlsx(vals,path = filename)
         if(is.null(show.sp.names.beta)){
           c.show.sp.names = (is.null(m$phyloTree) && m$ns<=30) 
         } else {
           c.show.sp.names = show.sp.names.beta[j]
         }
         c.plotTree = !is.null(m$phyloTree)
         if(!is.null(plotTree)){
           c.plotTree = c.plotTree & plotTree[j]
         }
         plotBeta(m, post=postBeta, supportLevel = support.level.beta, param="Sign",
                 plotTree = c.plotTree,
                 covNamesNumbers = c(TRUE,FALSE),
                 spNamesNumbers=c(c.show.sp.names,FALSE),
                 cex=c(0.6,0.6,0.8))
         mymain = paste0("BetaPlot, ",modelnames[j])
         if(!is.null(m$phyloTree)){
           mpost = convertToCodaObject(m)
           rhovals = unlist(poolMcmcChains(mpost$Rho))
           mymain = paste0(mymain,", E[rho] = ",round(mean(rhovals),2),", Pr[rho>0] = ",
                                                                      round(mean(rhovals>0),2))
         }
         title(main=mymain, line=2.5, cex.main=0.8)
       }
     }
     for(j in 1:nm){
       m = models[[j]]      
       if(m$nt>1 & m$nc>1){
         postGamma = getPostEstimate(m, parName="Gamma")
         plotGamma(m, post=postGamma, supportLevel = support.level.gamma, param="Sign",
                   covNamesNumbers = c(TRUE,FALSE),
                   trNamesNumbers=c(m$nt<21,FALSE),
                   cex=c(0.6,0.6,0.8))
         title(main=paste0("GammaPlot ",modelnames[j]), line=2.5,cex.main=0.8)
       }
     }
     for(j in 1:nm){
       m = models[[j]]
       if(m$nr>0 & m$ns>1){
         OmegaCor = computeAssociations(m)
         for (r in 1:m$nr){
           toPlot = ((OmegaCor[[r]]$support>support.level.omega) + 
                    (OmegaCor[[r]]$support<(1-support.level.omega))>0)*sign(OmegaCor[[r]]$mean)
           if(is.null(show.sp.names.omega)){
             c.show.sp.names = (m$ns<=30) 
           } else {
             c.show.sp.names = show.sp.names.omega[j]
           }
           if(!c.show.sp.names){
             colnames(toPlot)=rep("",m$ns)
             rownames(toPlot)=rep("",m$ns)
           }
           if(all(omega.order[[j]]==0)){
             plotOrder = 1:m$ns
           } else {
             if(all(omega.order[[j]]=="AOE")){
               plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
             } else {
               plotOrder = omega.order[[j]]
             }
           }
           cat(c("\n","omega.order","\n","\n"),file=text.file,sep="",append=TRUE)
           cat(m$spNames[plotOrder],file=text.file,sep="\n",append=TRUE)
           mymain = paste0("Associations, ",modelnames[j], ": ",names(m$ranLevels)[[r]])
           if(m$ranLevels[[r]]$sDim>0){
             mpost = convertToCodaObject(m)
             alphavals = unlist(poolMcmcChains(mpost$Alpha[[1]][,1]))
             mymain = paste0(mymain,", E[alpha1] = ",round(mean(alphavals),2),", Pr[alpha1>0] = ",
                                                                        round(mean(alphavals>0),2))
           }
           corrplot(toPlot[plotOrder,plotOrder], method = "color",
                   col=colorRampPalette(c("blue","white","red"))(3),
                   mar=c(0,0,1,0),
                   main=mymain,cex.main=0.8)
           me = as.data.frame(OmegaCor[[r]]$mean)
           me = cbind(m$spNames,me)
           colnames(me)[1] = ""
           po = as.data.frame(OmegaCor[[r]]$support)
           po = cbind(m$spNames,po)
           colnames(po)[1] = ""
           ne = as.data.frame(1-OmegaCor[[r]]$support)
           ne = cbind(m$spNames,ne)
           colnames(ne)[1] = ""
           vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
           filename = file.path(resultDir, paste("parameter_estimates_Omega_",modelnames[j],"_",
                                                                  names(m$ranLevels)[[r]],".xlsx"))
           writexl::write_xlsx(vals,path = filename)
         }
       }
     }
     dev.off()
   }


___________________________________________________

.. _make_predictions:

Make predictions
~~~~~~~~~~~~~~~~

Making predictions and generating plots. 

.. admonition:: input:
  
  Fitted models directory ``models``. 
  **Check the input path**, and typically other parts of the scripts do not need modifications.

.. admonition:: output:
  
  Predictions over environmental gradients (for highest RUN of :ref:`Fit models <fit_model>`) in the file *results/predictions.pdf*.


.. code-block:: R
   :caption: make predictions
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)
   library(ggplot2)

   # make the script reproducible
   set.seed(1)

   # set working directory 
   localDir = "."
   # specify 'models' dir (input dir) 
   modelDir = file.path(localDir, "models")
   # specify 'results' dir (output dir) 
   resultDir = file.path(localDir, "results")
   if (!dir.exists(resultDir)) dir.create(resultDir)

   # number of samples and thinning intervals for each iteration (AS ABOVE)
   samples_list = c(5, 250, 250, 250, 250, 250)   # number of iterations of the MCMC
   thin_list = c(1, 1, 10, 100, 1000, 10000)      # thinning interval;  keep every k-th sample 
   # number of MCMC chains
   nChains = 4

   species.list = NULL # one example species shown for each model,
   #                     selected as prevalence closest to 0.5 (probit models) 
   #                     or most abundant species (other models).
     # -> use species.list to select which species are used as examples for which 
     #                                                           predictions are shown
     # species.list should be a list of length the number of models. 
     # for each element provide either 0 (use default) or a vector of species indices
     # species.list = list()
     # species.list[[1]] = 0
     # species.list[[2]] = c(1,2)

   trait.list = NULL # community weighted mean shown for all traits.
     # -> use trait.list to select for which traits predictions for community weighted 
     #                                                                  mean traits are shown
     # trait.list should be a list of length the number of models. 
     # for each element provide either 0 (use default) or a vector of trait indices
     # see models[[j]]$trNames to see which trait each index corresponds to
     # trait.list = list()
     # trait.list[[1]] = c(2,10)
     # trait.list[[2]] = 0
   
   env.list = NULL # predictions constructed over all environmental gradients.
     # -> use env.list to select over which environmental gradients predictions are generated
     # env.list should be a list of length the number of models. 
     # for each element provide either 0 (use default) or a vector of environmental variables
     # env.list = list()
     # env.list[[1]] = 0
     # env.list[[2]] = c("tree","decay")
   #-------------------------------------------------------------------------------------------#

   for (Lst in nst:1) {
     thin = thin_list[Lst]
     samples = samples_list[Lst]
     filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
     if(file.exists(filename)){break}
   }
   if(file.exists(filename)){
     load(filename)
     nm = length(models)
     modelnames = names(models)
     if(is.null(species.list)){
       species.list = list()
       for(j in 1:nm) species.list[[j]] = 0
     }
     if(is.null(trait.list)){
       trait.list = list()
       for(j in 1:nm) trait.list[[j]] = 0
     }
     if(is.null(env.list)){
       env.list = list()
       for(j in 1:nm) env.list[[j]] = 0
     }
      
     pdf(file= file.path(resultDir,"predictions.pdf"))
     for(j in 1:nm){
       m = models[[j]]
       if(all(env.list[[j]]==0)){
         if(m$XFormula=="~."){
           covariates = colnames(m$XData)
         } else {
           covariates = all.vars(m$XFormula)
         }
       } else {
         covariates = env.list[[j]]
       }
       covariates = setdiff(covariates,"pi")
       ex.sp = which.max(colMeans(m$Y,na.rm = TRUE)) #most common species as example species
       if(m$distr[1,1]==2){
         ex.sp = which.min(abs(colMeans(m$Y,na.rm = TRUE)-0.5))
       }
       if(!all(species.list[[j]])==0){
         ex.sp = species.list[[j]]
       }
       if(length(covariates)>0){
         for(k in 1:(length(covariates))){
           covariate = covariates[[k]]
           Gradient = constructGradient(m,focalVariable = covariate,ngrid = 100)
           Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1,
                                                                                  ngrid = 100)
           predY = predict(m, Gradient=Gradient, expected = TRUE)  
           predY2 = predict(m, Gradient=Gradient2, expected = TRUE)  
           par(mfrow=c(2,1))
           pl = plotGradient(m, Gradient, pred=predY, yshow = 0, measure="S", showData = TRUE, 
                             main = paste0(modelnames[j],": summed response (total effect)"))
           if(inherits(pl, "ggplot")){
             print(pl + labs(title=paste0(modelnames[j],": summed response (total effect)")))
           }
           pl = plotGradient(m, Gradient2, pred=predY2, yshow = 0, measure="S", showData = TRUE, 
                             main = paste0(modelnames[j],": summed response (marginal effect)"))
           if(inherits(pl, "ggplot")){
             print(pl + labs(title=paste0(modelnames[j],": summed response (marginal effect)")))
           }
           for(l in 1:length(ex.sp)){
             par(mfrow=c(2,1))
             pl = plotGradient(m, Gradient, pred=predY, 
                              yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, 
                               measure="Y",index=ex.sp[l], 
                               showData = TRUE, 
                               main = paste0(modelnames[j],": example species (total effect)"))
             if(inherits(pl, "ggplot")){
               print(pl + labs(title=paste0(modelnames[j],": example species (total effect)")))
             }
             pl = plotGradient(m, Gradient2, pred=predY2, 
                              yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, 
                               measure="Y",index=ex.sp[l], showData = TRUE, 
                               main = paste0(modelnames[j],": example species (marginal effect)"))
             if(inherits(pl, "ggplot")){
               print(pl + labs(title=paste0(modelnames[j],": example species (marginal effect)")))
             }
           }
           if(m$nt>1){
             traitSelection = 2:m$nt
             if(!all(trait.list[[j]]==0)) traitSelection = trait.list[[j]]
             for(l in traitSelection){
               par(mfrow=c(2,1))
               pl = plotGradient(m, Gradient, pred=predY, 
                                 measure="T", index=l, showData = TRUE,yshow = 0,
                                 main = paste0(modelnames[j],": community weighted mean trait 
                                 (total effect)"))
               if(inherits(pl, "ggplot")){
                 print(pl + labs(title=paste0(modelnames[j],": community weighted mean trait 
                                                                          (total effect)")))
               }
               pl = plotGradient(m, Gradient2, pred=predY2, 
                                 measure="T", index=l, showData = TRUE, yshow = 0,
                                 main = paste0(modelnames[j],": community weighted mean trait 
                                 (marginal effect)"))
               if(inherits(pl, "ggplot")){
                 print(pl + labs(title=paste0(modelnames[j],": community weighted mean trait 
                                 (marginal effect)")))
               }
             }
           }
         }
       }
     }
     dev.off()
   }

___________________________________________________

|logo_BGE_small| |eufund| |chfund| |ukrifund|
