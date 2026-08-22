#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(glmnet)
})

set.seed(300821)

ROOT <- "/home/samuelajulo/MeteG/MMER_repo"
IN_DIR <- file.path(ROOT, "results/fourstudy_loso/input")
OUT <- file.path(ROOT, "results/fourstudy_loso/transport_validation")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------------
# REQUIRED INPUT CONTRACT
# -----------------------------------------------------------------------------
# baseline_family_long.csv
#   study, trajectory_id, cow_id, family, abundance
#   one row per baseline trajectory x family; abundance must be non-negative
#   relative abundance is preferred, but any within-sample proportional scale
#   is acceptable because the script renormalizes after feature restriction.
#
# outcomes.csv
#   study, trajectory_id, cow_id, ecological_resistance
#
# Allowed studies: Italy, Manitoba, Porcellato, VanBeeck
# ecological_resistance must be the frozen study-specific continuous phenotype.
# -----------------------------------------------------------------------------

FAMILY_FILE <- file.path(IN_DIR, "baseline_family_long.csv")
OUTCOME_FILE <- file.path(IN_DIR, "outcomes.csv")

if (!file.exists(FAMILY_FILE) || !file.exists(OUTCOME_FILE)) {
  stop(
    "Missing four-study LOSO standardized inputs. Expected:\n",
    FAMILY_FILE, "\n",
    OUTCOME_FILE, "\n",
    "See comments at top of scripts/30_fourstudy_loso_transport.R for schema."
  )
}

fam <- read_csv(FAMILY_FILE, show_col_types=FALSE) %>%
  transmute(
    study=as.character(study),
    trajectory_id=as.character(trajectory_id),
    cow_id=as.character(cow_id),
    family=trimws(as.character(family)),
    abundance=as.numeric(abundance)
  ) %>%
  filter(!is.na(family), family!="", is.finite(abundance), abundance>=0) %>%
  group_by(study, trajectory_id, cow_id, family) %>%
  summarise(abundance=sum(abundance), .groups="drop")

outcomes <- read_csv(OUTCOME_FILE, show_col_types=FALSE) %>%
  transmute(
    study=as.character(study),
    trajectory_id=as.character(trajectory_id),
    cow_id=as.character(cow_id),
    ecological_resistance=as.numeric(ecological_resistance)
  )

STUDIES <- c("Italy","Manitoba","Porcellato","VanBeeck")
stopifnot(setequal(sort(unique(outcomes$study)), sort(STUDIES)))
stopifnot(!anyNA(outcomes$ecological_resistance))
stopifnot(all(is.finite(outcomes$ecological_resistance)))
stopifnot(!anyDuplicated(outcomes[,c("study","trajectory_id")]))

# Ensure family and outcome tables describe exactly the same trajectories.
traj_f <- fam %>% distinct(study, trajectory_id, cow_id)
traj_o <- outcomes %>% distinct(study, trajectory_id, cow_id)
stopifnot(nrow(anti_join(traj_o, traj_f, by=c("study","trajectory_id","cow_id")))==0)
stopifnot(nrow(anti_join(traj_f, traj_o, by=c("study","trajectory_id","cow_id")))==0)

cat("\n====================================\n")
cat("FOUR-STUDY LOSO INPUT QC\n")
cat("====================================\n")
print(outcomes %>% group_by(study) %>% summarise(
  trajectories=n(), cows=n_distinct(cow_id),
  mean_R=mean(ecological_resistance), sd_R=sd(ecological_resistance),
  .groups="drop"
))

# Study-relative standardized resistance. This is the PRIMARY transport target.
# It removes cohort-specific mean/scale differences caused by different sampling
# intervals while preserving individual ordering within each study.
outcomes <- outcomes %>%
  group_by(study) %>%
  mutate(
    R_mean_study=mean(ecological_resistance),
    R_sd_study=sd(ecological_resistance),
    z_resistance=(ecological_resistance-R_mean_study)/R_sd_study,
    high_resistance=as.integer(ecological_resistance > median(ecological_resistance))
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
rmse <- function(y,p) sqrt(mean((y-p)^2))
mae <- function(y,p) mean(abs(y-p))
cv_r2 <- function(y,p) 1-sum((y-p)^2)/sum((y-mean(y))^2)
safe_cor <- function(y,p,method="pearson") {
  if (length(y)<3 || sd(y)==0 || sd(p)==0) return(NA_real_)
  suppressWarnings(cor(y,p,method=method))
}

auc_rank <- function(y,score) {
  y <- as.integer(y)
  if (length(unique(y))<2) return(NA_real_)
  n1 <- sum(y==1); n0 <- sum(y==0)
  r <- rank(score, ties.method="average")
  (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0)
}

c_index <- function(y,score) {
  n <- length(y); conc <- 0; ties <- 0; usable <- 0
  for (i in seq_len(n-1)) for (j in (i+1):n) {
    if (y[i]==y[j]) next
    usable <- usable+1
    dy <- sign(y[i]-y[j]); ds <- sign(score[i]-score[j])
    if (ds==0) ties <- ties+1 else if (dy==ds) conc <- conc+1
  }
  if (usable==0) return(NA_real_)
  (conc + 0.5*ties)/usable
}

study_equal_cow_weights <- function(meta) {
  meta %>%
    group_by(study, cow_id) %>% mutate(n_in_cow=n()) %>% ungroup() %>%
    group_by(study) %>% mutate(n_cows=n_distinct(cow_id)) %>% ungroup() %>%
    mutate(w=1/(n_cows*n_in_cow)) %>% pull(w)
}

clr_transform <- function(X,pseudo) {
  X <- as.matrix(X)
  X[X<=0] <- pseudo
  X <- X/rowSums(X)
  L <- log(X)
  L-rowMeans(L)
}

# Build a feature representation using ONLY training studies to choose families,
# pseudocount, PCA loadings and PC scaling. The held-out study contributes only
# its baseline feature distribution for unsupervised within-study centering.
build_transport_features <- function(train_studies, test_study, kmax=6,
                                     prevalence_cut=0.10,
                                     min_train_studies_present=2) {

  ids <- outcomes %>% select(study, trajectory_id, cow_id)
  all_families <- sort(unique(fam$family))

  wide <- fam %>%
    select(study, trajectory_id, cow_id, family, abundance) %>%
    complete(nesting(study, trajectory_id, cow_id), family=all_families,
             fill=list(abundance=0)) %>%
    pivot_wider(names_from=family, values_from=abundance, values_fill=0)

  meta <- wide %>% select(study, trajectory_id, cow_id)
  X <- as.matrix(wide %>% select(-study,-trajectory_id,-cow_id))
  rownames(X) <- paste(meta$study,meta$trajectory_id,sep="__")

  train_idx <- meta$study %in% train_studies
  test_idx <- meta$study==test_study
  Xtr0 <- X[train_idx,,drop=FALSE]
  Xte0 <- X[test_idx,,drop=FALSE]
  mtr <- meta[train_idx,,drop=FALSE]
  mte <- meta[test_idx,,drop=FALSE]

  # Renormalize prior to prevalence assessment.
  Xtr_ra <- Xtr0/rowSums(Xtr0)
  prev_pool <- colMeans(Xtr_ra>0)
  study_presence <- sapply(colnames(Xtr_ra), function(f) {
    sum(sapply(train_studies, function(s) {
      ii <- which(mtr$study==s)
      any(Xtr_ra[ii,f]>0)
    }))
  })
  keep <- names(prev_pool)[prev_pool>=prevalence_cut &
                            study_presence>=min_train_studies_present]
  if (length(keep)<kmax+1) stop("Too few transferable families in fold ",test_study)

  Xtr <- Xtr0[,keep,drop=FALSE]
  Xte <- Xte0[,keep,drop=FALSE]
  pos <- Xtr[Xtr>0]
  pseudo <- 0.5*min(pos)

  Ctr <- clr_transform(Xtr,pseudo)
  Cte <- clr_transform(Xte,pseudo)

  # Unsupervised study-centering. For the held-out study this uses baseline
  # composition only; no resistance outcomes are used.
  for (s in unique(mtr$study)) {
    ii <- which(mtr$study==s)
    Ctr[ii,] <- sweep(Ctr[ii,,drop=FALSE],2,colMeans(Ctr[ii,,drop=FALSE]),"-")
  }
  Cte <- sweep(Cte,2,colMeans(Cte),"-")

  # Remove zero-variance training features.
  sds <- apply(Ctr,2,sd)
  good <- is.finite(sds) & sds>0
  Ctr <- Ctr[,good,drop=FALSE]
  Cte <- Cte[,good,drop=FALSE]

  pca <- prcomp(Ctr, center=FALSE, scale.=FALSE)
  if (ncol(pca$x)<kmax) stop("PCA has fewer than required PCs")
  tr_score <- pca$x[,1:kmax,drop=FALSE]
  te_score <- Cte %*% pca$rotation[,1:kmax,drop=FALSE]

  mu <- colMeans(tr_score)
  sig <- apply(tr_score,2,sd)
  tr_z <- sweep(sweep(tr_score,2,mu,"-"),2,sig,"/")
  te_z <- sweep(sweep(te_score,2,mu,"-"),2,sig,"/")

  colnames(tr_z) <- paste0("PC",seq_len(kmax))
  colnames(te_z) <- paste0("PC",seq_len(kmax))

  ve <- pca$sdev^2/sum(pca$sdev^2)

  list(
    train=cbind(mtr,as.data.frame(tr_z)),
    test=cbind(mte,as.data.frame(te_z)),
    retained_families=keep,
    variance=ve,
    pseudo=pseudo
  )
}

# Inner leave-one-study-out tuning among the three outer-training studies.
# Feature engineering is itself re-fitted inside each inner fold.
lambda_grid <- 10^seq(-4,2,length.out=80)

tune_lambda <- function(outer_train_studies,k) {
  inner_cache <- lapply(outer_train_studies,function(val_study) {
    tr_studies <- setdiff(outer_train_studies,val_study)
    f <- build_transport_features(tr_studies,val_study,kmax=max(6,k))
    tr <- f$train %>% left_join(outcomes,by=c("study","trajectory_id","cow_id"))
    va <- f$test %>% left_join(outcomes,by=c("study","trajectory_id","cow_id"))
    list(train=tr,val=va)
  })

  scores <- lapply(lambda_grid,function(lam) {
    fold_rmse <- sapply(inner_cache,function(z) {
      Xtr <- as.matrix(z$train[,paste0("PC",seq_len(k)),drop=FALSE])
      Xva <- as.matrix(z$val[,paste0("PC",seq_len(k)),drop=FALSE])
      w <- study_equal_cow_weights(z$train)
      fit <- glmnet(Xtr,z$train$z_resistance,alpha=0,lambda=lam,
                    standardize=FALSE,intercept=TRUE,weights=w)
      pr <- as.numeric(predict(fit,Xva,s=lam))
      rmse(z$val$z_resistance,pr)
    })
    tibble(lambda=lam,mean_study_RMSE=mean(fold_rmse),sd_study_RMSE=sd(fold_rmse))
  }) %>% bind_rows()

  best <- scores %>% arrange(mean_study_RMSE,desc(lambda)) %>% slice(1)
  list(lambda=best$lambda,curve=scores)
}

fit_outer <- function(test_study,k) {
  train_studies <- setdiff(STUDIES,test_study)
  features <- build_transport_features(train_studies,test_study,kmax=6)
  tr <- features$train %>% left_join(outcomes,by=c("study","trajectory_id","cow_id"))
  te <- features$test %>% left_join(outcomes,by=c("study","trajectory_id","cow_id"))

  tuned <- tune_lambda(train_studies,k)
  Xtr <- as.matrix(tr[,paste0("PC",seq_len(k)),drop=FALSE])
  Xte <- as.matrix(te[,paste0("PC",seq_len(k)),drop=FALSE])
  w <- study_equal_cow_weights(tr)

  fit <- glmnet(Xtr,tr$z_resistance,alpha=0,lambda=tuned$lambda,
                standardize=FALSE,intercept=TRUE,weights=w)
  pred <- as.numeric(predict(fit,Xte,s=tuned$lambda))

  tibble(
    study=test_study,
    model=paste0("PC1-PC",k),
    trajectory_id=te$trajectory_id,
    cow_id=te$cow_id,
    observed_R=te$ecological_resistance,
    observed_z=te$z_resistance,
    predicted_score=pred,
    high_resistance=te$high_resistance,
    lambda=tuned$lambda,
    n_families=length(features$retained_families),
    cumulative_PCA_variance=100*sum(features$variance[seq_len(k)])
  )
}

# -----------------------------------------------------------------------------
# Run all 4 x 2 outer folds.
# -----------------------------------------------------------------------------
cat("\nRunning four-study LOSO transport validation...\n")
predictions <- bind_rows(lapply(c(2,6),function(k) {
  bind_rows(lapply(STUDIES,function(s) {
    cat("Model PC1-PC",k," | held-out study: ",s,"\n",sep="")
    fit_outer(s,k)
  }))
}))

metric_one <- function(d) {
  tibble(
    n=nrow(d),
    RMSE_z=rmse(d$observed_z,d$predicted_score),
    MAE_z=mae(d$observed_z,d$predicted_score),
    CV_R2_z=cv_r2(d$observed_z,d$predicted_score),
    Pearson_r=safe_cor(d$observed_z,d$predicted_score,"pearson"),
    Spearman_rho=safe_cor(d$observed_z,d$predicted_score,"spearman"),
    C_index=c_index(d$observed_z,d$predicted_score),
    AUC_high_vs_low=auc_rank(d$high_resistance,d$predicted_score)
  )
}

per_study <- predictions %>% group_by(model,study) %>% group_modify(~metric_one(.x)) %>% ungroup()
pooled <- predictions %>% group_by(model) %>% group_modify(~metric_one(.x)) %>% ungroup() %>% mutate(study="POOLED")
macro <- per_study %>% group_by(model) %>% summarise(
  n=sum(n),
  RMSE_z=mean(RMSE_z,na.rm=TRUE),
  MAE_z=mean(MAE_z,na.rm=TRUE),
  CV_R2_z=mean(CV_R2_z,na.rm=TRUE),
  Pearson_r=mean(Pearson_r,na.rm=TRUE),
  Spearman_rho=mean(Spearman_rho,na.rm=TRUE),
  C_index=mean(C_index,na.rm=TRUE),
  AUC_high_vs_low=mean(AUC_high_vs_low,na.rm=TRUE),
  .groups="drop"
) %>% mutate(study="MACRO_AVERAGE")
metrics <- bind_rows(per_study,pooled,macro) %>% select(model,study,everything())

# -----------------------------------------------------------------------------
# Cow-clustered held-out sign-flip permutation inference.
# Predictions remain fixed because each study's predictions were generated
# without using that study's outcomes. The null breaks prediction-outcome
# direction while preserving study membership, cow clustering and the magnitude
# of within-cow outcome deviations.
# -----------------------------------------------------------------------------
N_PERM <- as.integer(Sys.getenv("MMER_LOSO_PERM",unset="9999"))
cat("\nRunning ",N_PERM," cow-clustered sign-flip permutations...\n",sep="")

perm_results <- lapply(c("PC1-PC2","PC1-PC6"),function(m) {
  d <- predictions %>% filter(model==m)
  obs_macro_spear <- d %>% group_by(study) %>% summarise(r=safe_cor(observed_z,predicted_score,"spearman"),.groups="drop") %>% summarise(mean(r,na.rm=TRUE)) %>% pull()
  obs_pool_spear <- safe_cor(d$observed_z,d$predicted_score,"spearman")
  obs_pool_pear <- safe_cor(d$observed_z,d$predicted_score,"pearson")

  null <- replicate(N_PERM, {
    keys <- d %>% distinct(study,cow_id) %>% mutate(sign=sample(c(-1,1),n(),replace=TRUE))
    dp <- d %>% left_join(keys,by=c("study","cow_id")) %>% mutate(yperm=observed_z*sign)
    macro_s <- dp %>% group_by(study) %>% summarise(r=safe_cor(yperm,predicted_score,"spearman"),.groups="drop") %>% summarise(mean(r,na.rm=TRUE)) %>% pull()
    c(macro_s,
      safe_cor(dp$yperm,dp$predicted_score,"spearman"),
      safe_cor(dp$yperm,dp$predicted_score,"pearson"))
  })
  null <- t(null)
  tibble(
    model=m,
    statistic=c("macro_Spearman","pooled_Spearman","pooled_Pearson"),
    observed=c(obs_macro_spear,obs_pool_spear,obs_pool_pear),
    permutation_p=c(
      (1+sum(abs(null[,1])>=abs(obs_macro_spear),na.rm=TRUE))/(N_PERM+1),
      (1+sum(abs(null[,2])>=abs(obs_pool_spear),na.rm=TRUE))/(N_PERM+1),
      (1+sum(abs(null[,3])>=abs(obs_pool_pear),na.rm=TRUE))/(N_PERM+1)
    ),
    n_permutations=N_PERM
  )
}) %>% bind_rows()

cat("\n====================================\n")
cat("FOUR-STUDY LOSO TRANSPORT METRICS\n")
cat("====================================\n")
print(metrics,width=Inf)
cat("\n====================================\n")
cat("HELD-OUT PREDICTION PERMUTATION INFERENCE\n")
cat("====================================\n")
print(perm_results,width=Inf)

write_csv(predictions,file.path(OUT,"LOSO_predictions_PC2_PC6.csv"))
write_csv(metrics,file.path(OUT,"LOSO_metrics_PC2_PC6.csv"))
write_csv(perm_results,file.path(OUT,"LOSO_cluster_signflip_permutation.csv"))

# Simple publication-ready observed-vs-predicted plots if ggplot2 is available.
if (requireNamespace("ggplot2",quietly=TRUE)) {
  g <- ggplot2::ggplot(predictions,ggplot2::aes(x=observed_z,y=predicted_score)) +
    ggplot2::geom_point(alpha=0.8) +
    ggplot2::geom_smooth(method="lm",se=FALSE) +
    ggplot2::facet_grid(model~study,scales="free") +
    ggplot2::labs(
      x="Observed study-standardized ecological resistance",
      y="Predicted resistance score",
      title="Leave-one-study-out transport of baseline ecological architecture"
    ) + ggplot2::theme_bw(base_size=11)
  ggplot2::ggsave(file.path(OUT,"LOSO_observed_vs_predicted.pdf"),g,width=12,height=6.5)
  ggplot2::ggsave(file.path(OUT,"LOSO_observed_vs_predicted.png"),g,width=12,height=6.5,dpi=300)
}

capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"))
cat("\n[PASS] FOUR-STUDY LOSO TRANSPORT VALIDATION COMPLETE\n")
cat("Output:",OUT,"\n")
