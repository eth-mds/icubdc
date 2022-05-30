
train_predict <- function(trn_src, tst_srcs = trn_src, method = "rf",
                          cfg = config("features"), seed = 2021, n_folds = 5L,
                          train_loader = split_coh_train,
                          test_loader = split_coh_test, ...,
                          job_prefix = jobname(), run_name = unique_name()) {

  res_dir <- results_dir(paste(job_prefix, jobid(), sep = "_"))
  on.exit(if (!length(list.files(res_dir))) unlink(res_dir))

  ctrl <- new_controller(method, run_name = run_name, trn_src = trn_src,
                         n_folds = n_folds, seed = seed, cfg = cfg)

  msg_ts("writing results of run ", run_name, " to ", res_dir)

  trn_dat <- train_loader(ctrl, trn_src, ...)
  trn_fld <- split_folds(seed, trn_dat, "hypo", n_folds)

  keep <- !is.na(trn_fld)

  trn_dat <- trn_dat[keep, ]
  trn_fld <- trn_fld[keep]

  update_controller(ctrl,
    train_folds = lapply(split(id_col(trn_dat), trn_fld), unique)
  )

  trn_res <- trn_dat[["hypo"]]
  trn_dat <- trn_dat[, c(meta_vars(trn_dat), "hypo") := NULL]
  trn_dat <- hypo_prep(ctrl, trn_dat, cfg)
  model <- hypo_train(ctrl, trn_dat, trn_res, trn_fld, n_cores(TRUE), seed)

  res_fil <- file.path(res_dir, run_name)

  hypo_save(ctrl, model, res_fil)

  for (tst_src in tst_srcs) {

    tst_dat <- test_loader(ctrl, tst_src, ...)

    tst_res <- tst_dat[["hypo"]]
    tst_cnt <- tst_dat[["hypo_cnt_imp"]]
    tst_met <- tst_dat[, meta_vars(tst_dat), with = FALSE]
    tst_dat <- tst_dat[, c(meta_vars(tst_dat), "hypo") := NULL]
    tst_dat <- hypo_prep(ctrl, tst_dat, cfg)

    pred <- hypo_pred(ctrl, model, tst_dat)

    tst_met <- tst_met[, c("label", "prediction", "cnt") := 
                         list(tst_res, pred, tst_cnt)]

    cnts <- table(tst_met$label)
    aucs <- precrec::auc(
      precrec::evalmod(scores = tst_met$prediction,
                       labels = tst_met$label)
    )

    msg_ts("hypo prevalence: ", fmt4(cnts["TRUE"] / nrow(tst_met)), ", ",
           paste0("au", tolower(aucs$curvetypes), ": ", fmt4(aucs$aucs),
                  collapse = ", "))

    save_predictions(tst_met, res_fil, tst_src)
  }

  save_ctrl(ctrl, res_fil)
}

split_coh_train <- function(ctrl, src, ...) {

  res <- load_split_cohort(src, ctrl[["seed"]], "train", ...)

  update_controller(ctrl,
    lst = rename(get_settings(res), "train_ids", "pids")
  )
  
  res
}

# <<<<<<< HEAD
#   if (is.string(cohort)) {
#     if (identical(cohort, "default")) {
#       pids <- NULL
#     } else {
#       pids <- cohort(src, type = cohort)
#     }
#   } else {
#     pids <- cohort
#     cohort <- "custom"
#   }
# =======

split_coh_test <- function(ctrl, src, ...) {

  res <- load_split_cohort(src, ctrl[["seed"]], "test", ...)

  update_controller(ctrl,
    test_ids = c(ctrl$test_ids, rename(get_settings(res), src, "pids")[src])
  )

  res
}

load_split_cohort <- function(src, seed, set, cohort = "all",
                              test_prop = 0.2, train_prop = 1 - test_prop,
                              feat_set = "imp", pids = cohort(src, cohort),
                              ins_only = TRUE, ...) {

  props <- c(train_prop, test_prop)
  sets  <- c("train", "test")
  set   <- match.arg(set, sets)

  msg_ts("loading ", src, " data (\"", cohort, "\" cohort), splitting ",
         fmt_perc(props), " (train/test)")

  splt <- get_response(src, ..., pids = pids)
  splt <- split_folds(seed, splt, "hypo", 2L, props, ret_ids = TRUE)

  ind  <- match(set, sets)
  lens <- lengths(splt)

  msg_ts("using ", set, " (", fmt_perc(lens[ind] / sum(lens)), ")")

  splt <- splt[[ind]]
  data <- load_data(src, feat_regex(feat_set), splt, ...)
  # wild fix by Drago
  if (ins_only) {
    ins_sub <- slide(data, list(ins_on = any(insfx_imp > 0)), 
                    before = hours(12L))
    data <- merge(data, ins_sub, all.x = TRUE)
    data <- data[is_true(ins_on)]
    data[, ins_on := NULL]
  }
  set_settings(data, cohort = cohort, pids = splt, feat_set = feat_set,
               train_prop = train_prop, test_prop = test_prop)
}

feat_regex <- function(set) {
  switch(
    set,
    raw = "_raw$", imp = "_imp$", cnt = "_cnt$",
    ind = "_(imp|ind)$", mis = "_(imp|ind|cnt)$",
    lbs = "_(imp|4lbk)$", lbm = "_(imp|8lbk)$", lbl = "_(imp|16lbk)$",
    lbk = "_(imp|lbk)$", all = "_(imp|ind|cnt|lbk)$",
    stop("unknown feature set \"", set, "\"")
  )
}

split_folds <- function(seed, tbl, strat_var, n_folds, fold_props = 1 / n_folds,
                        ret_ids = FALSE) {

  assert_that(
    is_ts_tbl(tbl), is.count(n_folds), is.string(strat_var), is.count(seed),
    is.flag(ret_ids)
  )

  if (length(fold_props) == 1L) {
    fold_props <- rep(fold_props, n_folds)
  }

  assert_that(
    length(fold_props) == n_folds, sum(fold_props) <= 1
  )

  idvr <- id_var(tbl)
  flds <- tbl[, list(split_var = max(get(strat_var))), by = c(idvr)]

  pids <- split(flds, by = "split_var", keep.by = FALSE)
  pids <- lapply(pids, unlist, recursive = FALSE, use.names = FALSE)

  msg_ts("setting rng seed to ", seed)

  set.seed(seed)

  lens <- lengths(pids)
  pids <- lapply(pids, sample)

  msg_ts("splitting ", sum(lens), " values into ", n_folds, " folds (",
         fmt_perc(fold_props), "), stratified on ", strat_var, " (",
         fmt_perc(lens / sum(lens)), ")")

  inds <- lapply(lens, `*`, fold_props)
  inds <- lapply(inds, floor)
  inds <- Map(rep, list(seq_len(n_folds)), inds)

  if (sum(lengths(inds)) != sum(lengths(pids))) {
    msg_ts("not using ", sum(lengths(pids)) - sum(lengths(inds)),
           " ids due to rounding")
  }

  pids <- Map(head, pids, lengths(inds))
  pids <- do.call("c", unname(pids))
  inds <- do.call("c", inds)

  if (ret_ids) {
    unname(lapply(split(pids, inds), sample))
  } else {
    inds[match(id_col(tbl), pids)]
  }
}

save_ctrl <- function(ctrl, path) {
  qs::qsave(as.list.environment(ctrl), paste0(path, ".set"))
  invisible(NULL)
}

save_predictions <- function(pred, path, tst_src) {
  qs::qsave(pred, paste0(path, "-", tst_src, ".prd"))
  invisible(NULL)
}

new_controller <- function(method, ...) {

  msg_ts("initializing ", method, " controller")

  structure(list2env(c(list(method = method), list(...))),
            class = paste0("hypo_", method))
}

save_mod <- function(model, path) {
  qs::qsave(model, paste0(path, ".mod"))
  invisible(NULL)
}

hypo_prep <- function(ctrl, tbl, ...) UseMethod("hypo_prep", ctrl)

hypo_prep_rf <- function(ctrl, tbl, ...) tbl

.S3method("hypo_prep", "hypo_rf", hypo_prep_rf)

hypo_prep_lasso <- function(ctrl, tbl, cfg, ...) indicator_encoding(tbl, cfg)

.S3method("hypo_prep", "hypo_lasso", hypo_prep_lasso)

hypo_prep_lgbm <- function(ctrl, tbl, cfg, ...) as.matrix(tbl)

.S3method("hypo_prep", "hypo_lgbm", hypo_prep_lgbm)

hypo_train <- function(ctrl, x, y, folds, n_cores, ...) {

  assert_that(is.count(n_cores))

  UseMethod("hypo_train", ctrl)
}

hypo_train_rf <- function(ctrl, x, y, folds, n_cores, seed, ...) {

  # is this ok?
  folds <- as.integer(folds != 1)
  grid  <- list()

  for (mns in c(10, 30, 100, 200, 500, 1000)) {

    msg_ts("trying min node size: ", mns)

    score <- ranger::ranger(
      y = y, x = x, probability = TRUE, min.node.size = mns,
      num.threads = n_cores, case.weights = folds, holdout = TRUE,
      seed = seed, ...
    )$prediction.error

    grid <- c(grid,
      list(list(min_node_size = mns, score = score))
    )

    gc(verbose = FALSE)
  }

  opt <- grid[[which.min(vapply(grid, `[[`, numeric(1L), "score"))]]
  mns <- opt[["min_node_size"]]

  msg_ts("choosing min node size: ", mns)

  update_controller(ctrl, param_grid = grid)

  ranger::ranger(
    y = y, x = x, probability = TRUE, min.node.size = mns,
    importance = "impurity", num.threads = n_cores, seed = seed, ...
  )
}

.S3method("hypo_train", "hypo_rf", hypo_train_rf)

hypo_train_lasso <- function(ctrl, x, y, folds, n_cores, seed, ...) {

  res <- biglasso::cv.biglasso(
    x[["data"]], y, family = "binomial", ncores = n_cores, cv.ind = folds,
    nlambda = 50, verbose = TRUE, ...
  )

  update_controller(ctrl,
    param_grid = list(lambda = res$lambda, score = res$cve)
  )

  attr(res, "ind_encoding") <- x[["meta"]]

  res
}

.S3method("hypo_train", "hypo_lasso", hypo_train_lasso)

hypo_train_lgbm <- function(ctrl, x, y, folds, n_cores, seed, ...) {

  folds <- lapply(seq_along(unique(folds)), `==`, folds)
  folds <- lapply(folds, which)

  dtrain <- lightgbm::lgb.Dataset(x, label = y)

  grid <- list()

  for (num_leaf in c(31, 50, 100))  {
    for (num_tree in c(250, 500)) {
      for (lr in c(0.001, 0.01, 0.1, 1)) {

        msg_ts("trying num leaves: ", num_leaf, ", num trees: ", num_tree,
               ", learning ", "rate: ", lr)

        score <- lightgbm::lgb.cv(
          data = dtrain, num_leaves = num_leaf, nrounds = num_tree,
          learning_rate = lr, boosting = "gbdt", folds = folds,
          num_threads = n_cores, reset_data = TRUE, verbose = 0L,
          objective = "binary", force_col_wise = TRUE, ...
        )$best_score

        grid <- c(grid,
          list(list(num_leaves = num_leaf, num_trees = num_tree,
                    learn_rate = lr, score = score))
        )

        gc(verbose = FALSE)
      }
    }
  }

  opt <- grid[[which.min(vapply(grid, `[[`, numeric(1L), "score"))]]

  msg_ts("choosing: num leaves: ", opt[["num_leaves"]], ", num trees: ",
         opt[["num_trees"]], ", ", "learning rate: ", opt[["learn_rate"]])

  update_controller(ctrl, param_grid = grid)

  lightgbm::lgb.train(
    data = dtrain, num_leaves = opt[["num_leaves"]],
    nrounds = opt[["num_trees"]], learning_rate = opt[["learn_rate"]],
    boosting = "gbdt", num_threads = n_cores, verbose = 0L,
    objective = "binary", force_col_wise = TRUE, ...
  )
}

.S3method("hypo_train", "hypo_lgbm", hypo_train_lgbm)

hypo_save <- function(ctrl, model, path, ...) {
  UseMethod("hypo_save", ctrl)
}

hypo_save_rf <- function(ctrl, model, path, ...) save_mod(model, path)

.S3method("hypo_save", "hypo_rf", hypo_save_rf)

hypo_save_lasso <- function(ctrl, model, path, ...) save_mod(model, path)

.S3method("hypo_save", "hypo_lasso", hypo_save_lasso)

hypo_save_lgbm <- function(ctrl, model, path, ...) {

  set_na <- if (is.na(model$raw)) {
    model$save()
    TRUE
  }

  save_mod(model, path)

  if (isTRUE(set_na)) {
    model$raw <- NA
  }

  invisible(NULL)
}

.S3method("hypo_save", "hypo_lgbm", hypo_save_lgbm)

hypo_pred <- function(ctrl, model, new_data, ...) UseMethod("hypo_pred", ctrl)

hypo_pred_rf <- function(ctrl, model, new_data, ...) {

  res <- predict(model, new_data, type = "response")
  cls <- identical(res$treetype, "Probability estimation")

  res <- res$predictions

  if (cls) {
    res <- res[, 2L]
  }

  res
}

.S3method("hypo_pred", "hypo_rf", hypo_pred_rf)

hypo_pred_lasso <- function(ctrl, model, new_data, ...) {
  predict(model, new_data[["data"]], type = "response")[, 1L]
}

.S3method("hypo_pred", "hypo_lasso", hypo_pred_lasso)

hypo_pred_lgbm <- function(ctrl, model, new_data, ...) predict(model, new_data)

.S3method("hypo_pred", "hypo_lgbm", hypo_pred_lgbm)
