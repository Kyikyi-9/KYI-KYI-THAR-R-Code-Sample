# VAR & IRF (ggplot) + Cumulative IRF + (optional) SVAR
# Data: dddM.csv with columns ddGDP, ddM2, ddInf, ddNex, ddRex
# Author: KYI KYI THAR
############################################################

## ============ 0) Setup ============
pkgs <- c("urca","tseries","vars","svars","lmtest","forecast",
          "haven","tidyverse","ggplot2","reshape2","scales")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, Ncpus = 2)
invisible(lapply(pkgs, library, character.only = TRUE))
dir.create("var_outputs", showWarnings = FALSE)
set.seed(123)

## ============ 1) Load data ============
setwd("/Users/kyikyithar/Desktop/R programming")
Myan_var <- read.csv("dddM.csv")
vars_need <- c("ddGDP","ddM2","ddInf","ddNex","ddRex")
stopifnot(all(vars_need %in% names(Myan_var)))
Myan_var <- Myan_var |>
  dplyr::select(all_of(vars_need)) |>
  tidyr::drop_na()

cat("Sample size (T):", nrow(d), "\n")

## ============ 2) Unit-root tests (ADF & PP) ============

cat("\n--- ADF tests ---\n")
print(adf.test(Myan_var$ddGDP)); print(adf.test(Myan_var$ddM2)); print(adf.test(Myan_var$ddInf)); print(adf.test(Myan_var$ddNex)); print(adf.test(Myan_var$ddRex))

cat("\n--- PP tests ---\n")
print(pp.test(Myan_var$ddGDP)); print(pp.test(Myan_var$ddM2)); print(pp.test(Myan_var$ddInf)); print(pp.test(Myan_var$ddNex)); print(pp.test(Myan_var$ddRex))

## ============ 3) OLS sanity checks (optional) ============
cat("\n--- OLS: Inf eq ---\n")
print(summary(lm(ddInf ~ ddGDP + ddM2 + ddNex + ddRex, data = Myan_var)))
cat("\n--- OLS: GDP eq ---\n")
print(summary(lm(ddGDP ~ ddInf + ddM2 + ddNex + ddRex, data = Myan_var)))

## ============ 4) Build quarterly ts ============
m_df <- ts(Myan_var[, vars_need], start = c(2013,2), frequency = 4)
colnames(m_df) <- vars_need

png("var_outputs/plot_ts.png", width=1400, height=900, res=130)
plot(m_df, main = "Time Series Plot", ylab = "Level", xlab = "Year")
dev.off()

## ============ 5) Lag selection ============
sel <- VARselect(m_df, lag.max = 4, type = "const")
cat("\n--- VARselect ---\n"); print(sel$selection)
p_opt <- sel$selection[["AIC(n)"]]; if (is.null(p_opt)) p_opt <- 2

## ============ 6) Estimate VAR + diagnostics ============
var_model <- VAR(m_df, p = p_opt, type = "const", season = 4)
cat("\n--- VAR summary ---\n"); print(summary(var_model))

roots_mod <- roots(var_model)
cat("\n[INFO] roots(|λ|):", paste(round(roots_mod,3), collapse=", "), "\n")
if (any(Mod(roots_mod) >= 1)) warning("VAR may be unstable: some roots >= 1")

cat("\n--- Residual diagnostics ---\n")
print(serial.test(var_model, lags.pt = 12, type = "PT.asymptotic"))
print(arch.test(var_model, lags.multi = 5))
print(normality.test(var_model))

stab_obj <- stability(var_model, type = "OLS-CUSUM")
png("var_outputs/stability_OLS_CUSUM.png", width=1200, height=800, res=130)
plot(stab_obj)
dev.off()

## ============ 7) Granger causality ============
cat("\n--- Causality tests (VAR Wald) ---\n")
for (v in vars_need) {
  cat("\nCause =", v, "\n"); print(causality(var_model, cause = v))
}

## ============ 8) FEFD (Forecast Error Variance Decomposition) ============
cat("\n--- FEVD (h=1,4,8,12) ---\n")
fevd_h1  <- fevd(var_model, n.ahead = 1);  print(fevd_h1)
fevd_h4  <- fevd(var_model, n.ahead = 4);  print(fevd_h4)
fevd_h8  <- fevd(var_model, n.ahead = 8);  print(fevd_h8)
fevd_h12 <- fevd(var_model, n.ahead = 12); print(fevd_h12)

## ============ 9) IRFs (orthogonalized via Cholesky) ============
irf_ortho <- irf(
  var_model,
  impulse  = c("ddGDP","ddM2","ddNex","ddRex"),
  response = c("ddGDP","ddM2","ddInf","ddNex","ddRex"),
  n.ahead  = 20, ortho = TRUE, boot = TRUE, runs = 1000, ci = 0.95
)

# base描画
png("var_outputs/IRF_ortho_base.png", width=1400, height=900, res=130)
# (1) VAR の安定性確認
roots(var_model)

# (2) IRFオブジェクトの中身を確認
str(irf_ortho)

# (3) 各変数の欠損や定数列の有無
summary(m_df)
apply(m_df, 2, sd)

roots_mod <- roots(var_model)
if (any(Mod(roots_mod) >= 1)) {
  message("[WARN] VAR unstable → ラグ数を短くして再推定します。")
  var_model <- VAR(m_df, p = max(1, p_opt - 1), type = "const")
}

m_df <- m_df[, apply(m_df, 2, function(x) sd(x, na.rm = TRUE) > 0)]
m_df <- na.omit(m_df)

irf_ortho <- irf(var_model, n.ahead = 10, boot = TRUE, runs = 500, ci = 0.95)
plot(irf_ortho, main = "Orthogonalized IRFs (Cholesky)")
plot(irf_ortho, main = "Orthogonalized IRFs (Cholesky)")
dev.off()

## ============ 10) IRFをggplotで図化（CIリボン & ファセット） ============
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(reshape2); library(scales)
})

# --- 関数の定義 ---
irf_to_long <- function(irf_ortho) {
  stopifnot(!is.null(irf_ortho), !is.null(irf_ortho), !is.null(irf_ortho$Upper))
  imps <- names(irf_ortho$irf)
  out <- lapply(imps, function(imp) {
    base <- as.data.frame(irf_ortho$irf[[imp]])
    base$period <- seq_len(nrow(base))
    base_long <- reshape2::melt(base, id.vars = "period",
                                variable.name = "response", value.name = "irf")
    lo <- as.data.frame(irf_ortho$Lower[[imp]])
    hi <- as.data.frame(irf_ortho$Upper[[imp]])
    lo$period <- hi$period <- seq_len(nrow(lo))
    lo_long <- reshape2::melt(lo, id.vars = "period",
                              variable.name = "response", value.name = "low")
    hi_long <- reshape2::melt(hi, id.vars = "period",
                              variable.name = "response", value.name = "high")
    dplyr::left_join(base_long, lo_long, by = c("period","response")) |>
      dplyr::left_join(hi_long, by = c("period","response")) |>
      dplyr::mutate(impulse = imp)
  })
  dplyr::bind_rows(out)
}

plot_irf_gg <- function(irf_ortho, file_out,
                        facet_by = c("response","impulse"),
                        title    = "Impulse Responses (Cholesky Orthogonalized)") {
  facet_by <- match.arg(facet_by)
  df <- irf_to_long(irf_ortho) |>
    dplyr::filter(is.finite(irf) & is.finite(low) & is.finite(high))
  
  p <- ggplot(df, aes(x = period, y = irf)) +
    geom_ribbon(aes(ymin = low, ymax = high, group = interaction(impulse, response)), alpha = 0.18) +
    geom_line(aes(color = impulse), linewidth = 0.9) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(title = title, x = "Horizon", y = "Response (level/Δ)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  if (facet_by == "response") {
    p <- p + facet_wrap(~response, scales = "free_y")
  } else {
    p <- p + facet_wrap(~impulse, scales = "free_y")
  }
  
  ggsave(filename = file_out, plot = p, width = 11, height = 7, dpi = 300)
  message("[Saved] ", file_out)
}

exists("plot_irf_gg")

plot_irf_gg(
  irf_ortho = irf_ortho,
  file_out = "var_outputs/IRF_ortho_gg_byResponse.png",
  facet_by = "response",
  title    = "Orthogonalized IRFs (Cholesky) - Facet by Response"
)

# 応答別ファセット
plot_irf_gg(irf_ortho,
            file_out = "var_outputs/IRF_ortho_gg_byResponse.png",
            facet_by = "response",
            title    = "Orthogonalized IRFs (Cholesky) - Facet by Response")

# 衝撃別ファセット
plot_irf_gg(irf_ortho,
            file_out = "var_outputs/IRF_ortho_gg_byImpulse.png",
            facet_by = "impulse",
            title    = "Orthogonalized IRFs (Cholesky) - Facet by Impulse")


## ============ 11) 累積IRF（cumsum） ============
irf_cumsum <- function(irf_ortho){
  stopifnot(!is.null(irf_ortho$irf), !is.null(irf_ortho$Lower), !is.null(irf_ortho$Upper))
  imps <- names(irf_ortho$irf)
  out <- list(irf = list(), Lower = list(), Upper = list())
  for (imp in imps) {
    m  <- as.matrix(irf_ortho$irf[[imp]])
    lo <- as.matrix(irf_ortho$Lower[[imp]])
    hi <- as.matrix(irf_ortho$Upper[[imp]])
    out$irf[[imp]]   <- apply(m,  2, cumsum)
    out$Lower[[imp]] <- apply(lo, 2, cumsum)
    out$Upper[[imp]] <- apply(hi, 2, cumsum)
    colnames(out$irf[[imp]])   <- colnames(m)
    colnames(out$Lower[[imp]]) <- colnames(lo)
    colnames(out$Upper[[imp]]) <- colnames(hi)
  }
  out
}

# Cholesky-IRF の累積図（応答/衝撃で2種類）
irf_ortho_cum <- irf_cumsum(irf_ortho)
plot_irf_gg(
  irf_ortho_cum,
  file_out = "var_outputs/IRF_ortho_cumulative_byResponse.png",
  facet_by = "response",
  title    = "Cumulative IRFs (Cholesky) - Facet by Response"
)
plot_irf_gg(
  irf_ortho_cum,
  file_out = "var_outputs/IRF_ortho_cumulative_byImpulse.png",
  facet_by = "impulse",
  title    = "Cumulative IRFs (Cholesky) - Facet by Impulse"
)

## ============ 12) Structural VAR (AB) *optional — CLEAN & ROBUST ============

suppressPackageStartupMessages({
  library(vars); library(ggplot2); library(dplyr); library(reshape2)
})

## 出力フォルダ
dir.create("var_outputs", showWarnings = FALSE)

## ---- 前提チェック：m_df が存在＆列名が想定通り？
if (!exists("m_df")) stop("m_df が見つかりません。上流のVAR構築パートを先に実行してください。")
stopifnot(identical(colnames(m_df), c("ddGDP","ddM2","ddInf","ddNex","ddRex")))

## ---- 変数順の入れ替え（識別用）：(GDP, Inf, M2, Rex, Nex)
m_df_AB <- cbind(m_df[,"ddGDP"], m_df[,"ddInf"], m_df[,"ddM2"], m_df[,"ddRex"], m_df[,"ddNex"])
colnames(m_df_AB) <- c("ddGDP","ddInf","ddM2","ddRex","ddNex")

## ---- IRFユーティリティ --------------------------------
irf_to_long <- function(irf_list) {
  stopifnot(!is.null(irf_list$irf))
  imps <- names(irf_list$irf)
  has_ci <- !is.null(irf_list$Lower) && !is.null(irf_list$Upper)
  out <- lapply(imps, function(imp) {
    base <- as.data.frame(irf_list$irf[[imp]]); base$period <- seq_len(nrow(base))
    base_long <- reshape2::melt(base, id.vars = "period",
                                variable.name = "response", value.name = "irf")
    if (has_ci) {
      lo <- as.data.frame(irf_list$Lower[[imp]]); hi <- as.data.frame(irf_list$Upper[[imp]])
      lo$period <- hi$period <- seq_len(nrow(lo))
      lo_long <- reshape2::melt(lo, id.vars = "period",
                                variable.name = "response", value.name = "low")
      hi_long <- reshape2::melt(hi, id.vars = "period",
                                variable.name = "response", value.name = "high")
      dplyr::left_join(base_long, lo_long, by = c("period","response")) |>
        dplyr::left_join(hi_long, by = c("period","response")) |>
        dplyr::mutate(impulse = imp)
    } else {
      base_long |>
        dplyr::mutate(low = NA_real_, high = NA_real_, impulse = imp)
    }
  })
  dplyr::bind_rows(out)
}

plot_irf_gg <- function(irf_list, file_out,
                        facet_by = c("response","impulse"),
                        title    = "Impulse Responses") {
  facet_by <- match.arg(facet_by)
  df <- irf_to_long(irf_list)
  p <- ggplot(df, aes(x = period, y = irf)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line(aes(color = impulse), linewidth = 0.9) +
    labs(title = title, x = "Horizon", y = "Response") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  # CIがあればリボン
  if (any(is.finite(df$low) & is.finite(df$high))) {
    p <- p + geom_ribbon(aes(ymin = low, ymax = high,
                             group = interaction(impulse, response)),
                         alpha = 0.18)
  }
  if (facet_by == "response") p <- p + facet_wrap(~response, scales = "free_y")
  else                        p <- p + facet_wrap(~impulse,  scales = "free_y")
  ggsave(filename = file_out, plot = p, width = 11, height = 7, dpi = 300)
  message("[Saved] ", file_out)
}

irf_cumsum <- function(irf_obj) {
  out <- irf_obj
  out$irf   <- lapply(out$irf,   function(M) apply(M, 2, cumsum))
  if (!is.null(out$Lower)) out$Lower <- lapply(out$Lower, function(M) apply(M, 2, cumsum))
  if (!is.null(out$Upper)) out$Upper <- lapply(out$Upper, function(M) apply(M, 2, cumsum))
  out
}

## ---- ラグ選択とVAR推定（基準） --------------------------------------------
sel_AB <- VARselect(m_df_AB, lag.max = 4, type = "const")
p_AB   <- sel_AB$selection[["AIC(n)"]]; if (is.null(p_AB)) p_AB <- 2
var_AB <- VAR(m_df_AB, p = p_AB, type = "const", season = 4)

## ---- A/B制約（1=自己, 0=排除, NA=推定） -----------------------------------
K <- NCOL(m_df_AB)
A <- matrix(c(
  1,   0,   0,   0,   0,   # GDP eq: 同期で他の影響を受けない
  NA,  1,   0,   NA,  0,   # Inf eq: GDP, Rex から同期影響
  NA,  NA,  1,   0,   0,   # M2  eq: GDP, Inf から同期影響
  0,   0,   NA,  1,   0,   # Rex eq: M2 から同期影響
  NA,  0,   0,   NA,  1    # Nex eq: GDP, Rex から同期影響
), nrow = K, byrow = TRUE)
B <- diag(NA_real_, K)
storage.mode(A) <- "double"; storage.mode(B) <- "double"

## ---- SVARの堅牢推定：複数パターン自動試行 → 失敗時はCholeskyへ -------------
try_fit_svar <- function(Y, A, B, p_candidates = c(4,3,2,1),
                         use_season = TRUE,
                         methods = c("direct","scoring"),
                         standardize = c(FALSE, TRUE)) {
  errs <- list()
  for (std in standardize) {
    YY <- if (std) ts(scale(Y), start = start(Y), frequency = frequency(Y)) else Y
    for (p in p_candidates) {
      var_obj <- try(vars::VAR(YY, p = p, type = "const", season = if (use_season) 4 else NULL), silent = TRUE)
      if (inherits(var_obj, "try-error")) { errs <- append(errs, list(var_obj)); next }
      for (m in methods) {
        fit <- try(vars::SVAR(x = var_obj, Amat = A, Bmat = B,
                              estmethod = m, max.iter = 5000, conv.crit = 1e-8),
                   silent = TRUE)
        if (!inherits(fit, "try-error")) {
          attr(fit, "std") <- std; attr(fit, "p") <- p; attr(fit, "method") <- m; attr(fit, "var") <- var_obj
          return(fit)
        } else {
          errs <- append(errs, list(fit))
        }
      }
    }
  }
  attr(errs, "failed") <- TRUE
  errs
}

svar_ab <- try_fit_svar(m_df_AB, A, B,
                        p_candidates = c(4,3,2,1),
                        use_season   = TRUE,
                        methods      = c("direct","scoring"),
                        standardize  = c(FALSE, TRUE))

if (is.list(svar_ab) && isTRUE(attr(svar_ab, "failed"))) {
  message("[WARN] AB 識別は不成功。Cholesky（再帰型）にフォールバックします。")
  sel_AB <- vars::VARselect(m_df_AB, lag.max = 4, type = "const")
  p_AB   <- sel_AB$selection[["AIC(n)"]]; if (is.null(p_AB)) p_AB <- 2
  var_AB <- vars::VAR(m_df_AB, p = p_AB, type = "const", season = 4)
  A_chol <- diag(1, K); A_chol[lower.tri(A_chol)] <- NA   # 下三角NA：再帰構造
  B_diag <- diag(NA_real_, K)
  svar_ab <- vars::SVAR(x = var_AB, Amat = A_chol, Bmat = B_diag,
                        estmethod = "direct", max.iter = 3000, conv.crit = 1e-8)
  attr(svar_ab, "fallback") <- "Cholesky"
}

## ------------- 強制フォールバック付きの安全な SVAR / VAR-IRF 生成 -------------
make_struct_irf <- function(m_df_AB, svar_ab, var_AB = NULL,
                            n_ahead = 12, want_ci = TRUE, season = TRUE) {
  if (inherits(svar_ab, "SVAR")) {
    message("[INFO] Using AB/Cholesky SVAR for IRF.")
    if (want_ci) {
      return(vars::irf(svar_ab, n.ahead = n_ahead, boot = TRUE,  runs = 800, ci = 0.95))
    } else {
      return(vars::irf(svar_ab, n.ahead = n_ahead, boot = FALSE))
    }
  }
  message("[WARN] SVAR object not available. Falling back to reduced-form VAR IRF (Cholesky).")
  if (!inherits(var_AB, "varest")) {
    # ラグ選択（失敗時 p=1）
    sel_tmp <- try(vars::VARselect(m_df_AB, lag.max = 4, type = "const"), silent = TRUE)
    p_tmp <- if (!inherits(sel_tmp, "try-error")) sel_tmp$selection[["AIC(n)"]] else 1
    if (is.null(p_tmp) || is.na(p_tmp)) p_tmp <- 1
    
    var_AB <- try(
      vars::VAR(m_df_AB, p = p_tmp, type = "const", season = if (season) 4 else NULL),
      silent = TRUE
    )
    if (inherits(var_AB, "try-error")) {
      # さらに堅牢化：季節項なし + p=1 + 標準化
      message("[WARN] VAR re-fit failed. Retrying with p=1, no season, standardized series.")
      m_df_AB_std <- ts(scale(m_df_AB), start = start(m_df_AB), frequency = frequency(m_df_AB))
      var_AB <- vars::VAR(m_df_AB_std, p = 1, type = "const", season = NULL)
    }
  }
   # 直交化（Cholesky）IRFを返す
  if (want_ci) {
    vars::irf(var_AB, n.ahead = n_ahead, ortho = TRUE, boot = TRUE,  runs = 800, ci = 0.95)
  } else {
    vars::irf(var_AB, n.ahead = n_ahead, ortho = TRUE, boot = FALSE)
  }
}


# CIなし（安定版)
irf_struct <- make_struct_irf(
  m_df_AB  = m_df_AB,
  svar_ab  = if (exists("svar_ab")) svar_ab else NULL,
  var_AB   = if (exists("var_AB"))  var_AB  else NULL,
  n_ahead  = 12,
  want_ci  = FALSE,
  season   = FALSE
)

# baseプロット
png("var_outputs/IRF_Struct_base.png", width = 1400, height = 900, res = 130)
plot(irf_struct, main = "Structural / Orthogonalized IRFs (robust fallback)")
dev.off()

## ========= IRF baseプロット（irf_struct が有効なときだけ） =========
if (exists("irf_struct") && is.list(irf_struct) && !is.null(irf_struct$irf)) {
  dir.create("var_outputs", showWarnings = FALSE)
  png("var_outputs/IRF_Struct_base.png", width = 1400, height = 900, res = 130)
  plot(irf_struct, main = "Structural / Orthogonalized IRFs (robust fallback)")
  dev.off()
} else {
  message("[INFO] irf_struct が未定義 or 内容が空のため、base図はスキップします。")
}

## ========= CI付きIRFの生成を try で保護 =========
irf_struct_ci <- try(
  make_struct_irf(
    m_df_AB  = m_df_AB,
    svar_ab  = if (exists("svar_ab")) svar_ab else NULL,
    var_AB   = if (exists("var_AB"))  var_AB  else NULL,
    n_ahead  = 12,
    want_ci  = TRUE,
    season   = TRUE
  ),
  silent = TRUE
)

## ========= CI付きIRFの ggplot 保存（成功時のみ） =========
ok_ci <- ( !inherits(irf_struct_ci, "try-error") &&
             is.list(irf_struct_ci) &&
             !is.null(irf_struct_ci$irf) )

if (ok_ci) {
  plot_irf_gg(
    irf_list = irf_struct_ci,
    file_out = "var_outputs/IRF_Struct_gg_byResponse.png",
    facet_by = "response",
    title    = "Structural / Orthogonalized IRFs (with CI) - Facet by Response"
  )
  plot_irf_gg(
    irf_list = irf_cumsum(irf_struct_ci),
    file_out = "var_outputs/IRF_Struct_cumulative_byResponse.png",
    facet_by = "response",
    title    = "Cumulative Structural / Orthogonalized IRFs (with CI)"
  )
} else {
  message("[INFO] CI IRF skipped（bootstrap 失敗 or オブジェクト不備）。")
  
  ## --- 代替：CIなしのIRFで ggplot 図を生成（使える材料があれば） ---
  # irf_struct（base用）の中身が使えるなら、それを流用
  if (exists("irf_struct") && is.list(irf_struct) && !is.null(irf_struct$irf)) {
    plot_irf_gg(
      irf_list = irf_struct,
      file_out = "var_outputs/IRF_Struct_gg_byResponse_noCI.png",
      facet_by = "response",
      title    = "Structural / Orthogonalized IRFs (no CI) - Facet by Response"
    )
    plot_irf_gg(
      irf_list = irf_cumsum(irf_struct),
      file_out = "var_outputs/IRF_Struct_cumulative_byResponse_noCI.png",
      facet_by = "response",
      title    = "Cumulative Structural / Orthogonalized IRFs (no CI)"
    )
  } else {
    message("[INFO] no-CI IRF も材料不足でスキップしました。")
  }
}

# FEVDは、SVARが無い場合縮約形VARのFEVD（orthogonalized）で代用
fevd_out <- NULL
if (exists("svar_ab") && inherits(svar_ab, "SVAR")) {
  fevd_out <- vars::fevd(svar_ab, n.ahead = 12)
} else if (exists("var_AB") && inherits(var_AB, "varest")) {
  fevd_out <- vars::fevd(var_AB, n.ahead = 12)
}
if (!is.null(fevd_out)) {
  capture.output(fevd_out, file = "var_outputs/fevd_struct_or_var.txt")
  cat("\nFEVD saved to var_outputs/fevd_struct_or_var.txt\n")
} else {
  message("[INFO] FEVD not available.")
}

if (is.null(svar_ab$call$estmethod)) svar_ab$call$estmethod <- "direct"
cat("\n--- SVAR summary ---\n"); print(summary(svar_ab))

## ---- IRF作成：まず安定重視（boot=FALSE）、次にCI付き（可能なら） ----------
irf_svar_ab <- vars::irf(svar_ab, n.ahead = 12, boot = FALSE)

# baseプロット
png("var_outputs/IRF_SVAR_AB_base.png", width=1400, height=900, res=130)
plot(irf_svar_ab, main = ifelse(is.null(attr(svar_ab,"fallback")),
                                "Structural IRFs (AB)", "Structural IRFs (Cholesky fallback)"))
dev.off()


# CI付きIRF
set.seed(123)
irf_svar_ab_ci <- try(vars::irf(svar_ab, n.ahead = 12, boot = TRUE, runs = 800, ci = 0.95),
                      silent = FALSE)
if (!inherits(irf_svar_ab_ci, "try-error")) {
  plot_irf_gg(
    irf_list = irf_svar_ab_ci,
    file_out = "var_outputs/IRF_SVAR_AB_gg_byResponse.png",
    facet_by = "response",
    title    = ifelse(is.null(attr(svar_ab,"fallback")),
                      "Structural IRFs (AB) - Facet by Response",
                      "Structural IRFs (Cholesky) - Facet by Response")
  )
  # 累積IRFも保存
  plot_irf_gg(
    irf_list = irf_cumsum(irf_svar_ab_ci),
    file_out = "var_outputs/IRF_SVAR_AB_cumulative_byResponse.png",
    facet_by = "response",
    title    = ifelse(is.null(attr(svar_ab,"fallback")),
                      "Cumulative Structural IRFs (AB) - Facet by Response",
                      "Cumulative Structural IRFs (Cholesky) - Facet by Response")
  )
} else {
  message("[INFO] ブートCI付きIRFはスキップ（計算失敗/時間都合）")
}

## ---- FEVD（12期） ----------------------------------------------------------
fevd_svar_ab <- vars::fevd(svar_ab, n.ahead = 12)
capture.output(fevd_svar_ab, file = "var_outputs/fevd_svar_ab.txt")
cat("\nSVAR FEVD saved to var_outputs/fevd_svar_ab.txt\n")

## オブジェクト保存
saveRDS(list(
  svar_ab        = svar_ab,
  var_AB         = if (!is.null(attr(svar_ab, "var"))) attr(svar_ab, "var") else var_AB,
  irf_svar_ab    = irf_svar_ab,
  irf_svar_ab_ci = if (exists("irf_svar_ab_ci") && !inherits(irf_svar_ab_ci,"try-error")) irf_svar_ab_ci else NULL,
  fevd_svar_ab   = fevd_svar_ab
), file = "var_outputs/svar_results.rds")

cat("\n Done. See 'var_outputs/' for PNGs & RDS.\n")
