# Brace-aware block extraction, validation, and rewrites for
# custom-family-generated Stan code (Phase 1 + Phase 2 + Phase 3).

# -- Brace-aware helpers -------------------------------------------------------

extract_named_block <- function(stancode, block_name) {
  pattern <- paste0("\\b", block_name, "\\s*\\{")
  match <- regexpr(pattern, stancode)
  if (match == -1L)
    cli::cli_abort("Block {.val {block_name}} not found in Stan code.")

  brace_offset <- regexpr("\\{", substring(stancode, match))
  start_pos <- as.integer(match) + brace_offset - 1L

  chars <- strsplit(substring(stancode, start_pos), "")[[1]]
  depth <- 0L
  for (i in seq_along(chars)) {
    if (chars[i] == "{") depth <- depth + 1L
    if (chars[i] == "}") depth <- depth - 1L
    if (depth == 0L) {
      return(substring(stancode, start_pos + 1L, start_pos + i - 2L))
    }
  }
  cli::cli_abort("Unmatched braces in {.val {block_name}} block.")
}

replace_named_block <- function(stancode, block_name, new_content) {
  pattern <- paste0("\\b", block_name, "\\s*\\{")
  match <- regexpr(pattern, stancode)
  if (match == -1L)
    cli::cli_abort("Block {.val {block_name}} not found in Stan code.")

  brace_offset <- regexpr("\\{", substring(stancode, match))
  start_pos <- as.integer(match) + brace_offset - 1L

  chars <- strsplit(substring(stancode, start_pos), "")[[1]]
  depth <- 0L
  for (i in seq_along(chars)) {
    if (chars[i] == "{") depth <- depth + 1L
    if (chars[i] == "}") depth <- depth - 1L
    if (depth == 0L) {
      end_pos <- start_pos + i - 1L
      before <- substring(stancode, 1L, start_pos)
      after  <- substring(stancode, end_pos)
      return(paste0(before, new_content, after))
    }
  }
  cli::cli_abort("Unmatched braces in {.val {block_name}} block.")
}

# -- Validation ----------------------------------------------------------------

validate_subsampled_surface <- function(model_block, non_mu_dpars = character(0)) {
  if (!grepl("vector\\[N\\]\\s+mu\\s*=\\s*rep_vector\\(0\\.0,\\s*N\\)", model_block))
    cli::cli_abort("Unsupported model: mu initialization pattern not found.")

  has_fe     <- grepl("\\bXc\\s*\\*\\s*b\\b", model_block)
  has_sp_lin <- grepl("\\bXs\\s*\\*\\s*bs\\b", model_block)
  has_sp_bas <- grepl("\\bZs_\\d+_\\d+\\s*\\*\\s*s_\\d+_\\d+\\b", model_block)
  has_gp     <- grepl("\\bgp_pred_\\d+\\[Jgp_\\d+\\]", model_block)
  has_re     <- grepl("r_\\d+_\\d+\\[J_\\d+\\[", model_block)

  has_dpar_predictor <- FALSE
  for (dp in non_mu_dpars) {
    has_dp_fe <- grepl(paste0("\\bXc_", dp, "\\s*\\*\\s*b_", dp, "\\b"), model_block)
    has_dp_sp <- grepl(paste0("\\bXs_", dp, "\\s*\\*\\s*bs_", dp, "\\b"), model_block) ||
                 grepl(paste0("\\bZs_", dp, "_\\d+_\\d+\\s*\\*\\s*s_", dp, "_\\d+_\\d+"), model_block)
    has_dp_gp <- grepl(paste0("\\bgp_pred_", dp, "_\\d+\\[Jgp_", dp, "_\\d+\\]"), model_block)
    has_dp_re <- grepl(paste0("r_\\d+_", dp, "_\\d+\\[J_\\d+\\["), model_block)
    if (has_dp_fe || has_dp_sp || has_dp_gp || has_dp_re)
      has_dpar_predictor <- TRUE
  }

  if (!has_fe && !has_sp_lin && !has_sp_bas && !has_gp && !has_re && !has_dpar_predictor)
    cli::cli_abort("Unsupported model: no supported predictor accumulation found (expected Xc*b, Xs*bs, spline bases, GP terms, random effects, or distributional predictors).")

  if (has_re) validate_re_loop_shape(model_block)
  validate_dpar_re_loops(model_block, non_mu_dpars)

  unknown_dpars <- c("sigma", "shape", "hu", "zi", "disc", "quantile", "phi")
  unknown_dpars <- setdiff(unknown_dpars, non_mu_dpars)
  if (length(unknown_dpars) > 0) {
    unknown_pattern <- paste0("\\bXc_(", paste(unknown_dpars, collapse = "|"), ")\\b")
    if (grepl(unknown_pattern, model_block))
      cli::cli_abort("Unsupported model: distributional predictor matrices detected for unknown dpars.")
  }

  unsupported <- c(
    "\\bcens\\b"   = "censoring",
    "\\btrunc\\b"  = "truncation"
  )
  for (i in seq_along(unsupported)) {
    if (grepl(names(unsupported)[[i]], model_block))
      cli::cli_abort("Unsupported model: {unsupported[[i]]} detected.")
  }

  mu_matches <- gregexpr("vector\\[N\\]\\s+mu\\s*=", model_block)[[1]]
  mu_count <- if (mu_matches[1] == -1L) 0L else length(mu_matches)
  if (mu_count != 1L)
    cli::cli_abort("Unsupported model: expected exactly one mu initialization, found {mu_count}.")

  for (dp in non_mu_dpars) {
    dp_init_pattern <- paste0("vector\\[N\\]\\s+", dp, "\\s*=")
    dp_matches <- gregexpr(dp_init_pattern, model_block)[[1]]
    dp_count <- if (dp_matches[1] == -1L) 0L else length(dp_matches)
    if (dp_count > 1L)
      cli::cli_abort("Unsupported model: expected at most one {dp} initialization, found {dp_count}.")
  }

  invisible(NULL)
}

validate_re_loop_shape <- function(model_block) {
  if (!grepl("for \\(n in 1:N\\)", model_block))
    cli::cli_abort("Unsupported model: random effects present but no for(n in 1:N) loop found.")

  all_loops <- gregexpr("for \\(n in 1:N\\) \\{[^}]*\\}", model_block)[[1]]
  if (all_loops[1] == -1L)
    cli::cli_abort("Unsupported model: could not extract random effects loop.")

  for (j in seq_along(all_loops)) {
    loop_text <- substr(model_block, all_loops[j],
                        all_loops[j] + attr(all_loops, "match.length")[j] - 1L)
    if (!grepl("mu\\[n\\]\\s*\\+=", loop_text)) next
    if (!grepl("r_\\d+_\\d+\\[J_\\d+\\[n\\]\\]\\s*\\*\\s*Z_\\d+_\\d+\\[n\\]", loop_text))
      cli::cli_abort("Unsupported model: mu random effects loop does not match expected r_G_K[J_G[n]] * Z_G_K[n] pattern.")
  }

  invisible(NULL)
}

validate_dpar_re_loops <- function(model_block, non_mu_dpars) {
  all_loops <- gregexpr("for \\(n in 1:N\\) \\{[^}]*\\}", model_block)[[1]]
  if (all_loops[1] == -1L) return(invisible(NULL))

  for (j in seq_along(all_loops)) {
    loop_text <- substr(model_block, all_loops[j],
                        all_loops[j] + attr(all_loops, "match.length")[j] - 1L)
    for (dp in non_mu_dpars) {
      if (!grepl(paste0(dp, "\\[n\\]\\s*\\+="), loop_text)) next
      re_pattern <- paste0("r_\\d+_", dp, "_\\d+\\[J_\\d+\\[n\\]\\]\\s*\\*\\s*Z_\\d+_", dp, "_\\d+\\[n\\]")
      if (!grepl(re_pattern, loop_text))
        cli::cli_abort("Unsupported model: {dp} random effects loop does not match expected r_G_{dp}_K[J_G[n]] * Z_G_{dp}_K[n] pattern.")
    }
  }

  invisible(NULL)
}

# -- Phase 1 rewrites ---------------------------------------------------------

rewrite_mu_init <- function(code) {
  original <- code
  code <- sub(
    "([ \t]*)(vector\\[N\\] mu = rep_vector\\(0\\.0, N\\);)",
    "\\1int m_sub = pdmp_get_subsample_size();\n\\1vector[m_sub] mu = rep_vector(0.0, m_sub);",
    code
  )
  if (identical(code, original))
    cli::cli_abort("Unsupported model: mu initialization pattern not found.")
  code
}

rewrite_fe_accumulation <- function(code) {
  gsub("\\bXc\\b(\\s*\\*\\s*b\\b)", "get_subsampled_Xc(Xc)\\1", code)
}

rewrite_spline_matrices <- function(code) {
  code <- gsub("\\bXs\\b(\\s*\\*\\s*bs\\b)", "get_subsampled_Xc(Xs)\\1", code)
  gsub("\\b(Zs_\\d+_\\d+)\\b(\\s*\\*\\s*s_\\d+_\\d+\\b)", "get_subsampled_Xc(\\1)\\2", code)
}

rewrite_gp_indexing <- function(code) {
  gsub(
    "\\b(gp_pred_\\d+)\\[(Jgp_\\d+)\\]",
    "\\1[get_subsampled_int_array(\\2)]",
    code
  )
}

# -- Phase 2 rewrites (random effects) -----------------------------------------

rewrite_re_loops <- function(code, non_mu_dpars = character(0)) {
  if (!grepl("for \\(n in 1:N\\)", code)) return(code)

  code <- gsub(
    "([ \t]*)for \\(n in 1:N\\) \\{\n",
    "\\1for (i in 1:m_sub) {\n\\1  int n = pdmp_get_subsample_index(i);\n",
    code
  )
  code <- gsub("\\bmu\\[n\\]", "mu[i]", code)
  for (dp in non_mu_dpars) {
    code <- gsub(paste0("\\b", dp, "\\[n\\]"), paste0(dp, "[i]"), code)
  }
  code
}

# -- Phase 3 rewrites (distributional parameters) -----------------------------

rewrite_dpar_init <- function(code, dpar) {
  pattern <- paste0(
    "([ \t]*)(vector\\[N\\] ", dpar,
    " = rep_vector\\(0\\.0, N\\);)"
  )
  replacement <- paste0(
    "\\1vector[m_sub] ", dpar,
    " = rep_vector(0.0, m_sub);"
  )
  sub(pattern, replacement, code)
}

rewrite_dpar_fe <- function(code, dpar) {
  pattern <- paste0("\\bXc_", dpar, "\\b(\\s*\\*\\s*b_", dpar, "\\b)")
  replacement <- paste0("get_subsampled_Xc(Xc_", dpar, ")\\1")
  gsub(pattern, replacement, code)
}

rewrite_dpar_splines <- function(code, dpar) {
  pattern_xs <- paste0("\\bXs_", dpar, "\\b(\\s*\\*\\s*bs_", dpar, "\\b)")
  code <- gsub(pattern_xs, paste0("get_subsampled_Xc(Xs_", dpar, ")\\1"), code)
  pattern_zs <- paste0("\\b(Zs_", dpar, "_\\d+_\\d+)\\b(\\s*\\*\\s*s_", dpar, "_\\d+_\\d+\\b)")
  gsub(pattern_zs, "get_subsampled_Xc(\\1)\\2", code)
}

rewrite_dpar_gp_indexing <- function(code, dpar) {
  pattern <- paste0("\\b(gp_pred_", dpar, "_\\d+)\\[(Jgp_", dpar, "_\\d+)\\]")
  gsub(pattern, "\\1[get_subsampled_int_array(\\2)]", code)
}

# -- Main entry point ----------------------------------------------------------

validate_and_rewrite_subsampled_code <- function(stancode, non_mu_dpars = character(0)) {
  model_block <- extract_named_block(stancode, "model")
  validate_subsampled_surface(model_block, non_mu_dpars)
  model_block <- rewrite_mu_init(model_block)
  model_block <- rewrite_fe_accumulation(model_block)
  model_block <- rewrite_spline_matrices(model_block)
  model_block <- rewrite_gp_indexing(model_block)
  for (dp in non_mu_dpars) {
    dp_init_pattern <- paste0("vector\\[N\\]\\s+", dp, "\\s*=\\s*rep_vector\\(0\\.0,\\s*N\\)")
    if (grepl(dp_init_pattern, model_block)) {
      model_block <- rewrite_dpar_init(model_block, dp)
      model_block <- rewrite_dpar_fe(model_block, dp)
      model_block <- rewrite_dpar_splines(model_block, dp)
      model_block <- rewrite_dpar_gp_indexing(model_block, dp)
    }
  }
  model_block <- rewrite_re_loops(model_block, non_mu_dpars)
  replace_named_block(stancode, "model", model_block)
}
