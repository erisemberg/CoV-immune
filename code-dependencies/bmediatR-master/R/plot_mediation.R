#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# Conversion function taken from qtl2convert
map_df_to_list <- function (map, chr_column = "chr", pos_column = "cM", marker_column = "marker", 
          Xchr = c("x", "X")) {
  
  if (is.null(marker_column)) {
    marker_column <- "qtl2tmp_marker"
    map[, marker_column] <- rownames(map)
  }
  if (!(marker_column %in% colnames(map))) 
    stop("Column \"", marker_column, "\" not found.")
  if (!(chr_column %in% colnames(map))) 
    stop("Column \"", chr_column, "\" not found.")
  if (!(pos_column %in% colnames(map))) 
    stop("Column \"", pos_column, "\" not found.")
  marker <- map[, marker_column]
  chr <- map[, chr_column]
  uchr <- unique(chr)
  pos <- map[, pos_column]
  result <- split(as.numeric(pos), factor(chr, levels = uchr))
  marker <- split(marker, factor(chr, levels = uchr))
  for (i in seq(along = result)) names(result[[i]]) <- marker[[i]]
  is_x_chr <- rep(FALSE, length(result))
  names(is_x_chr) <- names(result)
  if (!is.null(Xchr)) {
    Xchr_used <- Xchr %in% names(is_x_chr)
    if (any(Xchr_used)) {
      Xchr <- Xchr[Xchr_used]
      is_x_chr[Xchr] <- TRUE
    }
  }
  attr(result, "is_x_chr") <- is_x_chr
  result
}

#' Barplot of posterior model probabilities
#'
#' This function takes posterior probability results from mediation_bf() and plots barplots of posterior model probabilities.
#'
#' @param bmediatR_object Output from bmediatR(). 
#' @param med_annot Annotation data for -omic mediators.
#' @param mediator_id Which mediator to plot posterior model probabilities for.
#' @param med_var DEFAULT: "protein.id". The column in med_annot to be used as a mediator id.
#' @param stack DEFAULT: FALSE. If TRUE, a stacked barplot is produced. If FALSE, a staggered
#' barplot is produced.
#' @param add_number_labels DEFAULT: FALSE. Add posterior probabilities as text above bars.
#' @param label_size DEFAULT: 5. Text size of probability labels, when included.
#' @param num_dig DEFAULT: 3. The number of digits after the decimal point if probabilities
#' are being included as labels.
#' @export
#' @examples plot_posterior_bar()
plot_posterior_bar <- function(bmediatR_object,
                               med_annot = NULL,
                               mediator_id,
                               med_var = "protein.id",
                               stack = FALSE,
                               bar_col = c("seagreen4", "seagreen1", "skyblue", "goldenrod1", "goldenrod4", "gray"),
                               relabel_x = NULL,
                               add_number_labels = FALSE,
                               label_size = 5,
                               num_dig = 3,
                               main = NULL) {
  
  ## Flag for reactive model
  post_mat <- bmediatR_object$ln_post_c[1,, drop = FALSE]
  model_flag <- is.finite(post_mat)
  names(model_flag) <- colnames(post_mat)
  
  ## Long names
  long_names <- c("other non-med", "other non-med", "other non-med", 
                  "complete med", "other non-med", "other non-med",
                  "co-local", "partial med", "other non-med",
                  "complete med (react)", "other non-med", "partial med (react)")
  names(long_names) <- colnames(post_mat)
  
  bar_col <- bar_col[c(model_flag[c("1,1,0", "1,1,1", "1,0,1", "1,*,1", "0,*,1")], TRUE)]

  posterior_dat <- exp(bmediatR_object$ln_post_c) %>%
    as.data.frame %>%
    tibble::rownames_to_column(med_var) %>%
    dplyr::rename("partial med" = `1,1,1`,
                  "complete med" = `1,1,0`,
                  "co-local" = `1,0,1`,
                  "partial med (react)" = `1,*,1`,
                  "complete med (react)" = `0,*,1`)
  
  # Using annotation information
  if (!is.null(med_annot)) {
    posterior_dat <- posterior_dat %>%
      dplyr::left_join(med_annot %>%
                         dplyr::select(tidyselect::all_of(med_var), symbol))
  } else {
    posterior_dat <- posterior_dat %>%
      dplyr::mutate(symbol = get(med_var))
  }
  
  ## Calculating non-mediation or co-local probability
  posterior_dat <- posterior_dat %>%
    dplyr::left_join(posterior_dat %>%
                       dplyr::select(tidyselect::all_of(med_var), contains(",")) %>%
                       dplyr::mutate("other non-med" = rowSums(.[-1]))) %>%
    dplyr::select(-contains(",")) %>%
    tidyr::gather(key = model, value = post_p, -c(tidyselect::all_of(med_var), symbol))

  ## Set factors
  models_use <- unique(long_names[model_flag])
  posterior_dat <- posterior_dat %>%
    dplyr::filter(model %in% models_use) %>%
    dplyr::mutate(model = factor(model, levels = c("complete med", "partial med", "co-local", "partial med (react)", "complete med (react)", "other non-med")))
  
  bar_theme <- ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                              panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(), 
                              axis.line = ggplot2::element_line(colour = "black"),
                              plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain"), 
                              axis.title.x = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_text(hjust = 0.5, size = 14, face = "plain"),
                              axis.title.y = ggplot2::element_text(size = 14, face = "plain"),
                              axis.text.y = ggplot2::element_text(size = 14, face = "plain"),
                              axis.ticks.y = ggplot2::element_blank(),
                              legend.title = ggplot2::element_text(size = 14),
                              legend.text = ggplot2::element_text(size = 14))
  
  if (!is.null(relabel_x)) {
    posterior_dat$symbol <- relabel_x
  }
  p <- ggplot2::ggplot(data = posterior_dat %>% 
                         dplyr::filter((!!as.symbol(med_var)) == mediator_id)) +
    ggplot2::scale_fill_manual(values = bar_col) +
    ggplot2::ylab("Posterior model probability") + 
    bar_theme
  if (!is.null(main)) {
    p <- p + ggplot2::ggtitle(main)
  }
  if (stack) {
    p <- p + ggplot2::geom_col(ggplot2::aes(x = symbol, y = post_p, fill = model), width = 0.5) 
  } else {
    p <- p + ggplot2::geom_bar(ggplot2::aes(x = symbol, y = post_p, fill = model), position = "dodge", stat = "summary", fun = "mean") +
      ggplot2::geom_hline(yintercept = c(0, 1), col = "gray", linetype = "dashed")
  }
  if (add_number_labels) {
    p <- p + geom_text(data = posterior_dat %>% 
                         dplyr::filter((!!as.symbol(med_var)) == mediator_id), 
                       aes(symbol, post_p + 0.04, group = model, label = round(post_p, digits = num_dig), col = model), 
                       position = position_dodge(width = 0.9),
                       size = label_size) +
      ggplot2::scale_color_manual(values = bar_col)
  }
  
  p
}

#' Barplot of prior model probabilities
#'
#' This function takes posterior probability results from mediation_bf() and plots barplots of prior model probabilities.
#'
#' @param bmediatR_object Output from bmediatR(). 
#' @param stack DEFAULT: FALSE. If TRUE, a stacked barplot is produced. If FALSE, a staggered
#' barplot is produced.
#' @param add_number_labels DEFAULT: FALSE. Add posterior probabilities as text above bars.
#' @param label_size DEFAULT: 5. Text size of probability labels, when included.
#' @param num_dig DEFAULT: 3. The number of digits after the decimal point if probabilities
#' are being included as labels.
#' @export
#' @examples plot_prior_bar()
plot_prior_bar <- function(bmediatR_object,
                           stack = FALSE,
                           bar_col = c("seagreen4", "seagreen1", "skyblue", "goldenrod1", "goldenrod4", "gray"),
                           relabel_x = NULL,
                           add_number_labels = FALSE,
                           label_size = 5,
                           num_dig = 3,
                           main = NULL) {
  
  ## Flag for reactive model
  prior_mat <- bmediatR_object$ln_prior_c
  model_flag <- is.finite(prior_mat)
  names(model_flag) <- colnames(prior_mat)
  
  ## Long names
  long_names <- c("other non-med", "other non-med", "other non-med", 
                  "complete med", "other non-med", "other non-med",
                  "co-local", "partial med", "other non-med",
                  "complete med (react)", "other non-med", "partial med (react)")
  names(long_names) <- colnames(prior_mat)
  
  bar_col <- bar_col[c(model_flag[c("1,1,0", "1,1,1", "1,0,1", "1,*,1", "0,*,1")], TRUE)]
  
  prior_dat <- exp(bmediatR_object$ln_prior_c) %>%
    as.data.frame %>%
    dplyr::rename("partial med" = `1,1,1`,
                  "complete med" = `1,1,0`,
                  "co-local" = `1,0,1`,
                  "partial med (react)" = `1,*,1`,
                  "complete med (react)" = `0,*,1`)
  
  ## Calculating non-mediation or co-local probability
  prior_dat <- prior_dat %>%
    dplyr::left_join(prior_dat %>%
                       dplyr::select(contains(",")) %>%
                       dplyr::mutate("other non-med" = rowSums(.[-1]))) %>%
    dplyr::select(-contains(",")) %>%
    tidyr::gather(key = model, value = prior_p)
  
  ## Set factors
  models_use <- unique(long_names[model_flag])
  prior_dat <- prior_dat %>%
    dplyr::filter(model %in% models_use) %>%
    dplyr::mutate(model = factor(model, levels = c("complete med", "partial med", "co-local", "partial med (react)", "complete med (react)", "other non-med")),
                  case = "Prior")
  
  bar_theme <- ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                              panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(), 
                              axis.line = ggplot2::element_line(colour = "black"),
                              plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "plain"), 
                              axis.title.x = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_text(hjust = 0.5, size = 14, face = "plain"),
                              axis.title.y = ggplot2::element_text(size = 14, face = "plain"),
                              axis.text.y = ggplot2::element_text(size = 14, face = "plain"),
                              axis.ticks.y = ggplot2::element_blank(),
                              legend.title = ggplot2::element_text(size = 14),
                              legend.text = ggplot2::element_text(size = 14))
  
  p <- ggplot2::ggplot(data = prior_dat) +
    ggplot2::scale_fill_manual(values = bar_col) +
    ggplot2::ylab("Prior model probability") + 
    bar_theme
  if (!is.null(main)) {
    p <- p + ggplot2::ggtitle(main)
  }
  if (stack) {
    p <- p + ggplot2::geom_col(ggplot2::aes(x = case, y = prior_p, fill = model), width = 0.5) 
  } else {
    p <- p + ggplot2::geom_bar(ggplot2::aes(x = case, y = prior_p, fill = model), position = "dodge", stat = "summary", fun = "mean") +
      ggplot2::geom_hline(yintercept = c(0, 1), col = "gray", linetype = "dashed")
  }
  if (add_number_labels) {
    p <- p + geom_text(data = prior_dat,
                       aes(symbol, prior_p + 0.04, group = model, label = round(post_p, digits = num_dig), col = model), 
                       position = position_dodge(width = 0.9),
                       size = label_size) +
      ggplot2::scale_color_manual(values = bar_col)
  }
  
  p
}

#' Posterior odds genome plot function
#'
#' This function takes the posterior odds results from bmediatR() and plots the genome-wide scan.
#'
#' @param bmediatR_object Output from bmediatR(). 
#' @param model_type DEFAULT: "mediation". Specifies which model(s)'s posterior probabilities are to be included in the numerator of the posterior odds and then displayed for
#' for genome-wide mediators. 
#' @param med_annot Annotation data for -omic mediators. This data provides the coordinate for plotting the mediation odds summary and is expected as a data.frame with the following variables:
#' \itemize{
#' \item{"med_var" identifier variable -- }{Mediator IDs (e.g., ENSEMBL gene IDs) that uniquely identify evaluated mediators. Colname should match the med_var argument and the IDs in bmediatR_object.}
#' \item{symbol -- }{Alternative IDs that are more readily meaningful to humans (e.g., gene symbols). DO not have to be unique.}
#' \item{chr -- }{Chromosome of the candidate mediator.}
#' \item{middle -- }{Genomic coordinate of the candidate mediator (e.g., midpoint of a coding gene).}
#' }
#' @param med_var Colname of unique identifier for candidate mediators (e.g., ENSEMBL gene IDs) in med_annot data.frame. Should match IDs in bmediatR_object. 
#' @param include_chr DEFAULT: c(1:19, "X"). Chromosomes to include in plot.
#' @param expland_lim_factor DEFAULT: 0.025. Scale to increase plot limits by.
#' @param label_thresh DEFAULT: NULL. Label mediators that surpass label_thresh. Default does not add labels.
#' @param label_thresh_greater_than DEFAULT: TRUE. If TRUE, passing mediators have log odds greater than the threshold.
#' If FALSE, passing mediators have log odds less than the threshold.  
#' @param label_only_chr DEFAULT: NULL. Only label mediators that pass label_thresh on the specified chromosome.
#' @param qtl_dat DEFAULT: NULL. QTL data that includes position of QTL and outcome. Adds ticks to the figure.
#' @export
#' @examples plot_posterior_odds()
plot_posterior_odds <- function(bmediatR_object, 
                                model_type = c("mediation", "partial", "complete", "colocal"),
                                med_annot, 
                                med_var = "protein.id",
                                include_chr = c(1:19, "X"), 
                                expand_lim_factor = 0.025, 
                                label_thresh = NULL, 
                                label_thresh_greater_than = TRUE,
                                label_only_chr = NULL,
                                bgcol = "white", altcol = "gray", altbgcol = "white", 
                                hlines_col = "gray80", col = "black", cex = 0.75,
                                qtl_dat = NULL,
                                outcome_symbol = NULL,
                                ymax = NULL,
                                ymin = NULL,
                                ...) {
  
  model_type <- model_type[1]
  
  post_odds <- matrix(bmediatR_object[["ln_post_odds"]][,model_type], ncol = 1)
  rownames(post_odds) <- rownames(bmediatR_object[["ln_post_odds"]])
  class(post_odds) <- "scan1"
  
  med_map_df <- med_annot %>%
    dplyr::select(tidyselect::all_of(med_var), symbol, chr, middle) %>%
    dplyr::filter(chr %in% include_chr) %>%
    dplyr::mutate(chr = factor(chr, levels = c(1:19, "X"))) %>%
    as.data.frame %>% 
    dplyr::arrange(chr)
  if (!is.null(qtl_dat)) {
    ## Add QTL to map for plotting
    med_map_df <- dplyr::bind_rows(med_map_df,
                                   qtl_dat %>%
                                     dplyr::mutate((!!as.symbol(med_var)) := "QTL",
                                                   symbol = "QTL") %>%
                                     dplyr::rename(middle = pos) %>%
                                     dplyr::select(tidyselect::all_of(med_var), symbol, chr, middle))
  }
  med_map <- map_df_to_list(map = med_map_df, marker_column = med_var, pos_column = "middle")
  
  gap <- sum(qtl2::chr_lengths(med_map))/100
  
  lim_shift <- (max(post_odds[,1]) - min(post_odds[,1])) * expand_lim_factor
  
  if (is.null(ymax)) { ymax <- max(post_odds[,1]) + lim_shift }
  if (is.null(ymin)) { ymin <- min(post_odds[,1]) - lim_shift }
  
  qtl2:::plot.scan1(post_odds, map = med_map, ylab = "Log posterior odds", type = "p", pch = 20, 
                    ylim = c(ymin, ymax),
                    bgcol = bgcol, altcol = altcol, altbgcol = altbgcol, hlines_col = hlines_col, col = col, 
                    cex = cex, gap = gap,
                    ...)
  
  xpos <- qtl2:::map_to_xpos(map = med_map, gap = gap)
  
  ## Mediator labels
  label_dat <- matrix(bmediatR_object[["ln_post_odds"]][,model_type], ncol = 1)
  colnames(label_dat) <- "post_odds"
  rownames(label_dat) <- rownames(bmediatR_object[["ln_post_odds"]])
  label_dat <- label_dat %>%
    as.data.frame %>%
    tibble::rownames_to_column(med_var) %>%
    dplyr::left_join(med_map_df)
  if (!is.null(label_only_chr)) {
    label_dat <- label_dat %>%
      dplyr::filter(chr == label_only_chr)
  } else {
    label_dat <- label_dat %>%
      dplyr::filter(chr %in% include_chr)
  }
  label_post_odds <- label_dat %>%
    dplyr::select(tidyselect::all_of(med_var), post_odds) %>%
    tibble::column_to_rownames(med_var) %>%
    as.matrix()
  
  
  if (!is.null(label_thresh)) {
    if (label_thresh_greater_than & any(label_post_odds > label_thresh)) {
      labels <- rownames(label_post_odds)[label_post_odds > label_thresh]
    }
    if (!label_thresh_greater_than & any(label_post_odds < label_thresh)) {
      labels <- rownames(label_post_odds)[label_post_odds < label_thresh]
    }
    
    if (!is.null(labels)) {
      label_map_df <- med_map_df %>%
        filter((!!as.symbol(med_var)) %in% labels) 
      
      for (i in 1:nrow(label_map_df)) {
        lab_pos <- xpos[label_map_df[i, med_var]]
        lab_post_odds <- post_odds[label_map_df[i, med_var],]
        
        text(x = lab_pos, y = lab_post_odds, label_map_df$symbol[i], font = 3)
      }
    }
  }
  if (!is.null(outcome_symbol)) {
    rug(x = xpos[med_annot %>% 
                   dplyr::filter(symbol == outcome_symbol) %>% 
                   pull(tidyselect::all_of(med_var))],
        lwd = 3,
        col = "black")
  }
  if (!is.null(qtl_dat)) {
    rug(x = xpos["QTL"],
        lwd = 3,
        col = "red")
  }
}




#' Barplot of posterior model probabilities
#'
#' This function takes posterior probability results from mediation_bf() and plots barplots of posterior model probabilities.
#'
#' @param bmediatR_object Output from bmediatR(). 
#' @param med_annot Annotation data for -omic mediators.
#' @param mediator_id Which mediator to plot posterior model probabilities for.
#' @param med_var DEFAULT: "protein.id". The column in med_annot to be used as a mediator id.
#' @param stack DEFAULT: FALSE. If TRUE, a stacked barplot is produced. If FALSE, a staggered
#' barplot is produced.
#' @param add_number_labels DEFAULT: FALSE. Add posterior probabilities as text above bars.
#' @param label_size DEFAULT: 5. Text size of probability labels, when included.
#' @param num_dig DEFAULT: 3. The number of digits after the decimal point if probabilities
#' are being included as labels.
#' @export
#' @examples plot_posterior_bar()
plot_posterior_bar_gxt <- function(bmediatR_object,
                                   med_annot = NULL,
                                   mediator_id,
                                   med_var = "protein.id",
                                   stack = FALSE,
                                   bar_col = c("seagreen4", "seagreen1", "skyblue", "gray70"),
                                   grp_labels = c(grp1 = "grp1", grp2 = "grp2"),
                                   grp_name = "grp", 
                                   relabel_x = NULL,
                                   add_number_labels = FALSE,
                                   label_size = 5,
                                   num_dig = 3,
                                   main = NULL) {
  
  odds_index <- return_preset_odds_index_gxt()
  
  post_mat <- exp(bmediatR_object$ln_post_c)
  post_mat <- post_mat / rowSums(post_mat)
  
  posterior_wide <- tibble::tibble(!!med_var := rownames(post_mat)) %>%
    dplyr::mutate(
      #complete_any       = rowSums(post_mat[, odds_index$complete_any,       drop = FALSE]),
      complete_grp1      = rowSums(post_mat[, odds_index$complete_grp1,      drop = FALSE]),
      complete_grp2      = rowSums(post_mat[, odds_index$complete_grp2,      drop = FALSE]),
      
      #partial_any        = rowSums(post_mat[, odds_index$partial_any,        drop = FALSE]),
      partial_grp1       = rowSums(post_mat[, odds_index$partial_grp1,       drop = FALSE]),
      partial_grp2       = rowSums(post_mat[, odds_index$partial_grp2,       drop = FALSE]),
      
      #colocal_any        = rowSums(post_mat[, odds_index$colocal_any,        drop = FALSE]),
      colocal_grp1       = rowSums(post_mat[, odds_index$colocal_grp1,       drop = FALSE]),
      colocal_grp2       = rowSums(post_mat[, odds_index$colocal_grp2,       drop = FALSE]),
      
      #nonmed_any         = rowSums(post_mat[, odds_index$other_nonmed_any,   drop = FALSE]),
      nonmed_grp1        = rowSums(post_mat[, odds_index$other_nonmed_grp1,  drop = FALSE]),
      nonmed_grp2        = rowSums(post_mat[, odds_index$other_nonmed_grp2,  drop = FALSE])
    )
  
  posterior_long <- posterior_wide %>%
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(med_var),
      names_to = c("class", "grp"),
      names_sep = "_",
      values_to = "post_p"
    ) %>%
    dplyr::mutate(
      class = dplyr::recode(class,
                            complete = "complete",
                            partial  = "partial",
                            colocal  = "colocal",
                            nonmed   = "other non-\nmediation"),
      class = factor(class, levels = c("complete", "partial", "colocal", "other non-\nmediation")),
      grp   = factor(grp, levels = c("grp1", "grp2"), labels = grp_labels)
    ) 
  
  # colors 
  base_cols <- c(
    "complete" = bar_col[1],
    "partial"  = bar_col[2],
    "colocal"  = bar_col[3],
    "other non-\nmediation" = bar_col[4]
  )
  
  plot_dat <- posterior_long %>%
    dplyr::filter((!!as.symbol(med_var)) == mediator_id) %>%
    dplyr::mutate(class = factor(class, levels = c("complete", "partial", "colocal", "other non-\nmediation")),
                  grp   = factor(grp, levels = grp_labels))
  
  bar_theme <- theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text.x = element_text(size = 11),
          legend.title = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 10))
  
  p <- ggplot(plot_dat, aes(x = class, y = post_p, fill = class, group = grp, alpha = grp)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_manual(name = "Model", values = base_cols) + 
    scale_alpha_manual(name = grp_name, values = setNames(c(0.7, 1.0), grp_labels)) +
    geom_hline(yintercept = c(0, 1), col = "gray", linetype = "dashed") +
    ylab("Posterior model probability") +
    xlab(NULL) + 
    bar_theme +
    guides(fill = guide_legend(order = 1),
           alpha = guide_legend(order = 2, override.aes = list(fill = base_cols["complete"])))

  p
}


plot_posterior_bar_gxt3 <- function(bmediatR_object,
                                    med_annot = NULL,
                                    mediator_id,
                                    med_var = "protein.id",
                                    stack = FALSE,
                                    bar_col = c("seagreen4", "seagreen1", "gray70"),
                                    grp_labels = c(grp1 = "grp1", grp2 = "grp2", grp3 = "grp3"),
                                    grp_name = "grp", 
                                    relabel_x = NULL,
                                    add_number_labels = FALSE,
                                    label_size = 5,
                                    num_dig = 3,
                                    main = NULL) {
  
  # infer group keys from names(grp_labels): c("grp1","grp2") or c("grp1","grp2","grp3")
  grp_keys <- names(grp_labels)
  if (is.null(grp_keys) || any(grp_keys == "")) {
    stop("grp_labels must be a named vector, e.g. c(grp1 = 'Group 1', grp2 = 'Group 2', grp3 = 'Group 3')")
  }
  n_grp <- length(grp_keys)
  
  post_mat <- exp(bmediatR_object$ln_post_c)
  post_mat <- post_mat / rowSums(post_mat)
  
  # helper: safe rowSums over a set of model indices
  sum_post <- function(idx) {
    if (length(idx) == 0) {
      return(rep(0, nrow(post_mat)))
    }
    if (max(idx) > ncol(post_mat)) {
      stop("Index exceeds number of models in post_mat (ncol(post_mat)). Check that indices match the 512-model ordering.")
    }
    rowSums(post_mat[, idx, drop = FALSE])
  }
  
  posterior_wide <- tibble::tibble(!!med_var := rownames(post_mat))
  
  for (gk in grp_keys) {
    g_num <- sub("^grp", "", gk)
    
    posterior_wide[[paste0("complete_", gk)]] <- sum_post(return_preset_odds_index_gxt3(paste0("complete_grp", g_num)))
    posterior_wide[[paste0("partial_", gk)]] <- sum_post(return_preset_odds_index_gxt3(paste0("partial_grp",  g_num)))
    # posterior_wide[[paste0("colocal_", gk)]] <- sum_post(return_preset_odds_index_gxt3(paste0("colocal_grp",  g_num)))
    posterior_wide[[paste0("nonmed_", gk)]] <- sum_post(return_preset_odds_index_gxt3(paste0("colocal_grp", g_num))) +
      sum_post(return_preset_odds_index_gxt3(paste0("other_nonmed_grp", g_num))) # remove colocal_grp from this if want to plot separately 
  }
  
  posterior_long <- posterior_wide %>%
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(med_var),
      names_to = c("class", "grp"),
      names_sep = "_",
      values_to = "post_p"
    ) %>%
    dplyr::mutate(
      class = dplyr::recode(class,
                            complete = "complete",
                            partial  = "partial",
                            nonmed   = "other non-\nmediation"), # add colocal to this if desired 
      class = factor(class, levels = c("complete", "partial", "other non-\nmediation")), # and this
      grp   = factor(grp, levels = grp_keys, labels = grp_labels[grp_keys])
    )
  
  
  # colors 
  base_cols <- c(
    "complete" = bar_col[1],
    "partial"  = bar_col[2],
    "other non-\nmediation" = bar_col[3]
  ) # "colocal"  = bar_col[3],
  
  plot_dat <- posterior_long %>%
    dplyr::filter((!!as.symbol(med_var)) == mediator_id) %>%
    dplyr::mutate(class = factor(class, levels = c("complete", "partial", "other non-\nmediation")),
                  grp   = factor(grp, levels = grp_labels)) # add colocal to class factor
  
  alpha_vals <- c(0.4, 0.7, 1.0)
  names(alpha_vals) <- grp_labels[grp_keys]
  
  bar_theme <- theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 11),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  p <- ggplot(plot_dat, aes(x = class, y = post_p, fill = class, group = grp, alpha = grp)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_manual(name = "Model", values = base_cols) +
    scale_alpha_manual(name = grp_name, values = alpha_vals) +
    geom_hline(yintercept = c(0, 1), col = "gray", linetype = "dashed") +
    ylab("Posterior model probability") +
    xlab(NULL) +
    bar_theme +
    guides(fill = guide_legend(order = 1),
           alpha = guide_legend(order = 2, override.aes = list(fill = base_cols["complete"])))
  
  p
}





# GxT grid
make_gxt_64_models <- function() {
  M_models <- expand.grid(a1 = 0:1, a2 = 0:1) # a1 fastest
  y_models <- expand.grid(b1 = 0:1, b2 = 0:1, c1 = 0:1, c2 = 0:1) # b1,b2,c1,c2 fastest-to-slowest
  out <- vector("list", length = nrow(M_models) * nrow(y_models))
  k <- 0
  for (m_idx in 1:nrow(M_models)) {
    for (y_idx in 1:nrow(y_models)) {
      k <- k + 1
      out[[k]] <- cbind(M_models[m_idx, , drop = FALSE], y_models[y_idx, , drop = FALSE]) |>
        mutate(model_id = k)
    }
  }
  bind_rows(out)
}

# axis levels for 8 causal models per group 
make_abc_levels <- function() {
  lev <- expand.grid(a = 0:1, b = 0:1, c = 0:1) |>
    mutate(code = paste0(a, b, c), label = paste0("{", a, ",", b, ",", c, "}"))
  list(codes = lev$code, labels = lev$label)
}

# make triangles for heatmap 
make_triangles <- function(df_xy, group_col = c("grp1","grp2")) {
  group_col <- match.arg(group_col)
  if (group_col == "grp1") {
    df_xy |>
      transmute(x, y, model_id,
                grp = "grp1",
                x0 = x - 0.5, y0 = y - 0.5,
                x1 = x + 0.5, y1 = y - 0.5,
                x2 = x - 0.5, y2 = y + 0.5)
  } else {
    df_xy |>
      transmute(x, y, model_id,
                grp = "grp2",
                x0 = x + 0.5, y0 = y + 0.5,
                x1 = x + 0.5, y1 = y - 0.5,
                x2 = x - 0.5, y2 = y + 0.5)
  }
}

# make grid 
make_gxt_med_plot <- function(grp_labels = c("Group 1","Group 2"),
                              complete_col = "seagreen4", partial_col = "seagreen1",
                              alpha_grp1 = 0.5, alpha_grp2 = 1,
                              border_col = "grey70", border_lwd = 0.3) {
  odds <- return_preset_odds_index_gxt()
  lev <- make_abc_levels()
  
  models <- make_gxt_64_models() |>
    mutate(g1_code = paste0(a1,b1,c1), g2_code = paste0(a2,b2,c2),
           g1_type = case_when(model_id %in% odds$complete_grp1 ~ "Complete",
                               model_id %in% odds$partial_grp1  ~ "Partial",
                               TRUE ~ "None"),
           g2_type = case_when(model_id %in% odds$complete_grp2 ~ "Complete",
                               model_id %in% odds$partial_grp2  ~ "Partial",
                               TRUE ~ "None"),
           g1_code = factor(g1_code, levels = lev$codes,labels = lev$labels),
           g2_code = factor(g2_code, levels = lev$codes,labels = lev$labels),
           x = as.integer(g1_code), y = as.integer(g2_code),
           g1_med = g1_type != "None", g2_med = g2_type != "None")
  
  # Single-fill tiles when mediation in only one group
  singles <- models |>
    mutate(tile_type = case_when(g1_med & !g2_med ~ g1_type,
                                !g1_med & g2_med ~ g2_type,
                                TRUE ~ NA_character_),
           tile_grp = case_when(g1_med & !g2_med ~ grp_labels[1],
                               !g1_med & g2_med ~ grp_labels[2],
                               TRUE ~ NA_character_)) |>
    filter(!is.na(tile_type)) |>
    mutate(tile_type = factor(tile_type, levels = c("Partial","Complete")),
           tile_grp = factor(tile_grp, levels = grp_labels))
  
  # Triangle tiles only when both groups have mediation
  both_xy <- models |>
    filter(g1_med & g2_med) |>
    select(model_id, x, y, g1_type, g2_type)
  
  tri1 <- make_triangles(both_xy, "grp1") |>
    left_join(both_xy |> select(model_id, g1_type), by = "model_id") |>
    mutate(type = g1_type)
  tri2 <- make_triangles(both_xy, "grp2") |>
    left_join(both_xy |> select(model_id, g2_type), by = "model_id") |>
    mutate(type = g2_type)
  
  tris <- bind_rows(tri1, tri2) |>
    mutate(grp = factor(grp, levels = c("grp1", "grp2"), labels = grp_labels),
           type = factor(type, levels = c("Partial", "Complete")))
  
  poly <- tris |>
    pivot_longer(cols = c(x0, y0, x1, y1, x2, y2),
                 names_to = c("coord", "pt"), names_pattern = "([xy])(\\d)",
                 values_to = "value") |>
    pivot_wider(names_from = coord, values_from = value, names_prefix = "p") |>
    mutate(order = as.integer(pt)) |>
    arrange(model_id, grp, order)
  
  alpha_legend_df <- tibble(grp = factor(grp_labels, levels = grp_labels),
                            alpha = c(alpha_grp1, alpha_grp2),
                            x = 1, y = 1)
  
  ggplot() +
    geom_tile(data = singles, aes(x = x, y = y, fill = tile_type, alpha = tile_grp), color = NA) +
    geom_polygon(data = poly, aes(x = px, y = py, group = interaction(model_id, grp),
                                  fill = type, alpha = grp), color = NA) +
    geom_tile(data = models, aes(x = x, y = y), fill = NA, color = border_col, linewidth = border_lwd) +
    scale_x_continuous(breaks = 1:8, labels = levels(models$g1_code), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:8, labels = levels(models$g2_code), expand = c(0, 0)) +
    scale_fill_manual(name = "Mediation", values = c("Partial" = partial_col, "Complete" = complete_col)) +
    scale_alpha_manual(name = NULL, values = setNames(c(alpha_grp1, alpha_grp2), grp_labels)) +
    guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2, override.aes = list(fill = "seagreen4"))) + 
    coord_equal() +
    labs(x = "Group 1", y = "Group 2") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
}



make_gxt_med_plot_ext <- function(grp_labels = c("Group 1","Group 2"),
                                            complete_col = "seagreen4", partial_col = "seagreen1",
                                            coloc_col = "skyblue", nonmed_col = "gray70",
                                            alpha_grp1 = 0.5, alpha_grp2 = 1,
                                            border_col = "grey70", border_lwd = 0.3) {
  odds <- return_preset_odds_index_gxt()
  lev <- make_abc_levels()
  
  models <- make_gxt_64_models() |>
    mutate(g1_code = paste0(a1,b1,c1), g2_code = paste0(a2,b2,c2),
           g1_type = case_when(model_id %in% odds$complete_grp1 ~ "Complete",
                               model_id %in% odds$partial_grp1 ~ "Partial",
                               model_id %in% odds$colocal_grp1 ~ "Colocalization",
                               model_id %in% odds$other_nonmed_grp1 ~ "Non-mediation",
                               TRUE ~ "None"),
           g2_type = case_when(model_id %in% odds$complete_grp2 ~ "Complete",
                               model_id %in% odds$partial_grp2 ~ "Partial",
                               model_id %in% odds$colocal_grp2 ~ "Colocalization",
                               model_id %in% odds$other_nonmed_grp2 ~ "Non-mediation",
                               TRUE ~ "None"),
           g1_code = factor(g1_code, levels = lev$codes, labels = lev$labels),
           g2_code = factor(g2_code, levels = lev$codes, labels = lev$labels),
           x = as.integer(g1_code), y = as.integer(g2_code))
  
  # Triangles for ALL tiles
  tri1 <- make_triangles(models |> select(model_id, x, y), "grp1") |>
    left_join(models |> select(model_id, g1_type), by = "model_id") |>
    mutate(type = g1_type)
  tri2 <- make_triangles(models |> select(model_id, x, y), "grp2") |>
    left_join(models |> select(model_id, g2_type), by = "model_id") |>
    mutate(type = g2_type)
  
  tris <- bind_rows(tri1, tri2) |>
    mutate(grp = factor(grp, levels = c("grp1","grp2"), labels = grp_labels),
           type = factor(type, levels = c("Non-mediation","Colocalization","Partial","Complete","None")))
  
  poly <- tris |>
    pivot_longer(cols = c(x0, y0, x1, y1, x2, y2),
                 names_to = c("coord","pt"), names_pattern = "([xy])(\\d)",
                 values_to = "value") |>
    pivot_wider(names_from = coord, values_from = value, names_prefix = "p") |>
    mutate(order = as.integer(pt)) |>
    arrange(model_id, grp, order)
  
  ggplot() +
    geom_polygon(data = poly, aes(x = px, y = py, group = interaction(model_id, grp),
                                  fill = type), color = NA) +
    geom_tile(data = models, aes(x = x, y = y), fill = NA, color = border_col, linewidth = border_lwd) +
    scale_x_continuous(breaks = 1:8, labels = levels(models$g1_code), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:8, labels = levels(models$g2_code), expand = c(0, 0)) +
    scale_fill_manual(name = "Evidence",
                      values = c("Non-mediation" = nonmed_col,
                                 "Colocalization" = coloc_col,
                                 "Partial" = partial_col,
                                 "Complete" = complete_col,
                                 "None" = "white")) +
    coord_equal() +
    labs(x = "Group 1", y = "Group 2") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.box = "vertical")
}


# old grid 
# plot_gxt_mediation_grid <- function(idx_grp1, idx_grp2, 
#                                     col_grp1 = "#1f77b4",   # default: blue
#                                     col_grp2 = "#ff7f0e",   # default: orange
#                                     col_both = "#9467bd",   # default: purple (visual blend),
#                                     col_none = "grey75") {
#   models <- make_gxt_64_models()
#   
#   # Build grp1/grp2 3-bit causal model IDs:
#   # grp1 uses (a1,b1,c1), grp2 uses (a2,b2,c2)
#   models <- models |>
#     mutate(g1_code = paste0(a1, b1, c1),
#            g2_code = paste0(a2, b2, c2),
#            med_grp1 = model_id %in% idx_grp1,
#            med_grp2 = model_id %in% idx_grp2,
#            state = case_when(
#              med_grp1 & med_grp2 ~ "both",
#              med_grp1 ~ "grp1",
#              med_grp2 ~ "grp2",
#              TRUE ~ "none"
#            ))
#   
#   # Axis ordering: consistent 8-level order (a varies fastest, then b, then c)
#   g1_levels <- make_group8_levels(prefix = "grp1")
#   g2_levels <- make_group8_levels(prefix = "grp2")
#   
#   models <- models |>
#     mutate(
#       g1_code = factor(g1_code, levels = g1_levels$level_codes, labels = g1_levels$level_labels),
#       g2_code = factor(g2_code, levels = g2_levels$level_codes, labels = g2_levels$level_labels),
#       state   = factor(state, levels = c("none", "grp1", "grp2", "both"))
#     )
#   
#   ggplot(models, aes(x = g1_code, y = g2_code)) +
#     geom_tile(aes(fill = state), color = "white", linewidth = 0.6) +
#     scale_fill_manual(
#       values = c(
#         "none" = col_none,
#         "grp1" = col_grp1,
#         "grp2" = col_grp2,
#         "both" = col_both
#       ),
#       name = "Mediation present"
#     ) +
#     coord_equal() +
#     labs(
#       x = "grp1 causal model (a=X→M, b=M→Y, c=X→Y)",
#       y = "grp2 causal model (a=X→M, b=M→Y, c=X→Y)"
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid = element_blank()
#     )
# }
