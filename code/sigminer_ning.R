count_components <- function(min, max, label, feature) {
  # Count all samples
  if (label == "point") {
    feature %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(.data$value == min, na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else if (label == "range") {
    feature %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(.data$value > min & .data$value <= max, na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else {
    stop("Bad labels for feature setting, can only be 'point' and 'range'!")
  }
}

count_components_wrapper <- function(feature_df, f_name, feature_setting) {
  samps <- unique(feature_df$ID)
  feature_df <- feature_df %>% dplyr::as_tibble()
  # Make sure sample order is consistent
  feature_df$ID <- factor(feature_df$ID, levels = samps)
  
  specific_f <- feature_setting[feature_setting$feature == f_name]
  
  # https://www.jianshu.com/p/24bbf44e4fa2
  specific_f %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$component) %>%
    dplyr::summarize(
      count = list(count_components(.data$min, .data$max, .data$label, feature_df)),
      sample = list(samps)
    ) %>%
    tidyr::unnest(cols = c("count", "sample")) %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "count", values_fill = list(count = 0L)) %>%
    # Should not have NA value, but take case of it with 0
    data.table::as.data.table()
  
  # Same result as above
  # But I don't know which is more efficient
  #
  # specific_f %>%
  #   dplyr::as_tibble() %>%
  #   dplyr::group_by(.data$component) %>%
  #   dplyr::summarize(
  #     count = paste(count_components(.data$min, .data$max, .data$label, feature_df), collapse = ","),
  #   ) %>%
  #   tidyr::separate(.data$count, samps, sep = ",", convert = TRUE) %>%
  #   data.table::as.data.table()
}


set_future_strategy <- function() {
  if (packageVersion("future") >= "1.20.0") {
    "multisession"
  } else {
    "multiprocess"
  }
}

getBPnum <- function(abs_profiles, chrlen) {
  res <- purrr::map_df(abs_profiles, function(x, chrlen) {
    calcBPnum <- function(df, c, chrlen) {
      intervals <-
        seq(1, chrlen[chrlen[, 1] == c, 2] + 10000000, 10000000)
      y <- tryCatch(
        hist(df$end[-nrow(df)],
             breaks = intervals,
             plot = FALSE
        )$counts,
        error = function(e) {
          stop(
            "Stop due to the following reason. Please check if your genome build is right.",
            "\n", e$message
          )
        }
      )
      y
    }
    
    x <- x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      tidyr::nest() %>%
      dplyr::mutate(value = purrr::map2(
        data, .data$chromosome, calcBPnum,
        chrlen = chrlen
      )) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("value"))
    data.table::data.table(
      value = purrr::reduce(x$value, c)
    )
  }, chrlen = chrlen, .id = "ID")
  
  return(res)
}

get_features_wang <- function(CN_data,
                              cores = 1,
                              genome_build = c("hg19", "hg38", "mm10", "mm9"),
                              feature_setting = sigminer::CN.features[1:8,]) {
  genome_build <- match.arg(genome_build)
  # get chromosome lengths and centromere locations
  chrlen <- get_genome_annotation(data_type = "chr_size", genome_build = genome_build)
  centromeres <- get_genome_annotation(data_type = "centro_loc", genome_build = genome_build)
  
  oplan <- future::plan()
  future::plan(set_future_strategy(), workers = cores)
  on.exit(future::plan(oplan), add = TRUE)
  
  features <- unique(feature_setting$feature)
  
  .get_feature <- function(i) {
    if (i == "SS") {
      send_info("Getting (log10 based) segment size...")
      zz <- getSegsize(CN_data)
      zz$value <- log10(zz$value)
      zz
    } else if (i == "BP10MB") {
      send_info("Getting breakpoint count per 10 Mb...")
      getBPnum(CN_data, chrlen)
    } else if (i == "OsCN") {
      send_info("Getting length of chains of oscillating copy number...")
      getOscilation(CN_data)
    } else if (i == "BPArm") {
      send_info("Getting breakpoint count per chromosome arm...")
      getCentromereDistCounts(CN_data, centromeres)
    } else if (i == "CNCP") {
      send_info("Getting change-point copy number change...")
      getChangepointCN(CN_data)
    } else if (i == "CN") {
      send_info("Getting copy number...")
      getCN(CN_data)
    } else if (i == "BoChr") {
      send_info("Getting burden of chromosome...")
      getBoChr(CN_data, genome_build)
    } else if (i == "NChrV") {
      send_info("Getting number of autosome with CNV...")
      getNChrV(CN_data, genome_build)
    } else if (i == "NC50") {
      send_info("Getting the minimal number of chromosome with 50% CNV...")
      getNC50(CN_data, genome_build)
    }
  }
  
  res <- furrr::future_map(features, .get_feature,
                           .progress = TRUE, .options = furrr::furrr_options(seed = TRUE)
  )
  res <- res %>% setNames(features)
  res
}


sig_tally.CopyNumber <- function(object,
                                 method = "Wang",
                                 ignore_chrs = NULL,
                                 indices = NULL,
                                 add_loh = FALSE,
                                 feature_setting = sigminer::CN.features[1:7,],
                                 cores = 1,
                                 keep_only_matrix = FALSE,
                                 ...) {
  method <- match.arg(method, choices = c( "W"))
  
  cn_list <- sigminer:::get_cnlist(object, ignore_chrs = ignore_chrs)
  
  if (startsWith(method, "W")) {
    # Method: Wang Shixiang, 'W'
    
    sigminer:::send_info("Step: getting copy number features.")
    cn_features <- sigminer:::get_features_wang(
      CN_data = cn_list, cores = cores,
      genome_build = object@genome_build,
      feature_setting = feature_setting
    )
    sigminer:::send_success("Gotten.")
    # Make order as unique(feature_setting)$feature
    # cn_features <- cn_features[unique(feature_setting$feature)]

    sigminer:::send_info("Step: generating copy number components.")
    # Check feature setting
    if (!inherits(feature_setting, "sigminer.features")) {
      feature_setting <- get_feature_components(feature_setting)
    }
    sigminer:::send_success("{.code feature_setting} checked.")
    
    sigminer:::send_info("Step: counting components.")
    cn_components <- purrr::map2(cn_features, names(cn_features),
                                 count_components_wrapper,
                                 feature_setting = feature_setting
    )
    sigminer:::send_success("Counted.")
    
    ## Remove BoChr value is 0 in features
    if ("BoChr" %in% names(cn_features)) {
      cn_features$BoChr <- cn_features$BoChr[cn_features$BoChr$value != 0]
    }
    
    sigminer:::send_info("Step: generating components by sample matrix.")
    cn_matrix <- data.table::rbindlist(cn_components, fill = TRUE, use.names = TRUE) %>%
      dplyr::as_tibble() %>%
      tibble::column_to_rownames(var = "component") %>%
      as.matrix()
    # Order the matrix as feature_setting
    cn_matrix <- cn_matrix[feature_setting$component, ] %>%
      t()
    
    if (any(is.na(cn_matrix))) {
      sigminer:::send_warning("{.code NA} detected. There may be an issue, please contact the developer!")
      sigminer:::send_warning("Data will still returned, but please take case of it.")
    }
    # cn_matrix[is.na(cn_matrix)] <- 0L
    feature_setting$n_obs <- colSums(cn_matrix, na.rm = TRUE)
  }  
  
  sigminer:::send_success("Matrix generated.")
}
a<-sig_tally.CopyNumber(cn,method = "W")
