#' Annotate OncoKB levels of evidence
#'
#' Adds OncoKB oncogenicity and actionability annotation to VEP-annotated MAF. See URLs below.
#'
#' @param api_token OncoKB API token. This can be retrieved from https://www.oncokb.org/account/settings after you have set up your account.
#' @param maf Input MAF.
#' @param cancer_types Data frame with samples mapped to cancer type for accurate levels of actionability.
#'
#' @param gene Gene.
#' @param protein_change Abbreviated form, e.g. "H1047R".
#' @param variant_type Lower case variant type, e.g. "missense".
#' @param start Genomic start position of variant.
#' @param end Genomic end position of variant.
#' @param cancer_type Oncotree code for cancer type. Can be left blank.
#' @param parallelize Boolean indicating whether to parallelize annotation, using \code{future} backend.
#'
#' @return Annotated MAF with columns indicating functionality of mutation and levels of actionability.
#'
#' @source \url{oncokb.org}
#' @source \url{github.com/oncokb/oncokb-annotator}
#'
#' @import purrr
#' @importFrom dplyr case_when bind_cols left_join select
#' @importFrom future plan
#' @importFrom furrr future_pmap_dfr
#' @importFrom plyr revalue
#' @importFrom httr modify_url GET content
#' @importFrom stringr str_replace str_extract
#' @importFrom tibble tibble add_column
#'
#' @examples
#' query_oncokb('PIK3CA', 'H1047R', 'missense')
#'
#' @name oncokb_annotate_maf
NULL

consequence_map = c('3\'Flank' = 'any',
                    '5\'Flank ' = 'any',
                    # 'Targeted_Region'= 'inframe_deletion', 'inframe_insertion',
                    'Frame_Shift_Del' = 'frameshift_variant',
                    'Frame_Shift_Ins' = 'frameshift_variant',
                    'In_Frame_Del' = 'inframe_deletion',
                    'In_Frame_Ins' = 'inframe_insertion',
                    'Missense_Mutation' = 'missense_variant',
                    'Nonsense_Mutation' = 'stop_gained',
                    'Nonstop_Mutation' = 'stop_lost',
                    'Splice_Site' = 'splice_region_variant',
                    'Splice_Region' = 'splice_region_variant',
                    'Translation_Start_Site' = 'start_lost')

coding_mutations = c('Frame_Shift_Del',
                     'Frame_Shift_Ins',
                     'In_Frame_Del',
                     'In_Frame_Ins',
                     'Missense_Mutation',
                     'Nonsense_Mutation',
                     'Nonstop_Mutation',
                     'Splice_Site',
                     'Targeted_Region',
                     'Translation_Start_Site')

# Allow parallellization
future::plan(future::multisession)

#' @export
#' @rdname oncokb_annotate_maf
#TODO: some results such as KIT D816Y, work correctly with mquery but not here.
query_oncokb = function(api_token, gene, protein_change, variant_type, cancer_type = 'CANCER') {

  if (variant_type != '') {

    base_url = 'https://www.oncokb.org/api/v1/'
    info_url = paste0(base_url, "info")
    query_base_url = paste0(base_url, "annotate/mutations/byProteinChange")
    if(!exists('oncokb_version')) {
      oncokb_version <<- httr::content(httr::GET(info_url))$dataVersion$version
    }
    tag = paste(gene, protein_change, cancer_type, sep = '-')

    evidence_types <- stringr::str_c(
      c("ONCOGENIC", "MUTATION_EFFECT",
        "STANDARD_THERAPEUTIC_IMPLICATIONS_FOR_DRUG_SENSITIVITY",
        "INVESTIGATIONAL_THERAPEUTIC_IMPLICATIONS_DRUG_SENSITIVITY",
        "INVESTIGATIONAL_THERAPEUTIC_IMPLICATIONS_DRUG_RESISTANCE"),
      collapse = ",")

    if (!exists('cached_entries')) cached_entries <<- vector(mode = 'list')

    if (!tag %in% names(cached_entries)) {
        query_url = httr::modify_url(query_base_url, query = list(
            hugoSymbol = gene,
            alteration = protein_change,
            consequence = variant_type,
            tumorType = cancer_type,
            evidenceType = evidence_types
        ))
        cat(query_url)

        headers <- httr::add_headers(Authorization = stringr::str_c("Bearer ", api_token),
                                     accept = "application/json")

        oncokb_response = httr::GET(query_url, headers)

        if(oncokb_response$headers$`content-type` == "application/problem+json") {
          json <- httr::content(oncokb_response, "text", "application/json")
          errmsg <- stringr::str_c("Problem querying oncokb:\n",
                                   jsonlite::prettify(json))
          rlang::abort(errmsg)
        }

        oncokb_response = httr::content(oncokb_response)

        cached_entries[[tag]] = oncokb_response
    } else {
        oncokb_response = cached_entries[[tag]]
    }

    drugs = purrr::map(oncokb_response$treatments, 'drugs') %>%
        purrr::map(., function(x) paste(unlist(x))) %>%
        purrr::simplify() %>%
        unique()

    tibble::tibble(oncogenic = as.character(oncokb_response$oncogenic),
                   oncokb_effect = ifelse(is.null(oncokb_response$mutationEffect$knownEffect), '',
                                          oncokb_response$mutationEffect$knownEffect),
                   oncokb_level = ifelse(is.null(oncokb_response$highestSensitiveLevel), '',
                                         oncokb_response$highestSensitiveLevel),
                   oncokb_resistance_level = ifelse(is.null(oncokb_response$highestResistanceLevel), '',
                                                    oncokb_response$highestResistanceLevel),
                   oncokb_drugs = ifelse(length(drugs) == 0, '',
                                         paste(unlist(unique(drugs)), collapse = ',')),
                   oncokb_version = oncokb_version)
  } else {
      tibble::tibble(oncogenic = '')
  }
}


#' @export
#' @rdname oncokb_annotate_maf
mquery_oncokb = function(api_token, gene, protein_change, variant_type, cancer_type = 'CANCER', chunk_size = 1000,
                         shorten_query = F, verbose = F) {

  if (any(variant_type != '')) {
    query_df <- data.frame(gene, protein_change, variant_type)
    if(shorten_query) {
      query_df <- dplyr::distinct(query_df)
    }
    query_inputs <- as.list(query_df)
    if(dplyr::n_distinct(lengths(query_inputs)) != 1) {
      rlang::abort("Must provide same number of query inputs to mquery_oncokb()")
    }

    is <- seq_along(gene)
    chunks <- split(is, ceiling(is/chunk_size))

    pb <- progress::progress_bar$new(format = "[:bar] :current/:total chunks eta: :eta", total = length(chunks), force = T)
    res <-
      purrr::map_dfr(chunks,
                     ~{
                       if(verbose) pb$tick()
                       mquery_oncokb_chunk(api_token, gene[.x], protein_change[.x],
                                           variant_type[.x], cancer_type,
                                           shortened_query = shorten_query)
                     })
    pb$terminate()
    res
  } else {
    tibble::tibble(oncogenic = '')
  }
}

mquery_oncokb_chunk = function(api_token, gene, protein_change, variant_type, cancer_type = 'CANCER',
                               shortened_query = F) {
  base_url = 'https://www.oncokb.org/api/v1/'
  info_url = paste0(base_url, "info")
  query_url = paste0(base_url, "annotate/mutations/byProteinChange")
  if(!exists('oncokb_version')) {
    oncokb_version <<- httr::content(httr::GET(info_url))$dataVersion$version
  }

  query_base <- list(tumor_type = cancer_type,
                     evidenceTypes = list("ONCOGENIC", "MUTATION_EFFECT",
                                          "STANDARD_THERAPEUTIC_IMPLICATIONS_FOR_DRUG_SENSITIVITY",
                                          "INVESTIGATIONAL_THERAPEUTIC_IMPLICATIONS_DRUG_SENSITIVITY",
                                          "INVESTIGATIONAL_THERAPEUTIC_IMPLICATIONS_DRUG_RESISTANCE"))
  query_inputs <- list(gene, protein_change, variant_type)
  body_structure <-  purrr::pmap(
    query_inputs,
    ~c(query_base,  list(gene = list(hugoSymbol = ..1),
                         alteration = ..2, consequence = ..3)))
  body <- jsonlite::toJSON(body_structure, auto_unbox = TRUE)

  headers <- c(
    httr::accept_json(),  httr::content_type_json(),
    httr::add_headers(Authorization = stringr::str_c("Bearer ", api_token)))

  oncokb_response = httr::POST(query_url, headers, body = body)

  if(oncokb_response$headers$`content-type` == "application/problem+json") {
    json <- httr::content(oncokb_response, "text", "application/json")
    errmsg <- stringr::str_c("Problem querying oncokb:\n",
                             jsonlite::prettify(json))
    rlang::abort(errmsg)
  }

  oncokb_response = httr::content(oncokb_response)


  purrr::map_dfr(seq_along(oncokb_response),
                 ~{
                   r <- oncokb_response[[.x]]
                   drugs = purrr::map(r$treatments, 'drugs') %>%
                     purrr::map(., function(x) paste(unlist(x))) %>%
                     purrr::simplify() %>%
                     unique()

                   res <-
                     tibble::tibble(
                       oncogenic = as.character(r$oncogenic),
                       oncokb_effect = ifelse(is.null(r$mutationEffect$knownEffect), '',
                                              r$mutationEffect$knownEffect),
                       oncokb_level = ifelse(is.null(r$highestSensitiveLevel), '',
                                             r$highestSensitiveLevel),
                       oncokb_resistance_level = ifelse(is.null(r$highestResistanceLevel), '',
                                                        r$highestResistanceLevel),
                       oncokb_drugs = ifelse(length(drugs) == 0, '',
                                             paste(unlist(unique(drugs)), collapse = ',')),
                       oncokb_version = oncokb_version)
                   if(shortened_query) {
                     res <- dplyr::bind_cols(
                       tibble::tibble(hugoSymbol = gene[.x],
                                      alteration = protein_change[.x],
                                      consequence = variant_type[.x]),
                       res
                     )
                   }
                   res
                 })
}


#' @export
#' @rdname oncokb_annotate_maf
oncokb_annotate_maf = function(api_token, maf, cancer_types = NULL, parallelize = TRUE)
{
    if (is.null(cancer_types) & !'cancer_type' %in% names(maf)) {
        message('No cancer types(s) specified, defaults to CANCER')
        maf$cancer_type = 'CANCER'
    } else if (is.character(cancer_types)) {
        maf = tibble::add_column(maf, cancer_type = cancer_types)
    } else {
        maf = dplyr::left_join(maf, cancer_types, by = 'Tumor_Sample_Barcode')
    }

    oncokb_cols = mutate(maf,
           api_token = api_token,
           gene = Hugo_Symbol,
           protein_change = stringr::str_replace(HGVSp_Short, 'p.', ''),
           variant_type = dplyr::case_when(
               (Variant_Classification %in% coding_mutations & HGVSp_Short != '') | # this is necessary to avoid poorly annotated but likely FP indel calls from Pindel
               (Variant_Classification %in% c('Splice_Site', 'Splice_Region') & HGVSc != '') |
                Hugo_Symbol == 'TERT' ~
                   plyr::revalue(Variant_Classification, consequence_map, warn_missing = F),
             TRUE ~ ''),
           start = stringr::str_extract(Protein_position, '^[0-9]+(?=\\/|\\-)'),
           end = stringr::str_extract(Protein_position, '(?<=-)[0-9]+(?=/)')
           ) %>%
        dplyr::select(gene, protein_change, variant_type, cancer_type)

    if (parallelize == TRUE) {
        oncokb_cols = furrr::future_pmap_dfr(oncokb_cols, query_oncokb)
    } else {
        oncokb_cols = purrr::pmap_dfr(oncokb_cols, query_oncokb)
    }

    dplyr::bind_cols(maf, oncokb_cols)

}
