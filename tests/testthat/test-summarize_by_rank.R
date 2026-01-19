test_tax <- tibble::tribble(
  ~seq_id, ~kingdom, ~phylum, ~class, ~order, ~family, ~genus, ~species,
  "seq1", "Fungi", "Ascomycota", "Dothideomycetes", "Pleosporales", "Pleosporaceae", "Alternaria", "Alternaria alternata",
  "seq2", "Fungi", "Ascomycota", "Dothideomycetes", "Pleosporales", "Pleosporaceae", "Alternaria", "Alternaria brassicae",
  "seq3", "Fungi", "Ascomycota", "Dothideomycetes", "Pleosporales", "Pleosporaceae", "Alternaria", "Alternaria alternata",
  "seq4", "Fungi", "Ascomycota", "Dothideomycetes", "Dothideales", "Dothideaceae", "Dothidea", "Dothidea umbrina",
  "seq5", "Fungi", "Basidiomycota", "Agaricomycetes", "Agaricales", "Agaricaceae", "Agaricus", "Agaricus bisporus",
  "seq6", "Fungi", "Basidiomycota", "Tremellomycetes", "Tremellales", "Tremellaceae", NA, NA,
  "seq7", "Fungi", "Basidiomycota", "Agaricomycetes", "Agaricales", "Amanitaceae", "Amanita", "Amanita muscaria",
  "seq8", "Fungi", "Basidiomycota", "Agaricomycetes", "Agaricales", "Amanitaceae", "Limacella", "Limacella guttata",
  "seq9", "Fungi", "Basidiomycota", "Agaricomycetes", "Agaricales", "Amanitaceae", "Amanita", "Amanita muscaria",
  "seq10", "Fungi", "Basidiomycota", "Agaricomycetes", "Agaricales", "Amanitaceae", "Amanita", "Amanita masasiensis"
)

to_int <- function(x) {
  as.integer(factor(x)) - 1L
}

test_summary <- with(test_tax,
  tibble::tribble(
    ~supertaxon, ~superrank, ~rank, ~n_taxa, ~n_seq, ~seq_id, ~true_partition,
    "Fungi", "kingdom", "phylum", 2L, 10L, seq_id, to_int(phylum),
    "Fungi", "kingdom", "class", 3L, 10L, seq_id, to_int(class),
    "Fungi", "kingdom", "order", 4L, 10L, seq_id, to_int(order),
    "Fungi", "kingdom", "family", 5L, 10L, seq_id, to_int(family),
    "Fungi", "kingdom", "genus", 5L, 9L, seq_id[!is.na(genus)], to_int(genus[!is.na(genus)]),
    "Fungi", "kingdom", "species", 7L, 9L, seq_id[!is.na(species)], to_int(species[!is.na(species)]),
    "Ascomycota", "phylum", "class", 1L, 4L, seq_id[phylum == "Ascomycota"], to_int(class[phylum == "Ascomycota"]),
    "Basidiomycota", "phylum", "class", 2L, 6L, seq_id[phylum == "Basidiomycota"], to_int(class[phylum == "Basidiomycota"]),
    "Ascomycota", "phylum", "order", 2L, 4L, seq_id[phylum == "Ascomycota"], to_int(order[phylum == "Ascomycota"]),
    "Basidiomycota", "phylum", "order", 2L, 6L, seq_id[phylum == "Basidiomycota"], to_int(order[phylum == "Basidiomycota"]),
    "Ascomycota", "phylum", "family", 2L, 4L, seq_id[phylum == "Ascomycota"], to_int(family[phylum == "Ascomycota"]),
    "Basidiomycota", "phylum", "family", 3L, 6L, seq_id[phylum == "Basidiomycota"], to_int(family[phylum == "Basidiomycota"]),
    "Ascomycota", "phylum", "genus", 2L, 4L, seq_id[phylum == "Ascomycota"], to_int(genus[phylum == "Ascomycota"]),
    "Basidiomycota", "phylum", "genus", 3L, 5L, seq_id[phylum == "Basidiomycota" & !is.na(genus)], to_int(genus[phylum == "Basidiomycota" & !is.na(genus)]),
    "Ascomycota", "phylum", "species", 3L, 4L, seq_id[phylum == "Ascomycota"], to_int(species[phylum == "Ascomycota"]),
    "Basidiomycota", "phylum", "species", 4L, 5L, seq_id[phylum == "Basidiomycota" & !is.na(species)], to_int(species[phylum == "Basidiomycota" & !is.na(species)]),
    "Agaricomycetes", "class", "order", 1L, 5L, seq_id[class == "Agaricomycetes"], to_int(order[class == "Agaricomycetes"]),
    "Dothideomycetes", "class", "order", 2L, 4L, seq_id[class == "Dothideomycetes"], to_int(order[class == "Dothideomycetes"]),
    "Tremellomycetes", "class", "order", 1L, 1L, seq_id[class == "Tremellomycetes"], to_int(order[class == "Tremellomycetes"]),
    "Agaricomycetes", "class", "family", 2L, 5L, seq_id[class == "Agaricomycetes"], to_int(family[class == "Agaricomycetes"]),
    "Dothideomycetes", "class", "family", 2L, 4L, seq_id[class == "Dothideomycetes"], to_int(family[class == "Dothideomycetes"]),
    "Tremellomycetes", "class", "family", 1L, 1L, seq_id[class == "Tremellomycetes"], to_int(family[class == "Tremellomycetes"]),
    "Agaricomycetes", "class", "genus", 3L, 5L, seq_id[class == "Agaricomycetes"], to_int(genus[class == "Agaricomycetes"]),
    "Dothideomycetes", "class", "genus", 2L, 4L, seq_id[class == "Dothideomycetes"], to_int(genus[class == "Dothideomycetes"]),
    "Agaricomycetes", "class", "species", 4L, 5L, seq_id[class == "Agaricomycetes"], to_int(species[class == "Agaricomycetes"]),
    "Dothideomycetes", "class", "species", 3L, 4L, seq_id[class == "Dothideomycetes"], to_int(species[class == "Dothideomycetes"]),
    "Agaricales", "order", "family", 2L, 5L, seq_id[order == "Agaricales"], to_int(family[order == "Agaricales"]),
    "Dothideales", "order", "family", 1L, 1L, seq_id[order == "Dothideales"], to_int(family[order == "Dothideales"]),
    "Pleosporales", "order", "family", 1L, 3L, seq_id[order == "Pleosporales"], to_int(family[order == "Pleosporales"]),
    "Tremellales", "order", "family", 1L, 1L, seq_id[order == "Tremellales"], to_int(family[order == "Tremellales"]),
    "Agaricales", "order", "genus", 3L, 5L, seq_id[order == "Agaricales"], to_int(genus[order == "Agaricales"]),
    "Dothideales", "order", "genus", 1L, 1L, seq_id[order == "Dothideales"], to_int(genus[order == "Dothideales"]),
    "Pleosporales", "order", "genus", 1L, 3L, seq_id[order == "Pleosporales"], to_int(genus[order == "Pleosporales"]),
    "Agaricales", "order", "species", 4L, 5L, seq_id[order == "Agaricales"], to_int(species[order == "Agaricales"]),
    "Dothideales", "order", "species", 1L, 1L, seq_id[order == "Dothideales"], to_int(species[order == "Dothideales"]),
    "Pleosporales", "order", "species", 2L, 3L, seq_id[order == "Pleosporales"], to_int(species[order == "Pleosporales"]),
    "Agaricaceae", "family", "genus", 1L, 1L, seq_id[family == "Agaricaceae"], to_int(genus[family == "Agaricaceae"]),
    "Amanitaceae", "family", "genus", 2L, 4L, seq_id[family == "Amanitaceae"], to_int(genus[family == "Amanitaceae"]),
    "Dothideaceae", "family", "genus", 1L, 1L, seq_id[family == "Dothideaceae"], to_int(genus[family == "Dothideaceae"]),
    "Pleosporaceae", "family", "genus", 1L, 3L, seq_id[family == "Pleosporaceae"], to_int(genus[family == "Pleosporaceae"]),
    "Agaricaceae", "family", "species", 1L, 1L, seq_id[family == "Agaricaceae"], to_int(species[family == "Agaricaceae"]),
    "Amanitaceae", "family", "species", 3L, 4L, seq_id[family == "Amanitaceae"], to_int(species[family == "Amanitaceae"]),
    "Dothideaceae", "family", "species", 1L, 1L, seq_id[family == "Dothideaceae"], to_int(species[family == "Dothideaceae"]),
    "Pleosporaceae", "family", "species", 2L, 3L, seq_id[family == "Pleosporaceae"], to_int(species[family == "Pleosporaceae"]),
    "Agaricus", "genus", "species", 1L, 1L, seq_id[genus == "Agaricus" & !is.na(genus)], to_int(species[genus == "Agaricus" & !is.na(genus)]),
    "Alternaria", "genus", "species", 2L, 3L, seq_id[genus == "Alternaria" & !is.na(genus)], to_int(species[genus == "Alternaria" & !is.na(genus)]),
    "Amanita", "genus", "species", 2L, 3L, seq_id[genus == "Amanita" & !is.na(genus)], to_int(species[genus == "Amanita" & !is.na(genus)]),
    "Dothidea", "genus", "species", 1L, 1L, seq_id[genus == "Dothidea" & !is.na(genus)], to_int(species[genus == "Dothidea" & !is.na(genus)]),
    "Limacella", "genus", "species", 1L, 1L, seq_id[genus == "Limacella" & !is.na(genus)], to_int(species[genus == "Limacella" & !is.na(genus)])
  )
)

test_that("summarize_by_rank() works", {
  expect_equal(
    optimotu:::summarize_by_rank(test_tax, c("kingdom", "phylum", "class", "order", "family", "genus", "species")),
    test_summary
  )
})

test_that("summarize_by_rank() accepts non-default ID column name", {
  expect_equal(
    optimotu:::summarize_by_rank(
      dplyr::rename(test_tax, id = seq_id),
      c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
      id_col = "id"
    ),
    test_summary
  )
})
