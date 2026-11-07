test_that("sample columns are detected dynamically", {
  expect_equal(
    plankannot:::plankannot_detect_sample_columns(c("ASV", "Sample_1", "Sample_2", "Sample_10", "Other")),
    c("Sample_1", "Sample_2", "Sample_10")
  )
})

test_that("only agreed traits are configured", {
  map <- plankannot:::plankannot_default_trait_map()
  expect_named(map$copepoda, c("Functional group", "Myelination", "Spawning mode", "Trophic group", "Feeding mode"))
  expect_named(map$mixoplankton, "MFT")
  expect_named(map$phytoplankton, c("Shape", "Elongation", "Group"))
  expect_named(map$habs, "Definition")
})

test_that("trait aliases resolve database-prefixed columns", {
  map <- plankannot:::plankannot_default_trait_map()
  cols <- c("ASV", "Sample_1", "mixoplankton_MFT (MIXOPLANKTON FUNCTIONAL CLASSIFICATION)")
  resolved <- plankannot:::plankannot_resolve_traits(cols, map$mixoplankton)
  expect_equal(unname(resolved[["MFT"]]), cols[[3]])
})

test_that("abundance long data ignores empty categories and non-positive abundances", {
  df <- data.frame(ASV = c("a", "b", "c"), Sample_1 = c(10, 20, 0), trait = c("A", "A", ""), check.names = FALSE)
  out <- plankannot:::plankannot_trait_long_data(df, "trait", "trait", "Sample_1")
  expect_equal(nrow(out), 2)
  expect_equal(sum(out$Abundance), 30)
})
