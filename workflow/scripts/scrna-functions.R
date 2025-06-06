#!/usr/bin/env Rscript



function_pca_dimensions <- function(Spatial_Data) {
  pct <- Stdev(object = Spatial_Data, reduction = "pca") / sum(Stdev(object = Spatial_Data, reduction = "pca")) * 100
  cum <- cumsum(pct)
  
  diffs <- pct[1:(length(pct)-1)] - pct[2:length(pct)]
  cat("\nDiffs between consecutive PCs (first 20):\n")
  print(round(diffs[1:20], 4))
  
  co1 <- which(cum > 90 & pct < 5)[1]
  co2 <- sort(which(diffs > 0.1), decreasing = TRUE)[1] + 1
  
  cat("\nco1:", co1, "\n")
  cat("co2:", co2, "\n")
  
  dimensionReduction <- min(co1, co2, na.rm = TRUE)
  return(dimensionReduction)
}


c(
  "#c15035",
  "#79ce5d",
  "#8449c1",
  "#cbbc50",
  "#c75291",
  "#83c8b7",
  "#4d2f48",
  "#c59581",
  "#8a90c4",
  "#506138"
) -> palette10

palette30 <- c(
  "#d2a33f",
  "#6844cb",
  "#81dc3c",
  "#c54cdb",
  "#6cda74",
  "#de43ac",
  "#d2d840",
  "#422679",
  "#c4da81",
  "#974096",
  "#659a3d",
  "#6677cd",
  "#db622a",
  "#66d6a9",
  "#d83e46",
  "#8dd4d6",
  "#ce4b7a",
  "#3e602b",
  "#ce90d4",
  "#8b7137",
  "#373053",
  "#d0c9a1",
  "#7d3d62",
  "#608e79",
  "#88352a",
  "#779ec5",
  "#d3866a",
  "#324a45",
  "#c898a6",
  "#503029"
)


palette50 <- c(
  "#42566f",
  "#66e741",
  "#7744d9",
  "#a0e340",
  "#d24cdd",
  "#5bb734",
  "#3d2897",
  "#e1e73b",
  "#9638ae",
  "#5fe277",
  "#e344a0",
  "#579836",
  "#6a7de4",
  "#a5bb36",
  "#5d50a4",
  "#c2e778",
  "#632267",
  "#84df99",
  "#a64194",
  "#e0bd42",
  "#31245c",
  "#d6d58a",
  "#c783de",
  "#4c9e66",
  "#e2462a",
  "#62dabf",
  "#db445e",
  "#77cfe3",
  "#d7862e",
  "#5c8ec9",
  "#ac482c",
  "#649aa4",
  "#9e345e",
  "#c4dfc5",
  "#251e33",
  "#948d34",
  "#e08bbb",
  "#416126",
  "#bcb4da",
  "#2f3922",
  "#de986c",
  "#4f2330",
  "#d8b6ac",
  "#6c2d23",
  "#417061",
  "#c5787b",
  "#909a6f",
  "#876690",
  "#825b26",
  "#84675c"
)

palette75 <- c(
  "#a95cdc",
  "#45c558",
  "#cd3ca6",
  "#72be3e",
  "#7570ec",
  "#aabf32",
  "#3d57c9",
  "#dcb434",
  "#5c86eb",
  "#e99827",
  "#6e4ca8",
  "#42972e",
  "#9b3ba5",
  "#38c481",
  "#d63b83",
  "#30853c",
  "#db74da",
  "#7c9828",
  "#9975d7",
  "#88bf5f",
  "#d864aa",
  "#38711a",
  "#c28bdb",
  "#bcb84e",
  "#3d67b2",
  "#e5682b",
  "#4cbee0",
  "#cd3c27",
  "#5eccb7",
  "#dc385a",
  "#80c687",
  "#8d5091",
  "#b59b31",
  "#6d65a9",
  "#c18b2b",
  "#4b96d2",
  "#b35714",
  "#9e9ee1",
  "#e68e40",
  "#5c6a9f",
  "#897b1c",
  "#ca85b8",
  "#66781d",
  "#ec88ae",
  "#347f4e",
  "#a84064",
  "#5ca267",
  "#c8515a",
  "#3da390",
  "#e76f57",
  "#2e8568",
  "#a1442f",
  "#1a6447",
  "#f1957e",
  "#296437",
  "#d2797d",
  "#4b600e",
  "#914b69",
  "#709145",
  "#bc6b45",
  "#507738",
  "#d0874b",
  "#476025",
  "#e9b176",
  "#616117",
  "#c5b672",
  "#814a28",
  "#a3b26a",
  "#765113",
  "#c18f60",
  "#62612c",
  "#9a6623",
  "#8d8444",
  "#91633b",
  "#705d18"
)

palette200 <- c(
  "#0060d6",
  "#b1019a",
  "#12e073",
  "#a42eb5",
  "#67db55",
  "#6c3fc4",
  "#84db4c",
  "#2b51d5",
  "#a7cb1c",
  "#5a5ae2",
  "#70ba19",
  "#4268f1",
  "#53c137",
  "#ae2db4",
  "#00b53d",
  "#b253da",
  "#519d00",
  "#af68f2",
  "#00a030",
  "#ec77fe",
  "#01d06f",
  "#cb38bb",
  "#008812",
  "#ca4fd3",
  "#699e00",
  "#9674ff",
  "#efc114",
  "#7a75ff",
  "#a2d650",
  "#8222a3",
  "#91d86e",
  "#eb4ecb",
  "#009538",
  "#de29a7",
  "#3bdf9e",
  "#f535a6",
  "#00791c",
  "#e784ff",
  "#518a00",
  "#327fff",
  "#aaa600",
  "#65b299",
  "#d1ca3f",
  "#7482ff",
  "#899400",
  "#0067d5",
  "#ff9d2e",
  "#304aaf",
  "#b8d155",
  "#5540aa",
  "#86da76",
  "#b10088",
  "#63dd8f",
  "#cd0088",
  "#00dfb2",
  "#f92564",
  "#02cdab",
  "#ff419d",
  "#01ae72",
  "#ff60c7",
  "#007e37",
  "#ff80f3",
  "#3c7000",
  "#8b218f",
  "#a9d46c",
  "#6f36a0",
  "#f2bf44",
  "#828eff",
  "#ff8d28",
  "#027ad7",
  "#ff7d2c",
  "#5f9aff",
  "#d46a00",
  "#008ce0",
  "#cd4b00",
  "#41adff",
  "#c22506",
  "#36d9ec",
  "#c80022",
  "#48dbc9",
  "#f42e4e",
  "#5fdcac",
  "#d70066",
  "#00975a",
  "#bd0070",
  "#02a97d",
  "#ca0063",
  "#006518",
  "#ff84e5",
  "#245f02",
  "#e296ff",
  "#a09000",
  "#015aac",
  "#c38d00",
  "#6ba6ff",
  "#ff6a35",
  "#01abf1",
  "#cc1f1b",
  "#02ace4",
  "#b63500",
  "#0175bc",
  "#c55b00",
  "#9ba6ff",
  "#a62300",
  "#01b5a0",
  "#cd0032",
  "#01a082",
  "#c80048",
  "#00844e",
  "#ff598e",
  "#8cd79f",
  "#9a1173",
  "#b0d276",
  "#8b2780",
  "#c1ce6e",
  "#384d9d",
  "#ffa951",
  "#654191",
  "#e3c363",
  "#c99cff",
  "#757800",
  "#f1afff",
  "#375c11",
  "#ff75bf",
  "#175f36",
  "#ff586d",
  "#7dd7ba",
  "#a90111",
  "#009abb",
  "#ff6a45",
  "#0086b4",
  "#ff8c52",
  "#624586",
  "#ffb560",
  "#80347c",
  "#a5d38d",
  "#9b1568",
  "#97d4ac",
  "#b00052",
  "#006751",
  "#ff6c63",
  "#2f5d29",
  "#ff96dd",
  "#48590b",
  "#ff72a5",
  "#625d00",
  "#d7bafd",
  "#ae6a00",
  "#546a9c",
  "#aa5c00",
  "#9c9ad0",
  "#956a00",
  "#f8afec",
  "#816e00",
  "#743e78",
  "#d6c76e",
  "#942858",
  "#c0cd88",
  "#75406d",
  "#ecbf74",
  "#7b3e63",
  "#fab974",
  "#845a85",
  "#ff9962",
  "#dbafe0",
  "#953004",
  "#688d5f",
  "#9c2521",
  "#849f6c",
  "#ff849c",
  "#694f06",
  "#ffa8cd",
  "#784c00",
  "#ff96ae",
  "#716d39",
  "#ff7a68",
  "#7b6535",
  "#cb8bac",
  "#8f4f00",
  "#823a5b",
  "#dac589",
  "#8e3248",
  "#ffa369",
  "#853c41",
  "#f3ba94",
  "#913235",
  "#bda370",
  "#903425",
  "#ea9eaa",
  "#804212",
  "#f2ac9c",
  "#873d12",
  "#ffad8e",
  "#724a26",
  "#ffa191",
  "#81412e",
  "#ff926b",
  "#a8626f",
  "#c68872",
  "#a56956"
)

function_gofunc <- function(df, algorithm = "weight01", statistics = "ks", mapping = "org.Hs.eg.db", ID = "symbol", ontology = "BP") {
  geneList <- df$p_val
  names(geneList) <- df$gene
  # Create topGOData object
  GOdata <- new("topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x,
    annot = annFUN.org, mapping = mapping, ID = ID, nodeSize = 10
  )

  resultsKS <- runTest(GOdata, algorithm = algorithm, statistic = statistics)


  tab <- GenTable(GOdata, raw.p.value = resultsKS, topNodes = length(resultsKS@score), numChar = 120)

  return(tab)
}


theme_cellsnake_classic <- function(base_size = 12, base_family = "Ubuntu") {
  theme_bw() %+replace%
    theme(
      panel.grid.major = element_line(color = "white"),
      # panel.background = element_rect(fill = "lightblue"),
      # panel.border = element_rect(color = "lightblue", fill = NA),
      # axis.line = element_line(color = "lightblue"),
      # axis.ticks = element_line(color = "lightblue"),
      # axis.text = element_text(color = "steelblue")
    )
}


function_color_palette <- function(n) {
  if (n > 75) {
    return(palette200)
  }
  if (n > 10 && n <= 30) {
    return(palette30)
  }
  if (n > 30 && n <= 50) {
    return(palette50)
  }

  if (n <= 10) {
    return(palette10)
  }
  if (n > 50 && n <= 75) {
    return(palette75)
  }

  if (!requireNamespace("randomcoloR", quietly = TRUE)) {
    install.packages("randomcoloR")
  }

  randomcoloR::distinctColorPalette(n)
}