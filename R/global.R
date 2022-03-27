## Set gloable variables to prevent no visible binding for global variable warning
utils::globalVariables(c(
  # checkPreloadKinaseTable
  "ID", "ABBREVIATION", "SPECIES", "kinasesPhosphatases",
  # confirmColumnNames
  "uniqueID", "seqWindowID", "CON", "REPs",
  # processData
  "Reverse", "Potential contaminant", "measure", "obs_tr", "Time",
  "intensity", "Localization Probability", "value", "ReplicateCount",
  "Amino Acid", "Protein", "UniqueID", "SampleID", "name", "set", "Amino acid",
  "condition", "obs_tp", "normValue", "m", "experiement", "AOV", "result",
  "argValue", "p.value", "term", "Multiplicity", "Reps",
  # findClusters
  "timepoint", "dist", "hclust", "Full", "Short", "Cluster", "cluster",
  # findPeptideIDs
  "Part1", "PeptideID",
  # networkAnlaysis
  "discreteState", "experiment", "potentialEdges", "regulator", "coReg",
  "target", "neighborCombinations", "searchSizeK", "regulatorisKinase",
  "searchID", "qi", "Nj", "Nj.count", "BDei", "regulatorCondition",
  "targetCondition", "Reg", "regulatorState", "coReg_1", "coReg_2",
  "regulator_1", "regulator_2", "percentChange", "regulatorChange",
  "targetChange", "FALSE", "TRUE", "timelapse",
  # plotClusterHeatmap
  "cutreeDynamic", ".", "foundClusterIDs",
  # plotDistribution
  "ConTime", "Value", "TimePoint",
  # regulationCheck
  "argValue", "argValue == 0", "tp", "Model_name", "conditionMedian", "mvalue",
  "change", "sigChange", "sigChangeSign", "# Fold change occurances",
  # validateKinaseTable
  "ABB", "Part1", "AA", "multiplicity", "FAMILY", "Total",
  # PlotMultiPeptides
  "ConDesign", "txtListNew",
  # plotNetPhorce
  "UniqueID_V2",
  # plotPCA
  "data.norm", "PC1", "PC2",
  #plotRegulation
  "avgValue", "signChangeText",
  #plotSInglePeptides
  "Normalized Value", "Description", "Replication",
  #plot_normalization
  "val", "var",
  # extractSummaryTable
  "Sample", "Sequence Window"
  ))



