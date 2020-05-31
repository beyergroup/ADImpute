### Example run (from the demo folder)

devtools::load_all("../")

network_path = # <absolute path to .rds or .txt coefficients file>

# Load demo data (counts) and normalize
data("demo_data", package = "ADImpute")
RPM <- NormalizeRPM(demo_data)
WriteTXT(RPM, "RPM.txt")
rpm_path <- paste0(getwd(),"/RPM.txt")

# Train method to obtain optimal per gene
methods_pergene <- EvaluateMethods(data = RPM,
                                   training.ratio = .7,
                                   mask.ratio = .2,
                                   training.only = T,
                                   split.seed = 12,
                                   mask.seed = 34,
                                   type = "count",
                                   cell.clusters = 2,
                                   cores = 4,
                                   cluster.type = "SOCK",
                                   network.path = network_path,
                                   drop.exclude = T)

# Impute full dataset
imputed <- Impute(do = "Ensemble",
                  method.choice = methods_pergene,
                  data = RPM,
                  count_path = rpm_path,
                  type = "count",
                  cell.clusters = 2,
                  cores = 4,
                  cluster.type = "SOCK",
                  network.path = network_path,
                  drop.exclude = T)
