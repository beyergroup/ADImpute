### Example run (from the demo folder)

devtools::load_all("../")

network_path = # <absolute path to .rds or .txt coefficients file>

# Load demo data (counts) and convert to TPM
data("demo_data", package = "ADImpute")
TPM <- NormalizeTPM(demo_data)
WriteTXT(TPM, "TPM.txt")
tpm_path <- paste0(getwd(),"/TPM.txt")

# Train method to obtain optimal per gene
methods_pergene <- EvaluateMethods(data = TPM,
                                   training.ratio = .7,
                                   mask.ratio = .2,
                                   training.only = T,
                                   split.seed = 12,
                                   mask.seed = 34,
                                   type = "TPM",
                                   cell.clusters = 2,
                                   cores = 4,
                                   cluster.type = "SOCK",
                                   network.path = network_path,
                                   drop.exclude = T)

# Impute full dataset
imputed <- Impute(method.choice = methods_pergene,
                  data = TPM,
                  count_path = tpm_path,
                  type = "TPM",
                  cell.clusters = 2,
                  cores = 4,
                  cluster.type = "SOCK",
                  network.path = network_path,
                  drop.exclude = T)
