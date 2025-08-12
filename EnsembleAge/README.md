# EnsembleAge

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D3.5.0-blue.svg)](https://www.r-project.org/)
[![Paper DOI](https://img.shields.io/badge/DOI-10.1007%2Fs11357--025--01808--1-blue)](https://doi.org/10.1007/s11357-025-01808-1)
[![Software DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16814673.svg)](https://doi.org/10.5281/zenodo.16814673)

**Multi-platform epigenetic age prediction using ensemble clock methods**

## 📖 Description

EnsembleAge provides tools for predicting epigenetic age using various clock methods including Universal Clocks, Elastic Epigenetic clocks, static ensemble methods, and dynamic ensemble clocks. The package **automatically detects your data platform** and prepares the data accordingly for accurate age prediction across various tissues and species.

**📄 This package accompanies the open access publication:**
> Haghani, A., et al. (2025). EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework. *GeroScience*. DOI: [10.1007/s11357-025-01808-1](https://doi.org/10.1007/s11357-025-01808-1)

The clocks and methods implemented in this package are freely available for research use under the MIT license.

## 🚀 Quick Start

**The easiest way to predict epigenetic age:**

```r
library(EnsembleAge)

# Ultimate one-liner - automatically detects platform and runs predictions
results <- predictEnsemble("path/to/methylation_data.RDS", "path/to/sample_sheet.csv")

# View predictions
head(results)
```

**That's it!** EnsembleAge automatically:
- ✅ Detects your platform (Human, Mammal40k, Mammal320k)
- ✅ Preprocesses your data correctly  
- ✅ Selects optimal clock methods
- ✅ Handles missing probes with imputation
- ✅ Generates comprehensive age predictions

## 📦 Installation

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install EnsembleAge
devtools::install_github("ahaghani/EnsembleAge")
```

### Dependencies

The package automatically installs these required R packages:
- `dplyr`, `tibble`, `tidyr`, `plyr`, `stringr`, `data.table`

## 💡 Simple Examples

### Example 1: Mammal40k Data

```r
library(EnsembleAge)

# Create simulated mammal40k methylation data
set.seed(123)
n_probes <- 1000  # Using subset for example
n_samples <- 10

# Methylation data: probes as rows, samples as columns, with CGid
methylation_data <- data.frame(
  CGid = paste0("cg", sprintf("%08d", 1:n_probes)),
  matrix(runif(n_probes * n_samples, 0.1, 0.9), nrow = n_probes)
)
names(methylation_data)[-1] <- paste0("Sample_", 1:n_samples)

# Sample sheet: required columns
sample_sheet <- data.frame(
  Basename = paste0("Sample_", 1:n_samples),
  Age = c(0.5, 1.2, 2.1, 0.8, 1.5, 2.8, 1.1, 0.9, 2.2, 1.7),  # Mouse ages in years
  Female = c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1),  # 1=female, 0=male
  Tissue = "Blood",
  SpeciesLatinName = "Mus musculus"
)

# Run predictions
results <- predictEnsemble(methylation_data, sample_sheet)

# View results
head(results)
#   Basename  Age epiClock   epiAge AgeAcceleration Female Tissue clockFamily
# 1 Sample_1  0.5   Static     0.8            0.12      1  Blood   Ensemble.Static
# 2 Sample_1  0.5 Dynamic1     0.7            0.05      1  Blood   EnsembleAge.Dynamic
# ...
```

### Example 2: Mammal320k Data (with file paths)

```r
# For mammal320k array data stored in files:
results <- predictEnsemble("mammal320k_data.RDS", "mouse_samples.csv")

# View results in wide format (traditional output)
results_wide <- predictEnsemble("mammal320k_data.RDS", "mouse_samples.csv", 
                                output_format = "wide")
head(results_wide)
#   Basename  Age  Ensemble.Static.clock.epiAge  EnsembleAge.Dynamic.X1.clock.epiAge
# 1 Mouse_01  1.2                          1.45                                  1.38
```

### Example 3: Human Data (with file paths)

```r
# For human EPIC array data stored in files:
results <- predictEnsemble("human_methylation.RData", "human_samples.csv")

# Get both long and wide format results
results_both <- predictEnsemble("human_methylation.RData", "human_samples.csv", 
                                output_format = "both")

long_format <- results_both$long    # One row per sample-clock combination
wide_format <- results_both$wide    # One row per sample, clocks as columns
```

### Example 4: Auto-generate Sample Sheet

```r
# If you don't have a sample sheet, let EnsembleAge create a template
results <- predictEnsemble(methylation_data, samples = NULL)

# This creates 'user_data_sample_sheet_template.csv' 
# Fill it out with your actual sample information and rerun
```

## 📊 Data Structure Guide

### Methylation Data Format
```r
# Required format: probes as rows, samples as columns
#        CGid   Sample_1  Sample_2  Sample_3
# 1  cg00001     0.234     0.567     0.123
# 2  cg00002     0.789     0.345     0.678
# 3  cg00003     0.456     0.234     0.567
```

### Sample Sheet Format
```r
# Required columns: Basename, Age
# Optional columns: Female, Tissue, SpeciesLatinName
#   Basename   Age  Female    Tissue  SpeciesLatinName
# 1 Sample_1   1.2       1     Blood      Mus musculus
# 2 Sample_2   2.1       0     Liver      Mus musculus  
# 3 Sample_3   0.8       1    Kidney      Mus musculus
```

### Output Format
```r
# Long format (default): one row per sample-clock combination
#   Basename  Age epiClock  epiAge AgeAcceleration  clockFamily
# 1 Sample_1  1.2   Static    1.45           0.25  Ensemble.Static
# 2 Sample_1  1.2 Dynamic1    1.38           0.18  EnsembleAge.Dynamic

# Wide format: one row per sample, clocks as columns  
#   Basename  Age  Ensemble.Static.clock.epiAge  EnsembleAge.Dynamic.X1.clock.epiAge
# 1 Sample_1  1.2                           1.45                                  1.38
```

## 🔧 Supported Platforms

EnsembleAge **automatically detects and handles** multiple methylation array platforms:

| Platform | Probe Count | File Types | Auto-Detection | File Path Support |
|----------|-------------|------------|----------------|-------------------|
| **Human EPIC/450k** | ~850k/450k | `.RData`, `.RDS`, `.csv`, `.txt` | ✅ Automatic | ✅ Yes |
| **Mammal320k** | ~320k | `.RData`, `.RDS`, `.csv`, `.txt` | ✅ Automatic | ✅ Yes |
| **Mammal40k** | ~40k | `.RData`, `.RDS`, `.csv`, `.txt` | ✅ Automatic | ✅ Yes |

**File Path Support**: All platforms work with both file paths (`"data.RDS"`) and pre-loaded R objects (`my_data_frame`).

### ⚡ What Happens Automatically

- 🔍 **Platform Detection**: Identifies your array type from data dimensions and patterns
- 🔄 **Data Orientation**: Detects if probes are rows or columns and fixes automatically  
- 🧬 **Missing Probe Handling**: Imputes missing CpG sites with neutral values (0.5)
- 📊 **Smart Method Selection**: Chooses best clock methods based on data coverage
- 💾 **Multiple Output Formats**: Long format (analysis-ready) or wide format (traditional)

## 📋 Data Requirements

### Methylation Data
- **Format**: Data frame or matrix with CpG sites and sample methylation values
- **Values**: Beta values (0-1 range) representing methylation levels
- **Orientation**: Probes as rows OR columns (auto-detected and fixed)
- **File Types**: `.RDS`, `.RData`, `.csv` supported

### Sample Sheet
- **Required**: `Basename` (sample IDs), `Age` (chronological age in years)
- **Optional**: `Female` (1=female, 0=male), `Tissue`, `SpeciesLatinName`
- **Format**: CSV file or R data frame

### 🎯 Quick Tips
- **Ages**: Use years (e.g., 1.5 for 18-month-old mouse, 25.3 for human)
- **Missing data**: Package handles missing probes automatically
- **No sample sheet?** Use `samples = NULL` to auto-generate a template
- **File paths**: Use full paths or ensure files are in working directory

## 🎯 Advanced Usage

### Multiple Output Formats

```r
# Get long format (default) - best for analysis
results_long <- predictEnsemble(methylation_data, sample_sheet, output_format = "long")

# Get wide format - one row per sample
results_wide <- predictEnsemble(methylation_data, sample_sheet, output_format = "wide") 

# Get both formats
results_both <- predictEnsemble(methylation_data, sample_sheet, output_format = "both")
```

### Method Selection

```r
# Automatic method selection (recommended)
results <- predictEnsemble(methylation_data, sample_sheet)  # Auto-selects best methods

# Available methods: "ensemble_static", "ensemble_dynamic", "ensemble_dual", "all"
```

### Working with Files vs Objects

```r
# Using file paths (works for ALL platforms: Human, Mammal320k, Mammal40k)
results <- predictEnsemble("methylation.RDS", "samples.csv")           # RDS files
results <- predictEnsemble("data.RData", "samples.csv")                # RData files (all platforms)
results <- predictEnsemble("mammal320k_data.RData", "mouse_samples.csv") # Mammal320k RData
results <- predictEnsemble("mammal40k_data.RData", "mouse_samples.csv")  # Mammal40k RData

# Using R objects directly  
results <- predictEnsemble(my_methylation_df, my_sample_df)

# Mixed approach (file + object, or object + file)
results <- predictEnsemble("methylation.RDS", my_sample_df)           # File + object
results <- predictEnsemble(my_methylation_df, "samples.csv")          # Object + file
```

### 💾 Result Saving Control

EnsembleAge gives you full control over where and how your results are saved:

```r
# Default: saves to current directory with automatic naming
results <- predictEnsemble("data.RDS", "samples.csv")
# Creates: EnsembleAge_data_long_20250811_143052.csv

# Custom directory
results <- predictEnsemble("data.RDS", "samples.csv", 
                          output_dir = "results/2024/study1/")

# Custom prefix instead of "EnsembleAge"
results <- predictEnsemble("data.RDS", "samples.csv", 
                          output_prefix = "MouseStudy")
# Creates: MouseStudy_data_long_20250811_143052.csv

# Completely custom filename
results <- predictEnsemble("data.RDS", "samples.csv", 
                          output_filename = "my_final_results.csv")

# Turn off saving entirely (just return results)
results <- predictEnsemble("data.RDS", "samples.csv", 
                          save_results = FALSE)

# Both formats with custom naming
results <- predictEnsemble("data.RDS", "samples.csv", 
                          output_format = "both",
                          output_dir = "final_results/",
                          output_prefix = "Project_Alpha")
# Creates: Project_Alpha_data_long_20250811_143052.csv
#          Project_Alpha_data_wide_20250811_143052.csv
```

**Automatic Features:**
- ✅ **Directory Creation**: Output directories are created automatically if they don't exist
- ✅ **Timestamp**: Automatic timestamps prevent file overwrites  
- ✅ **Smart Naming**: Custom filenames work with "both" format (adds `_long` and `_wide`)
- ✅ **File Size**: Shows saved file size for easy reference

### 📊 Clock Families Included

| Clock Family | Best For | # Clocks | Description | Reference |
|--------------|----------|----------|-------------|-----------|
| **EnsembleAge.Static** | 🐭 Mouse (primary) | 2 | Main mouse ensemble clocks | [GeroScience 2025](https://doi.org/10.1007/s11357-025-01808-1) |
| **EnsembleAge.Dynamic** | 🐭 Mouse (detailed) | 50 | Individual specialized mouse clocks | [GeroScience 2025](https://doi.org/10.1007/s11357-025-01808-1) |
| **EnsembleDualAge.Static** | 👨 Human (primary) | 1 | Human-mouse optimized clock | [GeroScience 2025](https://doi.org/10.1007/s11357-025-01808-1) |
| **UniClock2/3** | 🌐 Cross-species | 6 | Universal mammalian clocks | [Nature Aging 2023](https://doi.org/10.1038/s43587-023-00462-6) |
| **LifespanUberClock** | 🐭 Mouse variants | 12 | Lifespan-focused clocks | [bioRxiv 2022](https://doi.org/10.1101/2022.01.16.476530) |
| **DNAmAge*** | 🐭 Mouse categories | 39 | Development, Elastic, Intervention clocks | [eLife 2022](https://doi.org/10.7554/eLife.75244) |

**Total: 110+ individual clocks** covering development, aging, interventions, and cross-species predictions.

## 🛠️ Troubleshooting

### Common Issues

**"Cannot determine data orientation"**
- Ensure your data has a `CGid` column or CpG IDs as row names
- For human data, make sure methylation values are present (not just annotations)

**"No matching samples found"**  
- Check that `Basename` column in sample sheet matches column names in methylation data
- Use `head(colnames(methylation_data))` and `head(sample_sheet$Basename)` to compare

**"Low probe coverage warnings"**
- Normal for human data (relies on imputation)
- For mouse data, ensure you're using the correct platform (Mammal40k vs Mammal320k)

### Getting Help

```r
# Check what platform was detected
platform <- detect_platform(your_data)
cat("Detected platform:", platform)

# Examine your data structure  
str(your_methylation_data)
head(your_sample_sheet)

# Test with small subset first
small_test <- predictEnsemble(your_data[1:100, 1:3], your_samples[1:3, ])
```

## 📈 Understanding Results

### Output Columns

- **Basename**: Sample identifier  
- **Age**: Input chronological age
- **epiClock**: Clock name used for prediction
- **epiAge**: Predicted epigenetic age
- **AgeAcceleration**: Age acceleration (residualized difference)
- **Female**: Sex (1=female, 0=male, NA=unknown)
- **Tissue**: Tissue type
- **clockFamily**: Clock family group

### Interpreting Results

- **epiAge**: The biological/epigenetic age predicted by the clock
- **AgeAcceleration**: Positive = aging faster, Negative = aging slower
- **Multiple clocks**: Each row represents one clock prediction per sample

## 👥 Author & Citation

**Author:** Amin Haghani (Altos Labs)  
**Email:** dr.a.haghani@gmail.com

### Citation

If you use this package in your research, please cite:

**EnsembleAge Package:**
```
Haghani, A., et al. (2025). EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework. GeroScience.
DOI: 10.1007/s11357-025-01808-1
```

**Universal Clocks:**
```
Lu, A.T., Fei, Z., Haghani, A. et al. (2023). Universal DNA methylation age across mammalian tissues. 
Nature Aging 3, 1144–1166. DOI: 10.1038/s43587-023-00462-6
```

**Original Mammalian Clocks:**
```
Haghani, A., et al. (2022). Methylation-based epigenetic clocks for mammalian species. eLife 11:e75244.
DOI: 10.7554/eLife.75244
```

**Lifespan Uber Clock:**
```
Lu, A.T., et al. (2022). Universal methylation clocks enable analysis of aging effects in mammals.
bioRxiv. DOI: 10.1101/2022.01.16.476530
```

**BibTeX:**
```bibtex
@article{haghani2025ensembleage,
  title={EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework},
  author={Haghani, Amin and others},
  journal={GeroScience},
  year={2025},
  doi={10.1007/s11357-025-01808-1}
}

@article{lu2023universal,
  title={Universal DNA methylation age across mammalian tissues},
  author={Lu, Ake T and Fei, Zhe and Haghani, Amin and others},
  journal={Nature Aging},
  volume={3},
  pages={1144--1166},
  year={2023},
  doi={10.1038/s43587-023-00462-6}
}

@article{haghani2022methylation,
  title={Methylation-based epigenetic clocks for mammalian species},
  author={Haghani, Amin and others},
  journal={eLife},
  volume={11},
  pages={e75244},
  year={2022},
  doi={10.7554/eLife.75244}
}

@article{lu2022lifespan,
  title={Universal methylation clocks enable analysis of aging effects in mammals},
  author={Lu, Ake T and others},
  journal={bioRxiv},
  year={2022},
  doi={10.1101/2022.01.16.476530}
}
```

## 📄 License

MIT License - Free for research and commercial use.

## 🤝 Contributing & Issues

- **Contributing**: Contributions are welcome! Please submit a Pull Request.
- **Issues**: Report bugs or request features on our [GitHub issues page](https://github.com/ahaghani/EnsembleAge/issues).
- **Questions**: Contact Amin Haghani at dr.a.haghani@gmail.com

---

**🎯 TL;DR: Just run `predictEnsemble(your_data, your_samples)` - everything else is automatic!**