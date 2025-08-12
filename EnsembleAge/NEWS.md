# EnsembleAge NEWS

## Version 0.99.0

### NEW FEATURES
- Initial release of EnsembleAge package
- Support for multiple methylation array platforms (Mammal320k, Mammal40k, Human EPIC, Human 450k)
- Automatic platform detection and data preprocessing
- Multiple clock families: EnsembleAge.Static, EnsembleAge.Dynamic, EnsembleDualAge.Static, UniClock2/3, LifespanUberClock, DNAmAge*
- Cross-species epigenetic age prediction capabilities

### IMPROVEMENTS
- Efficient data loading for large datasets
- Automatic probe mapping and orientation detection
- Missing probe imputation with neutral methylation values
- Comprehensive error handling and validation

### BUG FIXES
- Fixed probe mapping issues for different platform formats
- Corrected age acceleration calculations
- Improved sample matching between data and sample sheets

### DOCUMENTATION
- Comprehensive README with usage examples
- Vignette with detailed workflow instructions
- Complete function documentation with examples

### CITATION
This package accompanies the publication: Haghani, A., et al. (2025). EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework. GeroScience. DOI: 10.1007/s11357-025-01808-1 