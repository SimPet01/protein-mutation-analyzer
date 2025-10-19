# Protein Mutation Analyzer - Usage Guide

## Quick Start

### 1. Installation
```bash
# Clone the repository
git clone https://github.com/SimPet01/protein-mutation-analyzer.git
cd protein-mutation-analyzer

# Install dependencies
pip install -r requirements.txt
```

### 2. Prepare Your Data

#### Protein sequence (FASTA format)
Create a file `data/my_protein.fasta`:
```
>sp|P12345|MY_PROTEIN My protein description
MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASED
```

#### Mutations (CSV format)
Create a file `data/my_mutations.csv`:
```csv
position,original_aa,mutant_aa,description
7,E,K,Glutamic acid to Lysine
22,I,V,Isoleucine to Valine
```

**Important:** 
- Positions use 1-based indexing (as in biology)
- Use single-letter amino acid codes
- Ensure the original amino acid matches the sequence

### 3. Run Analysis
```python
from scripts.analyzer import ProteinMutationAnalyzer

# Create analyzer instance
analyzer = ProteinMutationAnalyzer('data/my_protein.fasta')

# Load mutations
analyzer.load_mutations('data/my_mutations.csv')

# Analyze
results = analyzer.analyze_all_mutations()

# Generate outputs
analyzer.save_results('output/results.csv')
analyzer.visualize_mutations('output/visualization.png')
analyzer.generate_html_report('output/report.html')
```

## Understanding the Results

### Mutation Classification

- **Conservative**: Minor changes (e.g., L→I, similar hydrophobic residues)
  - Small MW change (< 5 Da)
  - No charge change
  - Small hydrophobicity change (< 0.5)

- **Moderate**: Intermediate changes
  - Medium impact on properties

- **Radical**: Major changes (e.g., E→K, charge reversal)
  - Large MW change (> 20 Da)
  - OR charge change
  - OR large hydrophobicity change (> 1.5)

### Amino Acid Groups

- **Hydrophobic**: A, V, I, L, M, F, W, P (water-avoiding)
- **Polar**: S, T, C, Y, N, Q (water-friendly, uncharged)
- **Charged+**: K, R, H (positively charged at pH 7)
- **Charged-**: D, E (negatively charged at pH 7)
- **Special**: G (smallest, most flexible)

## Output Files

### 1. CSV Report (`mutation_analysis.csv`)
Tab-separated table with all calculated properties for each mutation.

### 2. Visualization (`mutations_visualization.png`)
Four-panel figure showing:
- Molecular weight changes
- Charge changes
- Hydrophobicity changes
- Classification distribution

### 3. HTML Report (`mutation_report.html`)
Professional-looking report with color-coded mutations, viewable in any web browser.

## Troubleshooting

### Error: "Position X contains Y, not Z"
The amino acid you specified as original doesn't match the sequence at that position. Check:
- Position number (remember 1-based indexing)
- Amino acid code (must match sequence exactly)
- That you're using the correct protein sequence

### Error: "Position X out of range"
The position number exceeds the sequence length. Check:
- Sequence length matches your expectations
- Position numbers in CSV are correct

## Advanced Usage

### Custom Output Locations
```python
analyzer.save_results('custom/path/results.csv')
analyzer.visualize_mutations('custom/path/viz.png')
analyzer.generate_html_report('custom/path/report.html')
```

### Accessing Individual Results
```python
# Get results as DataFrame
results = analyzer.results_df

# Filter radical mutations
radical = results[results['classification'] == 'Radical']

# Get specific columns
important_cols = results[['position', 'mutation', 'classification']]
```

## Contact

For questions or issues, please open an issue on GitHub or contact me via email 'karelDOTcurinaATseznamDOTcz'.