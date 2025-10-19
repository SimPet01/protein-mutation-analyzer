# Installation and Testing Checklist

This document walks through installation and testing from scratch.
Use this to verify the project works for new users.

## Prerequisites Check

- [ ] Python 3.8+ installed (`python --version`)
- [ ] pip installed (`pip --version`)
- [ ] Git installed (`git --version`)

## Installation Steps

### 1. Clone Repository
```bash
git clone https://github.com/[your-username]/protein-mutation-analyzer.git
cd protein-mutation-analyzer
```

### 2. Verify Project Structure
Check that following directories and files exist:
- [ ] `data/` folder with FASTA and CSV files
- [ ] `scripts/` folder with Python files
- [ ] `output/` folder (may be empty)
- [ ] `README.md`
- [ ] `requirements.txt`

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

Expected output: All packages install successfully

Verify installation:
```bash
python -c "import pandas; import Bio; import matplotlib; print('All packages imported successfully')"
```

## Testing Steps

### Test 1: Run Main Analysis
```bash
cd scripts
python analyzer.py
```

Expected output:
- [ ] "Loaded sequence: 154 amino acids"
- [ ] "Protein ID: sp|P02144|MYO_HUMAN"
- [ ] "Loaded 5 mutations"
- [ ] Analysis results table displayed
- [ ] "Results saved to ../output/mutation_analysis.csv"
- [ ] "Visualization saved to ../output/mutations_visualization.png"
- [ ] "HTML report saved to ../output/mutation_report.html"
- [ ] "✓ Analysis complete!"

### Test 2: Verify Output Files
Check that following files exist in `output/` folder:
- [ ] `mutation_analysis.csv` (readable in text editor or Excel)
- [ ] `mutations_visualization.png` (4-panel graph, clear and readable)
- [ ] `mutation_report.html` (opens in browser, shows colored table)

### Test 3: Run Alternative Data Test
```bash
python test_analyzer.py
```

Expected output:
- [ ] "TESTING ANALYZER WITH DIFFERENT DATA"
- [ ] Test completes successfully
- [ ] "✓ TEST PASSED"
- [ ] Three new files created with "test_" prefix

### Test 4: Code Quality Check
Verify docstrings work:
```bash
python -c "from analyzer import ProteinMutationAnalyzer; help(ProteinMutationAnalyzer)"
```

Expected: Detailed help text with class description

## Common Issues

### Issue: "ModuleNotFoundError: No module named 'Bio'"
**Solution**: Run `pip install biopython`

### Issue: "FileNotFoundError: [Errno 2] No such file or directory"
**Solution**: Make sure you're running scripts from the `scripts/` directory

### Issue: Output folder doesn't exist
**Solution**: Create it manually: `mkdir output`

### Issue: Permission denied when saving files
**Solution**: Check folder permissions, may need to run with appropriate rights

## Success Criteria

All tests pass if:
✅ No error messages during execution
✅ All expected output files are created
✅ Visualizations display correctly
✅ HTML report opens and renders properly
✅ CSV files contain expected data

## Next Steps After Successful Installation

- Read `USAGE.md` for detailed usage instructions
- Try with your own protein sequences and mutations
- Review code in `analyzer.py` to understand implementation
- Check `PRESENTATION.md` for demo ideas