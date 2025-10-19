"""
Quick Demo Script for Protein Mutation Analyzer

This script provides a quick demonstration of the analyzer's capabilities.
Run this to see example output without needing to understand the code.

Usage:
    python quick_demo.py
"""

from analyzer import ProteinMutationAnalyzer
import os


def print_section(title):
    """Print formatted section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70 + "\n")


def main():
    """Run quick demonstration of the analyzer."""

    print_section("PROTEIN MUTATION ANALYZER - QUICK DEMO")

    print("This demo will:")
    print("  1. Load a protein sequence (Human Myoglobin)")
    print("  2. Analyze 5 example mutations")
    print("  3. Generate visualizations and reports")
    print("\nStarting analysis...\n")

    # Step 1: Load protein
    print_section("STEP 1: Loading Protein Sequence")
    try:
        analyzer = ProteinMutationAnalyzer('../data/example_protein.fasta')
        print("✓ Protein loaded successfully")
    except Exception as e:
        print(f"✗ Error loading protein: {e}")
        return

    # Step 2: Load mutations
    print_section("STEP 2: Loading Mutations")
    try:
        analyzer.load_mutations('../data/mutations_example.csv')
        print("✓ Mutations loaded successfully")
    except Exception as e:
        print(f"✗ Error loading mutations: {e}")
        return

    # Step 3: Analyze
    print_section("STEP 3: Analyzing Mutations")
    try:
        results = analyzer.analyze_all_mutations()
        print("✓ Analysis completed successfully")
        print("\nSample results:")
        print(results[['position', 'original_aa', 'mutant_aa',
                       'classification', 'mw_change', 'charge_change']].to_string())
    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        return

    # Step 4: Generate outputs
    print_section("STEP 4: Generating Output Files")

    try:
        # CSV
        analyzer.save_results()
        print("✓ CSV report generated")

        # Visualization
        analyzer.visualize_mutations()
        print("✓ Visualization created")

        # HTML report
        analyzer.generate_html_report()
        print("✓ HTML report generated")

    except Exception as e:
        print(f"✗ Error generating outputs: {e}")
        return

    # Summary
    print_section("DEMO COMPLETE")

    print("Output files created in '../output/' folder:")
    print("  • mutation_analysis.csv      - Detailed results table")
    print("  • mutations_visualization.png - Four-panel graph")
    print("  • mutation_report.html       - Interactive HTML report")

    print("\n" + "=" * 70)
    print("Summary of analyzed mutations:")
    print("=" * 70)

    # Classification summary
    class_counts = results['classification'].value_counts()
    for classification, count in class_counts.items():
        emoji = "🔴" if classification == "Radical" else "🟡" if classification == "Moderate" else "🟢"
        print(f"  {emoji} {classification}: {count} mutation(s)")

    print("\n" + "=" * 70)
    print("Next steps:")
    print("  • Open mutation_report.html in your browser")
    print("  • View mutations_visualization.png")
    print("  • Read USAGE.md for more details")
    print("  • Try with your own protein and mutations!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()