from analyzer import ProteinMutationAnalyzer


def test_with_different_data():
    """
    Testovací funkce pro ověření, že analyzátor funguje s různými daty.
    """
    print("=" * 60)
    print("TESTING ANALYZER WITH DIFFERENT DATA")
    print("=" * 60)

    try:
        # Načteme testovací data
        analyzer = ProteinMutationAnalyzer('../data/test_protein.fasta')
        analyzer.load_mutations('../data/test_mutations.csv')

        # Provedeme analýzu
        results = analyzer.analyze_all_mutations()

        # Zobrazíme výsledky
        print("\n=== TEST RESULTS ===")
        print(results[['position', 'original_aa', 'mutant_aa',
                       'classification', 'group_change']])

        # Uložíme výsledky s jiným názvem, aby se nepřepsaly původní
        analyzer.save_results('../output/test_mutation_analysis.csv')
        analyzer.visualize_mutations('../output/test_mutations_visualization.png')
        analyzer.generate_html_report('../output/test_mutation_report.html')

        print("\n✓ TEST PASSED - All outputs generated successfully!")
        print("- CSV: output/test_mutation_analysis.csv")
        print("- Visualization: output/test_mutations_visualization.png")
        print("- HTML: output/test_mutation_report.html")

        return True

    except Exception as e:
        print(f"\n✗ TEST FAILED - Error occurred: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_with_different_data()
    if success:
        print("\n" + "=" * 60)
        print("ALL TESTS COMPLETED SUCCESSFULLY")
        print("=" * 60)
    else:
        print("\n" + "=" * 60)
        print("TESTS FAILED - Please check errors above")
        print("=" * 60)