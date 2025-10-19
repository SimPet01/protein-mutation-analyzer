"""
Protein Mutation Analyzer

This module provides tools for analyzing the effects of point mutations
on protein properties. It calculates changes in molecular weight, charge,
and hydrophobicity, and classifies mutations based on their severity.

Main class:
    ProteinMutationAnalyzer: Analyzes protein mutations and generates reports

Author: Petr Šimandl
Date: 19.10.2025
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import numpy as np

# Settings for better viewing of Pandas DataFrames
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', None)  # Don't limit the width
pd.set_option('display.max_colwidth', 50)  # Max width of one column


class ProteinMutationAnalyzer:
    """
    Class for analyzing protein mutations and predicting their effects.

    This class provides a complete workflow for analyzing point mutations
    in proteins. It allows loading sequences, defining mutations, and calculating
    their effects on physicochemical properties.

    Main functionalities:
    - Loading protein sequences from FASTA format
    - Loading mutations from CSV files
    - Calculating changes in molecular weight, charge, and hydrophobicity
    - Classifying mutations by severity
    - Generating visualizations and HTML reports

    Attributes:
        sequence (SeqRecord): BioPython SeqRecord object with protein sequence
        sequence_str (str): Protein sequence as string
        mutations_df (DataFrame): Pandas DataFrame with loaded mutations
        results_df (DataFrame): Pandas DataFrame with analysis results

    Example:
        >>> analyzer = ProteinMutationAnalyzer('data/protein.fasta')
        >>> analyzer.load_mutations('data/mutations.csv')
        >>> results = analyzer.analyze_all_mutations()
        >>> analyzer.visualize_mutations()
        >>> analyzer.generate_html_report()
    """

    def __init__(self, fasta_file):
        """
        Initialize analyzer with protein sequence from FASTA file.

        Args:
            fasta_file: Path to FASTA file with protein sequence
        """
        self.sequence = self._load_fasta(fasta_file)
        # Store sequence twice: as SeqRecord (with metadata) and as string (for fast access)
        # String form is faster for indexing individual amino acids
        self.sequence_str = str(self.sequence.seq)
        print(f"Loaded sequence: {len(self.sequence_str)} amino acids")
        print(f"Protein ID: {self.sequence.id}")

    def _load_fasta(self, fasta_file):
        """
        Load FASTA file and return first sequence.

        Underscore prefix (_load_fasta) indicates this is a "private" method
        intended for internal use within the class.

        Args:
            fasta_file: Path to FASTA file

        Returns:
            SeqRecord object containing sequence and metadata
        """
        record = next(SeqIO.parse(fasta_file, "fasta"))
        return record

    def load_mutations(self, csv_file):
        """
        Load mutations from CSV file.

        CSV file must contain columns:
        - position: position number in sequence (1-based indexing)
        - original_aa: original amino acid (single letter code)
        - mutant_aa: mutant amino acid (single letter code)
        - description: optional mutation description

        Args:
            csv_file: Path to CSV file with mutations

        Returns:
            DataFrame with loaded mutations
        """
        self.mutations_df = pd.read_csv(csv_file)
        print(f"Loaded {len(self.mutations_df)} mutations")
        return self.mutations_df

    def get_aa_properties(self):
        """
        Return dictionary with physicochemical properties of amino acids.

        Hydrophobicity: Kyte-Doolittle scale
            - Positive values = hydrophobic (water-avoiding)
            - Negative values = hydrophilic (water-attracting)

        Molecular weight: in Daltons (Da)

        Charge: at pH 7.0 (physiological conditions)
            - +1 = positively charged amino acid
            - -1 = negatively charged amino acid
            - 0 = neutral

        Returns:
            dict: Dictionary containing three sub-dictionaries with properties
        """
        hydrophobicity = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }

        molecular_weight = {
            'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
            'Q': 146.1, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
            'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
            'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
        }

        charge = {
            'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
            'Q': 0, 'E': -1, 'G': 0, 'H': 0.1, 'I': 0,
            'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
            'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
        }

        return {
            'hydrophobicity': hydrophobicity,
            'molecular_weight': molecular_weight,
            'charge': charge
        }

    def get_aa_group(self, aa):
        """
        Classify amino acid into chemical group based on properties.

        Categorization based on physicochemical properties:
        - Hydrophobic: Non-polar amino acids with hydrophobic side chains
        - Polar: Polar uncharged amino acids
        - Charged+: Basic (positively charged at pH 7)
        - Charged-: Acidic (negatively charged at pH 7)
        - Special: Glycine (smallest, most flexible)

        Args:
            aa: Single letter amino acid code

        Returns:
            str: Name of the group the amino acid belongs to
        """
        groups = {
            'Hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'],
            'Polar': ['S', 'T', 'C', 'Y', 'N', 'Q'],
            'Charged+': ['K', 'R', 'H'],
            'Charged-': ['D', 'E'],
            'Special': ['G']
        }

        for group, members in groups.items():
            if aa in members:
                return group
        return 'Unknown'

    def classify_mutation(self, mw_change, charge_change, hydro_change):
        """
        Classify mutation as conservative, moderate, or radical.

        Classification is based on combination of property changes:

        Conservative:
            - Small MW change (< 5 Da)
            - No charge change
            - Small hydrophobicity change (< 0.5)
            Example: L→I (both hydrophobic, similar size)

        Radical:
            - Large MW change (> 20 Da)
            - OR charge change
            - OR large hydrophobicity change (> 1.5)
            Example: K→E (change from + to -)

        Moderate:
            - Everything between conservative and radical

        Classification rules are based on empirical biochemical knowledge:
        - Charge changes are almost always significant (affect electrostatic interactions)
        - Large MW changes can disrupt packing in protein core
        - Large hydrophobicity changes can cause misfolding

        Args:
            mw_change: Molecular weight change (Da)
            charge_change: Electric charge change
            hydro_change: Hydrophobicity change

        Returns:
            str: 'Conservative', 'Moderate', or 'Radical'
        """
        # Conservative = all changes are small
        if (abs(mw_change) < 5 and
                abs(charge_change) == 0 and
                abs(hydro_change) < 0.5):
            return 'Conservative'

        # Radical = at least one change is large
        elif (abs(mw_change) > 20 or
              abs(charge_change) != 0 or
              abs(hydro_change) > 1.5):
            return 'Radical'

        # Everything else is moderate
        else:
            return 'Moderate'

    def analyze_mutation(self, position, original_aa, mutant_aa):
        """
        Analyze single mutation and calculate property changes.

        Args:
            position: Position in sequence (1-based indexing, as in biology, 1 as first)
            original_aa: Original amino acid (single letter code)
            mutant_aa: Mutant amino acid (single letter code)

        Returns:
            dict: Dictionary with analysis results containing:
                - position: mutation position
                - original_aa: original amino acid
                - mutant_aa: new amino acid
                - mw_change: molecular weight change (Da)
                - charge_change: charge change
                - hydrophobicity_change: hydrophobicity change
                - valid: True if analysis is valid, otherwise error message
        """
        properties = self.get_aa_properties()

        # Check that position exists in sequence
        if position > len(self.sequence_str):
            return {'error': f'Position {position} out of range (sequence has {len(self.sequence_str)} amino acids)'}

        # IMPORTANT: Biology uses 1-based indexing (first AA is at position 1)
        # Python uses 0-based indexing (first element has index 0)
        # Therefore we subtract 1 everywhere when accessing the sequence
        seq_aa = self.sequence_str[position - 1]
        if seq_aa != original_aa:
            return {'error': f'Position {position} contains {seq_aa}, not {original_aa}'}

        # Calculate property changes (mutant - original)
        # .get() with default value 0 protects against errors if unknown AA code
        mw_change = (properties['molecular_weight'].get(mutant_aa, 0) -
                     properties['molecular_weight'].get(original_aa, 0))

        charge_change = (properties['charge'].get(mutant_aa, 0) -
                         properties['charge'].get(original_aa, 0))

        hydro_change = (properties['hydrophobicity'].get(mutant_aa, 0) -
                        properties['hydrophobicity'].get(original_aa, 0))

        return {
            'position': position,
            'original_aa': original_aa,
            'mutant_aa': mutant_aa,
            'mw_change': round(mw_change, 2),
            'charge_change': charge_change,
            'hydrophobicity_change': round(hydro_change, 2),
            'valid': True
        }

    def analyze_mutation_detailed(self, position, original_aa, mutant_aa):
        """
        Detailed mutation analysis including classification and group information.

        This method extends basic analyze_mutation with:
        - Mutation classification (conservative/moderate/radical)
        - Information about chemical group change
        - Information whether mutation stays in same group

        Args:
            position: Position in sequence
            original_aa: Original amino acid
            mutant_aa: Mutant amino acid

        Returns:
            dict: Extended dictionary with analysis results
        """
        # Call original method for basic analysis
        result = self.analyze_mutation(position, original_aa, mutant_aa)

        # If validation error occurred, return only error
        if 'error' in result:
            return result

        # Add mutation classification
        result['classification'] = self.classify_mutation(
            result['mw_change'],
            result['charge_change'],
            result['hydrophobicity_change']
        )

        # Add information about chemical groups
        orig_group = self.get_aa_group(original_aa)
        mut_group = self.get_aa_group(mutant_aa)
        result['group_change'] = f"{orig_group} → {mut_group}"
        result['same_group'] = orig_group == mut_group

        return result

    def analyze_all_mutations(self):
        """
        Analyze all mutations loaded from CSV file.

        Iterates through all rows of mutations DataFrame and calls
        analyze_mutation_detailed for each. Stores results in new DataFrame.

        Returns:
            DataFrame with results of all mutations analysis
        """
        results = []

        # iterrows() returns (index, row) tuple for each DataFrame row
        for idx, row in self.mutations_df.iterrows():
            # Now using detailed version of analysis
            result = self.analyze_mutation_detailed(
                row['position'],
                row['original_aa'],
                row['mutant_aa']
            )
            result['description'] = row.get('description', '')
            results.append(result)

        # Convert list of dictionaries to DataFrame
        self.results_df = pd.DataFrame(results)
        return self.results_df

    def save_results(self, output_file='../output/mutation_analysis.csv'):
        """
        Save analysis results to CSV file.

        Args:
            output_file: Path to output CSV file
        """
        self.results_df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")

    def visualize_mutations(self, output_file='../output/mutations_visualization.png'):
        """
        Create visualization of mutation effects using graphs.

        Creates 4 graphs in one image (2x2 grid):
        1. Molecular weight changes (bar chart)
        2. Charge changes (bar chart)
        3. Hydrophobicity changes (bar chart)
        4. Classification distribution (pie chart)

        Args:
            output_file: Path to output PNG file
        """
        if self.results_df is None or len(self.results_df) == 0:
            print("No results to visualize. Run analyze_all_mutations() first.")
            return

        # Create figure with 4 subplots (2 rows, 2 columns)
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Protein Mutation Analysis - Effects Overview',
                     fontsize=16, fontweight='bold')

        # Create labels for X axis (shortened mutation names)
        mutation_labels = [f"{row['original_aa']}{row['position']}{row['mutant_aa']}"
                           for _, row in self.results_df.iterrows()]
        x_positions = range(len(self.results_df))

        # ===== GRAPH 1: Molecular Weight Changes =====
        ax1 = axes[0, 0]
        # Color scheme: red = decrease (mass loss), green = increase
        # This helps quickly identify direction of change
        colors1 = ['red' if x < 0 else 'green'
                   for x in self.results_df['mw_change']]
        ax1.bar(x_positions, self.results_df['mw_change'], color=colors1, alpha=0.7)
        ax1.set_xlabel('Mutation', fontsize=11)
        ax1.set_ylabel('MW Change (Da)', fontsize=11)
        ax1.set_title('Molecular Weight Changes', fontsize=12, fontweight='bold')
        ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        ax1.set_xticks(x_positions)
        ax1.set_xticklabels(mutation_labels, rotation=45, ha='right')
        ax1.grid(axis='y', alpha=0.3)

        # ===== GRAPH 2: Charge Changes =====
        ax2 = axes[0, 1]
        # Colors: red for negative change, blue for positive, gray for none
        colors2 = ['red' if x < 0 else 'blue' if x > 0 else 'gray'
                   for x in self.results_df['charge_change']]
        ax2.bar(x_positions, self.results_df['charge_change'], color=colors2, alpha=0.7)
        ax2.set_xlabel('Mutation', fontsize=11)
        ax2.set_ylabel('Charge Change', fontsize=11)
        ax2.set_title('Charge Changes', fontsize=12, fontweight='bold')
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        ax2.set_xticks(x_positions)
        ax2.set_xticklabels(mutation_labels, rotation=45, ha='right')
        ax2.grid(axis='y', alpha=0.3)

        # ===== GRAPH 3: Hydrophobicity Changes =====
        ax3 = axes[1, 0]
        # Colors: blue for more hydrophilic, orange for more hydrophobic
        colors3 = ['blue' if x < 0 else 'orange'
                   for x in self.results_df['hydrophobicity_change']]
        ax3.bar(x_positions, self.results_df['hydrophobicity_change'],
                color=colors3, alpha=0.7)
        ax3.set_xlabel('Mutation', fontsize=11)
        ax3.set_ylabel('Hydrophobicity Change', fontsize=11)
        ax3.set_title('Hydrophobicity Changes', fontsize=12, fontweight='bold')
        ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        ax3.set_xticks(x_positions)
        ax3.set_xticklabels(mutation_labels, rotation=45, ha='right')
        ax3.grid(axis='y', alpha=0.3)

        # ===== GRAPH 4: Classification Distribution =====
        ax4 = axes[1, 1]
        classification_counts = self.results_df['classification'].value_counts()
        # Colors for different classifications
        colors_map = {'Conservative': 'green', 'Moderate': 'orange', 'Radical': 'red'}
        colors4 = [colors_map.get(label, 'gray') for label in classification_counts.index]

        wedges, texts, autotexts = ax4.pie(
            classification_counts.values,
            labels=classification_counts.index,
            autopct='%1.1f%%',
            colors=colors4,
            startangle=90
        )
        # Highlight texts in pie chart
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
        ax4.set_title('Mutation Classification Distribution',
                      fontsize=12, fontweight='bold')

        # Adjust layout so graphs don't overlap
        plt.tight_layout()

        # Save image
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()

    def generate_html_report(self, output_file='../output/mutation_report.html'):
        """
        Generate HTML report with analysis results.

        HTML report contains:
        - Basic information about protein
        - Table with all mutations and their effects
        - Color highlighting based on mutation classification

        Args:
            output_file: Path to output HTML file
        """
        # HTML template with CSS styles
        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="UTF-8">
            <title>Protein Mutation Analysis Report</title>
            <style>
                body {{ 
                    font-family: Arial, sans-serif; 
                    margin: 40px;
                    background-color: #f5f5f5;
                }}
                .container {{
                    background-color: white;
                    padding: 30px;
                    border-radius: 10px;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                    max-width: 1200px;
                    margin: 0 auto;
                }}
                h1 {{ 
                    color: #2c3e50;
                    border-bottom: 3px solid #3498db;
                    padding-bottom: 10px;
                }}
                h2 {{
                    color: #34495e;
                    margin-top: 30px;
                }}
                .info {{
                    background-color: #ecf0f1;
                    padding: 15px;
                    border-radius: 5px;
                    margin: 20px 0;
                }}
                .info p {{
                    margin: 8px 0;
                    font-size: 14px;
                }}
                table {{ 
                    border-collapse: collapse; 
                    width: 100%; 
                    margin-top: 20px;
                    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
                }}
                th, td {{ 
                    border: 1px solid #ddd; 
                    padding: 12px; 
                    text-align: left;
                    font-size: 13px;
                }}
                th {{ 
                    background-color: #3498db; 
                    color: white;
                    font-weight: bold;
                    text-transform: uppercase;
                    font-size: 12px;
                }}
                tr:nth-child(even) {{ 
                    background-color: #f9f9f9; 
                }}
                tr:hover {{
                    background-color: #e8f4f8;
                }}
                .Conservative {{ 
                    background-color: #d5f4e6 !important;
                }}
                .Moderate {{ 
                    background-color: #fff3cd !important;
                }}
                .Radical {{ 
                    background-color: #f8d7da !important;
                }}
                .legend {{
                    margin: 20px 0;
                    padding: 15px;
                    background-color: #f8f9fa;
                    border-radius: 5px;
                }}
                .legend-item {{
                    display: inline-block;
                    margin-right: 20px;
                    padding: 5px 15px;
                    border-radius: 3px;
                    font-size: 13px;
                }}
                .footer {{
                    margin-top: 30px;
                    padding-top: 20px;
                    border-top: 1px solid #ddd;
                    color: #7f8c8d;
                    font-size: 12px;
                    text-align: center;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>Protein Mutation Analysis Report</h1>

                <div class="info">
                    <p><strong>Protein ID:</strong> {protein_id}</p>
                    <p><strong>Sequence Length:</strong> {seq_length} amino acids</p>
                    <p><strong>Total Mutations Analyzed:</strong> {total_mutations}</p>
                    <p><strong>Analysis Date:</strong> {date}</p>
                </div>

                <div class="legend">
                    <strong>Classification Legend:</strong>
                    <span class="legend-item Conservative">Conservative</span>
                    <span class="legend-item Moderate">Moderate</span>
                    <span class="legend-item Radical">Radical</span>
                </div>

                <h2>Mutation Details</h2>
                <table>
                    <tr>
                        <th>Position</th>
                        <th>Mutation</th>
                        <th>MW Change (Da)</th>
                        <th>Charge Change</th>
                        <th>Hydrophobicity Change</th>
                        <th>Classification</th>
                        <th>AA Group Change</th>
                        <th>Description</th>
                    </tr>
                    {table_rows}
                </table>

                <div class="footer">
                    <p>Generated by "Protein Mutation Analyzer"</p>
                    <p>For questions or issues, please contact the Bioinformatics Team</p>
                </div>
            </div>
        </body>
        </html>
        """

        # Create table rows
        table_rows = ""
        for _, row in self.results_df.iterrows():
            mutation_str = f"{row['original_aa']}{row['position']}{row['mutant_aa']}"
            class_style = row.get('classification', 'Unknown')

            # Formate numeric values
            mw_change = f"{row.get('mw_change', 0):+.2f}"  # + displays the +/- sign even for positive numbers
            charge_change = f"{row.get('charge_change', 0):+.1f}"
            hydro_change = f"{row.get('hydrophobicity_change', 0):+.2f}"

            table_rows += f"""
            <tr class="{class_style}">
                <td><strong>{int(row['position'])}</strong></td>
                <td><strong>{mutation_str}</strong></td>
                <td>{mw_change}</td>
                <td>{charge_change}</td>
                <td>{hydro_change}</td>
                <td><strong>{row.get('classification', 'N/A')}</strong></td>
                <td>{row.get('group_change', 'N/A')}</td>
                <td>{row.get('description', '')}</td>
            </tr>
            """

        # Get the current time
        from datetime import datetime
        current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Add data into HTML template
        html_content = html_content.format(
            protein_id=self.sequence.id,
            seq_length=len(self.sequence_str),
            total_mutations=len(self.results_df),
            date=current_date,
            table_rows=table_rows
        )

        # Save the HTML file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"HTML report saved to {output_file}")


# Executed only if run directly (not if imported as a module)
if __name__ == "__main__":
    # Create the instance of Analyzer
    analyzer = ProteinMutationAnalyzer('../data/example_protein.fasta')

    # Load the mutations
    analyzer.load_mutations('../data/mutations_example.csv')

    # Run the analysis
    results = analyzer.analyze_all_mutations()

    print("\n=== ANALYSIS RESULTS ===")
    columns_to_show = ['position', 'original_aa', 'mutant_aa',
                       'mw_change', 'charge_change', 'hydrophobicity_change',
                       'classification', 'group_change']
    print(results[columns_to_show])

    # Save the results
    analyzer.save_results()

    # Create the visualisation
    analyzer.visualize_mutations()

    # Create the HTML report
    analyzer.generate_html_report()

    print("\n✓ Analysis complete!")
    print("- CSV report: output/mutation_analysis.csv")
    print("- Visualization: output/mutations_visualization.png")
    print("- HTML report: output/mutation_report.html")