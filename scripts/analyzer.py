import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import numpy as np

# Nastavení pro lepší zobrazení Pandas DataFrames
pd.set_option('display.max_columns', None)  # Zobraz všechny sloupce
pd.set_option('display.width', None)  # Neomezuj šířku
pd.set_option('display.max_colwidth', 50)  # Max šířka jednoho sloupce


class ProteinMutationAnalyzer:
    """
    Třída pro analýzu proteinových mutací a predikci jejich efektů.

    Tato třída umožňuje načíst proteinovou sekvenci, definovat mutace
    a vypočítat změny ve fyzikálně-chemických vlastnostech proteinu.
    """

    def __init__(self, fasta_file):
        """
        Inicializace analyzátoru s proteinovou sekvencí z FASTA souboru.

        Args:
            fasta_file: Cesta k FASTA souboru s proteinovou sekvencí
        """
        self.sequence = self._load_fasta(fasta_file)
        self.sequence_str = str(self.sequence.seq)
        print(f"Loaded sequence: {len(self.sequence_str)} amino acids")
        print(f"Protein ID: {self.sequence.id}")

    def _load_fasta(self, fasta_file):
        """
        Načte FASTA soubor a vrátí první sekvenci.

        Používáme podtržítko na začátku názvu metody (_load_fasta),
        což v Pythonu signalizuje, že jde o "privátní" metodu,
        která je určená pouze pro interní použití uvnitř třídy.

        Args:
            fasta_file: Cesta k FASTA souboru

        Returns:
            SeqRecord objekt obsahující sekvenci a metadata
        """
        record = next(SeqIO.parse(fasta_file, "fasta"))
        return record

    def load_mutations(self, csv_file):
        """
        Načte mutace z CSV souboru.

        CSV soubor musí obsahovat sloupce:
        - position: číslo pozice v sekvenci (1-based indexing)
        - original_aa: původní aminokyselina (jednopísmenný kód)
        - mutant_aa: mutantní aminokyselina (jednopísmenný kód)
        - description: volitelný popis mutace

        Args:
            csv_file: Cesta k CSV souboru s mutacemi

        Returns:
            DataFrame s načtenými mutacemi
        """
        self.mutations_df = pd.read_csv(csv_file)
        print(f"Loaded {len(self.mutations_df)} mutations")
        return self.mutations_df

    def get_aa_properties(self):
        """
        Vrací slovník s fyzikálně-chemickými vlastnostmi aminokyselin.

        Hydrofobicita: Kyte-Doolittle škála
            - Kladné hodnoty = hydrofobní (odpuzují vodu)
            - Záporné hodnoty = hydrofilní (přitahují vodu)

        Molekulární hmotnost: v daltonech (Da)

        Náboj: při pH 7.0 (fyziologické podmínky)
            - +1 = kladně nabitá aminokyselina
            - -1 = záporně nabitá aminokyselina
            - 0 = neutrální

        Returns:
            dict: Slovník obsahující tři pod-slovníky s vlastnostmi
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

    def classify_mutation(self, mw_change, charge_change, hydro_change):
        """
        Klasifikuje mutaci jako konzervativní, moderátní nebo radikální.

        Klasifikace je založená na kombinaci změn vlastností:

        Conservative (konzervativní):
            - Malá změna hmotnosti (< 5 Da)
            - Bez změny náboje
            - Malá změna hydrofobicity (< 0.5)
            Příklad: L→I (oba hydrofobní, podobná velikost)

        Radical (radikální):
            - Velká změna hmotnosti (> 20 Da)
            - NEBO změna náboje
            - NEBO velká změna hydrofobicity (> 1.5)
            Příklad: K→E (změna z + na -)

        Moderate (moderátní):
            - Vše mezi konzervativní a radikální

        Args:
            mw_change: Změna molekulární hmotnosti (Da)
            charge_change: Změna elektrického náboje
            hydro_change: Změna hydrofobicity

        Returns:
            str: 'Conservative', 'Moderate', nebo 'Radical'
        """
        # Konzervativní = všechny změny jsou malé
        if (abs(mw_change) < 5 and
                abs(charge_change) == 0 and
                abs(hydro_change) < 0.5):
            return 'Conservative'

        # Radikální = alespoň jedna změna je velká
        elif (abs(mw_change) > 20 or
              abs(charge_change) != 0 or
              abs(hydro_change) > 1.5):
            return 'Radical'

        # Všechno ostatní je moderátní
        else:
            return 'Moderate'

    def get_aa_group(self, aa):
        """
        Zařadí aminokyselinu do chemické skupiny podle vlastností.

        Kategorizace je založená na fyzikálně-chemických vlastnostech:
        - Hydrophobic: Nepolární aminokyseliny s hydrofobními postranními řetězci
        - Polar: Polární nenabitá aminokyseliny
        - Charged+: Bazické (kladně nabitá při pH 7)
        - Charged-: Kyselé (záporně nabitá při pH 7)
        - Special: Glycin (nejmenší, nejvíc flexibilní)

        Args:
            aa: Jednopísmenný kód aminokyseliny

        Returns:
            str: Název skupiny, do které aminokyselina patří
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

    def analyze_mutation(self, position, original_aa, mutant_aa):
        """
        Analyzuje jednotlivou mutaci a vypočítá změny vlastností.

        Args:
            position: Pozice v sekvenci (1-based indexing, jako v biologii)
            original_aa: Původní aminokyselina (jednopísmenný kód)
            mutant_aa: Mutantní aminokyselina (jednopísmenný kód)

        Returns:
            dict: Slovník s výsledky analýzy obsahující:
                - position: pozice mutace
                - original_aa: původní aminokyselina
                - mutant_aa: nová aminokyselina
                - mw_change: změna molekulární hmotnosti (Da)
                - charge_change: změna náboje
                - hydrophobicity_change: změna hydrofobicity
                - valid: True pokud je analýza validní, jinak error zpráva
        """
        properties = self.get_aa_properties()

        # Kontrola, že pozice existuje v sekvenci
        # V biologii se počítá od 1, v Pythonu od 0, proto -1
        if position > len(self.sequence_str):
            return {'error': f'Position {position} out of range (sequence has {len(self.sequence_str)} amino acids)'}

        # Ověření, že na dané pozici je skutečně původní aminokyselina
        # Toto je důležitá validace - pokud někdo zadá špatnou pozici nebo AA
        seq_aa = self.sequence_str[position - 1]  # -1 protože Python indexuje od 0
        if seq_aa != original_aa:
            return {'error': f'Position {position} contains {seq_aa}, not {original_aa}'}

        # Výpočet změn jednotlivých vlastností
        # Používáme .get() místo přímého přístupu [key], protože .get() vrací
        # defaultní hodnotu (0) pokud klíč neexistuje, místo vyvolání chyby
        mw_change = (properties['molecular_weight'].get(mutant_aa, 0) -
                     properties['molecular_weight'].get(original_aa, 0))

        charge_change = (properties['charge'].get(mutant_aa, 0) -
                         properties['charge'].get(original_aa, 0))

        hydro_change = (properties['hydrophobicity'].get(mutant_aa, 0) -
                        properties['hydrophobicity'].get(original_aa, 0))

        # Vrátíme slovník s výsledky
        # round() zaokrouhluje na 2 desetinná místa pro lepší čitelnost
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
        Detailní analýza mutace včetně klasifikace a informací o skupinách.

        Tato metoda rozšiřuje základní analyze_mutation o:
        - Klasifikaci mutace (conservative/moderate/radical)
        - Informaci o změně chemické skupiny
        - Informaci, jestli mutace zůstává ve stejné skupině

        Args:
            position: Pozice v sekvenci
            original_aa: Původní aminokyselina
            mutant_aa: Mutantní aminokyselina

        Returns:
            dict: Rozšířený slovník s výsledky analýzy
        """
        # Zavoláme původní metodu pro základní analýzu
        result = self.analyze_mutation(position, original_aa, mutant_aa)

        # Pokud byla chyba při validaci, vrátíme jen error
        if 'error' in result:
            return result

        # Přidáme klasifikaci mutace
        result['classification'] = self.classify_mutation(
            result['mw_change'],
            result['charge_change'],
            result['hydrophobicity_change']
        )

        # Přidáme informace o chemických skupinách
        orig_group = self.get_aa_group(original_aa)
        mut_group = self.get_aa_group(mutant_aa)
        result['group_change'] = f"{orig_group} → {mut_group}"
        result['same_group'] = orig_group == mut_group

        return result

    def analyze_all_mutations(self):
        """
        Analyzuje všechny mutace načtené z CSV souboru.

        Používá rozšířenou metodu analyze_mutation_detailed pro získání
        kompletních informací o každé mutaci včetně klasifikace.

        Returns:
            DataFrame s detailními výsledky analýzy všech mutací
        """
        results = []

        for idx, row in self.mutations_df.iterrows():
            # Teď používáme detailní verzi analýzy
            result = self.analyze_mutation_detailed(
                row['position'],
                row['original_aa'],
                row['mutant_aa']
            )
            result['description'] = row.get('description', '')
            results.append(result)

        self.results_df = pd.DataFrame(results)
        return self.results_df

    def save_results(self, output_file='../output/mutation_analysis.csv'):
        """
        Uloží výsledky analýzy do CSV souboru.

        Args:
            output_file: Cesta k výstupnímu CSV souboru
        """
        self.results_df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")

    def visualize_mutations(self, output_file='../output/mutations_visualization.png'):
        """
        Vytvoří vizualizaci efektů mutací pomocí grafů.

        Vytváří 4 grafy v jednom obrázku (2x2 grid):
        1. Změny molekulární hmotnosti (sloupcový graf)
        2. Změny náboje (sloupcový graf)
        3. Změny hydrofobicity (sloupcový graf)
        4. Distribuce klasifikací (koláčový graf)

        Args:
            output_file: Cesta k výstupnímu PNG souboru
        """
        if self.results_df is None or len(self.results_df) == 0:
            print("No results to visualize. Run analyze_all_mutations() first.")
            return

        # Vytvoříme figuru s 4 subploty (2 řádky, 2 sloupce)
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Protein Mutation Analysis - Effects Overview',
                     fontsize=16, fontweight='bold')

        # Vytvoříme popisky pro osu X (zkrácené názvy mutací)
        mutation_labels = [f"{row['original_aa']}{row['position']}{row['mutant_aa']}"
                           for _, row in self.results_df.iterrows()]
        x_positions = range(len(self.results_df))

        # ===== GRAF 1: Změny molekulární hmotnosti =====
        ax1 = axes[0, 0]
        # Barvy: červená pro pokles, zelená pro nárůst
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

        # ===== GRAF 2: Změny náboje =====
        ax2 = axes[0, 1]
        # Barvy: červená pro zápornou změnu, modrá pro kladnou, šedá pro žádnou
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

        # ===== GRAF 3: Změny hydrofobicity =====
        ax3 = axes[1, 0]
        # Barvy: modrá pro hydrofilnější, oranžová pro hydrofóbnější
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

        # ===== GRAF 4: Distribuce klasifikací =====
        ax4 = axes[1, 1]
        classification_counts = self.results_df['classification'].value_counts()
        # Barvy pro různé klasifikace
        colors_map = {'Conservative': 'green', 'Moderate': 'orange', 'Radical': 'red'}
        colors4 = [colors_map.get(label, 'gray') for label in classification_counts.index]

        wedges, texts, autotexts = ax4.pie(
            classification_counts.values,
            labels=classification_counts.index,
            autopct='%1.1f%%',
            colors=colors4,
            startangle=90
        )
        # Zvýraznění textů v koláčovém grafu
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
        ax4.set_title('Mutation Classification Distribution',
                      fontsize=12, fontweight='bold')

        # Upravíme rozložení, aby se grafy nepřekrývaly
        plt.tight_layout()

        # Uložíme obrázek
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()

    def generate_html_report(self, output_file='../output/mutation_report.html'):
        """
        Vytvoří HTML report s výsledky analýzy.

        HTML report obsahuje:
        - Základní informace o proteinu
        - Tabulku se všemi mutacemi a jejich efekty
        - Barevné zvýraznění podle klasifikace mutací

        Args:
            output_file: Cesta k výstupnímu HTML souboru
        """
        # HTML šablona s CSS styly
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

        # Vytvoříme řádky tabulky
        table_rows = ""
        for _, row in self.results_df.iterrows():
            mutation_str = f"{row['original_aa']}{row['position']}{row['mutant_aa']}"
            class_style = row.get('classification', 'Unknown')

            # Formátujeme číselné hodnoty
            mw_change = f"{row.get('mw_change', 0):+.2f}"  # + zobrazí znaménko i pro kladná čísla
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

        # Získáme aktuální datum
        from datetime import datetime
        current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Doplníme data do HTML šablony
        html_content = html_content.format(
            protein_id=self.sequence.id,
            seq_length=len(self.sequence_str),
            total_mutations=len(self.results_df),
            date=current_date,
            table_rows=table_rows
        )

        # Uložíme HTML soubor
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"HTML report saved to {output_file}")


# Tato část se spustí pouze pokud soubor spouštíme přímo
# (ne pokud ho importujeme jako modul)
if __name__ == "__main__":
    # Vytvoříme instanci analyzátoru
    analyzer = ProteinMutationAnalyzer('../data/example_protein.fasta')

    # Načteme mutace
    analyzer.load_mutations('../data/mutations_example.csv')

    # Provedeme analýzu
    results = analyzer.analyze_all_mutations()

    print("\n=== ANALYSIS RESULTS ===")
    columns_to_show = ['position', 'original_aa', 'mutant_aa',
                       'mw_change', 'charge_change', 'hydrophobicity_change',
                       'classification', 'group_change']
    print(results[columns_to_show])

    # Uložíme výsledky
    analyzer.save_results()

    # Vytvoříme vizualizaci
    analyzer.visualize_mutations()

    # NOVĚ: Vytvoříme HTML report
    analyzer.generate_html_report()

    print("\n✓ Analysis complete!")
    print("- CSV report: output/mutation_analysis.csv")
    print("- Visualization: output/mutations_visualization.png")
    print("- HTML report: output/mutation_report.html")