import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


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

    def analyze_all_mutations(self):
        """
        Analyzuje všechny mutace načtené z CSV souboru.

        Prochází všechny řádky DataFrame s mutacemi a pro každou
        volá metodu analyze_mutation. Výsledky ukládá do nového DataFrame.

        Returns:
            DataFrame s výsledky analýzy všech mutací
        """
        results = []

        # iterrows() vrací (index, row) tuple pro každý řádek DataFrame
        for idx, row in self.mutations_df.iterrows():
            result = self.analyze_mutation(
                row['position'],
                row['original_aa'],
                row['mutant_aa']
            )
            # Přidáme popis z CSV, pokud existuje
            result['description'] = row.get('description', '')
            results.append(result)

        # Převedeme list slovníků na DataFrame
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


# Tato část se spustí pouze pokud soubor spouštíme přímo
# (ne pokud ho importujeme jako modul)
if __name__ == "__main__":
    # Vytvoříme instanci analyzátoru s naší testovací sekvencí
    analyzer = ProteinMutationAnalyzer('../data/example_protein.fasta')

    # Načteme mutace
    analyzer.load_mutations('../data/mutations_example.csv')

    # Provedeme analýzu
    results = analyzer.analyze_all_mutations()

    print("\n=== ANALYSIS RESULTS ===")
    print(results)

    # Uložíme výsledky
    analyzer.save_results()

    print("\n✓ Analysis complete!")