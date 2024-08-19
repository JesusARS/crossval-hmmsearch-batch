import re
import pandas as pd


class HmmDataProcessor:
    """Class for the processing of HMMER data"""

    def parse_table(self, hmm_data):
        """
        Convert the output file of an HMMER search into a pandas DataFrame.

        Args:
            hmm_data (str): Path to the input hmm data file.

        Returns:
            pd.DataFrame: DataFrame containing the data from the file.
        """
        with open(hmm_data, 'r', encoding="utf-8") as file:
            lines = file.readlines()

        # Filter out lines that start with '#' and strip leading
        # and trailing spaces
        lines = [
            line.strip() for line in lines if not line.startswith('#') and line.strip()
        ]
        
        # Define column names
        columns = [
            "target_name", "accession", "query_name", "accession2", "E-value_fs",
            "score_fs", "bias_fs", "E-value_bd","score_bd", "bias_bd", "exp", "reg",
            "clu", "ov", "env", "dom", "rep", "inc", "description"
        ]

        # Regular expression pattern to extract data
        pattern = re.compile(
            r'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)'
        )

        data = []
        for line in lines:
            match = pattern.match(line)
            if match:
                data.append(match.groups())
            else:
                print(f"Warning: Line does not match the pattern: {line}")

        # Create DataFrame
        df = pd.DataFrame(data, columns=columns)
        
        return df

    def format_scientific_notation(self, value):
        """
        Formats a value in scientific notation with one decimal place.

        Args:
            value (str): Value to format.

        Returns:
            str: Formatted value or original value if it cannot be converted.
        """
        try:
            return '{:.1e}'.format(float(value))
        except ValueError:
            return value

    def parse_results_data(self, hmm_data, hmm_csv_path):
        """
        Extracts results from an HMM data file and saves them to a CSV file.

        Args:
            hmm_data (str): Path to the HMM data file.
            hmm_csv_path (str): Path to the output HMM CSV file.
        """
        df = self.parse_table(hmm_data)
        df['E-value_fs'] = df['E-value_fs'].apply(self.format_scientific_notation)
        df['E-value_bd'] = df['E-value_bd'].apply(self.format_scientific_notation)
        df.to_csv(hmm_csv_path, index=False)
