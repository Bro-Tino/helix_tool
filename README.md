#  HelixLab - DNA Sequence Analysis Tool

A comprehensive, interactive DNA sequence analysis tool built with Streamlit. HelixLab provides professional-grade bioinformatics analysis with an intuitive interface and beautiful visualizations.

![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/streamlit-1.28+-red.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

##  Features

###  Nucleotide Analysis
- **Base Composition**: Detailed counting of A, T, C, G nucleotides
- **GC Content Calculation**: Percentage of guanine and cytosine bases
- **AT/GC Ratios**: Comprehensive ratio analysis
- **Interactive Charts**: 
  - Bar chart showing nucleotide distribution
  - Pie chart for GC vs AT content visualization

###  Codon Analysis
- **Start/Stop Codon Detection**: Identifies ATG (start) and TAA, TAG, TGA (stop) codons
- **Codon Usage Statistics**: Frequency analysis of all codons
- **Full Amino Acid Names**: Displays both single-letter codes and full names (e.g., M → Methionine)
- **Interactive Visualization**: Top 20 most frequent codons chart
- **Comprehensive Table**: Sortable table with codon, amino acid, count, and frequency

###  Protein Translation & ORF Detection
- **Multi-frame Translation**: Translation in all three reading frames
- **Full Amino Acid Display**: Shows both abbreviated and full amino acid names
- **Open Reading Frame (ORF) Detection**: 
  - Configurable minimum length (30-300 bp)
  - Forward and reverse strand analysis
  - All six reading frames
  - Position and length information
- **ORF Visualization**: Scatter plot showing ORF distribution across the sequence

###  DNA Conversions
- **DNA to RNA**: Converts thymine (T) to uracil (U)
- **Reverse Complement**: Generates the reverse complement strand
- **Download Options**: Export conversions as text files

### Summary Report
- Complete sequence statistics
- Molecular weight estimation
- Downloadable text report

### User Interface
- **Dark/Light Mode Toggle**: Switch between themes for comfortable viewing
- **Responsive Design**: Works on desktop and mobile devices
- **Interactive Visualizations**: All charts are powered by Plotly for zoom, pan, and hover interactions

## Supported File Formats

### FASTA Variants
- `.fasta` - Standard FASTA format
- `.fa` - FASTA abbreviation
- `.fna` - FASTA nucleic acid
- `.ffn` - FASTA nucleotide coding regions
- `.faa` - FASTA amino acid
- `.frn` - FASTA non-coding RNA

### Other Formats
- `.gb`, `.gbk` - GenBank format
- `.txt` - Plain text sequences

## Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package manager)

### Step 1: Clone or Download
```bash
# If using git
git clone https://github.com/Bro-Tino/helixlab.git
cd helixlab

# Or download and extract the ZIP file
```

### Step 2: Install Dependencies
```bash
pip install -r requirements.txt
```

### Step 3: Run the Application
```bash
streamlit run helix_lab_enhanced.py
```

The application will open in your default web browser at `http://localhost:8501`

##  Usage Guide

### Method 1: Direct Input
1. Select **"Direct Input"** from the sidebar
2. Paste your DNA sequence in the text area
3. Click outside the text area to load the sequence
4. Explore the analysis tabs

### Method 2: File Upload
1. Select **"File Upload"** from the sidebar
2. Click **"Browse files"** and select your sequence file
3. Supported formats: FASTA, GenBank, or plain text
4. The sequence will load automatically
5. If multiple sequences exist in the file, select one from the dropdown

### Navigation
Use the tabs to access different analyses:
- ** Nucleotide Analysis** - Base composition and GC content
- ** Codon Analysis** - Codon usage and statistics
- ** Translation & ORF** - Protein translation and ORF detection
- ** DNA Conversions** - RNA and reverse complement
- ** Summary Report** - Comprehensive overview with download option

### Tips
- Use the **theme toggle** (☀️/🌙) in the top-right to switch between dark and light modes
- Adjust the **minimum ORF length** slider to filter ORFs by size
- **Download** converted sequences and reports using the download buttons
- **Expand** sequence view to see the full DNA sequence

##  Example Input

### Direct Input Example
```
ATGGCTAGCTAGCTAGCTAGTAGATGAAACCCGGGTTTAAATAG
```

### FASTA File Example
```
>Sample_Gene_1
ATGGCTAGCTAGCTAGCTAGTAG
ATGAAACCCGGGTTTAAATAG
>Sample_Gene_2
ATGCGATCGATCGATCGTAGCTAGCTAGCTAG
TAATAG
```

### GenBank File Example
```
LOCUS       SAMPLE_SEQ               45 bp    DNA     linear   01-JAN-2025
ORIGIN
        1 atggctagct agctagctag tagatgaaac ccgggtttaa atag
//
```

## 🧬 Scientific Background

### Genetic Code
HelixLab uses the standard genetic code for translation:
- **Start Codon**: ATG (Methionine)
- **Stop Codons**: TAA, TAG, TGA

### ORF Detection
Open Reading Frames are detected by:
1. Searching for start codons (ATG)
2. Translating in frame until a stop codon
3. Filtering by minimum length threshold
4. Analyzing both forward and reverse strands
5. Checking all three reading frames per strand

### GC Content
GC content is an important indicator of:
- DNA stability (higher GC = more stable)
- Gene density (gene-rich regions often have higher GC)
- Organism classification

##  Technical Details

### Dependencies
- **Streamlit**: Web application framework
- **Plotly**: Interactive visualization library
- **Pandas**: Data manipulation and analysis

### Architecture
- **Single-file application**: All code in one Python file for easy deployment
- **Session state management**: Theme preferences persist during session
- **Modular functions**: Clean, reusable code structure

##  Output Examples

### Nucleotide Metrics
```
Adenine (A): 150 (25.00%)
Thymine (T): 150 (25.00%)
Cytosine (C): 150 (25.00%)
Guanine (G): 150 (25.00%)
GC Content: 50.00%
```

### Protein Translation
```
Frame 1:
Single letter: MASLLLV
Full names: Methionine - Alanine - Serine - Leucine - Leucine - Leucine - Valine
```

### ORF Detection
```
Found 3 ORF(s)
ORF 1: Forward Strand, Frame 1 (150 bp)
Position: 1 - 150
Protein Length: 49 amino acids
```

##  Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Areas for Contribution
- Additional sequence analysis features
- Support for more file formats
- Performance optimizations
- UI/UX improvements
- Bug fixes

##  License

This project is licensed under the MIT License - see the LICENSE file for details.

##  Author

**T. Semwa**

##  Acknowledgments

- Built with [Streamlit](https://streamlit.io/)
- Visualizations powered by [Plotly](https://plotly.com/)
- Genetic code based on NCBI standards

##  Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Contact: [tinosemwa@gmail.com]

## 🔄 Version History

### v2.0 (Current)
- ✅ Added dark/light mode toggle
- ✅ Full amino acid names display
- ✅ Interactive Plotly visualizations
- ✅ Extended file format support (FASTA variants, GenBank)
- ✅ ORF distribution visualization
- ✅ Download options for conversions and reports

### v1.0
- Initial release
- Basic nucleotide and codon analysis
- ORF detection
- DNA-RNA conversion
- FASTA file support

---

Made with love and  by T. Semwa | HelixLab © 2026
