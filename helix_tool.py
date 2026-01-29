import streamlit as st
import re
from collections import Counter
from io import StringIO
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(page_title="HelixLab", page_icon="🧬", layout="wide")

# Genetic code dictionary for translation
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

# Full amino acid names
AMINO_ACID_NAMES = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
    'C': 'Cysteine', 'Q': 'Glutamine', 'E': 'Glutamic acid', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
    '*': 'Stop', 'X': 'Unknown'
}

START_CODONS = ['ATG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

# Theme configuration
if 'theme' not in st.session_state:
    st.session_state.theme = 'dark'

def toggle_theme():
    st.session_state.theme = 'light' if st.session_state.theme == 'dark' else 'dark'

def get_theme_colors():
    if st.session_state.theme == 'dark':
        return {
            'bg': '#0E1117',
            'secondary_bg': '#262730',
            'text': '#FAFAFA',
            'plot_template': 'plotly_dark'
        }
    else:
        return {
            'bg': '#FFFFFF',
            'secondary_bg': '#F0F2F6',
            'text': '#31333F',
            'plot_template': 'plotly_white'
        }

def parse_fasta(fasta_text):
    """Parse FASTA format files (.fasta, .fa, .fna, .ffn, .faa, .frn)"""
    sequences = {}
    current_header = None
    current_seq = []
    
    for line in fasta_text.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header:
                sequences[current_header] = ''.join(current_seq).upper()
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_header:
        sequences[current_header] = ''.join(current_seq).upper()
    
    return sequences

def parse_genbank(gb_text):
    """Basic GenBank format parser"""
    sequences = {}
    current_seq = []
    locus_name = "GenBank_Sequence"
    in_origin = False
    
    for line in gb_text.split('\n'):
        if line.startswith('LOCUS'):
            locus_name = line.split()[1]
        elif line.startswith('ORIGIN'):
            in_origin = True
        elif line.startswith('//'):
            in_origin = False
            if current_seq:
                sequences[locus_name] = ''.join(current_seq).upper()
                current_seq = []
        elif in_origin:
            # Remove numbers and spaces from sequence lines
            seq_part = re.sub(r'[^a-zA-Z]', '', line)
            current_seq.append(seq_part)
    
    return sequences

def validate_dna(seq):
    return bool(re.match('^[ATCG]+$', seq.upper()))

def count_nucleotides(seq):
    counter = Counter(seq.upper())
    return {
        'A': counter.get('A', 0),
        'T': counter.get('T', 0),
        'C': counter.get('C', 0),
        'G': counter.get('G', 0)
    }

def calculate_gc_content(seq):
    counts = count_nucleotides(seq)
    total = sum(counts.values())
    if total == 0:
        return 0
    gc = counts['G'] + counts['C']
    return (gc / total) * 100

def calculate_ratios(seq):
    counts = count_nucleotides(seq)
    at = counts['A'] + counts['T']
    gc = counts['G'] + counts['C']
    at_gc_ratio = at / gc if gc > 0 else float('inf')
    gc_at_ratio = gc / at if at > 0 else float('inf')
    return at_gc_ratio, gc_at_ratio

def dna_to_rna(seq):
    return seq.upper().replace('T', 'U')

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in seq.upper()[::-1]])

def translate_sequence(seq, frame=0):
    seq = seq.upper()[frame:]
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)

def find_orfs(seq, min_length=30):
    orfs = []
    seq = seq.upper()
    
    for frame in range(3):
        for strand, sequence in [('Forward', seq), ('Reverse', reverse_complement(seq))]:
            i = frame
            while i < len(sequence) - 2:
                codon = sequence[i:i+3]
                if codon in START_CODONS:
                    start = i
                    for j in range(i, len(sequence) - 2, 3):
                        stop_codon = sequence[j:j+3]
                        if stop_codon in STOP_CODONS:
                            length = j - start + 3
                            if length >= min_length:
                                orf_seq = sequence[start:j+3]
                                protein = translate_sequence(orf_seq)
                                orfs.append({
                                    'strand': strand,
                                    'frame': frame + 1,
                                    'start': start + 1,
                                    'end': j + 3,
                                    'length': length,
                                    'sequence': orf_seq,
                                    'protein': protein
                                })
                            i = j + 3
                            break
                    else:
                        i += 3
                else:
                    i += 3
    
    return orfs

def analyze_codons(seq):
    seq = seq.upper()
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    codon_counts = Counter(codons)
    start_count = sum(codon_counts[c] for c in START_CODONS)
    stop_count = sum(codon_counts[c] for c in STOP_CODONS)
    return codon_counts, start_count, stop_count

def create_nucleotide_chart(counts, theme):
    df = pd.DataFrame(list(counts.items()), columns=['Nucleotide', 'Count'])
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A']
    
    fig = px.bar(df, x='Nucleotide', y='Count', 
                 title='Nucleotide Distribution',
                 color='Nucleotide',
                 color_discrete_sequence=colors)
    fig.update_layout(
        template=theme['plot_template'],
        showlegend=False,
        height=400
    )
    return fig

def create_gc_pie_chart(counts, theme):
    gc = counts['G'] + counts['C']
    at = counts['A'] + counts['T']
    
    fig = go.Figure(data=[go.Pie(
        labels=['GC', 'AT'],
        values=[gc, at],
        hole=.3,
        marker_colors=['#4ECDC4', '#FF6B6B']
    )])
    fig.update_layout(
        title='GC vs AT Content',
        template=theme['plot_template'],
        height=400
    )
    return fig

def create_codon_usage_chart(codon_counts, theme, top_n=20):
    sorted_codons = sorted(codon_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    codons, counts = zip(*sorted_codons) if sorted_codons else ([], [])
    
    fig = px.bar(x=list(codons), y=list(counts),
                 title=f'Top {top_n} Codon Usage',
                 labels={'x': 'Codon', 'y': 'Count'})
    fig.update_layout(
        template=theme['plot_template'],
        height=400,
        xaxis_tickangle=-45
    )
    return fig

def create_orf_distribution(orfs, theme):
    if not orfs:
        return None
    
    df = pd.DataFrame(orfs)
    fig = px.scatter(df, x='start', y='length', color='strand',
                     title='ORF Distribution',
                     labels={'start': 'Position (bp)', 'length': 'Length (bp)'},
                     hover_data=['frame', 'end'])
    fig.update_layout(
        template=theme['plot_template'],
        height=400
    )
    return fig

# Main UI
col1, col2 = st.columns([6, 1])
with col1:
    st.title("🧬 Helix_Tool")
    st.subheader("DNA Sequence Analysis Tool")
with col2:
    st.write("")
    st.write("")
    theme_icon = "☀️" if st.session_state.theme == 'dark' else "🌙"
    if st.button(theme_icon, key="theme_toggle", help="Toggle theme"):
        toggle_theme()
        st.rerun()

theme = get_theme_colors()

# Sidebar for input
with st.sidebar:
    st.header("📥 Input Sequence")
    input_method = st.radio("Input Method", ["Direct Input", "File Upload"])
    
    if input_method == "Direct Input":
        sequence_input = st.text_area("Enter DNA Sequence", height=200, 
                                     placeholder="Enter DNA sequence (A, T, C, G)")
        sequences = {"User Input": sequence_input.strip().upper()} if sequence_input.strip() else {}
    else:
        uploaded_file = st.file_uploader(
            "Upload Sequence File", 
            type=['fasta', 'fa', 'fna', 'ffn', 'faa', 'frn', 'txt', 'gb', 'gbk', 'genbank']
        )
        
        if uploaded_file:
            file_content = StringIO(uploaded_file.getvalue().decode("utf-8")).read()
            file_extension = uploaded_file.name.split('.')[-1].lower()
            
            # Parse based on file type
            if file_extension in ['gb', 'gbk', 'genbank']:
                sequences = parse_genbank(file_content)
            else:  # FASTA and variants
                sequences = parse_fasta(file_content)
        else:
            sequences = {}

# Main content
if sequences:
    seq_name = st.selectbox("Select Sequence", list(sequences.keys()))
    sequence = sequences[seq_name]
    
    if not validate_dna(sequence):
        st.error("❌ Invalid DNA sequence! Please ensure the sequence contains only A, T, C, G characters.")
    else:
        st.success(f"✓ Valid DNA sequence loaded: {len(sequence):,} nucleotides")
        
        # Display sequence
        with st.expander("📄 View Sequence", expanded=False):
            st.code(sequence, language=None)
        
        # Analysis tabs
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            " Nucleotide Analysis",
            " Codon Analysis",
            " Translation & ORF",
            " DNA Conversions",
            " Summary Report"
        ])
        
        with tab1:
            st.header(" Nucleotide Analysis")
            
            counts = count_nucleotides(sequence)
            
            # Metrics
            col1, col2, col3, col4 = st.columns(4)
            total = len(sequence)
            with col1:
                st.metric("Adenine (A)", f"{counts['A']:,}", f"{(counts['A']/total*100):.2f}%")
            with col2:
                st.metric("Thymine (T)", f"{counts['T']:,}", f"{(counts['T']/total*100):.2f}%")
            with col3:
                st.metric("Cytosine (C)", f"{counts['C']:,}", f"{(counts['C']/total*100):.2f}%")
            with col4:
                st.metric("Guanine (G)", f"{counts['G']:,}", f"{(counts['G']/total*100):.2f}%")
            
            # Charts
            col1, col2 = st.columns(2)
            with col1:
                st.plotly_chart(create_nucleotide_chart(counts, theme), use_container_width=True)
            with col2:
                st.plotly_chart(create_gc_pie_chart(counts, theme), use_container_width=True)
            
            # GC Content & Ratios
            col1, col2, col3 = st.columns(3)
            gc_content = calculate_gc_content(sequence)
            at_gc, gc_at = calculate_ratios(sequence)
            
            with col1:
                st.metric("GC Content", f"{gc_content:.2f}%")
            with col2:
                st.metric("AT/GC Ratio", f"{at_gc:.3f}")
            with col3:
                st.metric("GC/AT Ratio", f"{gc_at:.3f}")
        
        with tab2:
            st.header(" Codon Analysis")
            
            codon_counts, start_count, stop_count = analyze_codons(sequence)
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Start Codons (ATG)", start_count)
            with col2:
                st.metric("Stop Codons", stop_count)
            with col3:
                st.metric("Total Codons", sum(codon_counts.values()))
            
            # Codon usage chart
            st.plotly_chart(create_codon_usage_chart(codon_counts, theme), use_container_width=True)
            
            # Codon usage table
            st.subheader("Codon Usage Table")
            if codon_counts:
                sorted_codons = sorted(codon_counts.items(), key=lambda x: x[1], reverse=True)
                codon_data = []
                for codon, count in sorted_codons:
                    aa = CODON_TABLE.get(codon, 'X')
                    aa_name = AMINO_ACID_NAMES.get(aa, 'Unknown')
                    percentage = (count / sum(codon_counts.values())) * 100
                    codon_type = "Start" if codon in START_CODONS else ("Stop" if codon in STOP_CODONS else "")
                    codon_data.append({
                        'Codon': codon,
                        'Amino Acid': aa,
                        'Full Name': aa_name,
                        'Count': count,
                        'Frequency (%)': f"{percentage:.2f}",
                        'Type': codon_type
                    })
                st.dataframe(codon_data, use_container_width=True)
        
        with tab3:
            st.header(" Protein Translation & ORF Detection")
            
            st.subheader("Standard Translation (All Frames)")
            for frame in range(3):
                protein = translate_sequence(sequence, frame)
                if protein:
                    # Convert to full names
                    full_names = [AMINO_ACID_NAMES.get(aa, 'Unknown') for aa in protein]
                    st.write(f"**Frame {frame + 1}:**")
                    st.write(f"*Single letter:* `{protein}`")
                    st.write(f"*Full names:* {' - '.join(full_names[:10])}{'...' if len(full_names) > 10 else ''}")
                else:
                    st.write(f"**Frame {frame + 1}:** No translation")
            
            st.subheader(" Open Reading Frame (ORF) Detection")
            min_orf_length = st.slider("Minimum ORF Length (bp)", 30, 300, 100, 10)
            
            orfs = find_orfs(sequence, min_orf_length)
            
            if orfs:
                st.write(f"Found {len(orfs)} ORF(s)")
                
                # ORF distribution chart
                orf_fig = create_orf_distribution(orfs, theme)
                if orf_fig:
                    st.plotly_chart(orf_fig, use_container_width=True)
                
                for idx, orf in enumerate(orfs, 1):
                    with st.expander(f"ORF {idx}: {orf['strand']} Strand, Frame {orf['frame']} ({orf['length']} bp)"):
                        st.write(f"**Position:** {orf['start']:,} - {orf['end']:,}")
                        st.write(f"**DNA Sequence:**")
                        st.code(orf['sequence'], language=None)
                        st.write(f"**Protein (Single Letter):** `{orf['protein']}`")
                        
                        # Show full names
                        full_names = [AMINO_ACID_NAMES.get(aa, 'Unknown') for aa in orf['protein']]
                        st.write(f"**Protein (Full Names):**")
                        st.write(' → '.join(full_names))
                        st.write(f"**Protein Length:** {len(orf['protein'])} amino acids")
            else:
                st.info(f"No ORFs found with minimum length {min_orf_length} bp")
        
        with tab4:
            st.header(" DNA Conversions")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("DNA → RNA")
                rna = dna_to_rna(sequence)
                st.code(rna, language=None)
                st.download_button(
                    label="Download RNA",
                    data=rna,
                    file_name=f"{seq_name}_RNA.txt",
                    mime="text/plain"
                )
            
            with col2:
                st.subheader("Reverse Complement")
                rev_comp = reverse_complement(sequence)
                st.code(rev_comp, language=None)
                st.download_button(
                    label="Download Reverse Complement",
                    data=rev_comp,
                    file_name=f"{seq_name}_RevComp.txt",
                    mime="text/plain"
                )
        
        with tab5:
            st.header(" Summary Report")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Sequence Information")
                st.write(f"**Name:** {seq_name}")
                st.write(f"**Length:** {len(sequence):,} nucleotides")
                st.write(f"**Molecular Weight:** ~{len(sequence) * 650:,} Da")
                
                st.subheader("Base Composition")
                counts = count_nucleotides(sequence)
                for base, count in counts.items():
                    percentage = (count / len(sequence)) * 100
                    st.write(f"**{base}:** {count:,} ({percentage:.2f}%)")
            
            with col2:
                st.subheader("Sequence Properties")
                gc_content = calculate_gc_content(sequence)
                st.write(f"**GC Content:** {gc_content:.2f}%")
                
                at_gc, gc_at = calculate_ratios(sequence)
                st.write(f"**AT/GC Ratio:** {at_gc:.3f}")
                st.write(f"**GC/AT Ratio:** {gc_at:.3f}")
                
                codon_counts, start_count, stop_count = analyze_codons(sequence)
                st.write(f"**Start Codons:** {start_count}")
                st.write(f"**Stop Codons:** {stop_count}")
                st.write(f"**Total Codons:** {sum(codon_counts.values())}")
                
                orfs = find_orfs(sequence, 100)
                st.write(f"**ORFs (≥100 bp):** {len(orfs)}")
            
            # Download full report
            report = f"""HelixLab Analysis Report
{'='*50}

Sequence Name: {seq_name}
Length: {len(sequence):,} nucleotides
GC Content: {gc_content:.2f}%

Base Composition:
- Adenine (A): {counts['A']:,} ({counts['A']/len(sequence)*100:.2f}%)
- Thymine (T): {counts['T']:,} ({counts['T']/len(sequence)*100:.2f}%)
- Cytosine (C): {counts['C']:,} ({counts['C']/len(sequence)*100:.2f}%)
- Guanine (G): {counts['G']:,} ({counts['G']/len(sequence)*100:.2f}%)

Ratios:
- AT/GC: {at_gc:.3f}
- GC/AT: {gc_at:.3f}

Codon Analysis:
- Start Codons: {start_count}
- Stop Codons: {stop_count}
- Total Codons: {sum(codon_counts.values())}

ORF Analysis:
- ORFs Found (≥100 bp): {len(orfs)}
"""
            
            st.download_button(
                label=" Download Full Report",
                data=report,
                file_name=f"{seq_name}_analysis_report.txt",
                mime="text/plain"
            )

else:
    st.info(" Please enter a DNA sequence or upload a file to begin analysis")
    
    with st.expander("ℹ How to Use Helix_Tool"):
        st.markdown("""
        ### Features:
        - **Nucleotide Analysis**: Base counting, GC content, AT/GC ratios with interactive charts
        - **Codon Analysis**: Start/stop codon detection, usage tables with full amino acid names
        - **Protein Translation**: Multi-frame translation with full amino acid names and ORF detection
        - **DNA Conversions**: DNA to RNA and reverse complement
        - **File Format Support**: FASTA (.fasta, .fa, .fna, .ffn, .faa, .frn), GenBank (.gb, .gbk)
        - **Visualizations**: Interactive charts for nucleotide distribution, GC content, codon usage, and ORF distribution
        - **Theme Toggle**: Switch between dark and light mode
        
        ### Instructions:
        1. Choose input method (Direct Input or File Upload)
        2. Enter or upload your DNA sequence
        3. Explore different analysis tabs with interactive visualizations
        4. Adjust parameters as needed (e.g., minimum ORF length)
        5. Download results and reports
        
        ### Supported File Formats:
        - **FASTA variants**: .fasta, .fa, .fna (nucleic acid), .ffn (nucleotide coding regions), 
          .faa (amino acid), .frn (non-coding RNA)
        - **GenBank**: .gb, .gbk, .genbank
        - **Plain text**: .txt
        
        ### Example FASTA Format:
        ```
        >Sample_Sequence_1
        ATGGCTAGCTAGCTAGCTAGTAG
        >Sample_Sequence_2
        ATGAAACCCGGGTTTAAATAG
        ```
        """)

st.sidebar.markdown("---")
st.sidebar.markdown("**HelixLab** v2.0 | By T.Semwa")

# Footer
st.markdown("---")
st.markdown(
    """
    <div style='text-align: center'>
        <p>HelixLab © 2026 | Developed by T.Semwa</p>
    </div>
    """,
    unsafe_allow_html=True
)