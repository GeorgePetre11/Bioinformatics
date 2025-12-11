"""
Genome Alignment Tool for Influenza and COVID-19
Uses Smith-Waterman local alignment algorithm with chunking for large sequences
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from urllib import request
import json
import time
from threading import Thread


class SmithWatermanAligner:
    """Implementation of Smith-Waterman local alignment algorithm."""
    
    def __init__(self, match_score=2, mismatch_score=-1, gap_score=-2):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        
    def align(self, seq1, seq2, max_length=5000):
        """
        Perform Smith-Waterman local alignment.
        For large sequences, truncate to max_length.
        """
        # Truncate if sequences are too long
        if len(seq1) > max_length:
            seq1 = seq1[:max_length]
        if len(seq2) > max_length:
            seq2 = seq2[:max_length]
            
        n = len(seq1)
        m = len(seq2)
        
        # Initialize matrices
        score_matrix = np.zeros((n + 1, m + 1))
        traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)
        
        max_score = 0
        max_pos = (0, 0)
        
        # Fill the score matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = score_matrix[i-1][j-1] + (
                    self.match_score if seq1[i-1] == seq2[j-1] 
                    else self.mismatch_score
                )
                delete = score_matrix[i-1][j] + self.gap_score
                insert = score_matrix[i][j-1] + self.gap_score
                
                score_matrix[i][j] = max(0, match, delete, insert)
                
                if score_matrix[i][j] == 0:
                    traceback_matrix[i][j] = -1
                elif score_matrix[i][j] == match:
                    traceback_matrix[i][j] = 0  # diagonal
                elif score_matrix[i][j] == delete:
                    traceback_matrix[i][j] = 1  # up
                else:
                    traceback_matrix[i][j] = 2  # left
                    
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)
        
        return self.traceback(seq1, seq2, score_matrix, traceback_matrix, max_pos)
    
    def traceback(self, seq1, seq2, score_matrix, traceback_matrix, max_pos):
        """Traceback from maximum score position."""
        aligned_seq1 = ""
        aligned_seq2 = ""
        
        i, j = max_pos
        start_i, start_j = i, j
        
        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            if traceback_matrix[i][j] == 0:  # diagonal
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                i -= 1
                j -= 1
            elif traceback_matrix[i][j] == 1:  # up
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                i -= 1
            else:  # left
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
                
        end_i, end_j = i, j
        
        # Calculate statistics
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                     if a == b and a != '-')
        length = len(aligned_seq1)
        identity = (matches / length * 100) if length > 0 else 0
        
        return {
            'seq1': aligned_seq1,
            'seq2': aligned_seq2,
            'matches': matches,
            'length': length,
            'identity': identity,
            'score': score_matrix[max_pos[0]][max_pos[1]],
            'start': (end_i, end_j),
            'end': max_pos,
            'score_matrix': score_matrix
        }
    
    def align_chunks(self, seq1, seq2, chunk_size=2000, overlap=200, progress_callback=None):
        """
        Align large sequences by dividing into overlapping chunks.
        Returns list of local alignments found in each chunk.
        """
        alignments = []
        n_chunks_seq1 = max(1, (len(seq1) - overlap) // (chunk_size - overlap))
        
        total_comparisons = n_chunks_seq1
        current = 0
        
        for i in range(n_chunks_seq1):
            start1 = i * (chunk_size - overlap)
            end1 = min(start1 + chunk_size, len(seq1))
            chunk1 = seq1[start1:end1]
            
            # Align chunk1 against the entire seq2 (or chunks of seq2 if too large)
            if len(seq2) > chunk_size * 2:
                # Also chunk seq2
                n_chunks_seq2 = max(1, (len(seq2) - overlap) // (chunk_size - overlap))
                total_comparisons = n_chunks_seq1 * n_chunks_seq2
                
                for j in range(n_chunks_seq2):
                    start2 = j * (chunk_size - overlap)
                    end2 = min(start2 + chunk_size, len(seq2))
                    chunk2 = seq2[start2:end2]
                    
                    result = self.align(chunk1, chunk2)
                    result['chunk1_offset'] = start1
                    result['chunk2_offset'] = start2
                    
                    if result['score'] > 50:  # Only keep significant alignments
                        alignments.append(result)
                    
                    current += 1
                    if progress_callback:
                        progress_callback(current, total_comparisons)
            else:
                result = self.align(chunk1, seq2)
                result['chunk1_offset'] = start1
                result['chunk2_offset'] = 0
                
                if result['score'] > 50:
                    alignments.append(result)
                
                current += 1
                if progress_callback:
                    progress_callback(current, total_comparisons)
        
        return alignments


class NCBIDownloader:
    """Download genome sequences from NCBI or load from local files."""
    
    @staticmethod
    def load_from_fasta(filepath):
        """Load sequence from local FASTA file."""
        try:
            with open(filepath, 'r') as f:
                data = f.read()
            
            # Parse FASTA
            lines = data.strip().split('\n')
            if not lines or not lines[0].startswith('>'):
                raise ValueError(f"Invalid FASTA format in {filepath}")
            
            header = lines[0]
            sequence = ''.join(lines[1:]).upper().replace(' ', '').replace('\n', '')
            
            # Extract accession from header if possible
            accession = header.split()[0].replace('>', '')
            
            return {
                'accession': accession,
                'header': header,
                'sequence': sequence,
                'length': len(sequence)
            }
        except Exception as e:
            raise Exception(f"Failed to load {filepath}: {str(e)}")
    
    @staticmethod
    def fetch_sequence(accession_id, retries=3):
        """
        Fetch sequence from NCBI using accession ID.
        Uses NCBI E-utilities API.
        """
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        for attempt in range(retries):
            try:
                # Fetch FASTA format
                url = f"{base_url}efetch.fcgi?db=nucleotide&id={accession_id}&rettype=fasta&retmode=text"
                
                with request.urlopen(url, timeout=30) as response:
                    data = response.read().decode('utf-8')
                    
                # Parse FASTA
                lines = data.strip().split('\n')
                if not lines or not lines[0].startswith('>'):
                    raise ValueError(f"Invalid FASTA format for {accession_id}")
                
                header = lines[0]
                sequence = ''.join(lines[1:]).upper().replace(' ', '').replace('\n', '')
                
                return {
                    'accession': accession_id,
                    'header': header,
                    'sequence': sequence,
                    'length': len(sequence)
                }
                
            except Exception as e:
                if attempt < retries - 1:
                    time.sleep(2)  # Wait before retry
                    continue
                else:
                    raise Exception(f"Failed to fetch {accession_id}: {str(e)}")
        
        return None
    
    @staticmethod
    def get_default_sequences():
        """Return default accession IDs for COVID-19 and Influenza."""
        return {
            'covid19': 'NC_045512.2',  # SARS-CoV-2 reference genome
            'influenza': 'NC_007373.1'  # Influenza A virus (H1N1) segment 1
        }


class GenomeAlignerGUI:
    """GUI application for genome alignment."""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Genome Aligner - COVID-19 vs Influenza")
        self.root.geometry("1400x900")
        
        self.aligner = SmithWatermanAligner()
        self.downloader = NCBIDownloader()
        
        self.covid_data = None
        self.influenza_data = None
        self.alignments = []
        
        self.create_widgets()
        
    def create_widgets(self):
        """Create all GUI widgets."""
        
        # Main container
        main_container = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel - Controls
        left_panel = ttk.Frame(main_container, width=350)
        main_container.add(left_panel, weight=0)
        
        # Right panel - Results
        right_panel = ttk.Frame(main_container)
        main_container.add(right_panel, weight=1)
        
        # === LEFT PANEL ===
        
        # Download section
        download_frame = ttk.LabelFrame(left_panel, text="Load Genomes", padding=10)
        download_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # COVID-19
        ttk.Label(download_frame, text="COVID-19:").grid(row=0, column=0, sticky=tk.W, pady=3)
        self.covid_entry = ttk.Entry(download_frame, width=20)
        self.covid_entry.grid(row=0, column=1, pady=3, padx=5)
        self.covid_entry.insert(0, "NC_045512.2")
        
        covid_btn_frame = ttk.Frame(download_frame)
        covid_btn_frame.grid(row=0, column=2, pady=3)
        self.covid_btn = ttk.Button(covid_btn_frame, text="NCBI", width=6,
                                    command=lambda: self.download_genome('covid'))
        self.covid_btn.pack(side=tk.LEFT, padx=2)
        ttk.Button(covid_btn_frame, text="File", width=6,
                  command=lambda: self.load_from_file('covid')).pack(side=tk.LEFT, padx=2)
        
        self.covid_status = ttk.Label(download_frame, text="Not loaded", foreground="gray")
        self.covid_status.grid(row=1, column=0, columnspan=3, sticky=tk.W, pady=2)
        
        # Influenza
        ttk.Label(download_frame, text="Influenza:").grid(row=2, column=0, sticky=tk.W, pady=3)
        self.influenza_entry = ttk.Entry(download_frame, width=20)
        self.influenza_entry.grid(row=2, column=1, pady=3, padx=5)
        self.influenza_entry.insert(0, "NC_007373.1")
        
        flu_btn_frame = ttk.Frame(download_frame)
        flu_btn_frame.grid(row=2, column=2, pady=3)
        self.influenza_btn = ttk.Button(flu_btn_frame, text="NCBI", width=6,
                                       command=lambda: self.download_genome('influenza'))
        self.influenza_btn.pack(side=tk.LEFT, padx=2)
        ttk.Button(flu_btn_frame, text="File", width=6,
                  command=lambda: self.load_from_file('influenza')).pack(side=tk.LEFT, padx=2)
        
        self.influenza_status = ttk.Label(download_frame, text="Not loaded", foreground="gray")
        self.influenza_status.grid(row=3, column=0, columnspan=3, sticky=tk.W, pady=2)
        
        # Parameters
        param_frame = ttk.LabelFrame(left_panel, text="Alignment Parameters", padding=10)
        param_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(param_frame, text="Match Score:").grid(row=0, column=0, sticky=tk.W)
        self.match_entry = ttk.Entry(param_frame, width=10)
        self.match_entry.grid(row=0, column=1, padx=5)
        self.match_entry.insert(0, "2")
        
        ttk.Label(param_frame, text="Mismatch Score:").grid(row=1, column=0, sticky=tk.W)
        self.mismatch_entry = ttk.Entry(param_frame, width=10)
        self.mismatch_entry.grid(row=1, column=1, padx=5)
        self.mismatch_entry.insert(0, "-1")
        
        ttk.Label(param_frame, text="Gap Penalty:").grid(row=2, column=0, sticky=tk.W)
        self.gap_entry = ttk.Entry(param_frame, width=10)
        self.gap_entry.grid(row=2, column=1, padx=5)
        self.gap_entry.insert(0, "-2")
        
        ttk.Label(param_frame, text="Chunk Size:").grid(row=3, column=0, sticky=tk.W)
        self.chunk_entry = ttk.Entry(param_frame, width=10)
        self.chunk_entry.grid(row=3, column=1, padx=5)
        self.chunk_entry.insert(0, "2000")
        
        # Options
        options_frame = ttk.LabelFrame(left_panel, text="Options", padding=10)
        options_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.add_n_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Add 'N' spacers between sequences",
                       variable=self.add_n_var).pack(anchor=tk.W)
        
        # Align button
        self.align_btn = ttk.Button(left_panel, text="Align Genomes", 
                                    command=self.start_alignment, state=tk.DISABLED)
        self.align_btn.pack(pady=10, ipadx=20, ipady=5)
        
        # Progress
        self.progress_var = tk.StringVar(value="Ready")
        ttk.Label(left_panel, textvariable=self.progress_var, 
                 wraplength=300).pack(pady=5)
        
        self.progress_bar = ttk.Progressbar(left_panel, mode='determinate')
        self.progress_bar.pack(fill=tk.X, padx=5, pady=5)
        
        # Info section
        info_frame = ttk.LabelFrame(left_panel, text="Sequence Info", padding=10)
        info_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.info_text = scrolledtext.ScrolledText(info_frame, height=10, width=40,
                                                   font=("Courier", 9))
        self.info_text.pack(fill=tk.BOTH, expand=True)
        
        # === RIGHT PANEL ===
        
        # Notebook for tabs
        self.notebook = ttk.Notebook(right_panel)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Visualization tab
        self.viz_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.viz_frame, text="Similarity Visualization")
        
        # Alignments tab
        self.align_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.align_frame, text="Local Alignments")
        
        self.align_text = scrolledtext.ScrolledText(self.align_frame, 
                                                    font=("Courier", 9))
        self.align_text.pack(fill=tk.BOTH, expand=True)
        
        # Statistics tab
        self.stats_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.stats_frame, text="Statistics")
        
        self.stats_text = scrolledtext.ScrolledText(self.stats_frame,
                                                    font=("Courier", 10))
        self.stats_text.pack(fill=tk.BOTH, expand=True)
        
    def load_from_file(self, genome_type):
        """Load genome from local FASTA file."""
        from tkinter import filedialog
        
        # Check for default files first
        default_file = "covid19_genome.fasta" if genome_type == 'covid' else "influenza_genome.fasta"
        import os
        
        if os.path.exists(default_file):
            filepath = default_file
        else:
            filepath = filedialog.askopenfilename(
                title=f"Select {genome_type.upper()} genome file",
                filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
            )
        
        if not filepath:
            return
        
        if genome_type == 'covid':
            status_label = self.covid_status
            button = self.covid_btn
        else:
            status_label = self.influenza_status
            button = self.influenza_btn
        
        status_label.config(text=f"Loading {os.path.basename(filepath)}...", foreground="blue")
        self.root.update()
        
        try:
            data = self.downloader.load_from_fasta(filepath)
            
            if genome_type == 'covid':
                self.covid_data = data
            else:
                self.influenza_data = data
            
            self.download_complete(genome_type, data, status_label, button)
            
        except Exception as e:
            self.download_error(genome_type, str(e), status_label, button)
    
    def download_genome(self, genome_type):
        """Download genome from NCBI."""
        if genome_type == 'covid':
            accession = self.covid_entry.get().strip()
            status_label = self.covid_status
            button = self.covid_btn
        else:
            accession = self.influenza_entry.get().strip()
            status_label = self.influenza_status
            button = self.influenza_btn
            
        if not accession:
            messagebox.showerror("Error", "Please enter an accession ID")
            return
        
        # Disable button and show status
        button.config(state=tk.DISABLED)
        status_label.config(text=f"Downloading {accession}...", foreground="blue")
        self.root.update()
        
        def download_thread():
            try:
                data = self.downloader.fetch_sequence(accession)
                
                if genome_type == 'covid':
                    self.covid_data = data
                else:
                    self.influenza_data = data
                
                # Update UI in main thread
                self.root.after(0, lambda: self.download_complete(genome_type, data, status_label, button))
                
            except Exception as e:
                self.root.after(0, lambda: self.download_error(genome_type, str(e), status_label, button))
        
        Thread(target=download_thread, daemon=True).start()
    
    def download_complete(self, genome_type, data, status_label, button):
        """Handle successful download."""
        status_label.config(
            text=f"✓ Downloaded: {data['length']} bp",
            foreground="green"
        )
        button.config(state=tk.NORMAL)
        
        self.update_info_text()
        self.check_ready_to_align()
    
    def download_error(self, genome_type, error_msg, status_label, button):
        """Handle download error."""
        status_label.config(text=f"✗ Error: {error_msg[:50]}", foreground="red")
        button.config(state=tk.NORMAL)
        messagebox.showerror("Download Error", f"Failed to download:\n{error_msg}")
    
    def update_info_text(self):
        """Update the info text with sequence details."""
        self.info_text.delete(1.0, tk.END)
        
        if self.covid_data:
            self.info_text.insert(tk.END, "COVID-19 Genome:\n")
            self.info_text.insert(tk.END, f"Accession: {self.covid_data['accession']}\n")
            self.info_text.insert(tk.END, f"Length: {self.covid_data['length']:,} bp\n")
            self.info_text.insert(tk.END, f"Preview: {self.covid_data['sequence'][:50]}...\n\n")
        
        if self.influenza_data:
            self.info_text.insert(tk.END, "Influenza Genome:\n")
            self.info_text.insert(tk.END, f"Accession: {self.influenza_data['accession']}\n")
            self.info_text.insert(tk.END, f"Length: {self.influenza_data['length']:,} bp\n")
            self.info_text.insert(tk.END, f"Preview: {self.influenza_data['sequence'][:50]}...\n")
    
    def check_ready_to_align(self):
        """Enable align button if both sequences are downloaded."""
        if self.covid_data and self.influenza_data:
            self.align_btn.config(state=tk.NORMAL)
    
    def start_alignment(self):
        """Start the alignment process in a separate thread."""
        self.align_btn.config(state=tk.DISABLED)
        self.progress_var.set("Preparing sequences...")
        self.progress_bar['value'] = 0
        
        Thread(target=self.perform_alignment, daemon=True).start()
    
    def perform_alignment(self):
        """Perform the alignment (runs in separate thread)."""
        try:
            # Get parameters
            match_score = int(self.match_entry.get())
            mismatch_score = int(self.mismatch_entry.get())
            gap_score = int(self.gap_entry.get())
            chunk_size = int(self.chunk_entry.get())
            
            self.aligner = SmithWatermanAligner(match_score, mismatch_score, gap_score)
            
            # Get sequences
            seq1 = self.covid_data['sequence']
            seq2 = self.influenza_data['sequence']
            
            # Add 'N' spacers if requested
            if self.add_n_var.get():
                seq1 = 'N' * 100 + seq1 + 'N' * 100
                seq2 = 'N' * 100 + seq2 + 'N' * 100
            
            # Perform chunked alignment
            def progress_callback(current, total):
                progress = (current / total) * 100
                self.root.after(0, lambda: self.update_progress(progress, current, total))
            
            self.root.after(0, lambda: self.progress_var.set("Aligning sequences in chunks..."))
            
            self.alignments = self.aligner.align_chunks(
                seq1, seq2, 
                chunk_size=chunk_size,
                overlap=200,
                progress_callback=progress_callback
            )
            
            # Update UI
            self.root.after(0, self.alignment_complete)
            
        except Exception as e:
            self.root.after(0, lambda: self.alignment_error(str(e)))
    
    def update_progress(self, progress, current, total):
        """Update progress bar and text."""
        self.progress_bar['value'] = progress
        self.progress_var.set(f"Processing chunk {current}/{total}...")
    
    def alignment_complete(self):
        """Handle completion of alignment."""
        self.align_btn.config(state=tk.NORMAL)
        self.progress_bar['value'] = 100
        self.progress_var.set(f"Complete! Found {len(self.alignments)} local alignments")
        
        # Display results
        self.display_alignments()
        self.display_statistics()
        self.visualize_similarities()
        
        messagebox.showinfo("Success", 
                          f"Alignment complete!\nFound {len(self.alignments)} significant local alignments")
    
    def alignment_error(self, error_msg):
        """Handle alignment error."""
        self.align_btn.config(state=tk.NORMAL)
        self.progress_var.set("Error during alignment")
        messagebox.showerror("Alignment Error", f"Error:\n{error_msg}")
    
    def display_alignments(self):
        """Display the local alignments found."""
        self.align_text.delete(1.0, tk.END)
        
        if not self.alignments:
            self.align_text.insert(tk.END, "No significant alignments found.\n")
            return
        
        # Sort by score
        sorted_aligns = sorted(self.alignments, key=lambda x: x['score'], reverse=True)
        
        self.align_text.insert(tk.END, f"Top Local Alignments (showing up to 10):\n")
        self.align_text.insert(tk.END, "=" * 80 + "\n\n")
        
        for idx, aln in enumerate(sorted_aligns[:10], 1):
            self.align_text.insert(tk.END, f"Alignment #{idx}\n")
            self.align_text.insert(tk.END, f"Score: {aln['score']:.1f}\n")
            self.align_text.insert(tk.END, f"Identity: {aln['identity']:.1f}%\n")
            self.align_text.insert(tk.END, f"Length: {aln['length']}\n")
            self.align_text.insert(tk.END, f"Matches: {aln['matches']}\n")
            
            # Show alignment (truncated if too long)
            seq1_display = aln['seq1'][:100]
            seq2_display = aln['seq2'][:100]
            if len(aln['seq1']) > 100:
                seq1_display += "..."
                seq2_display += "..."
            
            self.align_text.insert(tk.END, f"\nCOVID:  {seq1_display}\n")
            
            # Match line
            match_line = ""
            for a, b in zip(aln['seq1'][:100], aln['seq2'][:100]):
                if a == b and a != '-':
                    match_line += "|"
                elif a == '-' or b == '-':
                    match_line += " "
                else:
                    match_line += "."
            
            self.align_text.insert(tk.END, f"        {match_line}\n")
            self.align_text.insert(tk.END, f"FLU:    {seq2_display}\n")
            self.align_text.insert(tk.END, "\n" + "-" * 80 + "\n\n")
    
    def display_statistics(self):
        """Display alignment statistics."""
        self.stats_text.delete(1.0, tk.END)
        
        if not self.alignments:
            self.stats_text.insert(tk.END, "No alignments to analyze.\n")
            return
        
        # Calculate statistics
        total_alignments = len(self.alignments)
        avg_score = np.mean([a['score'] for a in self.alignments])
        avg_identity = np.mean([a['identity'] for a in self.alignments])
        avg_length = np.mean([a['length'] for a in self.alignments])
        max_score = max([a['score'] for a in self.alignments])
        max_identity = max([a['identity'] for a in self.alignments])
        
        self.stats_text.insert(tk.END, "GENOME ALIGNMENT STATISTICS\n")
        self.stats_text.insert(tk.END, "=" * 60 + "\n\n")
        
        self.stats_text.insert(tk.END, "Sequences:\n")
        self.stats_text.insert(tk.END, f"  COVID-19: {self.covid_data['accession']}\n")
        self.stats_text.insert(tk.END, f"    Length: {self.covid_data['length']:,} bp\n\n")
        self.stats_text.insert(tk.END, f"  Influenza: {self.influenza_data['accession']}\n")
        self.stats_text.insert(tk.END, f"    Length: {self.influenza_data['length']:,} bp\n\n")
        
        self.stats_text.insert(tk.END, "Local Alignments Found:\n")
        self.stats_text.insert(tk.END, f"  Total significant alignments: {total_alignments}\n")
        self.stats_text.insert(tk.END, f"  Average alignment score: {avg_score:.2f}\n")
        self.stats_text.insert(tk.END, f"  Average sequence identity: {avg_identity:.2f}%\n")
        self.stats_text.insert(tk.END, f"  Average alignment length: {avg_length:.0f} bp\n")
        self.stats_text.insert(tk.END, f"  Maximum score: {max_score:.2f}\n")
        self.stats_text.insert(tk.END, f"  Maximum identity: {max_identity:.2f}%\n\n")
        
        self.stats_text.insert(tk.END, "Score Distribution:\n")
        score_ranges = [(0, 100), (100, 200), (200, 500), (500, float('inf'))]
        for low, high in score_ranges:
            count = sum(1 for a in self.alignments if low <= a['score'] < high)
            label = f"{low}-{high}" if high != float('inf') else f"{low}+"
            self.stats_text.insert(tk.END, f"  {label}: {count} alignments\n")
        
        self.stats_text.insert(tk.END, "\n" + "=" * 60 + "\n")
        self.stats_text.insert(tk.END, "\nInterpretation:\n")
        self.stats_text.insert(tk.END, "These local alignments represent regions of similarity\n")
        self.stats_text.insert(tk.END, "between COVID-19 and Influenza genomes. Higher scores\n")
        self.stats_text.insert(tk.END, "indicate stronger local similarity, which may suggest\n")
        self.stats_text.insert(tk.END, "conserved functional regions or common viral elements.\n")
    
    def visualize_similarities(self):
        """Create visualization of similarities between genomes."""
        # Clear previous plot
        for widget in self.viz_frame.winfo_children():
            widget.destroy()
        
        if not self.alignments:
            ttk.Label(self.viz_frame, text="No alignments to visualize").pack()
            return
        
        fig = Figure(figsize=(12, 8))
        
        # Plot 1: Dot plot of alignments
        ax1 = fig.add_subplot(2, 2, 1)
        for aln in self.alignments:
            start1 = aln['chunk1_offset'] + aln['start'][0]
            end1 = aln['chunk1_offset'] + aln['end'][0]
            start2 = aln['chunk2_offset'] + aln['start'][1]
            end2 = aln['chunk2_offset'] + aln['end'][1]
            
            # Color by score
            color = plt.cm.hot(min(aln['score'] / 500, 1.0))
            ax1.plot([start2, end2], [start1, end1], 'o-', 
                    color=color, markersize=3, alpha=0.6)
        
        ax1.set_xlabel('Influenza Genome Position (bp)')
        ax1.set_ylabel('COVID-19 Genome Position (bp)')
        ax1.set_title('Alignment Positions Dot Plot')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Score distribution
        ax2 = fig.add_subplot(2, 2, 2)
        scores = [a['score'] for a in self.alignments]
        ax2.hist(scores, bins=30, color='steelblue', edgecolor='black')
        ax2.set_xlabel('Alignment Score')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Score Distribution')
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Identity vs Length scatter
        ax3 = fig.add_subplot(2, 2, 3)
        identities = [a['identity'] for a in self.alignments]
        lengths = [a['length'] for a in self.alignments]
        scatter = ax3.scatter(lengths, identities, c=scores, cmap='hot',
                             s=50, alpha=0.6, edgecolors='black', linewidth=0.5)
        ax3.set_xlabel('Alignment Length (bp)')
        ax3.set_ylabel('Sequence Identity (%)')
        ax3.set_title('Identity vs Length')
        ax3.grid(True, alpha=0.3)
        plt.colorbar(scatter, ax=ax3, label='Score')
        
        # Plot 4: Coverage along genomes
        ax4 = fig.add_subplot(2, 2, 4)
        
        # Create coverage arrays
        covid_len = self.covid_data['length']
        flu_len = self.influenza_data['length']
        
        bin_size = 1000
        covid_coverage = np.zeros(covid_len // bin_size + 1)
        flu_coverage = np.zeros(flu_len // bin_size + 1)
        
        for aln in self.alignments:
            start1 = aln['chunk1_offset'] + aln['start'][0]
            end1 = aln['chunk1_offset'] + aln['end'][0]
            start2 = aln['chunk2_offset'] + aln['start'][1]
            end2 = aln['chunk2_offset'] + aln['end'][1]
            
            covid_coverage[start1 // bin_size:end1 // bin_size + 1] += 1
            flu_coverage[start2 // bin_size:end2 // bin_size + 1] += 1
        
        x_covid = np.arange(len(covid_coverage)) * bin_size / 1000
        x_flu = np.arange(len(flu_coverage)) * bin_size / 1000
        
        ax4.plot(x_covid, covid_coverage, label='COVID-19', color='red', alpha=0.7)
        ax4.plot(x_flu, flu_coverage, label='Influenza', color='blue', alpha=0.7)
        ax4.set_xlabel('Genome Position (kbp)')
        ax4.set_ylabel('Alignment Coverage')
        ax4.set_title('Alignment Coverage Distribution')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        fig.suptitle('Genome Similarity Analysis: COVID-19 vs Influenza', 
                    fontsize=14, fontweight='bold')
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=self.viz_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


def main():
    root = tk.Tk()
    app = GenomeAlignerGUI(root)
    
    # Add instructions on startup
    instructions = """
GENOME ALIGNER - Instructions:

1. Download both genomes using default accession IDs or enter custom ones
2. Adjust alignment parameters if needed (defaults work well)
3. Click "Align Genomes" to start the analysis
4. View results in the tabs:
   - Similarity Visualization: Graphical analysis
   - Local Alignments: Detailed alignment regions
   - Statistics: Summary statistics

Note: Large genomes are processed in chunks for efficiency.
The 'N' spacers help prevent end-effects in chunked alignment.
"""
    app.info_text.insert(1.0, instructions)
    
    root.mainloop()


if __name__ == "__main__":
    main()
