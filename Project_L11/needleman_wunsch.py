import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


class NeedlemanWunschAligner:
    """Implements the Needleman-Wunsch algorithm for sequence alignment."""
    
    def __init__(self, seq1, seq2, match_score=1, mismatch_score=-1, gap_score=0):
        self.seq1 = seq1.upper()
        self.seq2 = seq2.upper()
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.score_matrix = None
        self.traceback_matrix = None
        
    def align(self):
        """Perform the alignment using Needleman-Wunsch algorithm."""
        n = len(self.seq1)
        m = len(self.seq2)
        
        # Initialize score matrix
        self.score_matrix = np.zeros((n + 1, m + 1))
        self.traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)
        
        # Initialize first row and column
        for i in range(n + 1):
            self.score_matrix[i][0] = i * self.gap_score
        for j in range(m + 1):
            self.score_matrix[0][j] = j * self.gap_score
            
        # Fill the score matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = self.score_matrix[i-1][j-1] + (
                    self.match_score if self.seq1[i-1] == self.seq2[j-1] 
                    else self.mismatch_score
                )
                delete = self.score_matrix[i-1][j] + self.gap_score
                insert = self.score_matrix[i][j-1] + self.gap_score
                
                self.score_matrix[i][j] = max(match, delete, insert)
                
                # Store traceback direction
                if self.score_matrix[i][j] == match:
                    self.traceback_matrix[i][j] = 0  # diagonal
                elif self.score_matrix[i][j] == delete:
                    self.traceback_matrix[i][j] = 1  # up
                else:
                    self.traceback_matrix[i][j] = 2  # left
                    
        return self.traceback()
    
    def traceback(self):
        """Traceback to find the optimal alignment."""
        aligned_seq1 = ""
        aligned_seq2 = ""
        match_line = ""
        
        i = len(self.seq1)
        j = len(self.seq2)
        
        traceback_path = [(i, j)]
        
        while i > 0 or j > 0:
            if i > 0 and j > 0 and self.traceback_matrix[i][j] == 0:
                # Diagonal
                aligned_seq1 = self.seq1[i-1] + aligned_seq1
                aligned_seq2 = self.seq2[j-1] + aligned_seq2
                if self.seq1[i-1] == self.seq2[j-1]:
                    match_line = "|" + match_line
                else:
                    match_line = " " + match_line
                i -= 1
                j -= 1
            elif i > 0 and (j == 0 or self.traceback_matrix[i][j] == 1):
                # Up (gap in seq2)
                aligned_seq1 = self.seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                match_line = " " + match_line
                i -= 1
            else:
                # Left (gap in seq1)
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = self.seq2[j-1] + aligned_seq2
                match_line = " " + match_line
                j -= 1
                
            traceback_path.append((i, j))
            
        # Calculate statistics
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        length = len(aligned_seq1)
        similarity = (matches / length * 100) if length > 0 else 0
        
        return {
            'seq1': aligned_seq1,
            'seq2': aligned_seq2,
            'match_line': match_line,
            'matches': matches,
            'length': length,
            'similarity': similarity,
            'traceback_path': traceback_path[::-1]
        }


class NeedlemanWunschGUI:
    """GUI application for Needleman-Wunsch sequence alignment."""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Needleman-Wunsch Sequence Aligner")
        self.root.geometry("1400x800")
        
        # Default sequences
        self.seq1_default = "ACCGTGAAGCCAATAC"
        self.seq2_default = "AGCGTGCAGCCAATAC"
        
        self.create_widgets()
        
    def create_widgets(self):
        """Create all GUI widgets."""
        
        # Left panel for input and controls
        left_frame = ttk.Frame(self.root, padding="10")
        left_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Sequences input
        seq_frame = ttk.LabelFrame(left_frame, text="Sequences", padding="10")
        seq_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(seq_frame, text="Sq 1 =").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.seq1_entry = ttk.Entry(seq_frame, width=30)
        self.seq1_entry.grid(row=0, column=1, pady=2)
        self.seq1_entry.insert(0, self.seq1_default)
        
        ttk.Label(seq_frame, text="Sq 2 =").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.seq2_entry = ttk.Entry(seq_frame, width=30)
        self.seq2_entry.grid(row=1, column=1, pady=2)
        self.seq2_entry.insert(0, self.seq2_default)
        
        # Parameters
        param_frame = ttk.LabelFrame(left_frame, text="Parameters", padding="10")
        param_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(param_frame, text="Gap =").grid(row=0, column=0, sticky=tk.W)
        self.gap_entry = ttk.Entry(param_frame, width=10)
        self.gap_entry.grid(row=0, column=1, padx=5)
        self.gap_entry.insert(0, "0")
        
        ttk.Label(param_frame, text="Match =").grid(row=1, column=0, sticky=tk.W)
        self.match_entry = ttk.Entry(param_frame, width=10)
        self.match_entry.grid(row=1, column=1, padx=5)
        self.match_entry.insert(0, "1")
        
        ttk.Label(param_frame, text="MMach =").grid(row=2, column=0, sticky=tk.W)
        self.mismatch_entry = ttk.Entry(param_frame, width=10)
        self.mismatch_entry.grid(row=2, column=1, padx=5)
        self.mismatch_entry.insert(0, "-1")
        
        # Options
        options_frame = ttk.LabelFrame(left_frame, text="Options", padding="10")
        options_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5)
        
        self.plot_traceback = tk.BooleanVar(value=True)
        self.plot_grid = tk.BooleanVar(value=True)
        
        ttk.Checkbutton(options_frame, text="Plot TraceBack", 
                       variable=self.plot_traceback).grid(row=0, column=0, sticky=tk.W)
        ttk.Checkbutton(options_frame, text="Plot grid", 
                       variable=self.plot_grid).grid(row=1, column=0, sticky=tk.W)
        
        # Align button
        ttk.Button(left_frame, text="Align", command=self.perform_alignment).grid(
            row=3, column=0, pady=10, ipadx=20, ipady=5)
        
        # Right panel for visualizations
        right_frame = ttk.Frame(self.root, padding="10")
        right_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(right_frame)
        self.notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Matrix visualization tab
        self.matrix_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.matrix_frame, text="Alignment Matrix")
        
        # Traceback visualization tab
        self.traceback_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.traceback_frame, text="Traceback Path")
        
        # Results text area
        result_frame = ttk.LabelFrame(right_frame, text="Show Alignment:", padding="10")
        result_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=10)
        
        self.result_text = tk.Text(result_frame, height=15, width=80, font=("Courier", 10))
        self.result_text.grid(row=0, column=0, sticky=(tk.W, tk.E))
        
        scrollbar = ttk.Scrollbar(result_frame, orient="vertical", command=self.result_text.yview)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.result_text.configure(yscrollcommand=scrollbar.set)
        
        # Configure grid weights
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)
        
    def perform_alignment(self):
        """Perform the sequence alignment and update visualizations."""
        # Get parameters
        seq1 = self.seq1_entry.get().strip()
        seq2 = self.seq2_entry.get().strip()
        
        try:
            gap_score = int(self.gap_entry.get())
            match_score = int(self.match_entry.get())
            mismatch_score = int(self.mismatch_entry.get())
        except ValueError:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(1.0, "Error: Invalid parameter values!")
            return
            
        if not seq1 or not seq2:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(1.0, "Error: Please enter both sequences!")
            return
        
        # Perform alignment
        aligner = NeedlemanWunschAligner(seq1, seq2, match_score, mismatch_score, gap_score)
        result = aligner.align()
        
        # Display results
        self.display_results(result)
        
        # Visualize matrix
        self.visualize_matrix(aligner, result)
        
        # Visualize traceback path
        self.visualize_traceback(aligner, result)
        
    def display_results(self, result):
        """Display alignment results in text area."""
        self.result_text.delete(1.0, tk.END)
        
        output = f"{result['seq1']}\n"
        output += f"{result['match_line']}\n"
        output += f"{result['seq2']}\n\n"
        output += f"Matches = {result['matches']}\n"
        output += f"Length = {result['length']}\n"
        output += f"Similarity = {result['similarity']:.0f} %\n\n"
        output += f"Tracing back: M[{len(result['seq1'].replace('-', ''))},{len(result['seq2'].replace('-', ''))}]\n"
        
        self.result_text.insert(1.0, output)
        
    def visualize_matrix(self, aligner, result):
        """Create heatmap visualization of the score matrix."""
        # Clear previous plot
        for widget in self.matrix_frame.winfo_children():
            widget.destroy()
            
        fig = Figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        
        # Create heatmap
        im = ax.imshow(aligner.score_matrix, cmap='hot', aspect='auto', origin='upper')
        
        # Add colorbar
        fig.colorbar(im, ax=ax)
        
        # Set labels
        ax.set_xlabel('Sequence 2', fontsize=10)
        ax.set_ylabel('Sequence 1', fontsize=10)
        ax.set_title('Graphic representation of the alignment matrix', fontsize=11)
        
        # Add grid if option is selected
        if self.plot_grid.get():
            ax.grid(True, which='both', color='gray', linewidth=0.5, alpha=0.3)
            
        # Plot traceback path if option is selected
        if self.plot_traceback.get() and result.get('traceback_path'):
            path = result['traceback_path']
            path_i = [p[0] for p in path]
            path_j = [p[1] for p in path]
            ax.plot(path_j, path_i, 'b-', linewidth=2, alpha=0.7, label='Traceback path')
            ax.legend()
        
        # Set tick labels
        ax.set_xticks(range(len(aligner.seq2) + 1))
        ax.set_yticks(range(len(aligner.seq1) + 1))
        ax.set_xticklabels(['-'] + list(aligner.seq2), fontsize=8)
        ax.set_yticklabels(['-'] + list(aligner.seq1), fontsize=8)
        
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=self.matrix_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def visualize_traceback(self, aligner, result):
        """Create visualization of traceback path deviation."""
        # Clear previous plot
        for widget in self.traceback_frame.winfo_children():
            widget.destroy()
            
        fig = Figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        
        n = len(aligner.seq1)
        m = len(aligner.seq2)
        
        # Create a matrix showing path deviation
        deviation_matrix = np.zeros((n + 1, m + 1))
        
        if result.get('traceback_path'):
            path = result['traceback_path']
            for i, j in path:
                deviation_matrix[i, j] = 1
                
        # Color the matrix
        colors = ['lightyellow', 'red']
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(colors)
        
        im = ax.imshow(deviation_matrix, cmap=cmap, aspect='auto', origin='upper')
        
        # Add grid
        ax.grid(True, which='both', color='gray', linewidth=0.5)
        ax.set_xticks(np.arange(-.5, m + 1, 1), minor=True)
        ax.set_yticks(np.arange(-.5, n + 1, 1), minor=True)
        ax.grid(which='minor', color='gray', linewidth=1)
        
        # Set labels
        ax.set_xlabel('Sequence 2', fontsize=10)
        ax.set_ylabel('Sequence 1', fontsize=10)
        ax.set_title('Traceback path deviation from optimal alignment (diagonal)', fontsize=11)
        
        # Set tick labels
        ax.set_xticks(range(m + 1))
        ax.set_yticks(range(n + 1))
        ax.set_xticklabels(['-'] + list(aligner.seq2), fontsize=8)
        ax.set_yticklabels(['-'] + list(aligner.seq1), fontsize=8)
        
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=self.traceback_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


def main():
    root = tk.Tk()
    app = NeedlemanWunschGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
