#!/usr/bin/env python3
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import csv
import os

# ---------- FASTA parsing ----------
def read_fasta(path):
    ids, seqs = [], []
    cur_id, cur_seq = None, []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    ids.append(cur_id)
                    seqs.append("".join(cur_seq))
                cur_id = line[1:].strip() or "unnamed"
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            ids.append(cur_id)
            seqs.append("".join(cur_seq))
    return list(zip(ids, seqs))

# ---------- core computations ----------
def apply_filters(seq, to_upper=True, ignore_gaps=False):
    if to_upper:
        seq = seq.upper()
    if ignore_gaps:
        seq = seq.replace("-", "")
    return seq

def unique_symbols(seq, to_upper=True, ignore_gaps=False):
    seq = apply_filters(seq, to_upper, ignore_gaps)
    return "".join(sorted(set(seq)))

def counts_and_percents(seq, to_upper=True, ignore_gaps=False):
    seq = apply_filters(seq, to_upper, ignore_gaps)
    total = len(seq)
    counts = {}
    for ch in seq:
        counts[ch] = counts.get(ch, 0) + 1
    rows = []
    for ch, c in counts.items():
        p = int((c / total) * 100) if total else 0
        rows.append((ch, c, p, total))
    rows.sort(key=lambda r: (-r[2], r[0]))  # by percent desc, then symbol
    return rows

# ---------- GUI ----------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Alphabet & Percentages — All Sequences")
        self.geometry("980x600")
        self.records = []
        self.current_path = ""

        # top bar
        top = tk.Frame(self); top.pack(fill=tk.X, padx=10, pady=8)
        tk.Button(top, text="Open FASTA…", command=self.open_file).pack(side=tk.LEFT)
        self.path_var = tk.StringVar()
        tk.Label(top, textvariable=self.path_var, anchor="w").pack(side=tk.LEFT, padx=8)

        # options + export
        opts = tk.Frame(self); opts.pack(fill=tk.X, padx=10, pady=(0,8))
        self.opt_upper = tk.BooleanVar(value=True)
        self.opt_gaps  = tk.BooleanVar(value=False)
        tk.Checkbutton(opts, text="Uppercase", variable=self.opt_upper, command=self.recompute).pack(side=tk.LEFT)
        tk.Checkbutton(opts, text="Ignore '-' gaps", variable=self.opt_gaps, command=self.recompute).pack(side=tk.LEFT, padx=(10,0))
        tk.Button(opts, text="Export CSV", command=self.export_csv).pack(side=tk.RIGHT)

        # big table
        cols = ("seq_idx","seq_id","length","alphabet","symbol","count","percent")
        self.tree = ttk.Treeview(self, columns=cols, show="headings", height=22)
        self.tree.heading("seq_idx", text="#")
        self.tree.heading("seq_id", text="Sequence ID")
        self.tree.heading("length", text="Length")
        self.tree.heading("alphabet", text="Alphabet")
        self.tree.heading("symbol", text="Symbol")
        self.tree.heading("count", text="Count")
        self.tree.heading("percent", text="Percent")

        self.tree.column("seq_idx", width=50, anchor="e")
        self.tree.column("seq_id", width=420, anchor="w")
        self.tree.column("length", width=80, anchor="e")
        self.tree.column("alphabet", width=160, anchor="w")
        self.tree.column("symbol", width=80, anchor="center")
        self.tree.column("count", width=90, anchor="e")
        self.tree.column("percent", width=90, anchor="e")
        self.tree.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0,10))

        # status
        self.status = tk.StringVar(value="Open a FASTA file to begin.")
        tk.Label(self, textvariable=self.status, anchor="w").pack(fill=tk.X, padx=10, pady=(0,10))

    def open_file(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA", "*.fa *.fasta *.faa *.fna *.fas"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            self.records = read_fasta(path)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file:\n{e}")
            return
        if not self.records:
            messagebox.showinfo("Empty", "No sequences found in the file.")
            return
        self.current_path = path
        self.path_var.set(path)
        self.recompute()

    def recompute(self):
        self.tree.delete(*self.tree.get_children())
        if not self.records:
            self.status.set("No sequences loaded.")
            return

        total_rows = 0
        for i, (rid, seq) in enumerate(self.records, 1):
            alphabet = unique_symbols(seq, self.opt_upper.get(), self.opt_gaps.get())
            rows = counts_and_percents(seq, self.opt_upper.get(), self.opt_gaps.get())
            # Insert one row per (sequence, symbol)
            for (ch, c, p, length) in rows:
                self.tree.insert(
                    "", tk.END,
                    values=(i, rid, length, alphabet, repr(ch), c, f"{p}%")
                )
                total_rows += 1

        self.status.set(f"Loaded {len(self.records)} sequence(s). Displaying {total_rows} rows.")

    def export_csv(self):
        if not self.records:
            messagebox.showinfo("Nothing to export", "Open a FASTA file first.")
            return
        default_name = os.path.splitext(os.path.basename(self.current_path or "results"))[0] + "_alphabet_stats.csv"
        save_path = filedialog.asksaveasfilename(
            title="Save CSV",
            defaultextension=".csv",
            initialfile=default_name,
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")]
        )
        if not save_path:
            return
        try:
            with open(save_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["seq_index","sequence_id","length","alphabet","symbol","count","percent"])
                for item in self.tree.get_children():
                    row = self.tree.item(item)["values"]
                    # values already formatted except percent has a '%' sign
                    w.writerow(row)
            messagebox.showinfo("Saved", f"CSV exported to:\n{save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save CSV:\n{e}")

def main():
    App().mainloop()

if __name__ == "__main__":
    main()