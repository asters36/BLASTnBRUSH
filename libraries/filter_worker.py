# This file is part of BLASTnBRUSH
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the GPL v3.0 License

from PyQt6.QtCore import QObject, pyqtSignal
import sqlite3
import re

class FilterWorker(QObject):
    finished = pyqtSignal(list)

    def __init__(self, rowid_range, filter_args):
        super().__init__()
        self.rowid_range = rowid_range 
        self.args = filter_args

    def run(self):
        results = []

        try:
            conn = sqlite3.connect("cleaned.db")
            c = conn.cursor()

            start_id, end_id = self.rowid_range
            query = "SELECT header, sequence FROM sequences WHERE rowid BETWEEN ? AND ?"
            c.execute(query, (start_id, end_id))

            for header, sequence in c:
                if not self.filters(header, sequence):
                    continue
                results.append((header, sequence))

            conn.close()

        except Exception as e:
            print("FilterWorker error:", e)

        self.finished.emit(results)

    def filters(self, header, sequence):
        args = self.args
        header_lower = header.lower()
        seq_lower = sequence.lower()

        if args.get("name_terms"):
            match = False
            for line in args["name_terms"]:
                terms = [t.strip().lower() for t in line.split(",") if t.strip()]
                if all(term in header_lower for term in terms):
                    match = True
                    break
            if not match:
                return False

        if args.get("seq_terms"):
            match_seq = False
            for term in args["seq_terms"]:
                term = term.lower()
                if not term:
                    continue
                if args.get("similarity_threshold", 100) == 100:
                    if term in seq_lower:
                        match_seq = True
                        break
                else:
                    max_mismatches = int((1 - args["similarity_threshold"] / 100) * len(term))
                    for i in range(len(seq_lower) - len(term) + 1):
                        window = seq_lower[i:i + len(term)]
                        mismatches = sum(1 for a, b in zip(term, window) if a != b)
                        if mismatches <= max_mismatches:
                            match_seq = True
                            break
                if match_seq:
                    break
            if not match_seq:
                return False

        # LENGTH
        if args.get("check_length", False):
            slen = len(sequence)
            if not (args["min_len"] <= slen <= args["max_len"]):
                return False

        # ATG/M
        if args.get("check_atg"):
            if args.get("check_m") and not sequence.startswith("M"):
                return False
            if args.get("check_atg_dna") and not sequence.startswith("ATG"):
                return False

        # Score
        if args.get("check_score"):
            m = re.search(r"score[:=]?\s*([\d.]+)", header_lower)
            if not m or not (args["min_score"] <= float(m.group(1)) <= args["max_score"]):
                return False

        # E-value
        if args.get("check_eval"):
            import math
            m = re.search(r"e[-\s]?value[:=]?\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", header_lower)
            if not m:
                return False

            try:
                e_val = float(m.group(1))
            except ValueError:
                return False

            zero = 1e-300 
            neglog = -math.log10(max(e_val, zero))
            if not (args["min_eval"] <= neglog <= args["max_eval"]):
                return False

        # Alignment length
        if args.get("check_alength"):
            m = re.search(r"alignment length[:=]?\s*(\d+)", header_lower)
            if not m or not (args["min_alength"] <= int(m.group(1)) <= args["max_alength"]):
                return False

        # Identities
        if args.get("check_identities"):
            m = re.search(r"identities[:=]?\s*(\d+)", header_lower)
            if not m or not (args["min_identities"] <= int(m.group(1)) <= args["max_identities"]):
                return False

        # Positives
        if args.get("check_positives"):
            m = re.search(r"positives[:=]?\s*(\d+)", header_lower)
            if not m or not (args["min_positives"] <= int(m.group(1)) <= args["max_positives"]):
                return False

        return True