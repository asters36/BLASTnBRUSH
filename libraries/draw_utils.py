# This file is part of BLASTnBRUSH
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the GPL v3.0 License

###################### LENGTH HISTOGRAM ##############################
def draw_length_histogram(chart_widget, lengths, bin_count=5):
    if not lengths or len(lengths) < 2:
        return

    min_len = min(lengths)
    max_len = max(lengths)
    bin_size = max(1, round((max_len - min_len) / bin_count))

    bins = list(range(min_len, max_len + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for length in lengths:
        for i in range(len(bins) - 1):
            if bins[i] <= length < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if length == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bar_container = ax.bar(labels, counts, picker=True)
    ax.set_ylim(0, max(counts) * 1.15)
    for rect, count in zip(bar_container, counts):
        height = rect.get_height()
        if height > 0:
            ax.text(rect.get_x() + rect.get_width() / 2, height + 1,
                    str(count), ha='center', va='bottom', fontsize=8)
    
    ax.set_title("Length")
    ax.set_xlabel("Length range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()

    parent = chart_widget.parentWidget()
    if parent:
        parent.length_bins = bins
        parent.length_bar_patches = list(bar_container.patches)

###################### SCORE HISTOGRAM ##############################
def draw_bitscore_histogram(chart_widget, headers):
    import re

    scores = []
    for header in headers:
        match = re.search(r"Score: ([\d\.]+)", header)
        if match:
            try:
                scores.append(float(match.group(1)))
            except ValueError:
                pass

    if not scores or len(scores) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Scores Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_score = min(scores)
    max_score = max(scores)
    bin_count = 5
    bin_size = max(1, round((max_score - min_score) / bin_count))
    bins = list(range(int(min_score), int(max_score) + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for score in scores:
        for i in range(len(bins) - 1):
            if bins[i] <= score < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if score == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bars = ax.bar(labels, counts)
    ax.set_ylim(0, max(counts) * 1.15)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    height + 1,
                    str(count),
                    ha='center',
                    va='bottom',
                    fontsize=8)
    
                    
    ax.set_title("Bit Score")
    ax.set_xlabel("Score range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()
    
###################### EVALUE HISTOGRAM ##############################   
def draw_evalue_histogram(chart_widget, headers):
    import re, math

    e_raw = []
    for header in headers:
        m = re.search(r"E-Value:\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", header)
        if m:
            try:
                e = float(m.group(1))
                if e >= 0:
                    e_raw.append(e)
            except ValueError:
                pass

    if len(e_raw) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No E-Values Info", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([]); ax.axis('off')
        chart_widget.canvas.draw()
        return

    zero_count = sum(1 for e in e_raw if e == 0)
    TINY = 1e-300
    values = [-math.log10(max(e, TINY)) for e in e_raw if e > 0]

    vmax = max(values)
    vmin = 0  
    bins_count = 5
    step = (vmax - vmin) / bins_count
    bins_edges = [vmin + i * step for i in range(bins_count + 1)]

    counts = [0] * bins_count
    for val in values:
        if val < bins_edges[0]:
            counts[0] += 1
        elif val >= bins_edges[-1]:
            counts[-1] += 1
        else:
            for i in range(bins_count):
                if bins_edges[i] <= val < bins_edges[i + 1]:
                    counts[i] += 1
                    break

    counts[-1] += zero_count

    def log_label(log_min, log_max, is_first=False):
        if is_first:
            return f"0–$10^{{-{int(abs(round(log_max)))}}}$"
        exp1 = int(abs(round(log_min)))
        exp2 = int(abs(round(log_max)))
        if exp1 == exp2:
            return f"$10^{{-{exp1}}}$"
        return f"$10^{{-{exp1}}}$–$10^{{-{exp2}}}$"

    labels = [
        log_label(bins_edges[i], bins_edges[i + 1], is_first=(i == 0))
        for i in range(bins_count)
    ]

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    bars = ax.bar(labels, counts)

    ymax = max(counts) if any(counts) else 1
    ax.set_ylim(0, ymax * 1.15)

    for bar, count in zip(bars, counts):
        if count > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (0.02 * ymax),
                str(count),
                ha='center', va='bottom', fontsize=8
            )

    ax.set_title("E-Value")
    ax.set_xlabel("E-Value(10⁻ⁿ)")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=9)
    chart_widget.figure.tight_layout()
    chart_widget.canvas.draw()
    
    
###################### ALIGNMENT HISTOGRAM ##############################    
def draw_alength_histogram(chart_widget, lengths, bin_count=5):
    if len(lengths) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Alignment Lengths Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_len = min(lengths)
    max_len = max(lengths)
    bin_size = max(1, round((max_len - min_len) / bin_count))

    bins = list(range(min_len, max_len + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)
    
    

    for length in lengths:
        for i in range(len(bins) - 1):
            if bins[i] <= length < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if length == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bar_container = ax.bar(labels, counts, picker=True)
    ax.set_ylim(0, max(counts) * 1.15)
    for rect, count in zip(bar_container, counts):
        height = rect.get_height()
        if height > 0:
            ax.text(rect.get_x() + rect.get_width() / 2, height + 1,
                    str(count), ha='center', va='bottom', fontsize=8)
    ax.set_title("Alignment Length")
    ax.set_xlabel("Length range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()

    parent = chart_widget.parentWidget()
    if parent:
        parent.length_bins = bins
        parent.length_bar_patches = list(bar_container.patches)
        
        
###################### IDENTITY HISTOGRAM ##############################        
def draw_identities_histogram(chart_widget, headers):
    import re

    scores = []
    for header in headers:
        match = re.search(r"Identities: ([\d\.]+)", header)
        if match:
            try:
                scores.append(float(match.group(1)))
            except ValueError:
                pass

    if not scores or len(scores) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Identity Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_score = min(scores)
    max_score = max(scores)
    bin_count = 5
    bin_size = max(1, round((max_score - min_score) / bin_count))
    bins = list(range(int(min_score), int(max_score) + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for score in scores:
        for i in range(len(bins) - 1):
            if bins[i] <= score < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if score == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bars = ax.bar(labels, counts)
    ax.set_ylim(0, max(counts) * 1.15)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    height + 1,
                    str(count),
                    ha='center',
                    va='bottom',
                    fontsize=8)
    
                    
    ax.set_title("Identity")
    ax.set_xlabel("Ident. range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()
    
    
###################### SIMILARITY HISTOGRAM ##############################
def draw_positives_histogram(chart_widget, headers):
    import re

    scores = []
    for header in headers:
        match = re.search(r"Positives: ([\d\.]+)", header)
        if match:
            try:
                scores.append(float(match.group(1)))
            except ValueError:
                pass

    if not scores or len(scores) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Similarity Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_score = min(scores)
    max_score = max(scores)
    bin_count = 5
    bin_size = max(1, round((max_score - min_score) / bin_count))
    bins = list(range(int(min_score), int(max_score) + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for score in scores:
        for i in range(len(bins) - 1):
            if bins[i] <= score < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if score == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bars = ax.bar(labels, counts)
    ax.set_ylim(0, max(counts) * 1.15)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    height + 1,
                    str(count),
                    ha='center',
                    va='bottom',
                    fontsize=8)
    
                    
    ax.set_title("Similarity")
    ax.set_xlabel("Sim. range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()