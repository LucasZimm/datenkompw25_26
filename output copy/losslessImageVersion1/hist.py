import matplotlib.pyplot as plt
import glob

# ---------------------------
# Funktion zum Einlesen der Histogramm-Dateien
# ---------------------------
def read_histogram(filename):
    xs = []
    ys = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            x, y = line.split()
            xs.append(int(x))
            ys.append(int(y))
    return xs, ys

# ---------------------------
# Plot-Funktion
# ---------------------------
def plot_histogram(xs, ys, title, xlabel, ylabel, output_file):
    plt.figure(figsize=(12,6))
    plt.bar(xs, ys, width=1.0, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Saved plot: {output_file}")

# ---------------------------
# Fehler-Histogramme (*.err.txt)
# ---------------------------
err_files = glob.glob("*_err.txt")
for fname in err_files:
    xs, ys = read_histogram(fname)
    plot_histogram(xs, ys, title=f"Error Histogram: {fname}",
                   xlabel="Error", ylabel="Count",
                   output_file=f"{fname[:-4]}.png")

# ---------------------------
# Kontext-Histogramme (*.ctx.txt)
# ---------------------------
ctx_files = glob.glob("*_ctx.txt")
for fname in ctx_files:
    xs, ys = read_histogram(fname)
    plot_histogram(xs, ys, title=f"Context Histogram: {fname}",
                   xlabel="Pixel value / Context", ylabel="Count",
                   output_file=f"{fname[:-4]}.png")
