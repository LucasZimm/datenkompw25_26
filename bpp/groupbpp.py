import re
import matplotlib.pyplot as plt
import pandas as pd

# Datei einlesen
with open("bitrates.txt", "r", encoding="utf-8") as f:
    lines = f.readlines()

gruppen = []
bpp_values = []

# Zeilen parsen
for line in lines:
    # Gruppengröße extrahieren
    match_group = re.search(r"Gruppengröße\s+(\d+)", line)
    if match_group:
        current_group = int(match_group.group(1))
        continue

    # bpp-Wert extrahieren
    match_bpp = re.search(r":\s*([\d\.]+)\s*bpp", line)
    if match_bpp:
        bpp = float(match_bpp.group(1))
        gruppen.append(current_group)
        bpp_values.append(bpp)

# DataFrame erzeugen
df = pd.DataFrame({
    "Gruppengröße": gruppen,
    "bpp": bpp_values
})

# Plot erstellen
plt.figure(figsize=(8, 5))
plt.plot(df["Gruppengröße"], df["bpp"], marker='o', linewidth=1)

plt.title("bpp vs Gruppengröße")
plt.xlabel("Gruppengröße")
plt.ylabel("bpp")
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

# Werte an Punkten anzeigen
for x, y in zip(df["Gruppengröße"], df["bpp"]):
    plt.annotate(f"{y:.3f}", xy=(x, y), xytext=(0, 6),
                 textcoords="offset points", ha='center')

plt.xticks(df["Gruppengröße"])
plt.tight_layout()

# Speichern
plt.savefig("bpp_plot.png", dpi=150)

# Anzeigen
plt.show()
