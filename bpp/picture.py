import pandas as pd
import matplotlib.pyplot as plt
import glob
import re

# Alle Dateien einlesen, die mit "bitrates" beginnen
files = glob.glob("bitrates*.txt")

data = {}

# Jede Datei parsen
for file in files:
    codec = re.findall(r'bitrates(.*)\.txt', file)[0]  # z. B. "AV1"
    values = {}
    with open(file, "r") as f:
        for line in f:
            match = re.match(r"kodim(\d+):\s*([\d.]+)", line)
            if match:
                img = f"kodim{match.group(1)}"
                bpp = float(match.group(2))
                values[img] = bpp
    data[codec] = values

# DataFrame bauen
df = pd.DataFrame(data)
df.index.name = "Image"

# 1. Plot: Bitrate pro Bild
plt.figure(figsize=(12,6))
df.plot(kind="bar", figsize=(12,6))
plt.ylabel("Bitrate (bpp)")
plt.title("Vergleich der Bitraten pro Bild und Codec")
plt.legend(title="Codec")
plt.tight_layout()
plt.show()

# Durchschnittliche Bitrate pro Codec
avg_bpp = df.mean()  # Mittelwert pro Codec

plt.figure(figsize=(8,5))
bars = plt.bar(avg_bpp.index, avg_bpp.values, color="skyblue")

# Werte über die Balken schreiben
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, height + 0.005, f"{height:.3f}", 
             ha='center', va='bottom', fontsize=10)

plt.ylabel("Average Bitrate (bpp)")
plt.title("Durchschnittliche Bitrate pro Codec")
plt.ylim(0, avg_bpp.max()*1.15)  # etwas Platz nach oben für die Zahlen
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()

