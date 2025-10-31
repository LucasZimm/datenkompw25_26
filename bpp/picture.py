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

# Plotten
plt.figure(figsize=(10,6))
df.plot(kind="bar", figsize=(12,6))
plt.ylabel("Bitrate (bpp)")
plt.title("Vergleich der Bitraten pro Bild und Codec")
plt.legend(title="Codec")
plt.tight_layout()
plt.show()
