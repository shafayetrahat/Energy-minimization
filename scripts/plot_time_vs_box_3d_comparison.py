import pandas as pd
import matplotlib.pyplot as plt


# !!!!!!!!!!!!! Relative to Energy-minimization path. Check if you are in <somthing>/Energy-minimization directory. !!!!!!!!!!!!!!!!!!! 
# Read .dat file

df = pd.read_csv(
    './run/time_vs_box_3d_linked_cell.dat', 
    sep="\s+",          # Matches any whitespace (space/tab)
    header=None,        # No header in file
    names=["box_size", "particle_number", "time"]  # Assign column names
)
df1 = pd.read_csv(
    './run/time_vs_box_3d.dat', 
    sep="\s+",          # Matches any whitespace (space/tab)
    header=None,        # No header in file
    names=["box_size", "particle_number", "time"]  # Assign column names
)

df.box_size = df.box_size**3
df1.box_size = df1.box_size**3

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(df['box_size'], df['time'], 'ro-', alpha=0.9, markersize=8)
plt.plot(df1['box_size'], df1['time'], 'go-', alpha=0.8, markersize=8)

# Add particle number labels to each point
# for i, row in df.iterrows():
#     plt.annotate(f'N={row["particle_number"]}', 
#                  (row['box_size'], row['time']),
#                  textcoords="offset points",
#                  xytext=(0, 10),
#                  ha='center',
#                  fontsize=9)

# Add labels and title
plt.xlabel('Box Size(Volume)', fontsize=12)
plt.ylabel('Time(sec)', fontsize=12)
plt.title('Simulation Time vs Box Size (Density=0.15)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(['Linked Cell', 'Naive'], fontsize=12)

# Show the plot
plt.tight_layout()
plt.savefig('./run/time_vs_box_3d_compare.png')
plt.show()
