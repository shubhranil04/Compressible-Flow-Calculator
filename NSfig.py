import matplotlib.pyplot as plt
import numpy as np

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 4))

# Plot normal shock wave (a vertical line)
shock_x = 0.5
ax.axvline(x=shock_x, color='black', linestyle='--', linewidth=2)
ax.text(shock_x + 0.01, 1.1, 'Shock Wave', fontsize=12, verticalalignment='center')

# Define flow directions and annotate flow properties
# Properties ahead of the shock (State 1)
ax.arrow(0.1, 0.8, 0.2, 0, head_width=0.05, head_length=0.02, fc='blue', ec='blue')
ax.text(0.1, 0.9, 'State 1 (Ahead of Shock)', fontsize=10, verticalalignment='center')
ax.text(0.1, 0.6, r'$\rho_1, T_1, p_1, M_1$', fontsize=12, color='blue')

# Properties behind the shock (State 2)
ax.arrow(0.6, 0.8, 0.2, 0, head_width=0.05, head_length=0.02, fc='red', ec='red')
ax.text(0.6, 0.9, 'State 2 (Behind Shock)', fontsize=10, verticalalignment='center')
ax.text(0.6, 0.6, r'$\rho_2, T_2, p_2, M_2$', fontsize=12, color='red')

# Set plot limits and labels
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.2)
ax.set_xlabel('Position along flow', fontsize=12)

# Title
ax.set_title('Normal Shock Wave', fontsize=14)

# Remove axes ticks
ax.set_xticks([])
ax.set_yticks([])

# Display the figure
plt.savefig("NormalShock.png",dpi=300)
plt.show()

