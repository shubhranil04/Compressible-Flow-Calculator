import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 4))

# Plot oblqiue shock wave (a vertical line)
shock_x = 0.05
beta = np.pi*60/180
delta = np.pi*30/180	

ax.axline(xy1=(shock_x,0), slope = np.tan(beta), color='black', linestyle='--', linewidth=2)

arc_beta = Arc((shock_x,0), 0.3, 0.4, theta1=0, theta2=60, color='black', lw=1.5)
ax.annotate(r'$\beta$', xy=(shock_x+0.125,0.125), fontsize=12, color='black')
ax.add_patch(arc_beta)

ax.axline(xy1=(shock_x,0), slope = np.tan(delta), color='black', linestyle='-', linewidth=2)
arc_delta = Arc((shock_x,0), 0.5, 0.6, theta1=0, theta2=30, color='black', lw=1.5)
ax.annotate(r'$\delta$', xy=(shock_x+0.25,0.08), fontsize=12, color='black')
ax.add_patch(arc_delta)

ax.text(0.5, 1.1, 'Shock Wave', fontsize=12, verticalalignment='center')

# Define flow directions and annotate flow properties
# Properties ahead of the shock (State 1)
ax.arrow(0.05, 0.6, 0.2, 0, head_width=0.05, head_length=0.02, fc='blue', ec='blue')
ax.text(0.1, 0.8, 'State 1 (Ahead of Shock)', fontsize=10, verticalalignment='center')
ax.text(0.1, 0.5, r'$\rho_1, T_1, p_1, M_1$', fontsize=12, color='blue')

# Properties behind the shock (State 2)
ax.arrow(0.4, 0.4, 0.2*np.cos(delta), 0.2*np.sin(delta), head_width=0.05, head_length=0.02, fc='red', ec='red')
ax.text(0.6, 0.8, 'State 2 (Behind Shock)', fontsize=10, verticalalignment='center')
ax.text(0.65, 0.5, r'$\rho_2, T_2, p_2, M_2$', fontsize=12, color='red')

# Set plot limits and labels
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.2)
ax.set_xlabel('Position along flow', fontsize=12)

# Title
ax.set_title('Oblique Shock Wave', fontsize=14)

# Remove axes ticks
ax.set_xticks([])
ax.set_yticks([])

# Display the figure
plt.savefig("ObliqueShock.png",dpi=300)
plt.show()

