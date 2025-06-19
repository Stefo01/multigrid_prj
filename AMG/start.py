import subprocess
import pyvista as pv
import numpy as np # Import numpy for potential scalar manipulation

# --- Configuration for plotting ---
# Define a custom color map (cmap). You can choose from many built-in colormaps
# listed here: https://matplotlib.org/stable/users/explain/colors/colormaps.html
# 'turbo', 'magma', 'plasma', 'cividis', 'hsv' are good alternatives.
CUSTOM_CMAP = "turbo" # or 'viridis', 'plasma', 'turbo', 'cividis', 'inferno', etc.

# Color for mesh edges (if shown)
EDGE_COLOR = "black"

# Background color
BACKGROUND_COLOR = "white"

# --- Main Script ---

# Step 1: Run your C++ program that outputs a VTU file
print("Running C++ program AMG")
try:
    result=subprocess.run(["./AMG"], check=True, capture_output=True, text=True)
    print("C++ program finished successfully.")
    print("\n--- Output from ./AMG (stdout) ---")
    print(result.stdout)
    print("-----------------------------------\n")

except subprocess.CalledProcessError as e:
    print(f"Error running C++ program: {e}")
    print(f"Stdout:\n{e.stdout}")
    print(f"Stderr:\n{e.stderr}")
    exit(1) # Exit if the C++ program fails

# Step 2: Load the generated VTU file
print("Loading output...")
try:
    mesh = pv.read("output.vtu")
    print("Mesh loaded successfully.")
except FileNotFoundError:
    print("Error: output.vtu not found. Make sure your C++ program creates it.")
    exit(1)
except Exception as e:
    print(f"Error loading VTU file: {e}")
    exit(1)

# Step 3: Use scalar field 'u' and ensure it's active
if 'u' in mesh.point_data:
    mesh.set_active_scalars('u')
    print(f"Active scalars set to 'u'. Min: {mesh.active_scalars.min():.2f}, Max: {mesh.active_scalars.max():.2f}")
else:
    print("Warning: Scalar field 'u' not found in mesh. Using default scalars if available.")
    # If 'u' is not found, you might want to create some dummy scalars
    # or handle this case based on your data.
    # For demonstration, let's create a dummy scalar if no 'u' for warping
    if mesh.n_points > 0 and mesh.active_scalars is None:
        print("Creating dummy scalars for demonstration.")
        mesh['dummy_scalars'] = np.arange(mesh.n_points)
        mesh.set_active_scalars('dummy_scalars')


# Step 4: Plot and save 2D PNG with enhanced colors (no rotation needed for default 'xy')
print("Generating 2D plot (output_2D.png)...")
plotter_2d = pv.Plotter(off_screen=True)
plotter_2d.add_mesh(mesh, show_edges=True, edge_color=EDGE_COLOR, cmap=CUSTOM_CMAP,
                    scalar_bar_args={'title': 'Scalar Field (u)', 'label_font_size': 12, 'title_font_size': 16})
plotter_2d.set_background(BACKGROUND_COLOR)
plotter_2d.camera_position = 'xy'
plotter_2d.show(screenshot="output_2D.png")
print("Saved: output_2D.png")


# --- Step 5: Create a 3D angled view with warp_by_scalar and custom rotation ---
print("Generating 3D warped plot (output_3D.png)...")
if mesh.active_scalars is not None:
    bounds = mesh.bounds
    max_dim = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
    warp_scale_factor = max_dim * 0.1 # Adjust this multiplier

    warped_mesh = mesh.warp_by_scalar(factor=warp_scale_factor)
    print(f"Mesh warped by scalar with factor: {warp_scale_factor:.2f}")

    plotter_3d = pv.Plotter(off_screen=True)
    plotter_3d.add_mesh(warped_mesh,
                        show_edges=True,
                        edge_color=EDGE_COLOR,
                        cmap=CUSTOM_CMAP,
                        scalar_bar_args={'title': 'Scalar Field (u)', 'label_font_size': 12, 'title_font_size': 16})
    plotter_3d.set_background(BACKGROUND_COLOR)

    # --- ROTATION PART STARTS HERE ---

    # Option A: Start from an isometric view and then rotate
    plotter_3d.view_isometric() # Default isometric view

    # Rotate the camera by specific angles
    plotter_3d.camera.azimuth += 30  # Rotate 30 degrees horizontally (around Z-axis by default)
    plotter_3d.camera.elevation += 15 # Rotate 15 degrees vertically (lift camera up)
    # plotter_3d.camera.roll += 5 # Tilt the camera slightly if needed

    # Option B: Set a specific camera position from scratch (more advanced)
    # (position_x, position_y, position_z), (focal_point_x, focal_point_y, focal_point_z), (view_up_x, view_up_y, view_up_z)
    # This gives you full control. You might need to experiment to find good values.
    # For a 2D mesh warped into 3D, if it's in the XY plane:
    # A good camera position might be something like:
    # cam_pos = (mesh.center[0] + max_dim, mesh.center[1] + max_dim, mesh.center[2] + max_dim)
    # cam_focal_point = mesh.center
    # cam_view_up = (0, 0, 1) # Z-axis is up
    # plotter_3d.camera_position = [cam_pos, cam_focal_point, cam_view_up]


    # --- ROTATION PART ENDS HERE ---

    plotter_3d.show(screenshot="output_3D.png")
    print("Saved: output_3D.png")
else:
    print("Skipping 3D warped plot: No active scalars found to warp by.")

print("\nAll plots generated successfully!")