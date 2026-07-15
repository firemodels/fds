'''
----------------------------------------------------------------------
This script visualizes the spatial coverage of radiation rays emitted 
from a spherical source located at the center of a cubic enclosure.

Randomly rotated FDS angular distributions are projected onto the cube 
surfaces to demonstrate how increasing the number of random rotations 
improves surface coverage and reduces directional bias (ray0effect).

The visualization is developed for illustration purposes in the 
FDS Technical Guide.
----------------------------------------------------------------------
'''
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path
from PIL import Image

# Use maximum lossless compression for raster layers embedded in PDF files.
plt.rcParams["pdf.compression"] = 9



# ============================================================
# USER PARAMETERS
# ============================================================

REQUESTED_NRA = 24


# ============================================================
# FDS-style angular discretization
# ============================================================

def final_number_radiation_angles(initial_nra):

    nra = initial_nra

    nrt = 2 * round(0.5 * 1.17 * nra ** (1.0 / 2.26))

    nrp = []
    total = 0

    for i in range(1, nrt + 1):

        theta_low = math.pi * (i - 1) / nrt
        theta_up  = math.pi * i / nrt

        nphi = round(
            0.5 * nra *
            (math.cos(theta_low) - math.cos(theta_up))
        )

        nphi = max(4, nphi)
        nphi = 4 * round(nphi / 4)

        nrp.append(int(nphi))
        total += int(nphi)

    return total, nrt, nrp



NRA, NRT, NRP = final_number_radiation_angles(
    REQUESTED_NRA
)

print("NRA =", NRA)
print("NRT =", NRT)
print("NRP =", NRP)



# ============================================================
# RANDOM BASIS (same as Fortran)
# ============================================================

def generate_random_unit_vector():

    z = 2*np.random.rand() - 1
    phi = 2*np.pi*np.random.rand()

    r = np.sqrt(1-z*z)

    return np.array([
        r*np.cos(phi),
        r*np.sin(phi),
        z
    ])



def generate_random_basis():

    # random axis
    axis = generate_random_unit_vector()

    # random spin angle
    psi = 2*np.pi*np.random.rand()


    # reference vector
    ref = np.array([0.0,0.0,1.0])


    if abs(np.dot(axis,ref)) > 1.0-1e-12:
        ref = np.array([1.0,0.0,0.0])


    # e1 = REF x AXIS
    e1 = np.cross(ref,axis)
    e1 /= np.linalg.norm(e1)


    # e2 = AXIS x e1
    e2 = np.cross(axis,e1)
    e2 /= np.linalg.norm(e2)


    # random rotation
    meridian = (
        np.cos(psi)*e1
        +
        np.sin(psi)*e2
    )

    meridian /= np.linalg.norm(meridian)


    # azimuth = axis x meridian
    azimuth = np.cross(axis,meridian)
    azimuth /= np.linalg.norm(azimuth)


    return axis, meridian, azimuth



# ============================================================
# Generate FDS directions
# ============================================================

def generate_fds_directions(axis, meridian, azimuth, NRT, NRP):

    directions=[]

    for i in range(NRT):

        theta_low=np.pi*i/NRT
        theta_high=np.pi*(i+1)/NRT

        nphi=NRP[i]

        # use midpoint theta of band
        theta=0.5*(theta_low+theta_high)

        for j in range(nphi):

            phi=2*np.pi*j/nphi

            d = (
                np.cos(theta)*axis
                +
                np.sin(theta)*
                (
                    np.cos(phi)*meridian
                    +
                    np.sin(phi)*azimuth
                )
            )

            directions.append(d)

    return np.array(directions)




def hit_cube(direction,L=20):

    half=L/2

    t=half/np.abs(direction)

    dist=np.min(t)

    return direction*dist


def generate_fds_hits(number_of_rotations):

    hits_by_rotation=[]

    for k in range(number_of_rotations):

        rotation_hits=[]

        # random FDS coordinate system
        axis,meridian,azimuth = generate_random_basis()


        dirs=generate_fds_directions(
            axis,
            meridian,
            azimuth,
            NRT,
            NRP
        )


        for d in dirs:

            rotation_hits.append(
                hit_cube(d)
            )

        hits_by_rotation.append(
            np.array(rotation_hits)
        )

    return hits_by_rotation


# ============================================================
# Draw source sphere at center
# ============================================================

def draw_source_sphere(ax, radius=1.0):

    u, v = np.mgrid[
        0:2*np.pi:40j,
        0:np.pi:20j
    ]

    x = radius*np.cos(u)*np.sin(v)
    y = radius*np.sin(u)*np.sin(v)
    z = radius*np.cos(v)

    ax.plot_surface(
        x, y, z,
        color="deepskyblue",
        alpha=0.25,
        linewidth=0
    )

    # optional sphere outline
    ax.plot_wireframe(
        x, y, z,
        color="blue",
        alpha=0.15,
        linewidth=0.3
    )

# ============================================================
# Draw cube enclosure
# ============================================================

def draw_cube(ax, L=20):

    h = L/2

    # 8 cube vertices
    vertices = np.array([
        [-h,-h,-h],
        [ h,-h,-h],
        [ h, h,-h],
        [-h, h,-h],
        [-h,-h, h],
        [ h,-h, h],
        [ h, h, h],
        [-h, h, h]
    ])


    # cube edges
    edges = [
        (0,1),(1,2),(2,3),(3,0),
        (4,5),(5,6),(6,7),(7,4),
        (0,4),(1,5),(2,6),(3,7)
    ]


    for e in edges:

        pts = vertices[list(e)]

        ax.plot(
            pts[:,0],
            pts[:,1],
            pts[:,2],
            color="black",
            linewidth=1.0
        )



# ============================================================
# 2x2 comparison plot
# ============================================================

rotation_cases = [1, 10, 100, 500]

fig = plt.figure(
    figsize=(14,12)
)


for k, nrot in enumerate(rotation_cases):

    ax = fig.add_subplot(
        2,2,k+1,
        projection="3d"
    )


    # ----------------------------
    # Draw cube
    # ----------------------------
    draw_cube(ax, L=20)


    # ----------------------------
    # Draw source sphere
    # ----------------------------
    draw_source_sphere(
        ax,
        radius=1.0
    )


    # ----------------------------
    # Generate FDS hits
    # ----------------------------
    hits_by_rotation = generate_fds_hits(nrot)


    # ----------------------------
    # Plot wall impacts
    # ----------------------------
    # Give every random rotation its own color within this snapshot.
    # Explicitly use red for rotation 1 and blue for rotation 2; obtain
    # additional distinct colors from a continuous colormap.
    rotation_colors = plt.cm.turbo(
        np.linspace(0.0, 1.0, nrot)
    )
    rotation_colors[0] = plt.matplotlib.colors.to_rgba("red")
    if nrot > 1:
        rotation_colors[1] = plt.matplotlib.colors.to_rgba("blue")

    # Combine every rotation into one scatter collection.  The color array
    # retains the identity of each rotation, while a single vector collection
    # avoids the overhead of creating as many as 500 separate PDF objects.
    all_hits = np.vstack(hits_by_rotation)
    point_colors = np.repeat(
        rotation_colors,
        [len(rotation_hits) for rotation_hits in hits_by_rotation],
        axis=0
    )

    ax.scatter(
        all_hits[:, 0],
        all_hits[:, 1],
        all_hits[:, 2],
        s=3,
        color=point_colors,
        alpha=0.7,
        rasterized=False
    )


    # ----------------------------
    # Formatting
    # ----------------------------
    if (nrot == 1):
        ax.set_title(rf"{REQUESTED_NRA} angles, {nrot} random rotation",fontsize=22)
    else:
        ax.set_title(rf"{REQUESTED_NRA} angles, {nrot} random rotations",fontsize=22)


    ax.set_xlim([-10,10])
    ax.set_ylim([-10,10])
    ax.set_zlim([-10,10])

    # Use 4 m tick spacing and include zero on every axis.
    major_ticks = np.arange(-8, 9, 4)
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    ax.set_zticks(major_ticks)

    ax.set_box_aspect(
        [1,1,1]
    )

    axis_label_style = {
        "fontsize": 22,
        "fontweight": "normal"
    }
    ax.set_xlabel("X (m)", **axis_label_style)
    ax.set_ylabel("Y (m)", **axis_label_style)
    ax.set_zlabel("Z (m)", **axis_label_style)

    # Larger numeric tick labels with normal font weight.
    ax.tick_params(axis="both", which="major", labelsize=18)


    ax.view_init(
        elev=25,
        azim=45
    )

fig.subplots_adjust(
    left=0.02,
    right=0.98,
    bottom=0.04,
    top=0.94,
    wspace=0.08,
    hspace=0.12
)

# 540 Kb PDF.
plt.savefig(
   "../../../Manuals/FDS_Technical_Reference_Guide/FIGURES//Rad_rand_rot_cube_hit.pdf",
   format="PDF"
)

# ## 220 Kb PDF.
# # Rasterize and JPEG-compress the complete figure inside a PDF.  This is the
# # practical way to keep this multicolor point-cloud figure below 150 kB.
# output_path = Path(
#     "../../../Manuals/FDS_Technical_Reference_Guide/FIGURES/"
#     "Rad_rand_rot_cube_hit.pdf"
# )
# temporary_png = output_path.with_name(
#     output_path.stem + "_temporary.png"
# )

# fig.savefig(
#     temporary_png,
#     format="png",
#     dpi=100,
#     facecolor="white"
# )

# with Image.open(temporary_png) as image:
#     image.convert("RGB").save(
#         output_path,
#         format="PDF",
#         resolution=150,
#         quality=70,
#         optimize=True,
#     )

# temporary_png.unlink()

plt.show()
