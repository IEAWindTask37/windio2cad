import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator as spline
import geometry_tools as geom
import openpyscad as ops

try:
    import ruamel_yaml as yaml

    loader = yaml.Loader
except:
    try:
        import ruamel.yaml as yaml

        loader = yaml.Loader
    except:
        import yaml

        loader = yaml.FullLoader


def load_yaml(fname_input):
    with open(fname_input, "r") as f:
        input_yaml = yaml.load(f, Loader=loader)
    outer_shape = input_yaml["components"]["blade"]["outer_shape_bem"]
    airfoils = input_yaml["airfoils"]
    return (outer_shape, airfoils)


def myinterp(xi, x, f):
    myspline = spline(x, f)
    return myspline(xi)


def generate_lofted(outer_shape, airfoils, n_span_min=100, n_xy=400):
    # Use yaml grid points and others that we add
    r_span = np.unique(
        np.r_[
            np.linspace(0.0, 1.0, n_span_min),
            outer_shape["chord"]["grid"],
            outer_shape["twist"]["grid"],
            outer_shape["pitch_axis"]["grid"],
            outer_shape["reference_axis"]["x"]["grid"],
            outer_shape["reference_axis"]["y"]["grid"],
            outer_shape["reference_axis"]["z"]["grid"],
        ]
    )
    n_span = len(r_span)

    # Read in blade spanwise geometry values and put on common grid
    chord = myinterp(
        r_span, outer_shape["chord"]["grid"], outer_shape["chord"]["values"]
    )
    twist = myinterp(
        r_span, outer_shape["twist"]["grid"], outer_shape["twist"]["values"]
    )
    pitch_axis = myinterp(
        r_span, outer_shape["pitch_axis"]["grid"], outer_shape["pitch_axis"]["values"]
    )
    ref_axis = np.c_[
        myinterp(
            r_span,
            outer_shape["reference_axis"]["x"]["grid"],
            outer_shape["reference_axis"]["x"]["values"],
        ),
        myinterp(
            r_span,
            outer_shape["reference_axis"]["y"]["grid"],
            outer_shape["reference_axis"]["y"]["values"],
        ),
        myinterp(
            r_span,
            outer_shape["reference_axis"]["z"]["grid"],
            outer_shape["reference_axis"]["z"]["values"],
        ),
    ]

    # Get airfoil names and thicknesses
    af_position = outer_shape["airfoil_position"]["grid"]
    af_used = outer_shape["airfoil_position"]["labels"]
    n_af_span = len(af_position)
    n_af = len(airfoils)
    name = n_af * [""]
    r_thick = np.zeros(n_af)
    for i in range(n_af):
        name[i] = airfoils[i]["name"]
        r_thick[i] = airfoils[i]["relative_thickness"]

    # Create common airfoil coordinates grid
    coord_xy = np.zeros((n_af, n_xy, 2))
    for i in range(n_af):
        points = np.c_[airfoils[i]["coordinates"]["x"], airfoils[i]["coordinates"]["y"]]

        # Check that airfoil points are declared from the TE suction side to TE pressure side
        idx_le = np.argmin(points[:, 0])
        if np.mean(points[:idx_le, 1]) > 0.0:
            points = np.flip(points, axis=0)

        # Remap points using class AirfoilShape
        af = geom.AirfoilShape(points=points)
        af.redistribute(n_xy, even=False, dLE=True)
        af_points = af.points

        # Add trailing edge point if not defined
        if [1, 0] not in af_points.tolist():
            af_points[:, 0] -= af_points[np.argmin(af_points[:, 0]), 0]
        c = max(af_points[:, 0]) - min(af_points[:, 0])
        af_points[:, :] /= c

        coord_xy[i, :, :] = af_points

    # Reconstruct the blade relative thickness along span with a pchip
    r_thick_used = np.zeros(n_af_span)
    coord_xy_used = np.zeros((n_af_span, n_xy, 2))
    coord_xy_interp = np.zeros((n_span, n_xy, 2))
    coord_xy_dim = np.zeros((n_span, n_xy, 2))

    for i in range(n_af_span):
        for j in range(n_af):
            if af_used[i] == name[j]:
                r_thick_used[i] = r_thick[j]
                coord_xy_used[i, :, :] = coord_xy[j, :, :]

    r_thick_interp = myinterp(r_span, af_position, r_thick_used)

    # Spanwise interpolation of the profile coordinates with a pchip
    r_thick_unique, indices = np.unique(r_thick_used, return_index=True)
    coord_xy_interp = np.flip(
        myinterp(np.flip(r_thick_interp), r_thick_unique, coord_xy_used[indices, :, :]),
        axis=0,
    )
    for i in range(n_span):
        # Correction to move the leading edge (min x point) to (0,0)
        af_le = coord_xy_interp[i, np.argmin(coord_xy_interp[i, :, 0]), :]
        coord_xy_interp[i, :, 0] -= af_le[0]
        coord_xy_interp[i, :, 1] -= af_le[1]
        c = max(coord_xy_interp[i, :, 0]) - min(coord_xy_interp[i, :, 0])
        coord_xy_interp[i, :, :] /= c
        # If the rel thickness is smaller than 0.4 apply a trailing ege smoothing step
        if r_thick_interp[i] < 0.4:
            coord_xy_interp[i, :, :] = geom.trailing_edge_smoothing(
                coord_xy_interp[i, :, :]
            )

    # Offset by pitch axis and scale for chord
    coord_xy_dim = coord_xy_interp.copy()
    coord_xy_dim[:, :, 0] -= pitch_axis[:, np.newaxis]
    coord_xy_dim = coord_xy_dim * chord[:, np.newaxis, np.newaxis]

    # Rotate to twist angle
    coord_xy_dim_twisted = np.zeros(coord_xy_interp.shape)
    for i in range(n_span):
        x = coord_xy_dim[i, :, 0]
        y = coord_xy_dim[i, :, 1]
        coord_xy_dim_twisted[i, :, 0] = x * np.cos(twist[i]) - y * np.sin(twist[i])
        coord_xy_dim_twisted[i, :, 1] = y * np.cos(twist[i]) + x * np.sin(twist[i])

    # Assemble lofted shape along reference axis
    lofted_shape = np.zeros((n_span, n_xy, 3))
    for i in range(n_span):
        for j in range(n_xy):
            lofted_shape[i, j, :] = (
                np.r_[coord_xy_dim_twisted[i, j, 1], coord_xy_dim_twisted[i, j, 0], 0.0]
                + ref_axis[i, :]
            )

    return lofted_shape


def set_axes_equal(ax):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().

    https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def write_openscad(lofted_shape):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    # ax.set_aspect('equal')

    n_span = lofted_shape.shape[0]

    extrusions = ops.Union()

    for i in range(n_span - 1):
        j = i + 1
        bottom = lofted_shape[i, 0, 2]
        top = lofted_shape[j, 0, 2]
        z_height = top - bottom

        # Create polygon / polyhedron of the airfoil shape in lofted_shape[k,:,:].
        # Use the hull(){} command to connect them

        points = [[row[0], row[1]] for row in lofted_shape[i, :, :]]
        extruded_section = ops.Polygon(points=points).linear_extrude(height=z_height).translate([0, 0, bottom])
        extrusions.append(extruded_section)

        # Graph version
        # ax.plot(
        #     lofted_shape[k, :, 0], lofted_shape[k, :, 1], lofted_shape[k, :, 2], "b"
        # )

    set_axes_equal(ax)
    plt.show()

    extrusions.write("foo.scad")


if __name__ == "__main__":
    fname = "IEA-15-240-RWT_FineGrid.yaml"
    outer_shape, airfoils = load_yaml(fname)
    lofted3d = generate_lofted(outer_shape, airfoils)
    write_openscad(lofted3d)
