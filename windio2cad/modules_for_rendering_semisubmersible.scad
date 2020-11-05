/*** Draw the various parts of the platform ***/

$fn = 100;

union() {
    z_column(main_column_joint1, main_column_joint2, main_column_outer_grid, main_column_outer_values);
    z_column(column1_joint1, column1_joint2, column1_outer_grid, column1_outer_values);
    z_column(column2_joint1, column2_joint2, column2_outer_grid, column2_outer_values);
    z_column(column3_joint1, column3_joint2, column3_outer_grid, column3_outer_values);

    pontoon(column1_joint1, column1_joint2, column2_joint1, column1_joint2, delta_pontoon_lower12_joint1, delta_pontoon_lower12_joint2, delta_pontoon_lower12_outer_grid, delta_pontoon_lower12_outer_values);
    pontoon(column1_joint1, column1_joint2, column3_joint1, column3_joint2, delta_pontoon_lower31_joint1, delta_pontoon_lower31_joint2, delta_pontoon_lower31_outer_grid, delta_pontoon_lower31_outer_values);
    pontoon(column2_joint1, column2_joint2, column3_joint1, column3_joint2, delta_pontoon_lower23_joint1, delta_pontoon_lower23_joint2, delta_pontoon_lower23_outer_grid, delta_pontoon_lower23_outer_values);

    pontoon(column1_joint1, column1_joint2, column2_joint1, column1_joint2, delta_pontoon_upper12_joint1, delta_pontoon_upper12_joint2, delta_pontoon_upper12_outer_grid, delta_pontoon_upper12_outer_values);
    pontoon(column1_joint1, column1_joint2, column3_joint1, column3_joint2, delta_pontoon_upper31_joint1, delta_pontoon_upper31_joint2, delta_pontoon_upper31_outer_grid, delta_pontoon_upper31_outer_values);
    pontoon(column2_joint1, column2_joint2, column3_joint1, column3_joint2, delta_pontoon_upper23_joint1, delta_pontoon_upper23_joint2, delta_pontoon_upper23_outer_grid, delta_pontoon_upper23_outer_values);

    pontoon(main_column_joint1, main_column_joint2, column1_joint1, column1_joint2, Y_pontoon_upper1_joint1, Y_pontoon_upper1_joint2, Y_pontoon_upper1_outer_grid, Y_pontoon_upper1_outer_values);
    pontoon(main_column_joint1, main_column_joint2, column2_joint1, column2_joint2, Y_pontoon_upper2_joint1, Y_pontoon_upper2_joint2, Y_pontoon_upper2_outer_grid, Y_pontoon_upper2_outer_values);
    pontoon(main_column_joint1, main_column_joint2, column3_joint1, column3_joint2, Y_pontoon_upper3_joint1, Y_pontoon_upper3_joint2, Y_pontoon_upper3_outer_grid, Y_pontoon_upper3_outer_values);

    pontoon(main_column_joint1, main_column_joint2, column1_joint1, column1_joint2, Y_pontoon_lower1_joint1, Y_pontoon_lower1_joint2, Y_pontoon_lower1_outer_grid, Y_pontoon_lower1_outer_values);
    pontoon(main_column_joint1, main_column_joint2, column2_joint1, column2_joint2, Y_pontoon_lower2_joint1, Y_pontoon_lower2_joint2, Y_pontoon_lower2_outer_grid, Y_pontoon_lower2_outer_values);
    pontoon(main_column_joint1, main_column_joint2, column3_joint1, column3_joint2, Y_pontoon_lower3_joint1, Y_pontoon_lower3_joint2, Y_pontoon_lower3_outer_grid, Y_pontoon_lower3_outer_values);

    pontoon(main_column_joint1, main_column_joint2, column1_joint1, column1_joint2, cross_pontoon1_joint1, cross_pontoon1_joint2, cross_pontoon1_outer_grid, cross_pontoon1_outer_values);
    pontoon(main_column_joint1, main_column_joint2, column2_joint1, column2_joint2, cross_pontoon2_joint1, cross_pontoon2_joint2, cross_pontoon2_outer_grid, cross_pontoon2_outer_values);
    pontoon(main_column_joint1, main_column_joint2, column3_joint1, column3_joint2, cross_pontoon3_joint1, cross_pontoon3_joint2, cross_pontoon3_outer_grid, cross_pontoon3_outer_values);
}

/*** 
 * Functions that draw the components for a semisubmersible
 * platform.
 ***/

/***
 * rod() creates a cyclinder between any two arbitrary points,
 * unlike the standalone cylinder() object, which creates cyclinders
 * exclusively along the z axis.
 * 
 * The original module is from:
 * http://forum.openscad.org/Rods-between-3D-points-td13104.html
 *
 * But it is modified here allow for cylinders that are parallel
 * to the z-axis buy do not lie on the z-axis. Also, this is
 * enhanced to allow differnt radii on the start and endpoints
 * of the rod.
 *
 * a: Vector of start coordinate
 * b: Vector of end coordinate
 * r1: Radius at start endpoint
 * r2: Radius at end endpoint
 ***/
 
module rod(a, b, r1, r2) {
    dir = b-a;
    h   = norm(dir);

    // If cylinder lies parallel to z axis
    if(dir[0] == 0 && dir[1] == 0) {
        translate([a[0], a[1], min(a[2], b[2])]) {
            cylinder(r1 = r1, r2 = r2, h=h);
        }
    }
    else {
        w  = dir / h;
        u0 = cross(w, [0,0,1]);
        u  = u0 / norm(u0);
        v0 = cross(w, u);
        v  = v0 / norm(v0);
        multmatrix(m=[[u[0], v[0], w[0], a[0]],
                      [u[1], v[1], w[1], a[1]],
                      [u[2], v[2], w[2], a[2]],
                      [0,    0,    0,    1]])
        cylinder(r1 = r1, r2 = r2, h=h);
    }
}

/***
 * The module z_column() creates a column that is parallel to the
 * z axis, but does not lie along the z axis.
 *
 * keel: The vector of coordinates below the water.
 * freeboard: The vector of coordinates above the water.
 * diam_grid: Relative positions along the z axis of the column,
 *            from 0 to 1
 * diam_values: The diameter at each corresponding grid point.
 ***/

module z_column(keel, freeboard, diam_grid, diam_values) {
    total_height = freeboard[2] - keel[2];

    for (i = [0:len(diam_grid) - 2]) {
        bottom_z = keel[2] + diam_grid[i] * total_height;
        top_z = keel[2] + diam_grid[i + 1] * total_height;
        bottom_radius = diam_values[i] / 2.0;
        top_radius = diam_values[i + 1] / 2.0;
        bottom = [keel[0], keel[1], bottom_z];
        top = [keel[0], keel[1], top_z];
        rod(bottom, top, bottom_radius, top_radius);
    }
}

/***
 * pontoon() creates a pontoon between two columns
 *
 * colA_keel, colA_freeboard: Keel and freeboard of the first column.
 * colB_keel, colB_freeboard: Keel and freeboard of the second column.
 * colA_pontoon, colB_pontoon: The pontoon positions on the first and
 *                             second columns.
 * diam_grid, diam_values: The grid and values for diameters along the
 *                         pontoon.
 ***/

module pontoon(colA_keel, colA_freeboard, colB_keel, colB_freeboard, colA_pontoon, colB_pontoon, diam_grid, diam_values) {
    colA_height = colA_freeboard[2] - colA_keel[2];
    colB_height = colB_freeboard[2] - colB_keel[2];
    
    for (i = [0:len(diam_grid) - 2]) {
        start = [colA_keel[0], colA_keel[1], colA_keel[2] + colA_height * colA_pontoon];
        end = [colB_keel[0], colB_keel[1], colB_keel[2] + colB_height * colB_pontoon];
        
        pontoon_start_x = start[0] + diam_grid[i] * (end[0] - start[0]);
        pontoon_start_y = start[1] + diam_grid[i] * (end[1] - start[1]);
        pontoon_start_z = start[2] + diam_grid[i] * (end[2] - start[2]);
        pontoon_end_x = start[0] + diam_grid[i + 1] * (end[0] - start[0]);
        pontoon_end_y = start[1] + diam_grid[i + 1] * (end[1] - start[1]);
        pontoon_end_z = start[2] + diam_grid[i + 1] * (end[2] - start[2]);
        pontoon_start_radius = diam_values[i] / 2;
        pontoon_end_radius = diam_values[i + 1] / 2;
        
        section_start = [pontoon_start_x, pontoon_start_y, pontoon_start_z];
        section_end = [pontoon_end_x, pontoon_end_y, pontoon_end_z];
        rod(section_start, section_end, pontoon_start_radius, pontoon_end_radius);
    }
}