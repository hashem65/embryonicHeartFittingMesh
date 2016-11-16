gfx read node UndeformedGeometry.part0.exnode time 0
gfx read node DeformedGeometry1.part0.exnode time 1
#gfx read node DeformedGeometry2.part0.exnode time 2

gfx read elem UndeformedGeometry.part0.exelem

gfx def faces egroup FittingRegion
gfx modify g_element "/" lines coordinate Coordinate tessellation default LOCAL native_discretization NONE select_on material default selected_material default_selected;
gfx modify g_element "/" surfaces coordinate Coordinate exterior tessellation default LOCAL native_discretization NONE select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element "/" node_points coordinate Coordinate LOCAL glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected;

gfx create spectrum distance_spectrum;
gfx modify spectrum distance_spectrum clear overwrite_colour;
gfx modify spectrum distance_spectrum linear range 0.0 300.0 extend_above extend_below white_to_red colour_range 0 1 component 1;

gfx read data DataPoints.part0.exdata
#gfx read data DataPoints.part1.exdata
gfx modify g_element DataPoints general clear;
gfx modify g_element DataPoints data_points coordinate data_coordinates LOCAL glyph point general size "0.1*0.1*0.1" centre 0,0,0 font default select_on material yellow selected_material default_selected;
gfx modify g_element DataPoints lines domain_mesh1d coordinate data_coordinates LOCAL line line_base_size 0 select_on material default data data_distance spectrum default selected_material default_selected render_shaded;
gfx modify g_element DataPoints points domain_datapoints coordinate data_coordinates tessellation default_points LOCAL glyph arrow_solid size "0.5*0.5*0.5" offset 0,0,0 font default orientation data_error scale_factors "1*0.05*0.05" select_on material gold data data_distance spectrum distance_spectrum selected_material default_selected render_shaded;

gfx cre wi
gfx cre time_editor
gfx edit scene
