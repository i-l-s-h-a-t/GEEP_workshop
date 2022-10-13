# case_3_sector
well_list = ['PROD005', 'PROD010', 'PROD012', 'PROD014', 'PROD021', 'INJ005', 'INJ006', 'INJ010', ]
well_type = ['PRD', 'PRD', 'PRD', 'PRD', 'PRD', 'INJ', 'INJ', 'INJ', ]
well_x = [2058.000000, 2197.000000, 2708.000000, 2919.000000, 1765.000000, 1955.000000, 2779.000000, 3134.000000, ]
well_y = [1677.000000, 2052.000000, 1936.000000, 1684.000000, 2034.000000, 1729.000000, 2514.000000, 2309.000000, ]
well_z = [1364.000000, 1373.000000, 1389.000000, 1458.000000, 1392.000000, 1360.000000, 1385.000000, 1396.000000, ]

# case_3


# base


from paraview.simple import *
well_cyl = Cylinder();
SetProperties(well_cyl, Height=1000, Radius=30);
for idx, val in enumerate(well_list):
	# wellbore
	t = Transform(well_cyl);
	t.Transform.Translate = [well_x[idx], well_y[idx], well_z[idx]];
	t.Transform.Rotate = [90, 0, 0];
	dp = GetDisplayProperties(t);
	if (well_type[idx] == 'PRD'):
		dp.DiffuseColor = [1, 0, 0];
	else:	
		dp.DiffuseColor = [0, 0, 1];
	Show(t);
	# well name
	title = a3DText();
	name = well_list[idx];
	SetProperties(title, Text=name);
	tt = Transform(title);
	tt.Transform.Translate = [well_x[idx], well_y[idx], well_z[idx] - 500];
	tt.Transform.Scale = [80, 80, 80];
	tt.Transform.Rotate = [60, 30, 0];
	dp = GetDisplayProperties(tt);
	dp.DiffuseColor = [0, 1, 0];
	Show(tt);


Render();

	
	





