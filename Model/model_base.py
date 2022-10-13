from darts.models.reservoirs.struct_reservoir import StructReservoir
from darts.models.darts_model import DartsModel
import numpy as np
from darts.tools.keyword_file_tools import load_single_keyword, save_few_keywords
import os


class BaseModel(DartsModel):
    def __init__(self, n_points=1000):
        # call base class constructor
        super().__init__()
        self.n_points = n_points
        # measure time spend on reading/initialization
        self.timer.node["initialization"].start()

        #case = 'case_3'
        case = 'case_3_sector'

        if case == 'case_3':
            self.nx = 82
            self.ny = 75
            self.nz = 22
            self.prop_filename = self.grid_filename = 'case_3.grdecl'
            self.width_filename = 'width_case_3.grdecl'
            self.well_perf_filename = 'wells_case_3.inc'
        elif case == 'case_3_sector':
            self.nx = 67
            self.ny = 51
            self.nz = 10
            self.prop_filename = self.grid_filename = 'case_3_sector.grdecl'
            self.width_filename = 'width_case_3_sector.grdecl'
            self.well_perf_filename = 'wells_case_3_sector.inc'
        else:  # original case
            self.nx = 81
            self.ny = 58
            self.nz = 20
            self.grid_filename = 'grid.grdecl'
            self.prop_filename = 'reservoir.in'
            self.width_filename = 'width.in'
            self.well_perf_filename = 'WELLS.INC'  # file with COMPDAT keyword (used columns: well_name, I, J, K1, K2)

        # create reservoir from UNISIM - 20 layers (81*58*20, Corner-point grid)
        self.permx = load_single_keyword(self.prop_filename, 'PERMX')
        self.permy = load_single_keyword(self.prop_filename, 'PERMY')
        self.permz = load_single_keyword(self.prop_filename, 'PERMZ')
        self.poro  = load_single_keyword(self.prop_filename, 'PORO')

        #self.permx *= 100
        #self.permy *= 100
        #self.permz *= 100

        self.is_CPG = True  # re-calculate dx, dy and dz from CPG grid and write them to width.in
        #self.is_CPG = False # read dx, dy and dz from width.in for faster mesh initialization

        if self.is_CPG is False:
            print('Reading dx, dy and dz specifications...')
            self.dx = load_single_keyword(self.width_filename, 'DX')
            self.dy = load_single_keyword(self.width_filename, 'DY')
            self.dz = load_single_keyword(self.width_filename, 'DZ')
            self.depth = load_single_keyword(self.width_filename, 'DEPTH')
        else: # read CPG (COORD and ZCORN arrays)
            self.dx = self.dy = self.dz = self.depth = 0

        self.coord = load_single_keyword(self.grid_filename, 'COORD')
        self.zcorn = load_single_keyword(self.grid_filename, 'ZCORN')

        # Import other properties from files
        self.actnum = load_single_keyword(self.grid_filename, 'ACTNUM')

        self.reservoir = StructReservoir(self.timer, nx=self.nx, ny=self.ny, nz=self.nz, dx=self.dx, dy=self.dy, dz=self.dz,
                                         permx=self.permx, permy=self.permy, permz=self.permz, poro=self.poro,
                                         depth=self.depth, actnum=self.actnum, coord=self.coord, zcorn=self.zcorn,
                                         is_cpg=self.is_CPG)

        poro = np.array(self.reservoir.mesh.poro, copy=False)
        poro[poro == 0.0] = 1.E-4

        self.reservoir.set_boundary_volume(yz_minus=1e15, yz_plus=1e15, xz_minus=1e15,
                                           xz_plus=1e15, xy_minus=-1, xy_plus=-1)

        if self.is_CPG:
            dx, dy, dz = self.reservoir.get_cell_cpg_widths()
            self.depth = self.reservoir.discretizer.cell_data['center'][:, :, :, 2].flatten(order='F')
            save_few_keywords(self.width_filename, ['DX', 'DY', 'DZ', 'DEPTH'], [dx, dy, dz, self.depth])

        self.read_and_add_perforations(self.well_perf_filename)

        self.timer.node["initialization"].stop()

    def set_initial_conditions(self):
        mesh = self.reservoir.mesh
        """ Uniform Initial conditions """
        # set initial pressure
        pressure = np.array(mesh.pressure, copy=False)
        pressure.fill(self.pressure_ini)

        nc = self.property_container.nc
        nb = mesh.n_blocks
        mesh.composition.resize(nb * (nc - 1))
        # set initial composition
        composition = np.array(mesh.composition, copy=False)
        for c in range(nc - 1):
            composition[c::(nc - 1)] = self.ini_comp[c]


    def read_and_add_perforations(self, filename):
        if filename is None:
            return

        well_dia = 0.152
        well_rad = well_dia / 2

        keep_reading = True
        prev_well_name = ''
        with open(filename) as f:
            while keep_reading:
                buff = f.readline()
                if 'COMPDAT' in buff:
                    while True:  # be careful here
                        buff = f.readline()
                        if len(buff) != 0:
                            CompDat = buff.split()

                            if len(CompDat) != 0 and '/' != CompDat[0]:  # skip the empty line and '/' line
                                # define well
                                if CompDat[0] == prev_well_name:
                                    pass
                                else:
                                    self.reservoir.add_well(CompDat[0], wellbore_diameter=well_dia)
                                    prev_well_name = CompDat[0]
                                # define perforation
                                for i in range(int(CompDat[3]), int(CompDat[4]) + 1):
                                    self.reservoir.add_perforation(self.reservoir.wells[-1],
                                                                   int(CompDat[1]), int(CompDat[2]), i,
                                                                   well_radius=well_rad,
                                                                   multi_segment=False, verbose=False)

                            if len(CompDat) != 0 and '/' == CompDat[0]:
                                keep_reading = False
                                break

    def wells4ParaView(self, filename):
        if not self.is_CPG:
            return
        x_array = self.reservoir.discretizer.cell_data['center'][:, :, :, 0]
        y_array = self.reservoir.discretizer.cell_data['center'][:, :, :, 1]
        z_array = self.reservoir.discretizer.cell_data['center'][:, :, :, 2]
        name = []
        type = []
        x = []
        y = []
        z = []
        keep_reading = True
        with open(self.well_perf_filename) as f:
            while keep_reading:
                buff = f.readline()
                if 'COMPDAT' in buff:
                    while True:  # be careful here
                        buff = f.readline()
                        if len(buff) != 0:
                            welspecs = buff.split()

                            if len(welspecs) != 0 and welspecs[0] != '/' and welspecs[0][:2] != '--':  # skip the empty line and '/' line
                                name += [welspecs[0]]
                                if 'INJ' in welspecs[0]:
                                    type += ['INJ']
                                else:
                                    type += ['PRD']
                                ix = int(welspecs[1]) - 1
                                iy = int(welspecs[2]) - 1
                                iz = 0
                                x += [x_array[ix, iy, iz]]
                                y += [y_array[ix, iy, iz]]
                                z += [z_array[ix, iy, iz]]

                            if len(welspecs) != 0 and welspecs[0] == '/':
                                keep_reading = False
                                break
        f.close()

        def str2file(fp, name_in, list_in):
            fp.write("%s = [" % name_in)
            for item in list_in:
                fp.write("\'%s\', " % item)
            fp.write("]\n")

        def num2file(fp, name_in, list_in):
            fp.write("%s = [" % name_in)
            for item in list_in:
                fp.write("%f, " % int(item))
            fp.write("]\n")

        f = open(filename, 'w')
        str2file(f, 'well_list', name)
        str2file(f, 'well_type', type)
        num2file(f, 'well_x', x)
        num2file(f, 'well_y', y)
        num2file(f, 'well_z', z)
        f.close()

    def save_to_grdecl(self, filename, cube, cubename):
        '''
        writes 1d-array of size nx*ny*nz from 1d-array of size n_active_cells, fills with zeros where actnum is zero
        might be useful to load into other software
        '''
        cube_full = np.zeros(self.nx * self.ny * self.nz)
        cube_full[self.reservoir.discretizer.local_to_global[:]] = cube[:]
        save_few_keywords(filename, [cubename], [cube_full])
