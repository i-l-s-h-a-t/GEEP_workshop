from darts.engines import value_vector, redirect_darts_output, sim_params
from model_deadoil import Model
import pandas as pd

import matplotlib.pyplot as plt

redirect_darts_output('out.log')

m = Model()
m.init()

# equilibrate (run simulation for 0.1 days)
m.run_python(0.1)
m.export_pro_vtk()

#exit()

m.set_wells()
# m.run_python(100, restart_dt=m.params.first_ts)
# m.export_pro_vtk()

#  simulation
n_steps = 10
time_step = 365  # days

for a in range(n_steps):
    m.run_python(time_step)
    m.export_pro_vtk()

m.print_timers()
m.print_stat()

m.wells4ParaView('well_gen.txt')

time_data = pd.DataFrame.from_dict(m.physics.engine.time_data)
time_data.to_pickle("darts_time_data.pkl")

writer = pd.ExcelWriter('time_data.xlsx')
time_data.to_excel(writer, 'Sheet1')
writer.save()
