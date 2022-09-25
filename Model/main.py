from darts.engines import value_vector, redirect_darts_output, sim_params
from model_deadoil import Model
import pandas as pd

import matplotlib.pyplot as plt

redirect_darts_output('out.log')
m = Model()
m.init()
m.run_python(0.1)
m.export_pro_vtk()

m.set_wells()
m.run_python(100, restart_dt=m.params.first_ts)
m.export_pro_vtk()

for a in range(10):
    m.run_python(100)
    m.export_pro_vtk()
m.print_timers()
m.print_stat()

m.wells4ParaView('well_gen.txt')

time_data = pd.DataFrame.from_dict(m.physics.engine.time_data)
time_data.to_pickle("darts_time_data.pkl")

writer = pd.ExcelWriter('time_data.xlsx')
time_data.to_excel(writer, 'Sheet1')
writer.save()

# rate plotting
from darts.tools.plot_darts import *

ax1 = plot_total_prod_oil_rate_darts(time_data, style='-', color='b')
ax1.set(xlabel="Days", ylabel="Produced gas rate, sm$^3$/day")

ax3 = plot_total_prod_water_rate_darts(time_data, style='-', color='b')
ax3.set(xlabel="Days", ylabel="Produced water rate, sm$^3$/day")

ax2 = plot_total_inj_gas_rate_darts(time_data, style='-', color='b')
ax2.set(xlabel="Days", ylabel="Injected gas rate, sm$^3$/day")

ax4 = plot_total_inj_water_rate_darts(time_data, style='-', color='b')
ax4.set(xlabel="Days", ylabel="Injected water rate, sm$^3$/day")

# ax5 = plot_temp_darts("PROD005", time_data, style='-', color='b')
# ax5.set(xlabel="Days", ylabel="Temperature, K")

plt.show()
