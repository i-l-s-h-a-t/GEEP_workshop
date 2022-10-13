import os
import pandas as pd
from darts.tools.plot_darts import *

plots_dir = 'plots/'

# WOPR = Well Oil Production Rate
# WOPT = Well Oil Production Total
# FWPR = Field Water Production Rate
# FGPT = Filed Gas Production Total

# plot PROD and INJ rates, summarized by wells
def plot_total_rates(time_data):
    try:
        ax = plot_total_prod_oil_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Produced gas rate, sm$^3$/day")
        plt.savefig(plots_dir + 'FOPR.png')
    except:
        pass

    try:
        ax = plot_total_prod_water_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Produced water rate, sm$^3$/day")
        plt.savefig(plots_dir + 'FWPR.png')
    except:
        pass

    try:
        ax = plot_total_prod_gas_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Produced gas rate, sm$^3$/day")
        plt.savefig(plots_dir + 'FGPR.png')
    except:
        pass

    try:
        ax = plot_total_inj_gas_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Injected gas rate, sm$^3$/day")
        plt.savefig(plots_dir + 'FGIR.png')
    except:
        pass

    try:
        ax = plot_total_inj_water_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Injected water rate, sm$^3$/day")
        plt.savefig(plots_dir + 'FWIR.png')
    except:
        pass

# plot gas rates summarized by wells
def plot_acc_rates(time_data):
    try:
        ax = plot_acc_prod_water_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Water production, sm3")
        plt.savefig(plots_dir + 'FWPT.png')
    except:
        pass

    try:
        ax = plot_acc_prod_oil_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Oil production, sm3")
        plt.savefig(plots_dir + 'FOPT.png')
    except:
        pass

    try:
        ax = plot_acc_prod_gas_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Gas production, sm3")
        plt.savefig(plots_dir + 'FGPT.png')
    except:
        pass

    try:
        ax = plot_total_inj_water_rate_darts(time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="Water injection, sm3")
        plt.savefig(plots_dir + 'FWIT.png')
    except:
        pass

# plot rates for each well
def plot_rates(time_data, well_fname):
    well_list = get_well_list(well_fname)
    for well_name in well_list:
        try:
            ax = plot_oil_rate_darts(well_name, time_data, style='-', color='b')
            ax.set(xlabel="Days", ylabel="Oil rate, sm3")
            plt.savefig(plots_dir + 'WOPR_' + well_name + '.png')
        except:
            pass

        try:
            ax = plot_water_rate_darts(well_name, time_data, style='-', color='b')
            ax.set(xlabel="Days", ylabel="Water rate, sm3")
            plt.savefig(plots_dir + 'WWPR_' + well_name + '.png')
        except:
            pass

        try:
            ax = plot_gas_rate_darts(well_name, time_data, style='-', color='b')
            ax.set(xlabel="Days", ylabel="Gas rate, sm3")
            plt.savefig(plots_dir + 'WGPR_' + well_name + '.png')
        except:
            pass

# plot temperature for each well
def plot_temperature(time_data, well_fname):
    try:
        well_list = get_well_list(well_fname)
        for well_name in well_list:
            ax = plot_temp_darts(well_name, time_data, style='-', color='b')
            ax.set(xlabel="Days", ylabel="Temperature, K")
            plt.savefig(plots_dir + 'TEMP_' + well_name + '.png')
    except:
        pass

# plot bottom hole pressure for each well
def plot_bhp(time_data, well_fname):
    well_list = get_well_list(well_fname)
    for well_name in well_list:
        ax = plot_bhp_darts(well_name, time_data, style='-', color='b')
        ax.set(xlabel="Days", ylabel="BHP, bar")
        # format the labels
        current_values = ax.get_yticks()
        ax.set_yticklabels(['{:.1f}'.format(x) for x in current_values])
        plt.savefig(plots_dir + 'BHP_' + well_name + '.png')


# read text file with COMPDAT and return the list with well names
def get_well_list(filename):
    if filename is None:
        return
    well_list = set()
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
                                prev_well_name = CompDat[0]
                            well_list.add(CompDat[0])
                        if len(CompDat) != 0 and '/' == CompDat[0]:
                            keep_reading = False
                            break
    return well_list


def plot_wells(well_fname, pkl_file='darts_time_data.pkl'):
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    time_data_orig = pd.read_pickle(pkl_file)
    time_data_report_orig = pd.read_pickle("darts_time_data_report.pkl")

    cut_first_steps = 10
    if cut_first_steps > 0:
        # cut few first timesteps to make the plots better
        time_data = pd.DataFrame()
        for col in time_data_orig.columns:
            time_data[col] = time_data_orig[col][10:]
        time_data_report = pd.DataFrame()
        for col in time_data_report_orig.columns:
            time_data_report[col] = time_data_report_orig[col][cut_first_steps:]
    else:
        time_data = time_data_orig
        time_data_report = time_data_report_orig

    plot_total_rates(time_data)

    plot_acc_rates(time_data_report)

    plot_temperature(time_data, well_fname)

    plot_bhp(time_data, well_fname)

    plot_rates(time_data, well_fname)

##########################################################

if __name__ == '__main__':
    plot_wells(well_fname='wells_case_3_sector.inc')