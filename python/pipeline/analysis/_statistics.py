import BioSimSpace as BSS

from ..utils import *
from ._network import *
from ._analysis import *
from ._plotting import *
from ._convergence import *
from ._dictionaries import *

import cinnabar
from cinnabar import wrangle,plotting,stats

class stats_engines(plotting_engines):


    def __init__(self, analysis_object=None, output_folder=None):
        # inherit the init from other protocol too
        super().__init__(analysis_object, output_folder)

        self._set_statistic_dicts()
    
    @staticmethod
    def available_statistics():

        available_statistics = ['RMSE', 'MUE', 'R2', 'rho','RAE','KTAU']

        return available_statistics
    
    def _set_statistic_dicts(self):

        self.statistics  = stats_engines.available_statistics()

        # for statistics compared to experimental value
        self.statistics_dict_exper = {}
        for statistic in self.statistics:
            self.statistics_dict_exper[statistic] ={}

    def _get_x_y(self, pert_val=None, data_x=None, data_y=None, x=None, y=None, xerr=None, yerr=None):

        # default is getting compared to experimental
        # this should be so can match any combo of a and b based on what is available in values dict

        # have dict to df in plotting and also have match engine to other results - dont wanna use dfs!
                    
        if data_x and data_y:
            
            pv = validate.pert_val(pert_val)

            data_x = self._validate_in_names_list(data_x)
            data_y = self._validate_in_names_list(data_y)
            df = self.freenrg_df_dict[data_x][data_y][pv]
            df = df.dropna()

            x = df[f"freenrg_{data_x}"]
            y = df[f"freenrg_calc"]
            xerr = df[f"err_{data_x}"]
            yerr = df[f"err_calc"]
  
        else:
            try:
                x = x
                y = y
                xerr = xerr
                yerr = yerr
            except:
                print("if not providing data_x and data_y (which should be a name in the self.names_list),\
                      please provide x,y,xerr,yerr values")
        
        return x,y,xerr,yerr

    @staticmethod
    def compute_stats(x=None, y=None, xerr=None, yerr=None, statistic=None):

        statistic = validate.string(statistic).upper()

        if statistic not in stats_engines.available_statistics():
            raise ValueError(f"please use one of the statistics in {stats_engines.available_statistics()}")
        
        # using cinnabar function
        s = stats.bootstrap_statistic(x, y, xerr, yerr, nbootstrap=10000, statistic=statistic)
        values = (s['mle'], s['stderr'])
        # string = f"{statistic}:   {s['mle']:.2f} [95%: {s['low']:.2f}, {s['high']:.2f}] " + "\n"
            
        return values

    def _compute_stats(self, pert_val=None, data_x=None, data_y=None, statistic=None, x=None, y=None, xerr=None, yerr=None):

        # get the x y values from the dictionaries, also validates pert val and engine
        x,y,xerr,yerr = self._get_x_y(pert_val, data_x, data_y, x, y, xerr, yerr)

        values = stats_engines.compute_stats(x, y, xerr, yerr, statistic)
            
        return values

    def compute_statistics(self):

        for pv in ["pert","val"]:
            self.compute_mue(pv, self.engines)
            # TODO all compute functions

    def _compute_base(self, pert_val=None, y=None, x=None, statistic=None):

        # validate from other that it is in names list
        x = self._validate_in_names_list(x)
        y = self._validate_in_names_list(y)
        pert_val = validate.pert_val(pert_val)

        if statistic not in self.statistics:
            raise ValueError(f"please use one of the statistics in {self.statistics}")
        
        values = self._compute_stats(pert_val, data_x=x, data_y=y, statistic=statistic)

        return values

    def compute_mue(self, pert_val=None, y=None, x="experimental"):

        return self._compute_base(pert_val=pert_val, y=y, x=x, statistic="MUE")

    def compute_rmse(self, pert_val=None, y=None, x="experimental"):

        return self._compute_base(pert_val=pert_val, y=y, x=x, statistic="RMSE")

    def compute_r2(self, pert_val=None, y=None, x="experimental"):

        return self._compute_base(pert_val=pert_val, y=y, x=x, statistic="R2")    

    def compute_rho(self, pert_val=None, y=None, x="experimental"):

        return self._compute_base(pert_val=pert_val, y=y, x=x, statistic="rho")

    def compute_rae(self, pert_val=None, y=None, x="experimental"):

        return self._compute_base(pert_val=pert_val, y=y, x=x, statistic="RAE")
    
    def compute_ktau(self, pert_val=None, y=None, x="experimental"):

        return self._compute_base(pert_val=pert_val, y=y, x=x, statistic="KTAU")

        # self.statistics_dict_exper["MUE"][y] = values # TODO sort out how all the dicts for stats work

