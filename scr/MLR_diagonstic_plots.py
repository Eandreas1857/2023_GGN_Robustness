from statsmodels.graphics.gofplots import ProbPlot
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.outliers_influence import variance_inflation_factor



def graph(formula, x_range, label=None):
    """
    Helper function for plotting cook's distance lines
    """
    x = x_range
    y = formula(x)
    plt.plot(x, y, label=label, lw=1, ls='--', color='grey')


def diagnostic_plots(X, y, keys, fs, labels, model_fit=None):
  """
  Function to reproduce the 4 base plots of an OLS model in R.

  ---
  Inputs:

  X: A numpy array or pandas dataframe of the features to use in building the linear regression model

  y: A numpy array or pandas series/dataframe of the target variable of the linear regression model

  fs: figure size tuple (n,m)

  labels: list of network labels  

  model_fit [optional]: a statsmodel.api.OLS model after regressing y on X. If not provided, will be
                        generated from X, y
  """
  plt.rc('font', size=12)
  plt.rc('figure', titlesize=18)
  plt.rc('axes', labelsize=15)
  plt.rc('axes', titlesize=18)
  plt.rc("figure", figsize=fs)

  df_x = pd.DataFrame(X, columns = keys)
  df_y = pd.DataFrame(y, columns = ['score'])

  if model_fit == None:
      model_fit = sm.OLS(df_y, sm.add_constant(df_x)).fit()

  
  # create dataframe from X, y for easier plot handling
  dataframe = pd.concat([df_x, df_y], axis=1)

  # model values
  model_fitted_y = model_fit.fittedvalues
  # model residuals
  model_residuals = model_fit.resid
  # normalized residuals
  model_norm_residuals = model_fit.get_influence().resid_studentized_internal
  # absolute squared normalized residuals
  model_norm_residuals_abs_sqrt = np.sqrt(np.abs(model_norm_residuals))
  # absolute residuals
  model_abs_resid = np.abs(model_residuals)
  # leverage, from statsmodels internals
  model_leverage = model_fit.get_influence().hat_matrix_diag
  # cook's distance, from statsmodels internals
  model_cooks = model_fit.get_influence().cooks_distance[0]


  plot_lm_1 = plt.figure()
  
  plot_lm_1.axes[0] = sns.residplot(model_fitted_y, dataframe.columns[-1], data=dataframe,
                            lowess=True,
                            scatter_kws={'color': 'b', 'alpha': 0.5},
                            line_kws={'color': 'red', 'lw': 1, 'alpha': 0.8})

  plot_lm_1.axes[0].set_title('Residuals vs Fitted')
  plot_lm_1.axes[0].set_xlabel('Fitted values')
  plot_lm_1.axes[0].set_ylabel('Residuals');


  with plt.rc_context():
    QQ = ProbPlot(model_norm_residuals)
    plot_lm_2 = QQ.qqplot(line='45', alpha=0.5, color='#4C72B0', lw=1)

  plot_lm_2.axes[0].set_title('Normal Q-Q')
  plot_lm_2.axes[0].set_xlabel('Theoretical Quantiles')
  plot_lm_2.axes[0].set_ylabel('Standardized Residuals');
  # annotations
  abs_norm_resid = np.flip(np.argsort(np.abs(model_norm_residuals)), 0)
  abs_norm_resid_top_3 = abs_norm_resid[:1]


  plot_lm_3 = plt.figure()
  plt.rcParams['text.usetex'] = True
  plt.scatter(model_fitted_y, model_norm_residuals_abs_sqrt, alpha=0.5, color = 'b');
  sns.regplot(model_fitted_y, model_norm_residuals_abs_sqrt,
              scatter=False,
              ci=False,
              lowess=True,
              line_kws={'color': 'red', 'lw': 1, 'alpha': 0.8});
  plot_lm_3.axes[0].set_title('Scale-Location')
  plot_lm_3.axes[0].set_xlabel('Fitted values')
  plot_lm_3.axes[0].set_ylabel(r'$\sqrt{|\textrm{Standardized Residuals}|}$');

  # annotations
  abs_sq_norm_resid = np.flip(np.argsort(model_norm_residuals_abs_sqrt), 0)
  abs_sq_norm_resid_top_3 = abs_sq_norm_resid[:1]
  

  plot_lm_4 = plt.figure();
  plt.scatter(model_leverage, model_norm_residuals, color = 'b', alpha=0.5);
  sns.regplot(model_leverage, model_norm_residuals,
              scatter=False,
              ci=False,
              lowess=True,
              line_kws={'color': 'red', 'lw': 1, 'alpha': 0.8});
  x_lim = max(model_leverage)+0.01
  plot_lm_4.axes[0].set_xlim(0, x_lim)

  plot_lm_4.axes[0].set_ylim(-10, 10)
  plot_lm_4.axes[0].set_title('Residuals vs Leverage')
  plot_lm_4.axes[0].set_xlabel('Leverage')
  plot_lm_4.axes[0].set_ylabel('Standardized Residuals');

  p = len(model_fit.params) # number of model parameters
  graph(lambda x: np.sqrt((0.5 * p * (1 - x)) / x),
        np.linspace(0.001, x_lim, 50),
        'Cook\'s distance') # 0.5 line
  graph(lambda x: np.sqrt((1 * p * (1 - x)) / x),
        np.linspace(0.001, x_lim, 50)) # 1 line
  plot_lm_4.legend(loc = 'upper left', bbox_to_anchor=(0.125, .88))
  graph(lambda x: -np.sqrt((0.5 * p * (1 - x)) / x),
        np.linspace(0.001, x_lim, 50),
        'Cook\'s distance') # 0.5 line
  graph(lambda x: -np.sqrt((1 * p * (1 - x)) / x),
        np.linspace(0.001, x_lim, 50)) # 1 line
  


def compute_vif(df, considered_features):
    
    X = df[considered_features]
    X['intercept'] = 1
    vif = pd.DataFrame()
    vif["Variable"] = X.columns
    vif["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    vif = vif[vif['Variable']!='intercept']
    return vif

  
