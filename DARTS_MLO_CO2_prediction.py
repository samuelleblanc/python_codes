#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Using DARTS to forecast CO2 at MLO
#     Based on <https://towardsdatascience.com/forecasting-atmospheric-co2-concentration-with-python-c4a99e4cf142>
#     Adaptation from <https://github.com/derevirn/co2-forecasting/blob/main/forecasting_co2.ipynb>
# 
# Input:
# 
#     arguments
# 
# Output:
# 
#     Figure and save files
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - DARTS
#     - pandas
#     - write_utils
#     - path_utils
#     - scipy
# 
# Needed Files:
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-03-23
#     Modified:
# 

# # Prepare python environment

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from darts import TimeSeries
from darts.models import *
from darts.metrics import *
from darts.dataprocessing.transformers import Scaler
import logging


# In[2]:


mpl.rcParams['figure.dpi'] = 150
logging.disable(logging.CRITICAL)


# In[3]:


from path_utils import getpath
get_ipython().magic(u'matplotlib notebook')


# In[4]:


#name = 'ORACLES'
vv = 'v1'
#fp = getpath(name)
name = 'DARTS'


# In[5]:


fp = '/data/sam/DARTS/'


# # Load files

# In[6]:


df = pd.read_csv(fp+'monthly_in_situ_co2_mlo.csv',
                  comment = '"', header = [0,1,2], na_values = '-99.99')

cols = [' '.join(col).replace(' ', '') for col in df.columns.values]
df.set_axis(cols, axis = 1, inplace = True)

# Converting Excel date format to datetime
# and setting as dataframe index
df['datetime'] = pd.to_datetime(df['DateExcel'], origin = '1899-12-30', unit = 'D')

df.set_index('datetime', inplace = True)

df = df[['CO2filled[ppm]']]
df.rename(columns = {'CO2filled[ppm]': 'CO2'}, inplace = True)
df.dropna(inplace = True)
df = df.resample('M').sum()


# # Plot out data

# In[7]:


df.plot(figsize=(8,5))
plt.title('Monthly CO2 Concentration (ppm)')

plt.show()


# In[8]:



mpl.rcParams['figure.figsize'] = (4, 3)

result = seasonal_decompose(df)
result.plot()

plt.show()


# In[9]:


fig, ax = plt.subplots(figsize = (4,2))

plot_acf(df, ax = ax)

plt.show()


# In[10]:



series = TimeSeries.from_dataframe(df)

start = pd.Timestamp('123115')
df_metrics = pd.DataFrame()

def plot_backtest(series, forecast, model_name):
    idx = -144
    series[idx:].plot(label='Actual Values')
    forecast[idx:].plot(label= 'Forecast')
    plt.title(model_name)
    plt.show()
    
def print_metrics(series, forecast, model_name):
    mae_ = mae(series, forecast)
    rmse_ = rmse(series, forecast)
    mape_ = mape(series, forecast)
    smape_ = smape(series, forecast)
    r2_score_ = r2_score(series, forecast)
    
    dict_ = {'MAE': mae_, 'RMSE': rmse_,
             'MAPE': mape_, 'SMAPE': smape_, 
             'R2': r2_score_}
    
    df = pd.DataFrame(dict_, index = [model_name])
    
    return(df.round(decimals = 2))      


# In[11]:



model = NaiveSeasonal(K = 12)
model_name = 'Naive Seasonal'

plt.figure(figsize = (8, 5))

forecast = model.historical_forecasts(series, start=start, forecast_horizon=12, verbose=True)
plot_backtest(series, forecast, model_name)
df_naive = print_metrics(series, forecast, model_name)
df_metrics = df_metrics.append(df_naive)

plt.show()
df_naive


# In[12]:


model = ExponentialSmoothing(seasonal_periods = 12)
model_name = 'Exponential Smoothing'

plt.figure(figsize = (8, 5))

forecast = model.historical_forecasts(series, start=start, forecast_horizon=12, verbose=True)
plot_backtest(series, forecast, model_name)
df_exp = print_metrics(series, forecast, model_name)
df_metrics = df_metrics.append(df_exp)

plt.show()
df_exp


# In[13]:



model = LinearRegressionModel(lags = 12)
model_name = 'Linear Regression'

plt.figure(figsize = (8, 5))

forecast = model.historical_forecasts(series, start=start, forecast_horizon=12, verbose=True)
plot_backtest(series, forecast, model_name)
df_lr = print_metrics(series, forecast, model_name)
df_metrics = df_metrics.append(df_lr)

plt.show()
df_lr


# In[14]:



model = TCNModel(
    input_chunk_length=24,
    output_chunk_length=12,
    n_epochs=100,
    dropout=0.1,
    dilation_base=3,
    weight_norm=True,
    kernel_size=5,
    num_filters=3,
    random_state=0,
)

model_name = 'TCN'

plt.figure(figsize = (8, 5))

scaler = Scaler()
scaled_series = scaler.fit_transform(series)
forecast = model.historical_forecasts(scaled_series, start=start,
                                      forecast_horizon=12, verbose=True)
plot_backtest(series, scaler.inverse_transform(forecast), model_name)
df_dl = print_metrics(series, scaler.inverse_transform(forecast), model_name)
df_metrics = df_metrics.append(df_dl)

plt.show()
df_dl


# In[15]:


model = ExponentialSmoothing(seasonal_periods = 12)
model_name = 'Exponential Smoothing'

model.fit(series)
forecast = model.predict(12)

plot_backtest(series, forecast, model_name)
print(forecast.pd_dataframe())

