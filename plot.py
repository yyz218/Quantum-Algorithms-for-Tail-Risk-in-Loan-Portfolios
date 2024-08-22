import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# plot for all three type of method 

mc = pd.read_csv('loss.csv')
mc['prob'] = mc['count'].cumsum() / mc['count'].sum()

cdf = pd.read_csv('cdf.csv')
prob = pd.read_csv('prob.csv')

plt.figure(figsize=(10, 6))
plt.plot(cdf['h'], cdf['F_h'], color='red', drawstyle='steps-post', linestyle='-', linewidth=2, label='Analytic F(h)')
plt.plot(prob['loss'], prob['prob'], color='orange', drawstyle='steps-post', linestyle='-', linewidth=2, label='K-bit simulation')
plt.plot(mc['loss'], mc['prob'], color='blue', drawstyle='steps-post', linestyle='-', linewidth=2, label='MC simulation')
plt.title('Cumulative Probability Distribution')
plt.xlabel('Loss')
plt.ylabel('Cumulative Probability')
plt.legend()
plt.show()

