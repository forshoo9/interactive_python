
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('boxplot_data.csv')
print(df)

fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 5))

ax.grid(axis='y', color='grey', ls=':')
ax.axhline(y=0, color='grey', lw=1)

ax.boxplot(df.EP.dropna().values, positions=[0], widths=0.5, showmeans=True, labels=['EP'])
ax.boxplot(df.LN.dropna().values, positions=[1], widths=0.5, showmeans=True, labels=['LN'])
ax.boxplot(df.PA.dropna().values, positions=[2], widths=0.5, showmeans=True, labels=['PA'])
ax.boxplot(df.NA.dropna().values, positions=[3], widths=0.5, showmeans=True, labels=['NA'])

ax.set_ylabel(r'Mean KR Temperature Anomalies [$\degree$C]'+'\n', fontsize=13)
ax.text(-0.1, 0.5, 'warming', rotation=90, transform=ax.transAxes)
ax.text(-0.1, 0.35, 'cooling', rotation=90, transform=ax.transAxes)

ax.tick_params(axis='x', which='major', labelsize=13)

plt.show()

fig.savefig('fig_06_boxplot.png', dpi=500, bbox_inches='tight')


##===============
##---for seaborn
##===============
#
#import seaborn as sns
#df = df.melt(value_vars=['EP','LN','PA','NA'], var_name='case')
#sns.set_theme(style="ticks")#, palette="pastel")
#ax = sns.boxplot(x="case", y="value", data=df, showmeans=True, linewidth=2.5, whis=1.5)

