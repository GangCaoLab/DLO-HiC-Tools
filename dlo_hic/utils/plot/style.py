import seaborn as sns

DEFAULT_SNS_STYLE = 'ticks'
sns.set_style(DEFAULT_SNS_STYLE)

def sns_despine():
    sns.despine(trim=True)
