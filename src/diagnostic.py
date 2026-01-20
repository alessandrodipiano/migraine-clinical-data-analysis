import pandas as pd 

def pre_stat(df):
    missing_fraction = df.isna().mean()
    missing_count = df.isna().sum()
    unique_value = df.nunique(dropna=True)

    result = pd.concat(
        [missing_fraction, missing_count, unique_value],
        axis=1,
        keys=["fraction_missing", "count_missing", "unique_value"]
    )

    return result




def missingness_matrix(df):
    cols = df.columns[df.isna().any(axis=0)]
    return df[cols].isna().corr().round(2)


def numerical_check(df, low_q=0.01, high_q=0.99, exclude= None , min_uniq=2):

    if exclude is None:
        exclude=set()

    

    issues={}
    for col in df.select_dtypes(include='number'):

        if col in exclude:
            continue
        
        s = df[col].dropna()
        if s.empty:
            continue
        
        if s.nunique() < min_uniq:
            continue

        

        q_low, q_high= s.quantile([low_q,high_q])
        extreme=df[(df[col] < q_low) | (df[col]>q_high)][col]

        if not extreme.empty:
            issues[col] = {
                "min": float(s.min()),
                "max": float(s.max()),
                "q01": float(q_low),
                "q99": float(q_high),
                "n_extreme": int(extreme.count()),
                "extreme_values": extreme 
            }
    return issues
