# Re-run this cell 
import pandas as pd

# Read in the data
schools = pd.read_csv("schools.csv")

# Preview the data
schools.head()

# Highest average math score
high_avg_score = schools[(schools["average_math"] >= 640)]
best_math_schools = high_avg_score.sort_values(["average_math"], ascending=False)
best_math_schools = best_math_schools[["school_name","average_math"]].reset_index(drop=True)
best_math_schools

# Top 10 SAT scores
top_10_schools = schools["total_SAT"] = schools[["average_math", "average_writing", "average_reading"]].sum(axis=1)
top_10_schools = schools.sort_values(["total_SAT"], ascending=False)
top_10_schools = top_10_schools[["school_name", "total_SAT"]].reset_index(drop=True).head(10)

top_10_schools

# Largest standard deviation in the combined SAT score
boroughs = schools.groupby('borough')['total_SAT'].agg(['count', 'mean', 'std']).round(2)
largest_std_dev = boroughs[boroughs['std'] == boroughs['std'].max()]
largest_std_dev = largest_std_dev.rename(columns={"count":"num_schools", "mean":"average_SAT", "std":"std_SAT"})
largest_std_dev.reset_index(inplace=True)

largest_std_dev
