---
sort: 2
---

# Preprocessing

We will now preprocess our raw data which will reduce noise from non-biological signal as well as filter out low quality cells. Part of the preprocessing stage is also tranforming the data such as log-normalizing or scaling data to a given distribution. This will help give more meaningful results and help remove outliers. These will all be explained in this section.

## Highest expressed genes

First let's which genes are the most expressed across our dataset:

<img style='border-radius:10px; box-shadow: 5px 5px 10px rgb(0 0 0 / 0.5);' alt='page screenshot' src='https://raw.githubusercontent.com/nuwa-genomics/Nuwa/main/docs/assets/images/screenshots/clustering_tutorial/highest_expr.png'>

```tip
Hover over the box plots to see their median, q1 & q3, lowest and highest counts.
```

We can see the most expressed gene by far with a median count of around 3.75, followed by some mitochondrial genes.

```note
## Malat1
Malat1 gene can be a result of unwanted technical noise so can often be removed from the dataset.
```

