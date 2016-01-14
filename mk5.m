SetDirectory["/Users/chris/G5/experiments/rhovec"];
m5=Import[".m5.csv"];
Export[".m5h.csv", Table[N[1 - CDF[HypergeometricDistribution[m5[[i]][[3]], m5[[i]][[2]], m5[[i]][[1]] ], m5[[i]][[4]] - 1]], {i, Length[m5]}]];
