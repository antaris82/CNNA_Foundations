# RESULTS: activity-gradient + patch-size star-closure scaling

## Gate

This diagnostic checks whether the strongest live DtN/Record drift is localized near the active growth/completion shell, and whether local patch operator-system closure improves with patch size/level. It remains pre-J and pre-C*.

## Last-shell incremental node-change gradient

This compares the same provenance nodes before and after the last growth shell. It is the cleaner test of your point: newest parent/frontier nodes receive the largest current update; the root receives the smallest shell-normalized marginal update.

```text
node_level count mean_abs_delta mean_rel_delta max_abs_delta
0          1     0.010722       0.0080205      0.010722      
1          3     0.01544        0.0079197      0.015665      
2          9     0.024125       0.010328       0.025001      
3          27    0.042889       0.016185       0.045746      
4          81    0.096499       0.033735       0.10685       
5          243   0.386          0.14122        0.44813       
```

```text
distance_to_frontier count mean_abs_delta mean_rel_delta max_abs_delta
0                    243   0.386          0.14122        0.44813       
1                    81    0.096499       0.033735       0.10685       
2                    27    0.042889       0.016185       0.045746      
3                    9     0.024125       0.010328       0.025001      
4                    3     0.01544        0.0079197      0.015665      
5                    1     0.010722       0.0080205      0.010722      
```

## Activity gradient by age

Age is measured from completion to the final growth level. Age 0 is the active/latest completed shell; root/interior cells have large age.

```text
age count aging_mean handoff_mean child_g_drift desc_load ancestor_env_drift
0   243   0.0013545   0.0063656    0             0         0.027093    
1   81    1.7048      1.7091       0.66903       14.867    0.15591     
2   27    2.0156      2.0195       0.79065       17.57     0.20595     
3   9     2.0422      2.0455       0.80066       17.792    0.21952     
4   3     1.9473      1.9498       0.76285       16.952    0.20294     
5   1     1.7756      1.7765       0.69464       15.437    0.14157     
```

## Activity gradient by cell level

```text
cell_level count aging_mean handoff_mean child_g_drift desc_load
0          1     1.7756      1.7765       0.69464       15.437   
1          3     1.9473      1.9498       0.76285       16.952   
2          9     2.0422      2.0455       0.80066       17.792   
3          27    2.0156      2.0195       0.79065       17.57    
4          81    1.7048      1.7091       0.66903       14.867   
5          243   0.0013545   0.0063656    0             0        
```

## Primary same-suffix patch closure scan

```text
L k count word_dim star_basis star_seed mult comm_rel
4  2  1     29.0     0.0705     0.0026    0.3597  0.2942  
4  3  1     77.0     0.9186     0.5632    0.5024  0.2453  
5  2  1     27.0     0.0091     0.0017    0.2792  0.2904  
5  3  1     73.0     0.8976     0.4862    0.5184  0.2517  
6  2  1     29.0     0.0800     0.0022    0.3770  0.2928  
6  3  1     68.0     0.0485     0.0026    0.4181  0.2921  
```

## Interpretation

The activity table should be read before the closure table: if the live changes concentrate near the active/latest shell while the root changes least, then a local operator-system candidate must be tested on active growth-defined patches, not on an isolated old/root subsystem. A positive local-net/*-limit signal would require star and multiplication residuals to decrease consistently as patch size and level grow. At this finite sampled scale, the scan is only diagnostic.