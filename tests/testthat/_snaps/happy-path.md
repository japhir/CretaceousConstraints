# bandpass_filter has a happy path

    Code
      bandpass_filter(dat, tibble::tibble(flow = 1, fhigh = 2, target = "group"), x = a,
      y = c, window = 0)
    Output
      # A tibble: 10 x 6
          flow fhigh target     a       c filter
         <dbl> <dbl> <chr>  <int>   <dbl>  <dbl>
       1     1     2 group      1 -0.560  0.0746
       2     1     2 group      2 -0.230  0.0746
       3     1     2 group      3  1.56   0.0746
       4     1     2 group      4  0.0705 0.0746
       5     1     2 group      5  0.129  0.0746
       6     1     2 group      6  1.72   0.0746
       7     1     2 group      7  0.461  0.0746
       8     1     2 group      8 -1.27   0.0746
       9     1     2 group      9 -0.687  0.0746
      10     1     2 group     10 -0.446  0.0746

