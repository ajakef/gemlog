# Missing-GPS Gem Data Pre-Processing Workflow

This sample dataset and notebook show a *non-typical* data conversion: dealing with Gem data where GPS information is missing or insufficient for obtaining precise sample timing. This may be caused by recording where GPS signal is weak (e.g., indoors or underground) or by instrument malfunction.

For ordinary data conversion, where you have no reason to suspect GPS data is missing, use the `demo_conversion` workflow.

Warning: Frequent GPS fixes are essential for precise timekeeping, so sample times in the data converted with this method should NOT be considered precise. Typical clock drifts can be around 10-30 ppm. gemconvert is restrictive when it comes to timing accuracy: it would rather refuse to convert a dataset with insufficient GPS data, than to convert the data and risk the timing being inaccurate. This workflow can bypass gemlog's restrictiveness, but *you should not do it if precise sample timing matters in your application* (e.g., array processing). 