Comparing methods to predict baseline mortality for excess mortality
calculations – unravelling ‘the German puzzle’ and its implications for
spline-regression
================
Tamás Ferenci[^1]

The study is available in [PDF
format](https://raw.githubusercontent.com/tamas-ferenci/MortalityPrediction/main/README.pdf).

## Abstract

Introduction: The World Health Organization presented global excess
mortality estimates for 2020 and 2021 on May 5, 2022, almost immediately
stirring controversy, one point of which was the suspiciously high
estimate for Germany. Later analysis revealed that the reason of this –
in addition to a data preparation issue – was the nature of the
spline-model underlying WHO’s method used for excess mortality
estimation. This paper aims to reproduce the problem using synthetic
datasets, thus allowing the investigation of its sensitivity to
parameters, both of the mortality curve and of the used method, thereby
shedding light on the conditions that gave rise to this error and its
possible remedies.

Material and Methods: A negative binomial model was used with constant
overdispersion, and a mean being composed of three terms: a long-term
change (modelled with a quadratic trend), deterministic seasonality,
modelled with a single harmonic term and random additional peaks during
the winter (flu season) and during the summer (heat waves). Simulated
mortality curves from this model were then analyzed with naive methods
(simple mean of the latest few years, simple linear trend projection
from the latest few years), with the WHO’s method and with the method of
Acosta and Irizarry. Four years of forecasting was compared with actual
data and mean squared error (MSE), mean absolute percentage error (MAPE)
and bias were calculated. Parameters of the simulation and parameters of
the methods were varied.

Results: Capturing only these three characteristics of the mortality
time series allowed the robust reproduction of the phenomenon underlying
WHO’s results. Using 2015 as the starting year – as in the WHO’s study –
results in very poor performance for the WHO’s method, clearly revealing
the problem as even simple linear extrapolation was better. However, the
Acosta-Irizarry method substantially outperformed WHO’s method despite
being also based on splines. In certain – but not all – scenarios,
errors were substantially affected by the parameters of the mortality
curve, but the ordering of the methods was very stable. Results were
highly dependent of the parameters of the estimation procedure: even
WHO’s method produced much better results if the starting year was
earlier, or if the basis dimension was lower. Conversely,
Acosta-Irizarry method can generate poor forecasts if the number of
knots is increased. Linear extrapolation could produce very good
results, but is highly dependent on the choice of the starting year,
while average was the worst in almost all cases.

Discussion and conclusion: WHO’s method is highly dependent on the
choice of parameters and is almost always dominated by the
Acosta-Irizarry method. Linear extrapolation could be better, but it is
highly dependent on the choice of the starting year; in contrast,
Acosta-Irizarry method exhibits a relatively stable performance
irrespectively of the starting year. Using the average method is almost
always the worst except for very special circumstances. This proves that
splines are not inherently unsuitable for predicting baseline mortality,
but care should be taken, in particular, these results suggest that the
key issue is that the structure of the splines should be rigid. No
matter what approach or parametrization is used, model diagnostics must
be performed before accepting the results, and used methods should be
preferably validated with extensive simulations on synthetic datasets.
Further research is warranted to understand how these results can be
generalized to other scenarios.

## Introduction

Excess mortality is the difference between the actual mortality (number
of deaths) of a time period in a given country or (sub- or
supranational) region and its “expected” mortality, defined as the
mortality statistically forecasted from the region’s historical data.
Calculation of excess mortality can be used to characterize the impact
of an event on mortality if the historical data is prior to the onset of
the event, therefore the prediction pertains to a counterfactual:
mortality that would have been observed without the event \[1\]. Thus,
the difference, i.e., excess mortality indeed measures the impact of the
event, assuming the prediction was correct.

Calculating excess mortality is particularly useful if the event’s
impact on mortality is hard to measure directly, for instance, one of
the typical applications is characterizing the mortality associated with
natural disasters \[2–4\], but is is also used for epidemics where
direct mortality registration is missing, incomplete or unreliable, a
prime example being the seasonal flu \[5,6\].

COVID-19 is no exception. While mortality is reported in developed
countries, typically daily or weekly, this suffers from two drawbacks,
first, the number of reported deaths is – while to much less extent than
the number of reported cases – but still contingent on testing activity,
which might be vastly different between countries or time periods, and
second, despite efforts at standardization the criteria for the
certification of deaths might be different between countries \[7\].
Excess mortality resolves both of these problems as it is completely
immune to differences in testing intensity and cause of death
certification procedure. This makes it especially suitable for
between-country comparisons, which is a critical issue to better
understand the pandemic, in particular, to evaluate different control
measures and responses \[8\].

This, however, comes at a price. First, and perhaps most importantly,
excess mortality is inherently a gross metric, measuring both the direct
and the indirect effects, the latter of which can be both positive
(e.g., COVID control measures also provided protection against flu) or
negative (e.g., the treatment of other diseases became less efficient)
\[9\]. Second, the excess mortality is the slowest indicator: the
necessary data, that is, the number of deaths usually becomes available
with a 4 week lag (and even that is typically revised to some extent
later) even in the developed countries. Finally, the whole calculation
depends on how accurate the forecast was.

The last of these issues will be the focus now: given the importance of
cross-country comparisons, it is crucial the results indeed reflect
differences and are not too sensitive to the used prediction method.

Only those methods will be considered now that use traditional
regression approach, i.e., methods using ARIMA models \[10–13\],
Holt-Winters method \[14\] or based on Gaussian process \[15\] won’t be
considered, just as ensemble methods \[16,17\]. Question concerning age-
or sex stratification or standardization \[18\], small area estimation
\[19,20\] and inclusion of covariates, such as temperature, to improve
modelling \[16,17,19\] will also not be considered here.

This leaves us with two questions, the handling of seasonality and the
handling of long-term trend. For the latter, these are the typical
solutions in the literature concerning COVID-19:

-   Using the last pre-pandemic year \[16,21\]. This is good – even if
    not perfect – considering the long-term trends, as it uses data
    closest to the investigated period, but it has a huge variance, due
    to the natural year-to-year variability of mortality.
-   Using the average of a few pre-pandemic years (typically 5)
    \[22–28\]. This is more reliable as averaging reduces variability,
    however, it is even more biased in case the mortality has a
    long-term trend (which it almost always has), for instance, if
    mortality is falling, this provides an overestimation, thus, excess
    mortality is underestimated.
-   Using a linear trend extrapolation \[29–31\]. This accounts for the
    potential trends in mortality, removing the bias of the above
    methods, at least as far as linearity is acceptable, but it depends
    on the selection of the starting year from which the linear curve is
    fitted to the data. It also has higher variance than the averaging
    approach, but is is usually less of concern, given the huge amount
    of data typically used (unless a small country, and/or age or sex
    strata is investigated).
-   Using splines \[32,33\]. The method of Acosta and Irizarry \[34,35\]
    is based on splines, just as many other custom implementation
    \[16,36\], which, crucially, includes the model used by the World
    Health Organization (WHO) \[37\].

The choice of this might have a highly relevant impact on the results of
the calculation, as evidenced by the case of the excess mortality
estimation of the WHO. On May 5, 2022, the WHO published its excess
mortality estimates \[38\], which immediately raised questions: the
estimates for Germany were surprisingly high, while that of Sweden were
low \[39\]. The case of Germany was especially intriguing, so much that
one paper termed it the “German puzzle” \[39\]. Figure
<a href="#fig:germanpuzzle">1</a> illustrates the “German puzzle” using
actual German data, with the curves fitted on 2015-2019 data and
extrapolated to 2020 and 2021 (as done by the WHO): while the dots
visually indicate rather clearly a simple upward trend (as shown by the
linear extrapolation), the spline prediction turns back.

``` r
fitSpline <- mgcv::gam(outcome ~ s(Year, k = 5) + s(WeekScaled, bs = "cc"),
                       data = RawData[Year>=2015&Year<=2019,],
                       family = mgcv::nb(), method = "REML")
fitLin <- mgcv::gam(outcome ~ Year + s(WeekScaled, bs = "cc"),
                    data = RawData[Year>=2015&Year<=2019,],
                    family = mgcv::nb(), method = "REML")

predgrid <- CJ(Year = seq(2015, 2022, length.out = 100), WeekScaled = 0.5)

predSpline <- predict(fitSpline, newdata = predgrid, newdata.guaranteed = TRUE, se.fit = TRUE,
                      type = "terms")
predSpline <- data.frame(fit = predSpline$fit[, 1] + attr(predSpline, "constant"),
                         se.fit = predSpline$se.fit[,1])
predLin <- predict(fitLin, newdata = predgrid, se.fit = TRUE, type = "terms")
predLin <- data.frame(fit = predLin$fit[, 1] + attr(predLin, "constant"), se.fit = predLin$se.fit[,1])

predgrid <- rbind(cbind(predgrid, Type = "Spline", predSpline),
                  cbind(predgrid, Type = "Linear", predLin))

ggplot(predgrid, aes(x = Year, y = exp(fit),  color = Type, fill = Type)) +
  geom_line() + geom_point(data = RawDataYearly[Year>=2015&Year<=2019,],
                           aes(x = Year, y = (outcome/52.25)), inherit.aes = FALSE) +
  labs(y = "Predicted mortality [/week]")
```

![Figure 1: A linear trend and a spline fitted on German mortality data
2015-2019 and extrapolated to 2020 and
2021.](README_files/figure-gfm/germanpuzzle-1.png)

Figure: A linear trend and a spline fitted on German mortality data
2015-2019 and extrapolated to 2020 and 2021.

The explanation later provided by the WHO \[37\] stated that the problem
was due to two issues, first, the WHO applied a rescaling method to the
raw data to compensate for underreporting (due to late registration, for
instance), but this was unnecessary in case of Germany, with excellent
death registration. Figure <a href="#fig:germanpuzzle">1</a> shows the
unadjusted German data avoiding this problem, so that focus can be
placed on the second issue that will be the subject of investigation
now: the usage of splines.

As described above, WHO’s method uses a spline to capture the long-term
trend, and it seems that the problem is that the lower data of 2015 had
a very high impact on the spline, with that single observation turning
the entire spline, despite earlier points showing an upward trend. It
seems to much weight is put on this – likely short-term, random,
noise-like – fluctuation, i.e., the extrapolation was too sensitive to
this. The culprit is quickly identified as spline-regression itself,
with one commentator saying “\[e\]xtrapolating a spline is a known bad
practice” \[39\].

(Interestingly enough, the problem only exists in this particular case
if year is used as a predictor. WHO’s paper suggests just this \[37\],
but this might be only a typo, as the long-term trend should be
represented not by the – abruptly changing – year indicator, but rather
a continuously changing indicator of time, such as days since a given
date.)

But really splines are to be blamed? Motivated by obtaining a better
understanding of the “German puzzle”, this paper aims to investigate the
following questions: 1) Really splines *per se* were the culprit? 2)
What were the particlar characteristics, both of the scenario and of the
used spline-regression, that gave rise to the problem? 3) Is there a
better way to predict baseline for excess mortality calculation avoiding
this problem?

To answer these questions, first a model will be devised that is able to
generate mortality curves that capture the relevant features exhibited
by the real-life German example. Thus, it’ll be possible to calculate
the accuracy of a forecast (as the ground truth is now known), and also
to investigate how parameters of the simulation influence it. With
averaging several simulations, the mean accuracy can be approximated,
allowing the comparison of the methods, and investigating its dependence
on the parameters – both of the mortality curve and of the parameters of
the method – thereby hopefully resolving the “German puzzle”.

## Material and Methods

### Data source

Weekly all-cause mortalities for Germany were obtained from the European
Statistical Office (Eurostat), database `demo_r_mwk_ts` \[40\]. No
additional preprocessing or correction was applied such as that for late
registration, i.e., the part of the problem with the WHO’s approach due
to upscaling was avoided, so that the focus is now solely on the
modeling aspect. A detailed comparison of the possible data sources can
be found in Supplementary Material 1.

Basic properties (raw weekly values, yearly trend, seasonal pattern) are
shown on Figure <a href="#fig:german-raw-plots">2</a>.

``` r
p1 <- ggplot(RawData[Year<=2019], aes(x = date, y = outcome)) + geom_line() +
  labs(x = "Date", y = "Mortality [/week]")
p2 <- ggplot(RawDataYearly[Year<=2019], aes(x = Year, y = outcome)) + geom_point() +
  geom_line() + geom_smooth(formula = y ~ x, method = "loess", se = FALSE) +
  labs(x = "Year", y = "Mortality [/year]")
p3 <- ggplot(RawData[Year<=2019], aes(x = Week, y = outcome, group = Year)) +
  geom_line(alpha = 0.2) + labs(x = "Week of year", y = "Mortality [/week]")

egg::ggarrange(p1, p2, p3, ncol = 1)
```

![Figure 2: Weekly mortalities (upper), yearly mortalities with
LOESS-smoother (middle) and seasonal pattern (bottom) of the German
mortality data,
2000-2019.](README_files/figure-gfm/german-raw-plots-1.png)

Figure: Weekly mortalities (upper), yearly mortalities with
LOESS-smoother (middle) and seasonal pattern (bottom) of the German
mortality data, 2000-2019.

### Simulation model

Based on the patterns that can be observed on Figure
<a href="#fig:german-raw-plots">2</a>, the following three components
will be used to create synthetic datasets:

-   Long-term change, modelled with quadratic trend; described by three
    parameters (constant, linear and quadratic term)
-   Deterministic seasonality, modelled with a single harmonic
    (sinusoidal) term; described by two parameters (amplitude and phase)
-   Random additional peaks during the winter (i.e., flu season) and
    during the summer (i.e., heat waves); described in each season by
    five parameters (probability of the peak, minimum and maximum value
    of the peak height, minimum and maximum value of the peak width)

These govern the expected value; the actual counts are obtained from a
negative binomial distribution with constant size. Detailed description
of how the above model is built, and what parameters are used can be
found in Supplementary Material 2.

### Baseline mortality prediction

Four methods will be used for predicting mortality, including the WHO’s
method, an advanced alternative method that also uses splines, developed
by Acosta and Irizarry in 2020 \[34\], and two naive methods as a
comparison. These cover the widely used, classical statistical methods
used for predicting baseline mortality in excess mortality studies.

-   Average: after accounting for seasonality with a single cyclic
    spline, the average of the preciding years will be used as the –
    constant – predicted value. Parameter: starting year (i.e., how many
    previous year is used for averaging). Some studies used the last
    pre-pandemic year (2019) as the predicted baseline mortality, this
    is just the special case of this method, with the starting year set
    to 2019.
-   Linear: after accounting for seasonality with a single cyclic
    spline, the long-term trend is modelled with a linear trend, that is
    extrapolated. Parameter: starting year (from which the model is
    fitted).
-   WHO’s method: the method is reconstructed from the description
    provided in \[37\]. In brief, seasonality is accounted with single a
    cyclic spline (as done in the previous cases), and the long-term
    trend is accounted with a thin plate regression spline. The only
    deviation compared to WHO’s paper is that – as noted above – the
    actual time (number of days since 1970-01-01) is used as the
    predictor for long-term trend, not the abruptly changing year. The
    model is estimated with restricted maximum likelihood (REML).
    Parameters: starting year (from which the model is fitted) and
    ![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k"),
    the dimension of the basis of the spline used for capturing the
    long-term trend.
-   Acosta-Irizarry (AI) method: the method described in \[34\] using
    their reference implementation. Parameters: starting year (from
    which the model is fitted) and
    ![tkpy](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;tkpy "tkpy"),
    the number of trend knots per year; other parameters are left on
    their default values (i.e., two harmonic term is used).

### Validation through simulation

First, a synthetic dataset is randomly generated from the model
described above using the investigated parameters (parameters of the
scenario). This dataset simulates mortalities from the beginning of 2020
to the end of 2023. Then, the investigated prediction method with the
investigated parametrization (parameters of the method) is applied, and
after fitting, a prediction is obtained for 2020 to 2023. The goodness
of prediction is quantified as mean squared error (MSE), mean absolute
percentage error (MAPE) and bias in those 4 years based on 1000
replications of this simulational procedure. This is repeated for all
parameters of the scenario, all prediction methods and all parameters of
the method.

Investigated parameters of the prediction methods were the following:

-   Average: starting year 2000, 2005, 2010, 2015, 2019
-   Linear: starting year 2000, 2005, 2010, 2015
-   WHO’s method: all possible combination of starting year 2000, 2005,
    2010, 2015 and
    ![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
    (basis dimension) 5, 10, 15, 20
-   Acosta-Irizarry method: all possible combination of starting year
    2000, 2005, 2010, 2015 and
    ![tkpy](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;tkpy "tkpy")
    (trend knots per year) 1/4, 1/5, 1/7, 1/9, 1/12

For the scenario, simulations were run with the optimal parameters
discussed in Supplementary Material 2 (base case scenario) and then all
parameters were varied in 5 steps from half of the base case value to
twice of that, with two exception: the constant term of the trend is
varied only from 90% to 110% (to avoid irrealistic values) and
probabilities are prevented from going above 100%.

Details are provided in Supplementary Material 3.

### Programs used

All calculations are carried out under the R statistical program package
version 42.0 \[41\] using packages `data.table` \[42\] (version 1.14.2)
and `ggplot2` \[43\] (version 3.3.6), as well as `excessmort` (version
0.6.1), `mgcv` (version `rpackageVersion("mgcv")`), `scorepeak` (version
0.1.2), `parallel` (version 4.2.0) `lubridate` (version 1.8.0),
`ISOweek` (version 0.6.2) and `eurostat` (version 3.7.10).

Full source code allowing complete reproducibility is openly available
at <https://github.com/tamas-ferenci/MortalityPrediction>.

## Results

Figure <a href="#fig:trajs">3</a> illustrates the predictions (for
2020-2023) for the base case scenario by showing the estimated yearly
deaths for 200 randomly selected simulation together with the ground
truth for all 4 method with all possible parameters.

``` r
p1 <- ggplot() + xlim(c(2018, 2023)) + ylim(c(0.5, 1.5)) +
  geom_line(data = predLongs$WHO[rep<=200&parsim==1],
            aes(x = Year, y = value/1e6, group = rep), alpha = 0.1) + 
  geom_line(data = predLongs$WHO[parsim==1, .(outcome = mean(outcome)), .(Year)],
            aes(x = Year, y = outcome/1e6), color = "red") +
  facet_grid(rows = vars(k), cols = vars(startyear)) +
  labs(y = "Outcome [million death / year]", title = "A) WHO's method")

p2 <- ggplot() + xlim(c(2018, 2023)) + ylim(c(0.5, 1.5)) +
  geom_line(data = predLongs$AI[rep<=200&parsim==1],
            aes(x = Year, y = value/1e6, group = rep), alpha = 0.1) + 
  geom_line(data = predLongs$AI[parsim==1, .(outcome = mean(outcome)), .(Year)],
            aes(x = Year, y = outcome/1e6), color = "red") +
  facet_grid(rows = vars(tkpy), cols = vars(startyear)) +
  labs(y = "Outcome [million death / year]", title = "B) Acosta-Irizarry method")

p3 <- ggplot() + xlim(c(2018, 2023)) + ylim(c(0.5, 1.5)) +
  geom_line(data = predLongs$Lin[rep<=200&parsim==1],
            aes(x = Year, y = value/1e6, group = rep), alpha = 0.1) + 
  geom_line(data = predLongs$Lin[parsim==1, .(outcome = mean(outcome)), .(Year)],
            aes(x = Year, y = outcome/1e6), color = "red") +
  facet_grid(cols = vars(startyear)) +
  labs(y = "Outcome [million death / year]", title = "C) Linear trend")

p4 <- ggplot() + xlim(c(2018, 2023)) + ylim(c(0.5, 1.5)) +
  geom_line(data = predLongs$Average[rep<=200&parsim==1],
            aes(x = Year, y = value/1e6, group = rep), alpha = 0.1) + 
  geom_line(data = predLongs$Average[parsim==1, .(outcome = mean(outcome)), .(Year)],
            aes(x = Year, y = outcome/1e6), color = "red") +
  facet_grid(cols = vars(startyear)) +
  labs(y = "Outcome [million death / year]", title = "D) Average")

egg::ggarrange(p1, p2, p3, p4, ncol = 1, heights = c(1, 1, 0.4, 0.4))
```

![Figure 3: Estimated yearly deaths (for 2020-2023) for 200 randomly
selected simulation together with the ground truth. A) WHO’s method, B)
Acosta-Irizarry method, C) Linear trend, D) Average. Parameters of the
methods are shown in column and row headers. Parameters of the scenario
are set to the base case values.](README_files/figure-gfm/trajs-1.png)

Figure: Estimated yearly deaths (for 2020-2023) for 200 randomly
selected simulation together with the ground truth. A) WHO’s method, B)
Acosta-Irizarry method, C) Linear trend, D) Average. Parameters of the
methods are shown in column and row headers.

Figure <a href="#fig:trajs">3</a> already strongly suggests some
tendencies, but to precisely evaluate it, error metrics have to be
calculated. Figure <a href="#fig:errorWHOAI">4</a>. shows all three
error metrics for WHO and Acosta-Irizarry methods, for all possible
parametrizations.

``` r
p1 <- ggplot(melt(predLongs$WHO[parsim==1,.(MSE = mean((outcome-value)^2)/1e6,
                                            MAPE = mean(abs(outcome-value)/value)*100,
                                      Bias = mean((outcome-value)/value)*100), .(startyear, k)],
                  id.vars = c("startyear", "k")),
       aes(x = startyear, y = value, color = factor(k))) + facet_wrap(~variable, scales = "free") +
  geom_line() + geom_point() + labs(color = "k", title = "A) WHO's method")
p2 <- ggplot(melt(predLongs$AI[parsim==1,.(MSE = mean((outcome-value)^2)/1e6,
                                           MAPE = mean(abs(outcome-value)/value)*100,
                                     Bias = mean((outcome-value)/value)*100) , .(startyear, tkpy)],
                  id.vars = c("startyear", "tkpy")),
       aes(x = startyear, y = value, color = factor(tkpy))) + facet_wrap(~variable, scales = "free") +
  geom_line() + geom_point() + labs(color = "tkpy", title = "B) Acosta-Irizarry method")
egg::ggarrange(p1, p2, ncol = 1)
```

![Figure 4: Different error metrics (MSE, MAPE, Bias) for the WHO’s
method (above) and the Acosta-Irizarry method (below) for all possible
parameter combinations of these two methods. Parameters of the scenario
are set to the base case
values.](README_files/figure-gfm/errorWHOAI-1.png)

Figure: Different error metrics (MSE, MAPE, Bias) for the WHO’s method
(above) and the Acosta-Irizarry method (below) for all possible
parameter combinations of these two methods. Parameters of the scenario
are set to the base case values.

As already suggested by Figure <a href="#fig:trajs">3</a>, it confirms
that
![k=5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%3D5 "k=5")
(WHO) and
![tkpy = 1/7](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;tkpy%20%3D%201%2F7 "tkpy = 1/7")
(Acosta-Irizarry) are the best parameters in this particular scenario.
It therefore worth comparing all four methods with different starting
years, but with the remaining parameters
(![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
for WHO’s method,
![tkpy](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;tkpy "tkpy")
for the Acosta-Irizarry method) set to the values that are optimal in
this particular scenario (Figure <a href="#fig:errorall">5</a>).

``` r
ggplot(melt(rbind(
  predLongs$WHO[parsim==1&k==5,.(Method = "WHO", MSE = mean((outcome-value)^2)/1e6,
                                 MAPE = mean(abs(outcome-value)/value)*100,
                                 Bias = mean((outcome-value)/value)*100), .(startyear)],
  predLongs$AI[parsim==1&invtkpy==7,.(Method = "AI", MSE = mean((outcome-value)^2)/1e6,
                                   MAPE = mean(abs(outcome-value)/value)*100,
                                   Bias = mean((outcome-value)/value)*100), .(startyear)],
  predLongs$Lin[parsim==1,.(Method = "Linear", MSE = mean((outcome-value)^2)/1e6,
                            MAPE = mean(abs(outcome-value)/value)*100,
                            Bias = mean((outcome-value)/value)*100), .(startyear)],
  predLongs$Average[parsim==1,.(Method = "Average", MSE = mean((outcome-value)^2)/1e6,
                                MAPE = mean(abs(outcome-value)/value)*100,
                                Bias = mean((outcome-value)/value)*100), .(startyear)]),
  id.vars = c("startyear", "Method")),
  aes(x = startyear, y = value, color = Method)) + facet_wrap(~variable, scales = "free") +
  geom_point() + geom_line() + labs(x = "Starting year")
```

![Figure 5: Different error metrics (MSE, MAPE, Bias) for all method,
with different starting years, but remaining parameters set to optimal
values for this particular scenario. Parameters of the scenario are set
to the base case values.](README_files/figure-gfm/errorall-1.png)

Figure: Different error metrics (MSE, MAPE, Bias) for all method, with
different starting years, but remaining parameters set to optimal values
for this particular scenario. Parameters of the scenario are set to the
base case values.

All the above investigations used the base case scenario for the
simulated mortality curve. Figure <a href="#fig:errorscenarios">6</a>.
shows the best error metrics achievable with each method in a given
scenario.

``` r
ggplot(melt(rbindlist(lapply(predLongs, function(pl)
  pl[,.(logMSE = log(mean((outcome-value)^2)/1e6), MAPE = mean(abs(outcome-value)/value)*100,
        Bias = mean((outcome-value)/value)*100),
     by = setdiff(names(pl), c("parmethod", "rep", "outcome", "Year", "value", "parsimName",
                               "parsimValue"))][, .(logMSE = min(logMSE), MAPE = min(MAPE),
                                                    Bias = min(Bias)), .(parsim, startyear)]),
  idcol = "Method"), id.vars = c("parsim", "Method", "startyear")),
  aes(x = parsim, y = value, color = Method)) +
  facet_grid(rows = vars(variable), cols = vars(startyear), scales = "free") + geom_point() +
  labs(x = "Scenario #")
```

![Figure 6: Best achievable error metrics with each method in each
simulational scenario of the mortality curve (#1 is the base case
scenario).](README_files/figure-gfm/errorscenarios-1.png)

Figure: Best achievable error metrics with each method in each
simulational scenario of the mortality curve (#1 is the base case
scenario).

Finally, note that as different methods were evaluated on the same
simulated dataset for each simulation, it is possible to compare not
only the averages, but directly compare the errors themselves. Figure
<a href="#fig:errordirect">7</a> shows direct comparison between the
best parametrization of the WHO’s method and the Acosta-Irizarry method
for 200 randomly selected simulations in each scenario. In the base case
scenario, Acosta-Irizarry performed better in 71.1% of the cases.

``` r
ggplot(merge(predLongs$WHO[startyear==2015&k==5, .(errorWHO = mean((value-outcome)^2)),
                           .(rep, parsimName, parsimValue)],
             predLongs$AI[startyear==2015&invtkpy==7, .(errorAI = mean((value-outcome)^2)),
                          .(rep, parsimName, parsimValue)])[parsimName!="Base"&rep<=200],
       aes(x = errorWHO, y = errorAI, color = factor(parsimValue))) + facet_wrap(~parsimName) +
  geom_point() + geom_abline(color = "red") +
  geom_point(data = merge(predLongs$WHO[startyear==2015&k==5&parsimName=="Base"&rep<=200,
                                        .(errorWHO = mean((value-outcome)^2)), .(rep)],
                          predLongs$AI[startyear==2015&invtkpy==7&parsimName=="Base"&rep<=200,
                                       .(errorAI = mean((value-outcome)^2)), .(rep)]),
             aes(x = errorWHO, y = errorAI), inherit.aes = FALSE) +
  scale_x_log10(labels = scales::label_log()) + scale_y_log10(labels = scales::label_log()) +
  annotation_logticks() +
  labs(x = "Mean squared error, WHO method (starting year: 2015, k = 5)",
       y = "Mean squared error, Acosta-Irizarry method (starting year: 2015, trend knots per year = 7)",
       color = "Scenario") + scale_color_discrete()
```

![Figure 7: Errors – squared distance of the predicted outcome from its
true value – of the WHO’s method and the Acosta-Irizarry method (under
best parametrization) on the same simulated datasets for 200 randomly
selected simulations with different scenarios. Black dots indicate the
base case scenario, scenario \#1 to \#5 represent varying the parameter
shown on the panel from half of its base case value to twice (with the
exception of the constant term where it is varied from 90% to 110%).
Probabilities are limited to be below
100%.](README_files/figure-gfm/errordirect-1.png)

Figure: Errors – squared distance of the predicted outcome from its true
value – of the WHO’s method and the Acosta-Irizarry method (under best
parametrization) on the same simulated datasets for 200 randomly
selected simulations with different scenarios. Black dots indicate the
base case scenario, scenario \#1 to \#5 represent varying the parameter
shown on the panel from half of its base case value to twice (with the
exception of the constant term where it is varied from 90% to 110%).
Probabilities are limited to be below 100%.

## Discussion

These results demonstrate that we were able to reliably reproduce the
“German puzzle” using synthetic datasets. This approach allowed a deep
investigation of how the results depend on the used method, its
parameters and on the parameters of the scenario.

As expected, prediction with average has the highest error, and is
highly biased (but is improved by shortening the fitting dataset). This
of course depends on the form of the historical mortality curve,
theoretically, for a more or less constant curve even this prediction
can work better.

Linear extrapolation seems to be a very promising alternative, the only
problem being that it is very sensitive to the appropriate choice of the
starting year: that phase should be covered where the change in
historical mortality is linear. This is prone to subjectivity and might
not work at all if the linear phase is too short (limiting the available
information), i.e., it depends on how wiggly is the historical curve.

Splines in contrast can work theoretically well even in those cases: it
can use all historical data, i.e., it is not abruptly cut off as with
linear extrapolation, but more weight is placed on the trends suggested
by the recent observations. At first glance, this seems to be the ideal
solution, but as this investigation reveals, what is meant by “more
weight” and “recent” is crucial, and certain choices can results in very
poor extrapolations, despite the tempting theoretical properties.

The overall picture in selecting the optimal parameters, confirmed by
the results of both spline-based methods, is that splines should be
quite rigid in baseline mortality prediction for excess mortality
calculation. This is the concordant conclusion from the experiences both
with the WHO method (increasing basis dimension decreased performance)
and the Acosta-Irizarry method (increasing trend knots per year
increased performance). The WHO method is only acceptable with
![k=5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%3D5 "k=5")
(but even that requires longer observation than starting from 2015, as
was done by the WHO), not higher. For the Acosta-Irizarry method, 1/4
trend knots per year was definitely too flexible, perhaps even 1/5 is
too high. Note that the default value in the reference implementation of
the Acosta-Irizarry method is 1/7, and authors in fact do not recommend
using a value much higher.

The likely explanation is that mortality curves exhibit only slow
changes, so high flexibility is not required, and – as with any
regression model with too high model capacity – can be downright
detrimental, as it allows the model to pick up noise, i.e., can result
in overfitting.

Note that data are presented using the ISO 8601 year defintion meaning
that “year” can be either 52 or 53 weeks long \[44\]; from 2015 to 2019
every year is 52 weeks long except 2015 which is one week longer. This
adds to the reasons why the value of 2015 is higher, increasing the
wiggliness of the German data and potentially contributing to the
problem.

Among the two spline-based methods, Acosta-Irizarry almost always
performed better than WHO’s method, this is especially true if the
fitting dataset was shorter (i.e., the starting year was later). This is
a crucial component in WHO’s experience.

We are aware of two previous works from the literature that are
comparable to the present investigation. Both Nepomuceno et al \[45\]
and Shöley \[46\] is similar to ours in a sense that they used – among
others – several models, partly overlapping those presented here, but,
importantly, neither of them considered splines for long-term trend at
all. Nepomuceno et al did not try to evaluate the methods, it compared
them to each other without having ground truth, i.e., the aim was to
investigate concordance. In contrast, Shöley did try to give an
objective evaluation, but in contrast to the synthetic dataset
simulation approach used here, it applied time series cross-validation
with historical data to measure accuracy. Cross-validation has the
advantage that is guaranteed to be realistic – as opposed to a
simulation – but there is less freedom, as investigators are bound to
empirical data, with limited possibility in varying the parameters.

Thus, we believe that this is the first systematical study to
investigate the application of splines for the long-term trend component
in mortality prediction for excess mortality calculation, the first to
investigate the impact of their parametrization, and the first to use
synthetic dataset validation in general for that end.

The most important limitation of the present study is that every
simulation is committed to the same model structure, for instance,
long-term trend is limited to be quadratic. It would be important to
extend the present investigation to other long-term trends, and to other
models in general (such as those with different seasonality, or
interaction between seasonality and trend etc.).

We did not investigate the impact of using the population and modelling
death rates versus modelling death counts directly, nor did we examine
the potential impact of the frequency of the data. Weekly data was used
throughout this study: this might not be available for developing
countries, but it is almost universally available for developed
countries, in which case, the use of the most frequent data seems to be
logical (with the appropriate handling of seasonality). In the same
vein, adjustment for late registration and imputation of missing data,
which might be needed where full data is not available, is not
considered here, as the focus is on developed countries.

## Conclusion

WHO’s method is highly dependent on the choice of parameters and is
almost always dominated by the Acosta-Irizarry method. Linear
extrapolation could be better, but it is highly dependent on the choice
of the starting year; in contrast, Acosta-Irizarry method exhibits a
relatively stable performance irrespectively of the starting year. Using
the average method is almost always the worst except for very special
circumstances.

This proves that splines are not inherently unsuitable for predicting
baseline mortality, but care should be taken, in particular, these
results suggest that the key issue is that the structure of the splines
should be rigid. No matter what approach or parametrization is used,
model diagnostics must be performed before accepting the results, and
used methods should be preferably validated with extensive simulations
on synthetic datasets. Further research is warranted to understand how
these results can be generalized to other scenarios.

## Supplementary Material 1: Comparison of data sources

For a country like Germany, four data sources come into consideration
for weekly mortality data: Eurostat \[40\], the Short-Term Mortality
Fluctuations (STMF) dataset of the Human Mortality Database (WMD)
\[47\], the World Mortality Database \[1\] and the national data
provider (in this case, the Federal Statistical Office of Germany). The
last is usually more complicated, limits extension to other countries
and is unnecessary for developed countries, so it’ll be avoided in this
case. Also, for Germany, WMD simply copies the data of the STMF (“We
collect the weekly STMF data for the following countries: \[…\] Germany,
\[…\].”) leaving us with two options.

We shall compare whether these two report identical data (Figure
<a href="#fig:datasources">8</a>).

``` r
RawDataEurostat <- as.data.table(eurostat::get_eurostat("demo_r_mwk_ts", time_format = "raw"))
RawDataEurostat <- RawDataEurostat[geo=="DE"&sex=="T"]
RawDataEurostat$Year <- as.numeric(substring(RawDataEurostat$time, 1, 4))
RawDataEurostat$Week <- as.numeric(substring(RawDataEurostat$time, 6, 7))

RawDataSTMF <- fread("https://www.mortality.org/File/GetDocument/Public/STMF/Outputs/stmf.csv")
RawDataSTMF <- RawDataSTMF[CountryCode=="DEUTNP"&Sex=="b"]

RawDataEurostatSTMF <- merge(RawDataEurostat[, .(Year, Week, ES = values)],
                             RawDataSTMF[, .(Year, Week, STMF = DTotal)])

ggplot(RawDataEurostatSTMF, aes(x = ES, y = STMF)) + geom_point() + geom_abline(color = "red") +
  labs(x = "Eurostat value [/week]", y = "STMF value [/week]")
```

![Figure 8: Supplementary Figure 1: Weekly number of deaths according to
the Eurostat (horizontal axis) and the STMF database (vertical axis) in
Germany. Red line indicates the line of
equality.](README_files/figure-gfm/datasources-1.png)

Supplementary Figure 1: Weekly number of deaths according to the
Eurostat (horizontal axis) and the STMF database (vertical axis) in
Germany. Red line indicates the line of equality.

The two are almost identical (with a correlation of 0.9999996), with
differences only occuring for the latest data and of minimal magnitude,
so we can safely use the Eurostat database.

## Supplementary Material 2: Generating realistic synthetic datasets

First, Figure <a href="#fig:german-raw-plots">2</a> should be inspected,
as it already gives important clues on the setup of a realistic model
from which synthetic datasets could be generated. Figure
<a href="#fig:yearsonpanels">9</a> gives further insight by plotting
each year separately.

``` r
RawData2019 <- RawData[Year<=2019]

ggplot(RawData2019, aes(x = Week, y = outcome)) + geom_line() + facet_wrap(~Year)
```

![Figure 9: Weekly number of deaths in Germany, separated according to
year.](README_files/figure-gfm/yearsonpanels-1.png)

Figure: Weekly number of deaths in Germany, separated according to year.

The following observations can be made:

-   There is a long-term trend, seemingly quadratic.
-   There is a strong seasonality with winter peak and summer trough.
-   There are peaks – in addition to the seasonality – in the winter and
    also during the summer (although the shape seems to be different,
    with winter peaks seeming to be broader and higher).

To investigate these, first a spline-smoothing – with thin plate
regression spline \[32\] – is applied to obtain the long-term trend, and
a single harmonic term is included as a covariate to remove seasonality.
No interaction is assumed between the two, i.e., it is assumed that the
seasonal pattern is the same every year. (The peaks are not accounted
for at this stage which means that the curve is above the true one, but
the difference is likely has minimal due to the rarity and short
duration of the peaks. This will be later corrected, after the peaks
were identified.) All analysis will be carried out on the log scale
(meaning the effect of covariates is multiplicative) using negative
binomial response distribution to allow for potential overdispersion
\[48\].

Figure <a href="#fig:longterm">10</a> shows the results, overplotted
with the model where the long-term trend is a completely parametric
quadratic trend. One can observe very good fit between the two, so all
subsequent investigation will use the quadratic trend which is much
easier to handle. This is only meaningful for short-term extrapolation,
but this is what will be needed now (two years of extrapolation will be
used in the present study); also it is not possible to better
differentiate between functional forms at this sample size.

``` r
fitSpline <- mgcv::gam(outcome ~ s(NumTrend) + cos(2*pi*WeekScaled) + sin(2*pi*WeekScaled),
                       data = RawData2019, family = mgcv::nb(), method = "REML")
fit <-  mgcv::gam(outcome ~ poly(NumTrend, 2, raw = TRUE) + cos(2*pi*WeekScaled) + sin(2*pi*WeekScaled),
                  data = RawData2019, family = mgcv::nb(), method = "REML")

fitTrendMinPoint <- -coef(fit)["poly(NumTrend, 2, raw = TRUE)1"]/
  (2*coef(fit)["poly(NumTrend, 2, raw = TRUE)2"])
fitTrendMinValue <- coef(fit)["(Intercept)"]-
  coef(fit)["poly(NumTrend, 2, raw = TRUE)1"]^2/(4*coef(fit)["poly(NumTrend, 2, raw = TRUE)2"])
fitTrend2020End <- coef(fit)["poly(NumTrend, 2, raw = TRUE)2"]*18624^2 +
  coef(fit)["poly(NumTrend, 2, raw = TRUE)1"]*18624 + coef(fit)["(Intercept)"]
fitSeasonAmplitude <- sqrt(coef(fit)["sin(2 * pi * WeekScaled)"]^2 +
                             coef(fit)["cos(2 * pi * WeekScaled)"]^2)
fitSeasonPhase <- atan(-coef(fit)["sin(2 * pi * WeekScaled)"]/coef(fit)["cos(2 * pi * WeekScaled)"])

predgrid <- data.frame(NumTrend = seq(min(RawData2019$NumTrend),
                                      max(RawData2019$NumTrend), length.out = 100),
                       WeekScaled = rep(0.5, 100))
predgrid <- rbind(cbind(predgrid, Type = "Spline",
                        with(predict(fitSpline, newdata = predgrid, newdata.guaranteed = TRUE,
                                     se.fit = TRUE), data.frame(fit, se.fit))),
                  cbind(predgrid, Type = "Quadratic",
                        with(predict(fit, newdata = predgrid, newdata.guaranteed = TRUE, se.fit = TRUE),
                                                           data.frame(fit, se.fit))))

ggplot(predgrid, aes(x = lubridate::as_date(NumTrend), y = exp(fit), ymin = exp(fit - 1.96*se.fit),
                     ymax = exp(fit + 1.96*se.fit), color = Type, fill = Type)) +
  geom_line() + geom_ribbon(alpha = 0.2, linetype = 0) +
  labs(x = "Date", y = "Predicted number of weekly deaths (adjusted to June-30)")
```

![Figure 10: Fitting long-term trend as spline (black) and as quadratic
trend (red); shaded area indicated 95% confidence interval. Seasonality
is removed by including a single harmonic term in the regression in both
cases.](README_files/figure-gfm/longterm-1.png)

Figure: Fitting long-term trend as spline (black) and as quadratic trend
(red). Seasonality is removed by including a single harmonic term in the
regression in both cases.

The coefficients can be transformed to equivalent forms that are more
meaningful. Thus, the three parameters of the quadratic trend can be
expressed as a minimum point (2003-07-15), value at the minimum
(15918.54) and value at the end of 2020 (18825.83). (This differs from
the value seen on Figure <a href="#fig:longterm">10</a>, as that also
includes the effect of the harmonic term.) The two parameters of the
harmonic regression can be expressed as an amplitude, a multiplier
(9.4%) and a phase shift (-0.7, i.e., minimum at week 32 of the year).

Figure <a href="#fig:longtermfit">11</a> shows the predictions of the
above model.

``` r
RawData2019$pred <- predict(fit)

ggplot(RawData2019, aes(x = Week, y = outcome)) +
  geom_line() + geom_line(aes(y = exp(pred)), color = "red") + facet_wrap(~Year)
```

![Figure 11: Weekly number of deaths in Germany, separated according to
year, showing the predictions from the model with quadratic long-term
trend and a single, fixed harmonic
term.](README_files/figure-gfm/longtermfit-1.png)

Figure: Weekly number of deaths in Germany, separated according to year,
showing the predictions from the model with quadratic long-term trend
and a single, fixed harmonic term.

A good fit can be observed, apart from summer and winter peaks. Thus, to
capture them, the predictions are subtracted; with the results shown on
Figure <a href="#fig:longtermfitresid">12</a>.

``` r
RawData2019$resid <- residuals(fit, type = "working")

peakdet <- scorepeak::detect_localmaxima(RawData2019$resid, 35) &
  scorepeak::score_type1(RawData2019$resid, 35)>0.1315
peaks <- data.table(x = RawData2019$NumTrend[peakdet], y = RawData2019$resid[peakdet])
peaks$peakDate <- lubridate::as_date(peaks$x)
peaks$Year <- lubridate::isoyear(peaks$peakDate)
peaks$peakWeek <- lubridate::isoweek(peaks$peakDate)
peaks$peakID <- 1:nrow(peaks)

ggplot(RawData2019, aes(x = Week, y = resid)) +  geom_line() + facet_wrap(~Year) +
  geom_point(data = peaks, aes(x = peakWeek, y = y)) + labs(y = "Working residual (log scale)")
```

![Figure 12: Residuals of the fitted model with quadratic long-term
trend and a single, fixed harmonic term. Dots indicate identified
peaks.](README_files/figure-gfm/longtermfitresid-1.png)

Figure: Residuals of the fitted model with quadratic long-term trend and
a single, fixed harmonic term. Dots indicate identified peaks.

Peaks in the residuals were identified with the peak detector of
Palshikar \[49\] using parameters that were empirically tuned to
identify to visually clear peaks. Results are shown on Figure
<a href="#fig:longtermfitresid">12</a> with black dots; an indeed good
identification of the unequivocal peaks can be seen.

Figure <a href="#fig:peakneigh">13</a> shows the peaks themselves with a
![\\pm](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpm "\pm")
100 days neighbourhood. This reinforces the idea that summer and winter
peaks are somewhat different, but more importantly, suggests that the
rescaled probability density function of the Cauchy distribution, i.e.,
![\\frac{a}{\\pi s}\\frac{1}{1+\\left(\\frac{x-x_0}{s}\\right)^2}+b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7Ba%7D%7B%5Cpi%20s%7D%5Cfrac%7B1%7D%7B1%2B%5Cleft%28%5Cfrac%7Bx-x_0%7D%7Bs%7D%5Cright%29%5E2%7D%2Bb "\frac{a}{\pi s}\frac{1}{1+\left(\frac{x-x_0}{s}\right)^2}+b")
might be a good – and parsimonious – function form to capture the shape
of the peaks.

``` r
RawData2019$peakID <- sapply(1:nrow(RawData2019),
                             function(i) which.min(abs(RawData2019$NumTrend[i]-peaks$x)))
RawData2019$peakdist <- sapply(1:nrow(RawData2019),
                               function(i) RawData2019$NumTrend[i]-peaks[RawData2019$peakID[i],]$x)

RawData2019 <- merge(RawData2019, peaks[, .(peakID, peakWeek)])
RawData2019$peakText <- paste0(RawData2019$peakID, " (Week: ", RawData2019$peakWeek, ")")
RawData2019$peakText <- factor(RawData2019$peakText, levels = unique(RawData2019$peakText))

minfun <- function(x, data) sum((data$resid-(x["a"]*dcauchy(data$peakdist, x["x0"], exp(x["s"]))+x["b"]))^2)

RawData2019[, c("fitpeak", "x0", "s", "a", "b") :=
              with(optim(c(a = 10, x0 = 0, s = 0, b = 0), minfun,
                         data = .SD[abs(peakdist)<100, .(resid = resid, peakdist)]),
                   list(dcauchy(peakdist, par["x0"], exp(par["s"]))*par["a"],
                        par["x0"], exp(par["s"]), par["a"],par["b"])), .(peakID)]

ggplot(RawData2019[abs(peakdist)<100], aes(x = peakdist, y = resid)) + geom_line() +
  geom_point() + facet_wrap(~peakText) + geom_line(aes(y = fitpeak + b), color = "red") +
  labs(x = "Distance from the peak [day]", y = "Working residual (log scale)")
```

![Figure 13: 100-day width neighbourhood of the identified peaks. Red
line indicates the best fitting rescaled Cauchy
density.](README_files/figure-gfm/peakneigh-1.png)

Figure: 100-day width neighbourhood of the identified peaks. Red line
indicates the best fitting rescaled Cauchy density.

To check this theory, the best fitting function was found for each peak
individually using the Nelder-Mead method \[50\] with mean squared error
objective function. Results are shown on <a href="#fig:peakneigh">13</a>
as red lines; an almost perfect fit can be observed for all peaks
confirming the initial idea of using Cauchy density.

This now puts us in a position to investigate the distribution of the
parameters (i.e.,
![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a"),
![b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b "b"),
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
and
![x_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0 "x_0")),
which is shown on Figure <a href="#fig:peakparamdist">14</a> separated
according to whether the peak is during the summer or not. Peak height
is also calculated, defined as height at zero (which is
![\\frac{a}{\\pi s \\left(1+\\frac{x_0^2}{s^2}\\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7Ba%7D%7B%5Cpi%20s%20%5Cleft%281%2B%5Cfrac%7Bx_0%5E2%7D%7Bs%5E2%7D%5Cright%29%7D "\frac{a}{\pi s \left(1+\frac{x_0^2}{s^2}\right)}"))
not the actual maximum height (which is at
![x_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0 "x_0"))
to avoid extremely large heights – which were never actually observed –
due to peaks with small
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s"),
i.e., very narrow peaks.

``` r
peaks <- merge(peaks, unique(RawData2019[, .(peakID, x0, s, a, b)]))
peaks$Summer <- peaks$peakWeek>10 & peaks$peakWeek<35
peaks$amplitude <- peaks$a/(pi*peaks$s*(1+peaks$x0^2/peaks$s^2))

fitWinterAmplitudeMin <- min(peaks[Summer==FALSE]$amplitude)
fitWinterAmplitudeMax <- max(peaks[Summer==FALSE]$amplitude)
fitWinterSigmaMin <- min(peaks[Summer==FALSE]$s)
fitWinterSigmaMax <- max(peaks[Summer==FALSE]$s)
fitWinterProb <- sum(!peaks$Summer)/(diff(range(RawData2019$Year)) + 1)

fitSummerAmplitudeMin <- min(peaks[Summer==TRUE]$amplitude)
fitSummerAmplitudeMax <- max(peaks[Summer==TRUE]$amplitude)
fitSummerSigmaMin <- min(peaks[Summer==TRUE]$s)
fitSummerSigmaMax <- max(peaks[Summer==TRUE]$s)
fitSummerProb <- sum(peaks$Summer)/(diff(range(RawData2019$Year)) + 1)

ggplot(melt(peaks[, .(peakID, Summer, s, a, peakWeek, amplitude)], id.vars = c("peakID", "Summer")),
       aes(x = value, y = Summer)) + facet_wrap(~variable, scales = "free") + geom_point()
```

![Figure 14: Distribution of the parameters of the best fitting rescaled
Cauchy densities for each peak, separated according to whether the peak
is during the summer.](README_files/figure-gfm/peakparamdist-1.png)

Figure: Distribution of the parameters of the best fitting rescaled
Cauchy densities for each peak, separated according to whether the peak
is during the summer.

This verifies that the width is indeed different, with the
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
of the summer peaks being below 10, and the winter peaks being above,
i.e., summer peaks are shorter in duration, raise and fall faster.
Interestingly, the peak heights are not substantially different between
winter and summer. Also note that the probability of having a peak at
all is different: there are 8 summer peaks and 9 winter peaks (from 20
years). Winter peaks occur between weeks 4 and 10, summer peaks occur
from weeks 25 to 35.

This allows the removal of the peaks (Figure
<a href="#fig:peaksremoved">15</a>), and, after these peaks are removed,
it is possible to re-estimate trend and seasonality, now without the
biasing effect of the peaks. This “bootstrap” procedure is adequate
after this second iteration, as no further peaks can be seen after the
removal of the re-estimated trend and seasonality.

``` r
ggplot(RawData2019, aes(x = Week, y = log(outcome))) + geom_line() + facet_wrap(~Year) +
  geom_line(aes(y = log(outcome) - fitpeak), color = "red") + labs(y = "Outcome (log scale)")
```

![Figure 15: Weekly number of deaths in Germany, separated according to
year (black) and with peaks removed
(red).](README_files/figure-gfm/peaksremoved-1.png)

Figure: Weekly number of deaths in Germany, separated according to year
(black) and with peaks removed (red).

Creating the appropriate model to simulate such peaks is not
straightforward: there is stochasticity in the position of the peaks and
in its shape, i.e., height and width. (Actually, even the presence of a
peak is stochastic.) The following procedure will be used: the presence
is generated as a Bernoulli random variate (with the probabilities
described above, different for summer and winter), the onset date is
uniformly distributed (from 0 to 0.2 in scaled weeks for the winter peak
and from 0.5 to 0.7 for the summer peak), i.e., the position itself is
random, but the parameters of the underyling distribution are fixed. The
![b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b "b")
parameter is set to zero (irrespectively of its estimated value, to
really capture only the peak, locally – a non-zero
![b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b "b")
would mean a non-local effect), while
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
and
![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a")
are randomly generated for each peak, again, separately for summer and
winter peaks. Given the high correlation between
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
and
![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a"),
not these, but rather
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
(width) and the amplitude will be generated as a random variate from –
independent – uniform distributions. The parameters (minimum and
maximum) of the uniform distributions both for
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
and the amplitude will be considered as a parameter (hyperparameter) of
the simulational procedure, just as the
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
probability of the Bernoulli distribution, all different for summer and
winter.

The following table summarizes the parameters:

``` r
fit <- mgcv::gam(I(log(outcome)-fitpeak) ~ poly(NumTrend, 2, raw = TRUE) + cos(2*pi*WeekScaled) +
                   sin(2*pi*WeekScaled), data = RawData2019,  method = "REML")

fitSeasonAmplitude <- sqrt(coef(fit)["sin(2 * pi * WeekScaled)"]^2 +
                             coef(fit)["cos(2 * pi * WeekScaled)"]^2)
fitSeasonPhase <- atan(-coef(fit)["sin(2 * pi * WeekScaled)"]/coef(fit)["cos(2 * pi * WeekScaled)"])

fittedpars <- setNames(c(coef(fit)[1:3], fitSeasonAmplitude, fitSeasonPhase,
                         fitWinterAmplitudeMin, fitWinterAmplitudeMax, fitWinterSigmaMin,
                         fitWinterSigmaMax, fitWinterProb,
                         fitSummerAmplitudeMin, fitSummerAmplitudeMax, fitSummerSigmaMin,
                         fitSummerSigmaMax, fitSummerProb),
                       c("TrendConst", "TrendLin", "TrendQuadr",
                         "SeasonAmplitude", "SeasonPhase",
                         "WinterAmplitudeMin", "WinterAmplitudeMax", "WinterSigmaMin",
                         "WinterSigmaMax", "WinterProb",
                         "SummerAmplitudeMin", "SummerAmplitudeMax", "SummerSigmaMin",
                         "SummerSigmaMax", "SummerProb"))
knitr::kable(data.table(
  Parameter = c("Linear term of the trend", "Constant term of the trend",
                "Quadratic term of the trend", "Amplitude of seasonality (log scale)",
                "Phase of seasonality", "Minimum of winter peak amplitude (log scale)",
                "Maximum of winter peak amplitude (log scale)", "Minimum of winter peak width",
                "Maximum of winter peak width", "Probability of winter peak",
                "Minimum of summer peak amplitude (log scale)",
                "Maximum of summer peak amplitude (log scale)",
                "Minimum of summer peak width", "Maximum of summer peak width",
                "Probability of summer peak"),
  Value = fittedpars), digits = 2)
```

| Parameter                                    | Value |
|:---------------------------------------------|------:|
| Linear term of the trend                     | 10.11 |
| Constant term of the trend                   |  0.00 |
| Quadratic term of the trend                  |  0.00 |
| Amplitude of seasonality (log scale)         |  0.07 |
| Phase of seasonality                         | -0.61 |
| Minimum of winter peak amplitude (log scale) |  0.11 |
| Maximum of winter peak amplitude (log scale) |  0.33 |
| Minimum of winter peak width                 |  8.41 |
| Maximum of winter peak width                 | 35.46 |
| Probability of winter peak                   |  0.45 |
| Minimum of summer peak amplitude (log scale) |  0.10 |
| Maximum of summer peak amplitude (log scale) |  0.24 |
| Minimum of summer peak width                 |  0.86 |
| Maximum of summer peak width                 |  9.24 |
| Probability of summer peak                   |  0.40 |

``` r
saveRDS(fittedpars, "fittedpars.rds")
```

Given the mechanism described above, the procedure to simulate synthetic
datasets can easily be created. Of course, as every calculation is
carried out on the log scale, the mean should be exponentiated at the
last step.

``` r
simdat <- function(TrendConst, TrendLin, TrendQuadr, SeasonAmplitude, SeasonPhase,
                   WinterAmplitudeMin, WinterAmplitudeMax, WinterSigmaMin, WinterSigmaMax, WinterProb,
                   SummerAmplitudeMin, SummerAmplitudeMax, SummerSigmaMin, SummerSigmaMax, SummerProb) {
  
  SimData <- data.frame(date = seq(as.Date("2000-01-03"), as.Date("2023-12-25"), by = 7))
  SimData$NumTrend <- as.numeric(SimData$date)
  SimData$WeekScaled <- lubridate::isoweek(SimData$date)/
    lubridate::isoweek(paste0(lubridate::isoyear(SimData$date), "-12-28"))
  
  SimData$logmu <- TrendConst + TrendLin*SimData$NumTrend + TrendQuadr*SimData$NumTrend^2 +
    SeasonAmplitude*cos(SimData$WeekScaled*2*pi + SeasonPhase)
  
  SimData$logmu <- SimData$logmu + rowSums(sapply(2000:2023, function(y) {
    if(rbinom(1, 1, WinterProb)==0) return(rep(0, nrow(SimData))) else {
      amplitude <- runif(1, WinterAmplitudeMin, WinterAmplitudeMax)
      sigm <- runif(1, WinterSigmaMin, WinterSigmaMax)
      dcauchy(SimData$NumTrend,
              as.numeric((as.Date(paste0(y, "-01-01")) + runif(1, 0, 0.2)*7*52.25)), 
              sigm)*(pi*sigm*amplitude)
    }
  }))
  
  SimData$logmu <- SimData$logmu + rowSums(sapply(2000:2023, function(y) {
    if(rbinom(1, 1, SummerProb)==0) return(rep(0, nrow(SimData))) else {
      amplitude <- runif(1, SummerAmplitudeMin, SummerAmplitudeMax)
      sigm <- runif(1, SummerSigmaMin, SummerSigmaMax)
      dcauchy(SimData$NumTrend,
              as.numeric((as.Date(paste0(y, "-01-01")) + runif(1, 0.5, 0.7)*7*52.25)),
              sigm)*(pi*sigm*amplitude)
    }
  }))
  
  SimData$outcome <- rnbinom(nrow(SimData), mu = exp(SimData$logmu), size = 1000)
  
  SimData
}
```

Figure <a href="#fig:simillustration">16</a> illustrates the synthetic
data set creation with a single simulation. (Of course, to assess the
correctness of the simulation, several realizations have to be
inspected.) In addition to the plots of
<a href="#fig:german-raw-plots">2</a>, it also gives the autocorrelation
function so that it can also be compared. (The simulated outcomes are
themselves independent – meaning that effects like the increased
probability of a flu season if the previous year did not have one are
neglected –, but the trend, seasonality and the peaks induce a
correlation structure.)

``` r
set.seed(1)

SimData <- as.data.table(do.call(simdat, as.list(fittedpars)))
SimData$Type <- "Simulated"
SimData$Week <- lubridate::isoweek(SimData$date)
SimData$Year <- lubridate::isoyear(SimData$date)

SimDataYearly <- SimData[, .(outcome = sum(outcome)), .(Year, Type)]

p1 <- ggplot(rbind(SimData[Year<=2019], RawData[Year<=2019]),
             aes(x = date, y = outcome, group = Type, color = Type)) + geom_line()

p2 <- ggplot(rbind(SimDataYearly[Year<=2019], RawDataYearly[Year<=2019]),
             aes(x = Year, y = outcome, group = Type, color = Type)) + geom_point() + geom_line()

p3 <- ggplot(rbind(SimData[Year<=2019], RawData[Year<=2019]),
             aes(x = Week, y = outcome, group = Year)) + facet_wrap(~Type) + geom_line(alpha = 0.2)

p4 <- ggplot(rbind(with(acf(RawData$outcome, plot = FALSE),
                        data.table(Type = "Actual", acf = acf[, 1, 1], lag = lag[, 1, 1])),
             with(acf(SimData$outcome, plot = FALSE),
                  data.table(Type = "Simulated", acf = acf[, 1, 1], lag = lag[, 1, 1]))),
       aes(x = lag - 1/4 + as.numeric(Type=="Simulated")/2,
           xend = lag - 1/4 + as.numeric(Type=="Simulated")/2, y = acf, yend = 0, color = Type)) +
  geom_line() + geom_point() + geom_hline(yintercept = 0, col = "black") +
  labs(x = "Lag", y = "Autocorrelation")

egg::ggarrange(p1, p2, p3, p4, ncol = 1)
```

![Figure 16: From top to bottom: weekly mortalities, yearly mortalities,
seasonal pattern and autocorrelation function of the actual German
mortality data and a single simulated dataset,
2000-2019.](README_files/figure-gfm/simillustration-1.png)

Figure: From top to bottom: weekly mortalities, yearly mortalities,
seasonal pattern and autocorrelation function of the actual German
mortality data and a single simulated dataset, 2000-2019.

Several simulations confirm an overall good fit between the actual data
set and the simulated ones. Thus, it is now possible to investigate the
properties of the mortality prediction on algorithms using simulated
datasets, where the actual outcome is known, and the parameters can be
varied.

## Supplementary Material 3: Validation through simulation

Two sets of parameters have to be set up: parameters of the simulation
(i.e., parameters of the scenario, on which the methods will be run) and
parameters of the methods. They’re set up as given in the main text.

``` r
pargridSim <- rbind(as.data.table(t(fittedpars)),
                    rbindlist(lapply(1:length(fittedpars), function(i) {
                      temp <- as.data.table(t(fittedpars))[rep(1, 5)]
                      temp[[names(fittedpars)[i]]] <- seq(fittedpars[i]*if(i==1) 0.9 else 0.5,
                                                          fittedpars[i]*if(i==1) 1.1 else 2,
                                                          length.out = 5)
                      temp
                    })))
pargridSim$SummerProb[pargridSim$SummerProb>1] <- 1
pargridSim$WinterProb[pargridSim$WinterProb>1] <- 1

pargridWHO <- expand.grid(startyear = c(2000, 2005, 2010, 2015), k = c(5, 10, 15, 20))
pargridWHO$parmethod <- paste0("WHO", seq_len(nrow(pargridWHO)))
pargridAI <- expand.grid(startyear = c(2000, 2005, 2010, 2015), invtkpy = c(4, 5, 7, 12))
pargridAI <- merge(pargridAI, data.table(invtkpy = c(4, 5, 7, 12),
                                         tkpy = c("1/4", "1/5", "1/7", "1/12")))
pargridAI$parmethod <- paste0("AI", seq_len(nrow(pargridAI)))
pargridAI$tkpy <- factor(pargridAI$tkpy, levels = c("1/12", "1/7", "1/5", "1/4"))
pargridAverage <- data.frame(startyear = c(2000, 2005, 2010, 2015, 2019))
pargridAverage$parmethod <- paste0("Average", seq_len(nrow(pargridAverage)))
pargridLin <- data.frame(startyear = c(2000, 2005, 2010, 2015))
pargridLin$parmethod <- paste0("Lin", seq_len(nrow(pargridLin)))
```

One thousand simulation will be run for each parameter of the scenario,
and for each of the 1000 simulated data set all 4 methods with all
possible parameters of the methods will be evaluated. To increase the
speed, simulations will be run in parallel. (The problem is
embarrassingly parallel, as different simulations are completely
independent of each other \[51\].)

``` r
if(!file.exists("predLongs.rds")) {
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c("simdat", "pargridSim", "pargridWHO", "pargridAI",
                                "pargridAverage", "pargridLin"))
  
  for(r in 1:10) {
    pred <- do.call(rbind, parallel::parLapply(cl, 1:100, function(j) {
      do.call(rbind, lapply(1:nrow(pargridSim), function(k) {
        SimData <- do.call(simdat, as.list(pargridSim[k, ]))
        SimData$Year <- lubridate::isoyear(SimData$date)
        
        predWHO <- sapply(1:nrow(pargridWHO), function(i)
          predict(mgcv::gam(outcome ~ s(NumTrend, k = pargridWHO$k[i]) + s(WeekScaled, bs = "cc"),
                            data = SimData[SimData$Year>=pargridWHO$startyear[i]&SimData$Year<=2019,],
                            family = mgcv::nb(), method = "REML"),
                  newdata = SimData[SimData$Year>=2020,], type = "response"))
        
        predAI <- sapply(1:nrow(pargridAI), function(i)
          with(excessmort::compute_expected(
            cbind(SimData[SimData$Year>=pargridAI$startyear[i],], population = 1),
            exclude = seq(as.Date("2020-01-01"), max(SimData$date), by = "day"),
            frequency = nrow(SimData)/(as.numeric(diff(range(SimData$date)))/365.25),
            trend.knots.per.year = 1/pargridAI$invtkpy[i], verbose = FALSE),
            expected[date>=as.Date("2019-12-30")]))
        
        predAverage <- sapply(1:nrow(pargridAverage), function(i)
          predict(mgcv::gam(outcome ~ s(WeekScaled, bs = "cc"),
                            data = SimData[SimData$Year>=pargridAverage$startyear[i]&SimData$Year<=2019,],
                            family = mgcv::nb(), method = "REML"),
                  newdata = SimData[SimData$Year>=2020,], type = "response"))
        
        predLin <- sapply(1:nrow(pargridLin), function(i)
          predict(mgcv::gam(outcome ~ NumTrend + s(WeekScaled, bs = "cc"),
                            data = SimData[SimData$Year>=pargridLin$startyear[i]&SimData$Year<=2019,],
                            family = mgcv::nb(), method = "REML"),
                  newdata = SimData[SimData$Year>=2020,], type = "response"))
        
        setNames(data.frame(j, k, SimData[SimData$Year>=2020, c("date", "outcome")], predWHO,
                            predAI, predAverage, predLin, row.names = NULL),
                 c("rep", "parsim", "date", "outcome", pargridWHO$variable, pargridAI$variable,
                   pargridAverage$variable, pargridLin$variable))
      }))
    }))
    
    pred$rep <- (r-1)*100 + pred$rep
    saveRDS(pred, paste0("pred_", r, ".rds"))
  }
  
  parallel::stopCluster(cl)
  
  pred <- rbindlist(lapply(1:10, function(r) readRDS(paste0("pred_", r, ".rds"))))
  saveRDS(pred, "pred.rds")

  pred$Year <- lubridate::isoyear(pred$date)
  pred$date <- NULL
  
  predYearly <- pred[, lapply(.SD, sum), .(rep, parsim, Year)]
  
  predLongs <- lapply(
    list(WHO = pargridWHO, AI = pargridAI, Average = pargridAverage,
         Lin = pargridLin),
    function(pg)
      merge(merge(melt(predYearly[, c("rep", "parsim", "outcome", "Year", pg$parmethod), with = FALSE],
                       id.vars = c("rep", "parsim", "outcome", "Year"), variable.name = "parmethod"), pg),
            data.table(parsim = 1:76,
                       parsimName = factor(c("Base", rep(names(fittedpars), each = 5)),
                                           levels = c("Base", names(fittedpars))),
                       parsimValue = c("Base", rep(paste0("#", 1:5), length(fittedpars)))), by = "parsim"))
  
  saveRDS(predLongs, "predLongs.rds")
} else predLongs <- readRDS("predLongs.rds")
```

## Acknowledgement

On behalf of Project KOMPLEXEPI we thank for the usage of ELKH Cloud
(<https://science-cloud.hu/>) that significantly helped us achieving the
results published in this paper.

## References

<div id="refs" class="references csl-bib-body">

<div id="ref-karlinskyTrackingExcessMortality2021b" class="csl-entry">

1\. Karlinsky A, Kobak D. Tracking excess mortality across countries
during the COVID-19 pandemic with the World Mortality Dataset. eLife
\[Internet\]. 2021 \[cited 2022 Jun 28\];10:e69336. Available from:
<https://elifesciences.org/articles/69336>

</div>

<div id="ref-santos-burgoaDifferentialPersistentRisk2018"
class="csl-entry">

2\. Santos-Burgoa C, Sandberg J, Suárez E, Goldman-Hawes A, Zeger S,
Garcia-Meza A, et al. Differential and persistent risk of excess
mortality from Hurricane Maria in Puerto Rico: A time-series analysis.
The Lancet Planetary Health \[Internet\]. 2018 \[cited 2022 Jul
13\];2:e478–88. Available from:
<https://linkinghub.elsevier.com/retrieve/pii/S2542519618302092>

</div>

<div id="ref-riveraModelingExcessDeaths2019" class="csl-entry">

3\. Rivera R, Rolke W. Modeling excess deaths after a natural disaster
with application to Hurricane Maria. Statistics in Medicine
\[Internet\]. 2019 \[cited 2022 Jul 13\];38:4545–54. Available from:
<https://onlinelibrary.wiley.com/doi/10.1002/sim.8314>

</div>

<div id="ref-moritaExcessMortalityDue2017" class="csl-entry">

4\. Morita T, Nomura S, Tsubokura M, Leppold C, Gilmour S, Ochi S, et
al. Excess mortality due to indirect health effects of the 2011 triple
disaster in Fukushima, Japan: A retrospective observational study. J
Epidemiol Community Health \[Internet\]. 2017 \[cited 2022 Jul
13\];71:974–80. Available from:
<https://jech.bmj.com/lookup/doi/10.1136/jech-2016-208652>

</div>

<div id="ref-simonsenImpactInfluenzaEpidemics1997" class="csl-entry">

5\. Simonsen L, Clarke MJ, Williamson GD, Stroup DF, Arden NH,
Schonberger LB. The impact of influenza epidemics on mortality:
Introducing a severity index. Am J Public Health \[Internet\]. 1997
\[cited 2022 Jul 13\];87:1944–50. Available from:
<https://ajph.aphapublications.org/doi/full/10.2105/AJPH.87.12.1944>

</div>

<div id="ref-rosanoInvestigatingImpactInfluenza2019" class="csl-entry">

6\. Rosano A, Bella A, Gesualdo F, Acampora A, Pezzotti P, Marchetti S,
et al. Investigating the impact of influenza on excess mortality in all
ages in Italy during recent seasons (2013/14–2016/17 seasons).
International Journal of Infectious Diseases \[Internet\]. 2019 \[cited
2022 Jul 13\];88:127–34. Available from:
<https://linkinghub.elsevier.com/retrieve/pii/S1201971219303285>

</div>

<div id="ref-leonCOVID19NeedRealtime2020a" class="csl-entry">

7\. Leon DA, Shkolnikov VM, Smeeth L, Magnus P, Pechholdová M, Jarvis
CI. COVID-19: A need for real-time monitoring of weekly excess deaths.
The Lancet \[Internet\]. 2020 \[cited 2022 Jul 13\];395:e81. Available
from: <https://linkinghub.elsevier.com/retrieve/pii/S0140673620309338>

</div>

<div id="ref-pearceComparisonsCountriesAre2020" class="csl-entry">

8\. Pearce N, Lawlor DA, Brickley EB. Comparisons between countries are
essential for the control of COVID-19. International Journal of
Epidemiology \[Internet\]. 2020 \[cited 2022 Jul 13\];49:1059–62.
Available from: <https://academic.oup.com/ije/article/49/4/1059/5864919>

</div>

<div id="ref-beaneyExcessMortalityGold2020" class="csl-entry">

9\. Beaney T, Clarke JM, Jain V, Golestaneh AK, Lyons G, Salman D, et
al. Excess mortality: The gold standard in measuring the impact of
COVID-19 worldwide? J R Soc Med \[Internet\]. 2020 \[cited 2022 Jul
13\];113:329–34. Available from:
<http://journals.sagepub.com/doi/10.1177/0141076820956802>

</div>

<div id="ref-faustAllCauseExcessMortality2021" class="csl-entry">

10\. Faust JS, Krumholz HM, Du C, Mayes KD, Lin Z, Gilman C, et al.
All-Cause Excess Mortality and COVID-19–Related Mortality Among US
Adults Aged 25-44 Years, March-July 2020. JAMA \[Internet\]. 2021
\[cited 2022 Jul 14\];325:785. Available from:
<https://jamanetwork.com/journals/jama/fullarticle/2774445>

</div>

<div id="ref-kirpichExcessMortalityBelarus2022" class="csl-entry">

11\. Kirpich A, Shishkin A, Weppelmann TA, Tchernov AP, Skums P, Gankin
Y. Excess mortality in Belarus during the COVID-19 pandemic as the case
study of a country with limited non-pharmaceutical interventions and
limited reporting. Sci Rep \[Internet\]. 2022 \[cited 2022 Jul
14\];12:5475. Available from:
<https://www.nature.com/articles/s41598-022-09345-z>

</div>

<div id="ref-faustExcessMortalityMassachusetts2022" class="csl-entry">

12\. Faust JS, Du C, Liang C, Mayes KD, Renton B, Panthagani K, et al.
Excess Mortality in Massachusetts During the Delta and Omicron Waves of
COVID-19. JAMA \[Internet\]. 2022 \[cited 2022 Jul 14\];328:74.
Available from:
<https://jamanetwork.com/journals/jama/fullarticle/2792738>

</div>

<div id="ref-rossenDisparitiesExcessMortality2021" class="csl-entry">

13\. Rossen LM, Ahmad FB, Anderson RN, Branum AM, Du C, Krumholz HM, et
al. Disparities in Excess Mortality Associated with COVID-19 — United
States, 2020. MMWR Morb Mortal Wkly Rep \[Internet\]. 2021 \[cited 2022
Jul 14\];70:1114–9. Available from:
<http://www.cdc.gov/mmwr/volumes/70/wr/mm7033a2.htm?s_cid=mm7033a2_w>

</div>

<div id="ref-bradshawTrackingMortalityReal2021" class="csl-entry">

14\. Bradshaw D, Dorrington RE, Laubscher R, Moultrie TA, Groenewald P.
Tracking mortality in near to real time provides essential information
about the impact of the COVID-19 pandemic in South Africa in 2020. S Afr
Med J \[Internet\]. 2021 \[cited 2022 Jul 14\];111:732. Available from:
<http://www.samj.org.za/index.php/samj/article/view/13304>

</div>

<div id="ref-modiEstimatingCOVID19Mortality2021" class="csl-entry">

15\. Modi C, Böhm V, Ferraro S, Stein G, Seljak U. Estimating COVID-19
mortality in Italy early in the COVID-19 pandemic. Nat Commun
\[Internet\]. 2021 \[cited 2022 Jul 14\];12:2729. Available from:
<http://www.nature.com/articles/s41467-021-22944-0>

</div>

<div id="ref-wangEstimatingExcessMortality2022a" class="csl-entry">

16\. Wang H, Paulson KR, Pease SA, Watson S, Comfort H, Zheng P, et al.
Estimating excess mortality due to the COVID-19 pandemic: A systematic
analysis of <span class="nocase">COVID-19-related</span> mortality,
2020–21. The Lancet \[Internet\]. 2022 \[cited 2022 Jul
14\];399:1513–36. Available from:
<https://linkinghub.elsevier.com/retrieve/pii/S0140673621027963>

</div>

<div id="ref-kontisMagnitudeDemographicsDynamics2020a"
class="csl-entry">

17\. Kontis V, Bennett JE, Rashid T, Parks RM, Pearson-Stuttard J,
Guillot M, et al. Magnitude, demographics and dynamics of the effect of
the first wave of the COVID-19 pandemic on all-cause mortality in 21
industrialized countries. Nat Med \[Internet\]. 2020 \[cited 2022 Jul
14\];26:1919–28. Available from:
<https://www.nature.com/articles/s41591-020-1112-0>

</div>

<div id="ref-stangExcessMortalityDue2020" class="csl-entry">

18\. Stang A, Standl F, Kowall B, Brune B, Böttcher J, Brinkmann M, et
al. Excess mortality due to COVID-19 in Germany. Journal of Infection
\[Internet\]. 2020 \[cited 2022 Jul 16\];81:797–801. Available from:
<https://linkinghub.elsevier.com/retrieve/pii/S016344532030596X>

</div>

<div id="ref-konstantinoudisRegionalExcessMortality2022"
class="csl-entry">

19\. Konstantinoudis G, Cameletti M, Gómez-Rubio V, Gómez IL, Pirani M,
Baio G, et al. Regional excess mortality during the 2020 COVID-19
pandemic in five European countries. Nat Commun \[Internet\]. 2022
\[cited 2022 Jul 14\];13:482. Available from:
<https://www.nature.com/articles/s41467-022-28157-3>

</div>

<div id="ref-daviesCommunityFactorsExcess2021" class="csl-entry">

20\. Davies B, Parkes BL, Bennett J, Fecht D, Blangiardo M, Ezzati M, et
al. Community factors and excess mortality in first wave of the COVID-19
pandemic in England. Nat Commun \[Internet\]. 2021 \[cited 2022 Jul
14\];12:3755. Available from:
<http://www.nature.com/articles/s41467-021-23935-x>

</div>

<div id="ref-joyExcessMortalityFirst2020" class="csl-entry">

21\. Joy M, Hobbs FR, Bernal JL, Sherlock J, Amirthalingam G, McGagh D,
et al. Excess mortality in the first COVID pandemic peak:
Cross-sectional analyses of the impact of age, sex, ethnicity, household
size, and long-term conditions in people of known SARS-CoV-2 status in
England. Br J Gen Pract \[Internet\]. 2020 \[cited 2022 Jul
14\];70:e890–8. Available from:
<https://bjgp.org/lookup/doi/10.3399/bjgp20X713393>

</div>

<div id="ref-alicandroItalyFirstWave2020" class="csl-entry">

22\. Alicandro G, Remuzzi G, La Vecchia C. Italy’s first wave of the
COVID-19 pandemic has ended: No excess mortality in May, 2020. The
Lancet \[Internet\]. 2020 \[cited 2022 Jul 13\];396:e27–8. Available
from: <https://linkinghub.elsevier.com/retrieve/pii/S0140673620318651>

</div>

<div id="ref-haklaiExcessMortalityCOVID192021" class="csl-entry">

23\. Haklai Z, Aburbeh M, Goldberger N, Gordon E-S. Excess mortality
during the COVID-19 pandemic in Israel, March–November 2020: When,
where, and for whom? Isr J Health Policy Res \[Internet\]. 2021 \[cited
2022 Jul 14\];10:17. Available from:
<https://ijhpr.biomedcentral.com/articles/10.1186/s13584-021-00450-4>

</div>

<div id="ref-modigExcessMortalityCOVID192021" class="csl-entry">

24\. Modig K, Ahlbom A, Ebeling M. Excess mortality from COVID-19:
Weekly excess death rates by age and sex for Sweden and its most
affected region. European Journal of Public Health \[Internet\]. 2021
\[cited 2022 Jul 14\];31:17–22. Available from:
<https://academic.oup.com/eurpub/article/31/1/17/5968985>

</div>

<div id="ref-kriegerExcessMortalityMen2020" class="csl-entry">

25\. Krieger N, Chen JT, Waterman PD. Excess mortality in men and women
in Massachusetts during the COVID-19 pandemic. The Lancet \[Internet\].
2020 \[cited 2022 Jul 14\];395:1829. Available from:
<https://linkinghub.elsevier.com/retrieve/pii/S0140673620312344>

</div>

<div id="ref-michelozziTemporalDynamicsTotal2020" class="csl-entry">

26\. Michelozzi P, de’Donato F, Scortichini M, Pezzotti P, Stafoggia M,
De Sario M, et al. Temporal dynamics in total excess mortality and
COVID-19 deaths in Italian cities. BMC Public Health \[Internet\]. 2020
\[cited 2022 Jul 14\];20:1238. Available from:
<https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-020-09335-8>

</div>

<div id="ref-vieiraRapidEstimationExcess2020" class="csl-entry">

27\. Vieira A, Peixoto VR, Aguiar P, Abrantes A. Rapid Estimation of
Excess Mortality during the COVID-19 Pandemic in Portugal -Beyond
Reported Deaths: JEGH \[Internet\]. 2020 \[cited 2022 Jul 14\];10:209.
Available from: <https://www.atlantis-press.com/article/125941684>

</div>

<div id="ref-StatisticsEurostat" class="csl-entry">

28\. Statistics \| Eurostat \[Internet\]. \[cited 2022 Jul 14\].
Available from:
<https://ec.europa.eu/eurostat/databrowser/view/demo_mexrt/default/table?lang=en>

</div>

<div id="ref-ghafariExcessDeathsAssociated2021" class="csl-entry">

29\. Ghafari M, Kadivar A, Katzourakis A. Excess deaths associated with
the Iranian COVID-19 epidemic: A province-level analysis. International
Journal of Infectious Diseases \[Internet\]. 2021 \[cited 2022 Jul
13\];107:101–15. Available from:
<https://linkinghub.elsevier.com/retrieve/pii/S120197122100326X>

</div>

<div id="ref-kobakExcessMortalityReveals2021" class="csl-entry">

30\. Kobak D. Excess mortality reveals Covid’s true toll in Russia.
Significance \[Internet\]. 2021 \[cited 2022 Jul 13\];18:16–9. Available
from: <https://onlinelibrary.wiley.com/doi/10.1111/1740-9713.01486>

</div>

<div id="ref-scortichiniExcessMortalityCOVID192021" class="csl-entry">

31\. Scortichini M, Schneider dos Santos R, De’ Donato F, De Sario M,
Michelozzi P, Davoli M, et al. Excess mortality during the COVID-19
outbreak in Italy: A two-stage interrupted time-series analysis.
International Journal of Epidemiology \[Internet\]. 2021 \[cited 2022
Jul 14\];49:1909–17. Available from:
<https://academic.oup.com/ije/article/49/6/1909/5923437>

</div>

<div id="ref-woodGeneralizedAdditiveModels2017b" class="csl-entry">

32\. Wood SN. Generalized additive models: An introduction with R.
Second edition. Boca Raton: CRC Press/Taylor & Francis Group; 2017.

</div>

<div id="ref-jaroslawSemiparametricRegression2018" class="csl-entry">

33\. Jaroslaw H, David R, Matt W. Semiparametric regression with R. New
York, NY: Springer Science+Business Media; 2018.

</div>

<div id="ref-acostaFlexibleStatisticalFramework2022" class="csl-entry">

34\. Acosta RJ, Irizarry RA. A Flexible Statistical Framework for
Estimating Excess Mortality. Epidemiology \[Internet\]. 2022 \[cited
2022 Jul 13\];33:346–53. Available from:
<https://journals.lww.com/10.1097/EDE.0000000000001445>

</div>

<div id="ref-islamExcessDeathsAssociated2021a" class="csl-entry">

35\. Islam N, Shkolnikov VM, Acosta RJ, Klimkin I, Kawachi I, Irizarry
RA, et al. Excess deaths associated with covid-19 pandemic in 2020: Age
and sex disaggregated time series analysis in 29 high income countries.
BMJ \[Internet\]. 2021 \[cited 2022 Jul 13\];n1137. Available from:
<https://www.bmj.com/lookup/doi/10.1136/bmj.n1137>

</div>

<div id="ref-riveraExcessMortalityUnited2020" class="csl-entry">

36\. Rivera R, Rosenbaum JE, Quispe W. Excess mortality in the United
States during the first three months of the COVID-19 pandemic. Epidemiol
Infect \[Internet\]. 2020 \[cited 2022 Jul 14\];148:e264. Available
from:
<https://www.cambridge.org/core/product/identifier/S0950268820002617/type/journal_article>

</div>

<div id="ref-knutsonEstimatingGlobalCountrySpecific2022"
class="csl-entry">

37\. Knutson V, Aleshin-Guendel S, Karlinsky A, Msemburi W, Wakefield J.
Estimating Global and Country-Specific Excess Mortality During the
COVID-19 Pandemic. arXiv; 2022 \[cited 2022 Jul 14\]; Available from:
<https://arxiv.org/abs/2205.09081>

</div>

<div id="ref-adam15MillionPeople2022" class="csl-entry">

38\. Adam D. 15 million people have died in the pandemic, WHO says.
Nature \[Internet\]. 2022 \[cited 2022 Jul 15\];605:206–6. Available
from: <https://www.nature.com/articles/d41586-022-01245-6>

</div>

<div id="ref-vannoordenCOVIDDeathTolls2022" class="csl-entry">

39\. Van Noorden R. COVID death tolls: Scientists acknowledge errors in
WHO estimates. Nature \[Internet\]. 2022 \[cited 2022 Jul
15\];606:242–4. Available from:
<https://www.nature.com/articles/d41586-022-01526-0>

</div>

<div id="ref-EurostatDataExplorer" class="csl-entry">

40\. Eurostat - Data Explorer \[Internet\]. \[cited 2022 Jun 28\].
Available from:
<https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=demo_r_mwk_ts&lang=en>

</div>

<div id="ref-rcoreteamLanguageEnvironmentStatistical2022"
class="csl-entry">

41\. R Core Team. R: A Language and Environment for Statistical
Computing \[Internet\]. Vienna, Austria: R Foundation for Statistical
Computing; 2022. Available from: <https://www.R-project.org/>

</div>

<div id="ref-dowleDataTableExtension2021" class="csl-entry">

42\. Dowle M, Srinivasan A. Data.table: Extension of ‘data.frame‘
\[Internet\]. 2021. Available from:
<https://CRAN.R-project.org/package=data.table>

</div>

<div id="ref-wickhamGgplot2ElegantGraphics2016" class="csl-entry">

43\. Wickham H. Ggplot2: Elegant Graphics for Data Analysis
\[Internet\]. Springer-Verlag New York; 2016. Available from:
<https://ggplot2.tidyverse.org>

</div>

<div id="ref-internationalorganizationforstandardizationISO86012019"
class="csl-entry">

44\. International Organization for Standardization. ISO 8601. 2019.

</div>

<div id="ref-nepomucenoSensitivityAnalysisExcess2022" class="csl-entry">

45\. Nepomuceno MR, Klimkin I, Jdanov DA, Alustiza‐Galarza A, Shkolnikov
VM. Sensitivity Analysis of Excess Mortality due to the COVID‐19
Pandemic. Population & Development Rev \[Internet\]. 2022 \[cited 2022
Jul 17\];48:279–302. Available from:
<https://onlinelibrary.wiley.com/doi/10.1111/padr.12475>

</div>

<div id="ref-scholeyRobustnessBiasEuropean2021" class="csl-entry">

46\. Schöley J. Robustness and bias of European excess death estimates
in 2020 under varying model specifications \[Internet\]. Epidemiology;
2021 Jun. Available from:
<http://medrxiv.org/lookup/doi/10.1101/2021.06.04.21258353>

</div>

<div id="ref-jdanovShorttermMortalityFluctuation2021a"
class="csl-entry">

47\. Jdanov DA, Galarza AA, Shkolnikov VM, Jasilionis D, Németh L, Leon
DA, et al. The short-term mortality fluctuation data series, monitoring
mortality shocks across time and space. Sci Data \[Internet\]. 2021
\[cited 2022 Jun 28\];8:235. Available from:
<https://www.nature.com/articles/s41597-021-01019-1>

</div>

<div id="ref-hilbeNegativeBinomialRegression2011" class="csl-entry">

48\. Hilbe JM. Negative Binomial Regression \[Internet\]. 2nd ed.
Cambridge University Press; 2011 \[cited 2022 Jul 15\]. Available from:
<https://www.cambridge.org/core/product/identifier/9780511973420/type/book>

</div>

<div id="ref-palshikarSimpleAlgorithmsPeak2009" class="csl-entry">

49\. Palshikar G. Simple algorithms for peak detection in time-series.
Proc 1st Int Conf Advanced Data Analysis, Business Analytics and
Intelligence. 2009.

</div>

<div id="ref-nelderSimplexMethodFunction1965" class="csl-entry">

50\. Nelder JA, Mead R. A Simplex Method for Function Minimization. The
Computer Journal \[Internet\]. 1965 \[cited 2022 Jul 15\];7:308–13.
Available from:
<https://academic.oup.com/comjnl/article-lookup/doi/10.1093/comjnl/7.4.308>

</div>

<div id="ref-matloffArtProgrammingTour2011" class="csl-entry">

51\. Matloff NS. The art of R programming: Tour of statistical software
design. San Francisco: No Starch Press; 2011.

</div>

</div>

[^1]: Physiological Controls Research Center, Óbuda University and
    Department of Statistics, Corvinus University of Budapest,
    <ferenci.tamas@nik.uni-obuda.hu>
